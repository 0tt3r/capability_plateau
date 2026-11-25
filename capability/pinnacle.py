"""Pinnacle fault tolerance (draft Sec. 3.6) for QUOPS, in the form of the surface-code model of surface.py.

Architecture [paper: Webster et al., arXiv:2602.11457v2]: one processing unit of beta = ceil(w/k) bridged generalised-
bicycle (GB) blocks holds the w working logical qubits; there is no separate memory. Every allocated block stays live.
Only the all-to-all hardware classes H2 and H4 are modeled: the paper's block-scale reconfiguration has not been
compiled onto the repository's nearest-neighbour grid, so H1 and H3 return no design (not zero capability).

Schedule, in code cycles t_c, by the same pauli-based compilation as the surface-code model: each synthesized T gate is
one magic-state-assisted logical measurement, the QUOPS Cliffords update a classical frame, and the w transformed final
measurements run after the N_T consuming measurements. A slot lasts D = max(d_t, reaction time, supply slot) code cycles:
  T_run = max[N_T D, ceil(N_T / F) t_supply] + (w + 1) D                 (cycles; 1 startup slot, w readouts)
where F identical supplies each deliver a T state every t_supply code cycles on average, rejection included.

Loss (minus the log no-fault probability), the form of draft Eq. (surface-code-failure):
  blocks      eta beta (T_run / d_t) Lambda(P_PB),  P_PB = 6.2 (p / 0.0158)^(d/2 + 0.47) per whole block per d_t cycles
  magic       N_T [Lambda(P_supply) + Lambda(P_port)]        (supply output error, plus the interface to the unit)
  synthesis   R eps^2 or R eps, shared with the stochastic loss exactly as in surface.py
with Lambda(x) = -ln(1 - x). Rejections stall the unit; the live blocks are charged for the mean stall. This substitutes
the mean for the geometric stall distribution (a relative difference below 1e-9 in the loss for the designs checked); the engine blocks
themselves (d = 24) contribute below 1e-13 per period and are not charged, as in the paper's black-box treatment.

Search: exhaustive over the GB block code, the supply (capability/pinnacle_factories.py), the T count per rotation, and
the least number F of supplies that meets the run-time limit, for the fewest physical qubits
  N = K (beta n_pb + F N_supply + N_port).
The qubits needed never decrease with the width, so best_at bisects exactly as in surface.py.
"""
from dataclasses import dataclass, replace
from math import ceil, expm1, log1p

import numpy as np

from . import pinnacle_factories as pf
from .config import N_MAX, WEEK_SECONDS
from .quops import ENSEMBLES, layer_counts, size
from .surface import EPS, FIT, NT


@dataclass(frozen=True)
class Code:
    n: int
    k: int
    d: int
    dt: int               # native logical period in code cycles
    qubits: int           # n_pb: check qubits, four gadgets, and four bridges included
    error: float = None   # failure per block per native period; None uses the original-Pinnacle fit


# [paper] Pinnacle v2, Table 1. Do not add the checks, gadgets, or bridges again.
CODES = (Code(30, 8, 4, 6, 140), Code(62, 10, 6, 8, 244), Code(126, 12, 10, 12, 452),
         Code(254, 14, 16, 18, 860), Code(510, 16, 24, 26, 1620))

# [paper] Frozen Pinnacle (arXiv:2609.21249v1): the [[510,16,24]] block with eight processing logical qubits, a 30-cycle
# logical period, the full 2522-qubit block including its unused memory-access gadget, and the summed Frozen Shift and
# Bridged Measurement fits per block per period (draft Eq. frozen-error). Used for Frozen with adapter-delivered T states.
P_FROZEN = 4.01*(1e-3/0.0340)**12 + 5.49*(1e-3/0.0284)**12
FROZEN_CODES = (Code(510, 8, 24, 30, 2522, P_FROZEN),)


@dataclass(frozen=True)
class Choices:
    eta: float = 2.0                   # [choice] extend the X-observable fit of the paper to both Pauli sectors
    reaction: int = 10                 # [paper] t_r / t_c
    codes: tuple = (4, 6, 10, 16, 24)  # GB code distances allowed
    families: tuple = ('surface', 'surface_jws', 'reference', 'gb_cult', 'gb_inject', 'gb_inject2')
    frozen: bool = False               # True: Frozen Pinnacle block (FROZEN_CODES) instead of the original GB blocks
    adapter_distances: tuple = pf.ADAPTER_DISTANCES
    supply_filter: tuple = ()          # if given, keep only supplies whose name contains one of these strings


BASELINE = Choices()


def runs_per_copy(M):
    """Executions per copy worth trying: ceil(M/K) for K = 1..M copies."""
    return np.unique(np.ceil(M/np.arange(1, M + 1)))


def synthesis_loss(R, ch):
    """Loss charged to the R synthesized rotations: eps^2 (infidelity) or eps per rotation, times the shared scale."""
    return ch.synthesis_scale*R*(EPS**2 if ch.synthesis_charge == 'squared' else EPS)


def vary(**changes):
    return replace(BASELINE, **changes)


def support_reason(hw, pc=None):
    """None if modeled; otherwise the missing mapping."""
    if hw.connectivity != 'all_to_all':
        return '2D nearest-neighbor placement, routing, and error costs are not modeled'
    if hw.p != pf.P:
        return 'the published block fits and engine parameters are for p = 1e-3 only'
    return None


def p_block(code, hw):
    """Failure of one whole block per native logical period. Both noise classes use the p = 1e-3 depolarizing fit."""
    return code.error if code.error is not None else pf.p_block_d(code.d, hw.p)


_CACHE = {}


def supply_arrays(pc):
    key = (pc.families, pc.adapter_distances, pc.supply_filter)
    if key not in _CACHE:
        S = pf.all_supplies(pc.families, pc.adapter_distances)
        if pc.supply_filter:
            S = [x for x in S if any(f in x.name for f in pc.supply_filter)]
        a = lambda f, dt=float: np.array([f(s) for s in S], dtype=dt)
        _CACHE[key] = dict(S=S, qubits=a(lambda s: s.qubits), cycles=a(lambda s: s.cycles),
                           port=a(lambda s: s.port_qubits), slot=a(lambda s: s.min_slot),
                           loss=a(lambda s: -log1p(-s.error) - log1p(-s.port_error)))
    return _CACHE[key]


def min_qubits(w, ensemble, hw, ch, pc=BASELINE):
    """Fewest physical qubits for one QUOPS width and ensemble; None if no design passes or hw is unsupported."""
    if support_reason(hw):
        return None
    n_rot, _ = layer_counts(w, ensemble)
    R = w*n_rot
    synthesis = synthesis_loss(R, ch)
    if ch.synthesis_budget is None:
        stochastic = ch.threshold - synthesis
    else:
        stochastic = np.where(synthesis <= ch.synthesis_budget, ch.threshold, -np.inf)
    week = WEEK_SECONDS/ch.cycle_seconds
    runs = runs_per_copy(ch.executions)
    sa = supply_arrays(pc)
    NTR = (R*NT).astype(float)                                   # T states consumed, per T count
    magic = NTR[None, :]*sa['loss'][:, None]                     # (supplies, T counts)
    best = None
    for code in (FROZEN_CODES if pc.frozen else CODES):
        if not pc.frozen and code.d not in pc.codes:
            continue
        beta = ceil(w/code.k)
        per_cycle = pc.eta*beta*(-log1p(-p_block(code, hw)))/code.dt          # loss per code cycle of one block unit
        D = np.maximum(max(code.dt, pc.reaction), sa['slot'])[:, None]         # slot per supply
        tail = (w + 1)*D                                                       # startup plus final readouts
        room = (stochastic - magic)/per_cycle                                  # cycles the loss budget allows
        for r in runs:
            limit = np.minimum(room, week/r)
            avail = limit - tail
            k_f = np.floor(avail/sa['cycles'][:, None])                        # states one supply delivers in time
            ok = (avail >= NTR[None, :]*D) & (k_f >= 1)
            if not ok.any():
                continue
            F = np.ceil(NTR[None, :]/np.where(ok, k_f, 1))
            cycles = np.maximum(NTR[None, :]*D, np.ceil(NTR[None, :]/F)*sa['cycles'][:, None]) + tail
            per_copy = np.minimum(np.floor(week/cycles), ch.executions)
            copies = np.ceil(ch.executions/np.maximum(per_copy, 1))
            qubits = np.where(ok, copies*(beta*code.qubits + F*sa['qubits'][:, None] + sa['port'][:, None]), np.inf)
            i, j = np.unravel_index(np.argmin(qubits), qubits.shape)
            if best is None or qubits[i, j] < best['qubits']:
                s = sa['S'][i]
                c = cycles[i, j]
                best = dict(qubits=float(qubits[i, j]), distance=code.d, code_n=code.n, code_k=code.k, blocks=beta,
                            factory=s.name, family=s.family, factories=int(F[i, j]), eps=float(EPS[j]),
                            T_per_rotation=int(NT[j]), T_total=float(NTR[j]), slot=float(D[i, 0]), cycles=float(c),
                            run_seconds=float(c*ch.cycle_seconds), copies=int(copies[i, j]),
                            runs_per_copy=int(per_copy[i, j]),
                            loss_blocks=float(per_cycle*c), loss_magic=float(magic[i, j]),
                            loss_synthesis=float(synthesis[j]), supply_qubits=s.qubits, supply_cycles=s.cycles,
                            supply_error=s.error)
    return best


def widths(w_max=200000):
    small = list(range(2, 101, 2))
    large = sorted({int(2*round(v/2)) for v in np.geomspace(102, w_max, 120)})
    return small + [v for v in large if v > 100]


def sweep(hw, ch, pc=BASELINE):
    """Designs at the sampled widths up to N_MAX qubits; stops each ensemble at the first width that fails."""
    if support_reason(hw):
        return []
    rows = []
    for ensemble in ENSEMBLES:
        for w in widths():
            design = min_qubits(w, ensemble, hw, ch, pc)
            if design is None or design['qubits'] > N_MAX:
                break
            rows.append(dict(ensemble=ensemble, w=w, s=size(w, ensemble), passes=True, **design))
    return rows


def best_at(N, hw, ch, rows, pc=BASELINE):
    """Largest-s design with at most N qubits, exact over every even width (bisection, as in surface.best_at).

    For a fixed code, supply, and T count, the blocks, T demand, loss, run time, copies, and area all grow with w, so a
    design that passes at w + 2 also passes at w with no more qubits. Minima over a fixed finite menu keep this."""
    if support_reason(hw):
        return None
    top = None
    for ensemble in ENSEMBLES:
        fit = sorted((r['w'], r) for r in rows if r['ensemble'] == ensemble and r['qubits'] <= N)
        if not fit:
            continue
        lo, row = fit[-1]
        hi = next((w for w in widths() if w > lo), lo + 2)
        while hi - lo > 2:
            mid = lo + 2*((hi - lo)//4)
            x = min_qubits(mid, ensemble, hw, ch, pc)
            if x is not None and x['qubits'] <= N:
                lo, row = mid, dict(ensemble=ensemble, w=mid, s=size(mid, ensemble), passes=True, **x)
            else:
                hi = mid
        if top is None or row['s'] > top['s']:
            top = row
    return top


def capability(N, hw, ch, rows, pc=BASELINE):
    row = best_at(N, hw, ch, rows, pc)
    return row['s'] if row else 0


# Used by capability/frozen_pinnacle.py (Frozen Pinnacle, from the earlier serial estimate; not part of the main model).
def loss_for_retries(tau, fixed, h, rejection):
    """Minus log mean data survival for tau geometric waits and fixed slots.

    With per-slot survival s=exp(-h), E[s**G]=(1-r)*s/(1-r*s)
    for G=1,2,... . This expression remains stable when h is tiny.
    """
    if h < 0 or not 0 <= rejection < 1:
        raise ValueError('invalid slot hazard or engine rejection probability')
    return fixed*h + tau*(h + log1p(rejection*(-expm1(-h))/(1 - rejection)))


# Frozen Pinnacle with adapter-delivered T states: its fixed degree-eight graph cannot host the reconfigurable GB engines.
FROZEN_T = Choices(frozen=True, families=('surface', 'surface_jws'))
