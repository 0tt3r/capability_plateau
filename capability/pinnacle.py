"""Original Pinnacle implementation for QUOPS (arXiv:2602.11457v2).

Layout [paper, Table 1 and Eqs. 9--10, 18]: one processing unit of
ceil(w/k) bridged GB processing blocks and one 15-to-1 magic engine. All
allocated blocks remain live; there is no separate memory. Only all-to-all
hardware is supported here: the paper's nonlocal block connections have
not been compiled onto the repository's nearest-neighbour grid.

Schedule [paper, Sec. IV.2.1]: commute the Clifford gates into a classical
frame, then perform one Pauli measurement per T gate and w final readouts.
A slot lasts max(data d_t, engine d_t, engine protocol, reaction) code
cycles. An additional slot primes the engine and prepares the data
[choice]. Rejected engine outputs stall the data; no whole-circuit
restart is charged. Runtime uses the paper's *expected* retry count.

Error model [paper, Eq. 12 and Table 2]: use the logical-measurement fit
per GB block per native d_t cycles, including padded logical qubits.
Stretched slots incur proportionally more exposure [choice]. The default
factor of two accounts for X and Z errors by a symmetry approximation
[choice]; the source simulations used ideal most-likely decoding of X-type
detectors, not a demonstrated real-time decoder. Both hardware
noise classes use the same p=1e-3 fit; no bimodal threshold gain is assumed.
Engine output infidelity is a published black-box target, charged once per
consumed state. Mean data survival is averaged over geometric retry waits.
Synthesis uses the repository's fitted T counts and shared error-budget
choices; its fidelity charge stays separate from q_no_fault.

Search is over the five published codes, two published p=1e-3 engines,
and 5--160 T gates per rotation. Independent copies supply the requested
executions within the week at the *mean-time* level, not a probabilistic
deadline guarantee. See docs/pinnacle.md for formulas and limitations.
"""
import json
from dataclasses import dataclass, replace
from math import expm1, isclose, isfinite, log1p
from operator import index
from pathlib import Path

import numpy as np

from .config import BASELINE as SHARED_BASELINE, N_MAX, WEEK_SECONDS
from .quops import ENSEMBLES, layer_counts, size

PAPER = 'https://arxiv.org/html/2602.11457v2'
FIT = json.loads((Path(__file__).resolve().parents[1]/'results/synthesis_fit.json').read_text())
NT = np.arange(5, 161, dtype=np.int64)          # [choice] same finite synthesis menu as surface.py
EPS = 2.**(-(NT - FIT['b'])/FIT['a'])


@dataclass(frozen=True)
class Code:
    n: int
    k: int
    d: int
    code_cycles: int
    processing_qubits: int


# [paper] Table 1: processing footprints already include check qubits,
# four measurement gadgets, and four bridges. Do not add them again.
CODES = (
    Code(30, 8, 4, 6, 140),
    Code(62, 10, 6, 8, 244),
    Code(126, 12, 10, 12, 452),
    Code(254, 14, 16, 18, 860),
    Code(510, 16, 24, 26, 1620),
)


@dataclass(frozen=True)
class Engine:
    name: str
    error: float
    qubits: int
    rejection: float
    code_distance: int
    ancilla_distance: int
    postselection_rounds: int


# [paper] Eqs. 9--10 and Sec. V.2.2. These are p=1e-3 designs, not
# interpolated engines at arbitrary physical error rates. Eq. 10's label
# repeats 1e-9 in the paper; its surrounding text specifies 1e-11.
ENGINES = (
    Engine('15-to-1_1e-9', 1e-9, 4410, 0.10, 24, 7, 2),
    Engine('15-to-1_1e-11', 1e-11, 5430, 0.10, 24, 9, 2),
)


@dataclass(frozen=True)
class Choices:
    fit_a: float = 6.2                       # [paper] Table 2, logical measurements
    fit_b: float = 0.0158                    # [paper]
    fit_c: float = 0.47                      # [paper]
    logical_error_scale: float = 2.0         # [choice] extend X-observable fit to X and Z
    reaction_cycles: float = 10.0            # [paper] t_r/t_c; includes classical reaction
    startup_logical_cycles: int = 1          # [choice] data preparation and engine priming overlap
    code_distances: tuple = (4, 6, 10, 16, 24)
    engines: tuple = ('15-to-1_1e-9', '15-to-1_1e-11')


BASELINE = Choices()


def vary(**changes):
    return replace(BASELINE, **changes)


def support_reason(hw, pc=BASELINE):
    """None if modeled; otherwise explain the missing mapping or parameter fit."""
    if hw.connectivity != 'all_to_all':
        return '2D nearest-neighbor placement, routing, and error costs are not modeled'
    if not isclose(hw.p, 1e-3, rel_tol=1e-12, abs_tol=1e-15):
        return 'the implemented magic-engine menu is published only for p = 1e-3'
    if hw.noise not in ('uniform', 'bimodal'):
        return 'no noise mapping is declared for this hardware class'
    return None


def _validate(w, ensemble, ch, pc):
    if index(w) < 2 or w % 2:
        raise ValueError('QUOPS width must be an even integer of at least two')
    if ensemble not in ENSEMBLES:
        raise ValueError(f'unknown QUOPS ensemble: {ensemble}')
    if index(ch.executions) < 1:
        raise ValueError('executions must be a positive integer')
    if not isfinite(ch.cycle_seconds) or ch.cycle_seconds <= 0:
        raise ValueError('cycle_seconds must be finite and positive')
    if not isfinite(ch.threshold) or ch.threshold < 0:
        raise ValueError('threshold must be finite and nonnegative')
    if ch.synthesis_charge not in ('squared', 'linear'):
        raise ValueError('synthesis_charge must be squared or linear')
    if not isfinite(ch.synthesis_scale) or ch.synthesis_scale < 0:
        raise ValueError('synthesis_scale must be finite and nonnegative')
    if ch.synthesis_budget is not None and (not isfinite(ch.synthesis_budget) or ch.synthesis_budget < 0):
        raise ValueError('synthesis_budget must be nonnegative or None')
    if not isfinite(pc.logical_error_scale) or pc.logical_error_scale <= 0:
        raise ValueError('logical_error_scale must be finite and positive')
    if not isfinite(pc.reaction_cycles) or pc.reaction_cycles < 0:
        raise ValueError('reaction_cycles must be finite and nonnegative')
    if index(pc.startup_logical_cycles) < 1:
        raise ValueError('at least one startup slot is needed for the empty engine')
    if not pc.code_distances or not set(pc.code_distances) <= {c.d for c in CODES}:
        raise ValueError('code_distances must select from the published code menu')
    if not pc.engines or not set(pc.engines) <= {e.name for e in ENGINES}:
        raise ValueError('engines must select from the published engine menu')
    if not (isfinite(pc.fit_a) and pc.fit_a > 0 and isfinite(pc.fit_b) and pc.fit_b > 0
            and isfinite(pc.fit_c)):
        raise ValueError('invalid logical-error fit parameters')


def p_block(code, hw, pc=BASELINE):
    """Eq. 12: X-observable failure per whole block per native logical cycle."""
    p = pc.fit_a*(hw.p/pc.fit_b)**(code.d/2 + pc.fit_c)
    if not 0 <= p < 1:
        raise ValueError('logical-error fit is outside the probability range')
    return p


def engine_period(engine, pc=BASELINE):
    """Code cycles per output attempt, including Eq. 11 and engine error correction."""
    da, r, reaction = engine.ancilla_distance, engine.postselection_rounds, pc.reaction_cycles
    protocol = max(2*da + 4*r, reaction + 4*r, da + reaction + 3*r)
    native = next(c.code_cycles for c in CODES if c.d == engine.code_distance)
    return max(native, protocol, reaction)


def loss_for_retries(tau, fixed, h, rejection):
    """Minus log mean data survival for tau geometric waits and fixed slots.

    With per-slot survival s=exp(-h), E[s**G]=(1-r)*s/(1-r*s)
    for G=1,2,... . This expression remains stable when h is tiny.
    """
    if h < 0 or not 0 <= rejection < 1:
        raise ValueError('invalid slot hazard or engine rejection probability')
    return fixed*h + tau*(h + log1p(rejection*(-expm1(-h))/(1 - rejection)))


def _designs(w, ensemble, hw, ch, code, engine, counts, pc):
    """Vectorized candidate arithmetic. Integer products use explicit 64-bit counts."""
    rotations = int(w)*int(layer_counts(w, ensemble)[0])
    counts = np.asarray(counts, dtype=np.int64)
    if rotations*int(counts.max()) > np.iinfo(np.int64).max:
        raise OverflowError('T count exceeds the supported int64 range')
    tau = counts*rotations
    eps = 2.**(-(counts - FIT['b'])/FIT['a'])
    blocks = (int(w) + code.k - 1)//code.k
    period = max(code.code_cycles, engine_period(engine, pc))
    fixed = int(w) + pc.startup_logical_cycles
    slots = tau.astype(float)/(1 - engine.rejection) + fixed
    cycles = slots*period
    seconds = cycles*ch.cycle_seconds
    h = blocks*(period/code.code_cycles)*pc.logical_error_scale*(-log1p(-p_block(code, hw, pc)))
    loss_data = loss_for_retries(tau, fixed, h, engine.rejection)
    loss_magic = tau*(-log1p(-engine.error))
    synthesis = ch.synthesis_scale*rotations*(eps**2 if ch.synthesis_charge == 'squared' else eps)
    stochastic = loss_data + loss_magic
    loss = stochastic + synthesis
    if ch.synthesis_budget is None:
        passes = loss <= ch.threshold
    else:
        passes = (stochastic <= ch.threshold) & (synthesis <= ch.synthesis_budget)
    per_copy = np.minimum(np.floor(WEEK_SECONDS/seconds), ch.executions)
    copies = np.ceil(ch.executions/np.maximum(per_copy, 1))
    per_machine = blocks*code.processing_qubits + engine.qubits
    qubits = np.where(passes & (per_copy >= 1), copies*per_machine, np.inf)
    return dict(qubits=qubits, blocks=blocks, period=period, slots=slots, cycles=cycles,
                seconds=seconds, per_copy=per_copy, copies=copies, counts=counts, tau=tau,
                eps=eps, data=loss_data, magic=loss_magic, synthesis=synthesis,
                stochastic=stochastic, loss=loss, per_machine=per_machine)


def _row(a, i, w, ensemble, code, engine):
    return dict(ensemble=ensemble, w=int(w), s=size(int(w), ensemble), passes=True,
                qubits=int(a['qubits'][i]), qubits_per_copy=a['per_machine'],
                distance=code.d, code_n=code.n, code_k=code.k, blocks=a['blocks'],
                logical_capacity=a['blocks']*code.k, engine=engine.name,
                engine_qubits=engine.qubits, engine_error=engine.error,
                engine_rejection=engine.rejection, eps=float(a['eps'][i]),
                T_per_rotation=int(a['counts'][i]), T_total=int(a['tau'][i]),
                logical_cycle_cycles=float(a['period']), logical_cycles=float(a['slots'][i]),
                cycles=float(a['cycles'][i]), run_seconds=float(a['seconds'][i]),
                copies=int(a['copies'][i]), runs_per_copy=int(a['per_copy'][i]),
                loss_data=float(a['data'][i]), loss_magic=float(a['magic'][i]),
                loss_synthesis=float(a['synthesis'][i]), loss_stochastic=float(a['stochastic'][i]),
                loss_total=float(a['loss'][i]), q_no_fault=float(np.exp(-a['stochastic'][i])),
                quality=float(np.exp(-a['loss'][i])))


def evaluate(w, ensemble, hw, ch, code, engine, nT, pc=BASELINE):
    """One finite-menu design; None if unsupported, over budget, or too slow."""
    _validate(w, ensemble, ch, pc)
    if support_reason(hw, pc) is not None:
        return None
    if index(nT) < int(NT[0]) or nT > int(NT[-1]):
        raise ValueError('T gates per rotation must be between 5 and 160')
    if code not in CODES or engine not in ENGINES:
        raise ValueError('evaluate requires a published code and engine')
    a = _designs(w, ensemble, hw, ch, code, engine, [nT], pc)
    return _row(a, 0, w, ensemble, code, engine) if np.isfinite(a['qubits'][0]) else None


def min_qubits(w, ensemble, hw, ch=SHARED_BASELINE, pc=BASELINE):
    """Fewest qubits in the finite serial model; break area ties by mean runtime."""
    _validate(w, ensemble, ch, pc)
    if support_reason(hw, pc) is not None:
        return None
    best = None
    for code in CODES:
        if code.d not in pc.code_distances:
            continue
        for engine in ENGINES:
            if engine.name not in pc.engines:
                continue
            a = _designs(w, ensemble, hw, ch, code, engine, NT, pc)
            i = np.lexsort((a['loss'], a['seconds'], a['qubits']))[0]
            if not np.isfinite(a['qubits'][i]):
                continue
            row = _row(a, i, w, ensemble, code, engine)
            if best is None or (row['qubits'], row['run_seconds'], row['loss_total']) < (
                    best['qubits'], best['run_seconds'], best['loss_total']):
                best = row
    return best


def widths(w_max=200000):
    """Same plotting resolution as surface.py: dense to 100, then about 6.5% apart."""
    small = list(range(2, min(100, w_max) + 1, 2))
    large = sorted({int(2*round(v/2)) for v in np.geomspace(102, max(102, w_max), 120)})
    return small + [w for w in large if 100 < w <= w_max]


def sweep(hw, ch=SHARED_BASELINE, pc=BASELINE):
    """Passing sampled widths, stopping at the first infeasible or over-N_MAX width."""
    _validate(2, 'shallow', ch, pc)
    if support_reason(hw, pc) is not None:
        return []
    rows = []
    for ensemble in ENSEMBLES:
        for w in widths():
            row = min_qubits(w, ensemble, hw, ch, pc)
            if row is None or row['qubits'] > N_MAX:
                break
            rows.append(row)
    return rows


def best_at(N, hw, ch, rows, pc=BASELINE):
    """Largest size at N, refining sampled gaps over every even width.

    For each fixed code/engine/T count, blocks, T demand, slot hazard,
    synthesis charge, mean runtime, copies, and area are nondecreasing
    in w. Feasible designs therefore form a prefix; minima over the same
    finite candidate menu preserve monotonicity. This justifies bisection
    and sweep's early termination even at block-size jumps.
    """
    _validate(2, 'shallow', ch, pc)
    if not isfinite(N) or N < 0 or N > N_MAX:
        raise ValueError(f'N must lie between zero and {N_MAX:g}')
    if support_reason(hw, pc) is not None:
        return None
    top = None
    for ensemble in ENSEMBLES:
        fit = sorted((r['w'], r) for r in rows if r['ensemble'] == ensemble
                     and r['passes'] and r['qubits'] <= N)
        if not fit:
            continue
        lo, row = fit[-1]
        hi = next((w for w in widths() if w > lo), lo + 2)
        while hi - lo > 2:
            mid = lo + 2*((hi - lo)//4)
            candidate = min_qubits(mid, ensemble, hw, ch, pc)
            if candidate is not None and candidate['qubits'] <= N:
                lo, row = mid, candidate
            else:
                hi = mid
        if top is None or row['s'] > top['s']:
            top = row
    return top


def capability(N, hw, ch, rows, pc=BASELINE):
    row = best_at(N, hw, ch, rows, pc)
    return row['s'] if row else 0
