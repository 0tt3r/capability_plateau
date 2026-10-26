"""Surface-code implementation (draft Sec. 3.5) for QUOPS.

Layout [choice, not compiled]: w data patches, each with one routing patch beside it (tiles = 2w by
default), all live for the whole run, plus F identical factories that deliver T states into the routing
region. A patch holds 2 d^2 physical qubits [draft].

Schedule, in surface-code cycles [choice, not compiled]:
  preparation and final measurement   2d
  CNOT layer                          two lattice-surgery steps (Z_c Z_a, then X_a X_t), d cycles each; all w/2
                                      CNOTs at once ('parallel', assumes no routing congestion) or one at a time
                                      ('serial', which any connected routing region can realize)
  rotation layer                      R_P(theta) is synthesized to n_T T gates. Each T consumes one factory
                                      state in one surgery step (d cycles), measured in the rotation's own Pauli
                                      basis, so no logical H is needed; Clifford corrections are tracked in
                                      software. Rotations on different qubits run in parallel, so the layer takes
                                      max(n_T d, ceil(w n_T / F) * factory period) cycles.

Loss (minus the log no-fault probability) [draft Sec. 3.5, plus the synthesis rule]:
  every live patch-cycle fails with P_SC(d) = A (p/p_th)^((d+1)/2); every consumed T state fails with the
  factory's output error. Each rotation's synthesis error eps counts as eps^2 ('squared') or eps ('linear'),
  either inside the same budget ('synthesis_budget' None) or against a separate budget of its own.

A design passes if the loss is within budget and the requested executions complete in the week: K copies,
each running floor(week / run time) executions back to back, must reach the requested number. For each width
the search is exhaustive over the discrete design variables (distance, factory type, T count per rotation,
factory count), so min_qubits is the exact minimum within this model.
"""
import json
from math import ceil, log1p
from pathlib import Path

import numpy as np

from .config import N_MAX, WEEK_SECONDS
from .factories import menu
from .quops import ENSEMBLES, layer_counts, size

FIT = json.loads((Path(__file__).resolve().parents[1]/'results/synthesis_fit.json').read_text())
NT = np.arange(5, 161)                          # T gates per rotation
EPS = 2.**(-(NT - FIT['b'])/FIT['a'])           # loosest tolerance the fit gives exactly NT T gates (0.42 to 3e-16)


def t_count(eps):
    """T gates per rotation from the gridsynth fit (results/synthesis_fit.json), rounded up."""
    return np.ceil(FIT['a']*np.log2(1/eps) + FIT['b'])


def p_sc(d, hw, ch):
    p_th = ch.sc_threshold_uniform if hw.noise == 'uniform' else ch.sc_threshold_bimodal
    return ch.sc_prefactor*(hw.p/p_th)**((d + 1)/2)


def synthesis_loss(R, ch):
    return ch.synthesis_scale*R*(EPS**2 if ch.synthesis_charge == 'squared' else EPS)


def runs_per_copy(M):
    """Executions per copy worth trying: ceil(M/K) for K = 1..M copies. Any other count needs as many copies
    as the next one below it while leaving each run less time."""
    return np.unique(np.ceil(M/np.arange(1, M + 1)))


def min_qubits(w, ensemble, hw, ch):
    """Fewest physical qubits for one QUOPS width and ensemble; None if no design passes."""
    n_rot, n_cnot = layer_counts(w, ensemble)
    R = w*n_rot
    tiles = ceil(ch.tiles_per_data_patch*w)
    synthesis = synthesis_loss(R, ch)
    if ch.synthesis_budget is None:
        stochastic = ch.threshold - synthesis
    else:
        stochastic = np.where(synthesis <= ch.synthesis_budget, ch.threshold, -np.inf)
    week = WEEK_SECONDS/ch.cycle_seconds                  # cycles in the week
    M = ch.executions
    runs = runs_per_copy(M)
    demand = (w*NT)[:, None]                              # T states per rotation layer
    best = None
    for f in menu(ch.factories):
        tau = f.cycles_per_output
        magic = R*NT*(-log1p(-f.error))
        for d in range(3, 101, 2):
            patches = tiles*2*d*d
            if best is not None and patches >= best['qubits']:
                break                                       # area only grows with d
            lam = -log1p(-p_sc(d, hw, ch))
            t_cnot = 2*d if ch.cnot_layer == 'parallel' else (w//2)*2*d
            inject = (ch.t_per_injection_steps*d*NT)[:, None]  # one qubit's T gates in sequence
            fixed = n_cnot*t_cnot + 2*d
            room = (stochastic - magic)/(tiles*lam)         # cycles the loss budget allows one run
            limit = np.minimum(room[:, None], week/runs[None, :])
            per_layer = (limit - fixed)/n_rot
            k = np.floor(per_layer/tau)                     # T states one factory can deliver per layer
            ok = (per_layer >= inject) & (k >= 1)
            if not ok.any():
                continue
            F = np.ceil(demand/np.where(ok, k, 1))          # fewest factories that keep the run within limit
            cycles = n_rot*np.maximum(inject, np.ceil(demand/F)*tau) + fixed
            per_copy = np.minimum(np.floor(week/cycles), M)
            copies = np.ceil(M/np.maximum(per_copy, 1))
            qubits = np.where(ok, copies*(patches + F*f.qubits), np.inf)
            i, j = np.unravel_index(np.argmin(qubits), qubits.shape)
            if best is None or qubits[i, j] < best['qubits']:
                c = cycles[i, j]
                best = dict(qubits=float(qubits[i, j]), distance=d, factory=f.name, factories=int(F[i, j]),
                            eps=float(EPS[i]), T_per_rotation=int(NT[i]), T_total=float(R*NT[i]),
                            cycles=float(c), run_seconds=float(c*ch.cycle_seconds), copies=int(copies[i, j]),
                            runs_per_copy=int(per_copy[i, j]), tiles=tiles,
                            loss_patches=float(tiles*c*lam), loss_magic=float(magic[i]),
                            loss_synthesis=float(synthesis[i]))
    return best


def widths(w_max=200000):
    """Sampled widths: every even width to 100, then about 6.5% apart."""
    small = list(range(2, 101, 2))
    large = sorted({int(2*round(v/2)) for v in np.geomspace(102, w_max, 120)})
    return small + [v for v in large if v > 100]


def sweep(hw, ch):
    """Designs at the sampled widths up to N_MAX qubits. The fewest qubits never decrease with w (see best_at),
    so the first width that fails or needs more than N_MAX ends each ensemble."""
    rows = []
    for ensemble in ENSEMBLES:
        for w in widths():
            design = min_qubits(w, ensemble, hw, ch)
            if design is None or design['qubits'] > N_MAX:
                break
            rows.append(dict(ensemble=ensemble, w=w, s=size(w, ensemble), passes=True, **design))
    return rows


def best_at(N, hw, ch, rows):
    """The largest-s design with at most N qubits, exact over every even width.

    The fewest qubits needed never decreases with w: a design (distance, factory, T count, factory count) that
    passes at w + 2 also passes at w, with no more qubits. So it suffices to bisect the gap between the largest
    sampled width that fits and the next sampled width."""
    top = None
    for ensemble in ENSEMBLES:
        fit = sorted((r['w'], r) for r in rows if r['ensemble'] == ensemble and r['qubits'] <= N)
        if not fit:
            continue
        lo, row = fit[-1]
        hi = next((w for w in widths() if w > lo), lo + 2)     # first sampled width that does not fit
        while hi - lo > 2:
            mid = lo + 2*((hi - lo)//4)
            x = min_qubits(mid, ensemble, hw, ch)
            if x is not None and x['qubits'] <= N:
                lo, row = mid, dict(ensemble=ensemble, w=mid, s=size(mid, ensemble), passes=True, **x)
            else:
                hi = mid
        if top is None or row['s'] > top['s']:
            top = row
    return top


def capability(N, hw, ch, rows):
    row = best_at(N, hw, ch, rows)
    return row['s'] if row else 0
