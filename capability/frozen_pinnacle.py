"""Serial Frozen Pinnacle QUOPS estimate with a declared phase-catalysis adapter.

Published hardware [arXiv:2609.21249v1, Secs. III--V]: distance-24 GB blocks,
eight processing logicals per 2,522-qubit block, 30 code cycles per logical
measurement, and a 7,855-qubit CCZ engine (including its magic port). One
Toffoli uses five logical cycles; output error is 2e-10, rejection is 1%.
Its three dedicated injection ancillas occupy three distinct data blocks.
The native fixed-degree-eight graph is not a nearest-neighbour square grid.

QUOPS adaptation [choice, not a published Frozen benchmark]: use Fig. 14 of
Gidney--Fowler, arXiv:1812.01238v3. One Toffoli, a reusable T catalyst, an AND
scratch, and Clifford frame updates implement two T gates. Pair gates on
distinct qubits within each rotation synthesis round. Reserve five logical
ancillas: catalyst, scratch, and three magic-port ancillas. Charge scratch
preparation and readout as two additional serial slots per pair.

Bootstrap [choice]: an extra 454-qubit cultivated-state patch and 29-qubit
bridge feed an existing magic-port gadget. Six slots transfer the seed to a
dedicated catalyst logical and reset the port; three prime the CCZ pipeline.
The seed is freshly prepared per execution, with error 2.7e-6 and mean wait
18/0.35 code cycles (Frozen Appendix A). These source input-state parameters
do not establish the proposed transfer circuit; its fixed layout is not
compiled. A surface-patch storage-error allowance is charged during transfer.

Errors [paper, Eq. 1/Table 2]: add the Shift and Bridged Measurement fits per
whole block per full Frozen logical cycle. The 30-cycle normalization already
includes automorphisms. A factor of two for full-Pauli error is a choice. Both
noise classes use the same fixed-p=1e-3 fit, with no claimed bimodal gain.
Mean data survival averages geometric CCZ and seed retries. CCZ faults are
charged once per catalytic pair and seed faults once per execution: poisoning
correlates subsequent T errors, so they are not independent T-state faults.

Only the published distance-24 configuration is searched, with 5--160 fitted
T gates per rotation. Copies and the one-week budget use mean runtime, not a
deadline guarantee. See docs/frozen_pinnacle.md for assumptions and formulas.
"""
from dataclasses import dataclass, replace
from math import isclose, isfinite, log1p
from operator import index

import numpy as np

from .config import BASELINE as SHARED_BASELINE, N_MAX, WEEK_SECONDS
from .pinnacle import FIT, NT, loss_for_retries, widths
from .quops import ENSEMBLES, layer_counts, size

PAPER = 'https://arxiv.org/html/2609.21249v1'
CATALYSIS_PAPER = 'https://arxiv.org/html/1812.01238v3'
ANCILLA_LOGICALS = 5


@dataclass(frozen=True)
class Code:
    n: int = 510
    k: int = 16
    d: int = 24
    processing_k: int = 8
    code_cycles: int = 30


@dataclass(frozen=True)
class Engine:
    name: str = '8T-to-CCZ'
    qubits: int = 7855
    error: float = 2e-10
    rejection: float = 0.01
    toffoli_slots: int = 5


CODE, ENGINE = Code(), Engine()


@dataclass(frozen=True)
class Choices:
    processing_qubits: int = 2522             # [paper] includes optional-memory gadget, retained here
    logical_error_scale: float = 2.0          # [choice] full-Pauli multiplier on published observable fit
    reaction_cycles: float = 10.0             # [paper] reaction time / code-cycle time
    extra_batch_slots: int = 2                # [choice] scratch preparation and X readout
    seed_patch_qubits: int = 454              # [paper input patch] proposed dedicated bootstrap adapter
    seed_bridge_qubits: int = 29              # [paper input bridge] 2*15 - 1
    seed_error: float = 2.7e-6                 # [paper] cultivated input-state error, charged once/run
    seed_success: float = 0.35                # [paper] approximate success per escaped-state attempt
    seed_attempt_cycles: int = 18             # [paper] cultivation rounds, equated to code cycles
    seed_transfer_slots: int = 6              # [choice] seed transfer, internal teleport, port reset
    engine_startup_slots: int = 3             # [choice] first Toffoli completes at slot eight, then every five
    seed_patch_error_per_cycle: float = 3e-10  # [choice] uniform d=15 surface-law storage estimate
    final_measurements: bool = True           # [choice] w readout slots; False gives N_M = 0, as for QUOPS in the paper


BASELINE = Choices()


def vary(**changes):
    return replace(BASELINE, **changes)


def support_reason(hw, fc=BASELINE):
    if hw.connectivity != 'all_to_all':
        return '2D nearest-neighbor placement, routing, and error costs are not modeled'
    if not isclose(hw.p, 1e-3, rel_tol=1e-12, abs_tol=1e-15):
        return 'the Frozen error fits and CCZ engine are specified only at p = 1e-3'
    if hw.noise not in ('uniform', 'bimodal'):
        return 'no noise mapping is declared for this hardware class'
    return None


def p_block():
    """Eq. 1/Table 2: per-block Shift plus Bridged Measurement, already at 30 cycles."""
    return 4.01*(1e-3/0.0340)**12 + 5.49*(1e-3/0.0284)**12


def _validate(w, ensemble, ch, fc):
    if index(w) < 2 or w % 2:
        raise ValueError('QUOPS width must be an even integer of at least two')
    if ensemble not in ENSEMBLES:
        raise ValueError(f'unknown QUOPS ensemble: {ensemble}')
    if index(ch.executions) < 1 or not isfinite(ch.cycle_seconds) or ch.cycle_seconds <= 0:
        raise ValueError('executions and cycle_seconds must be positive')
    if not isfinite(ch.threshold) or ch.threshold < 0:
        raise ValueError('threshold must be finite and nonnegative')
    if ch.synthesis_charge not in ('squared', 'linear'):
        raise ValueError('synthesis_charge must be squared or linear')
    if not isfinite(ch.synthesis_scale) or ch.synthesis_scale < 0:
        raise ValueError('synthesis_scale must be finite and nonnegative')
    if ch.synthesis_budget is not None and (not isfinite(ch.synthesis_budget) or ch.synthesis_budget < 0):
        raise ValueError('synthesis_budget must be nonnegative or None')
    if fc.processing_qubits not in (2522, 2372):
        raise ValueError('choose the published 2522-qubit block or derived 2372-qubit no-memory block')
    if not isfinite(fc.logical_error_scale) or fc.logical_error_scale <= 0:
        raise ValueError('logical_error_scale must be finite and positive')
    if not isfinite(fc.reaction_cycles) or fc.reaction_cycles < 0:
        raise ValueError('reaction_cycles must be finite and nonnegative')
    for name, minimum in (('extra_batch_slots', 2), ('seed_patch_qubits', 454),
                          ('seed_bridge_qubits', 29), ('seed_attempt_cycles', 1),
                          ('seed_transfer_slots', 6), ('engine_startup_slots', 3)):
        if index(getattr(fc, name)) < minimum:
            raise ValueError(f'{name} must be at least {minimum}')
    if not 0 < fc.seed_success <= 1:
        raise ValueError('seed_success must lie in (0, 1]')
    if not 0 <= fc.seed_error < 1 or not 0 <= fc.seed_patch_error_per_cycle < 1:
        raise ValueError('seed error probabilities must lie in [0, 1)')


def _designs(w, ensemble, hw, ch, counts, fc):
    rotations = int(w)*int(layer_counts(w, ensemble)[0])
    counts = np.asarray(counts, dtype=np.int64)
    if rotations*int(counts.max()) > np.iinfo(np.int64).max:
        raise OverflowError('T count exceeds the supported int64 range')
    tau = rotations*counts
    batches = tau//2                         # even width permits pairs in each synthesis round
    eps = 2.**(-(counts - FIT['b'])/FIT['a'])
    blocks = max(3, (int(w) + ANCILLA_LOGICALS + CODE.processing_k - 1)//CODE.processing_k)
    period = max(CODE.code_cycles, fc.reaction_cycles)
    fixed = (int(w) if fc.final_measurements else 0) + fc.seed_transfer_slots + fc.engine_startup_slots + fc.extra_batch_slots*batches
    slots = fixed + ENGINE.toffoli_slots*batches.astype(float)/(1 - ENGINE.rejection)
    seed_wait = fc.seed_attempt_cycles/fc.seed_success
    cycles = period*slots + seed_wait
    seconds = cycles*ch.cycle_seconds
    lam = -log1p(-p_block())
    h = blocks*(period/CODE.code_cycles)*fc.logical_error_scale*lam
    seed_wait_loss = loss_for_retries(1, 0, h*fc.seed_attempt_cycles/period, 1 - fc.seed_success)
    data = (fixed*h + loss_for_retries(batches, 0, ENGINE.toffoli_slots*h, ENGINE.rejection)
            + seed_wait_loss)
    magic = batches*(-log1p(-ENGINE.error))
    # Pausing the engine during catalytic scratch operations adds storage
    # beyond its published five-slot output estimate. Use the block fit as
    # an explicit storage proxy for that one engine block, without recharging production.
    engine_idle = fc.extra_batch_slots*batches*(period/CODE.code_cycles)*fc.logical_error_scale*lam
    seed = (-log1p(-fc.seed_error)
            + fc.seed_transfer_slots*period*(-log1p(-fc.seed_patch_error_per_cycle)))
    stochastic = data + magic + seed + engine_idle
    synthesis = ch.synthesis_scale*rotations*(eps**2 if ch.synthesis_charge == 'squared' else eps)
    loss = stochastic + synthesis
    if ch.synthesis_budget is None:
        passes = loss <= ch.threshold
    else:
        passes = (stochastic <= ch.threshold) & (synthesis <= ch.synthesis_budget)
    per_copy = np.minimum(np.floor(WEEK_SECONDS/seconds), ch.executions)
    copies = np.ceil(ch.executions/np.maximum(per_copy, 1))
    per_machine = (blocks*fc.processing_qubits + ENGINE.qubits
                   + fc.seed_patch_qubits + fc.seed_bridge_qubits)
    qubits = np.where(passes & (per_copy >= 1), copies*per_machine, np.inf)
    return dict(qubits=qubits, per_machine=per_machine, blocks=blocks, period=period,
                slots=slots + seed_wait/period, cycles=cycles, seconds=seconds,
                seed_wait=seed_wait, per_copy=per_copy, copies=copies, counts=counts,
                tau=tau, batches=batches, eps=eps, data=data, magic=magic, seed=seed,
                engine_idle=engine_idle, stochastic=stochastic, synthesis=synthesis, loss=loss)


def _row(a, i, w, ensemble, fc):
    return dict(model='frozen_phase_catalysis', ensemble=ensemble, w=int(w),
                s=size(int(w), ensemble), passes=True, qubits=int(a['qubits'][i]),
                qubits_per_copy=a['per_machine'], distance=CODE.d, code_n=CODE.n,
                code_k=CODE.k, processing_k=CODE.processing_k, blocks=a['blocks'],
                logical_capacity=a['blocks']*CODE.processing_k,
                workload_capacity=a['blocks']*CODE.processing_k - ANCILLA_LOGICALS,
                ancilla_logicals=ANCILLA_LOGICALS, processing_qubits=fc.processing_qubits,
                engine=ENGINE.name, engine_qubits=ENGINE.qubits, engine_error=ENGINE.error,
                engine_rejection=ENGINE.rejection, seed_adapter_qubits=fc.seed_patch_qubits+fc.seed_bridge_qubits,
                seed_wait_cycles=float(a['seed_wait']), eps=float(a['eps'][i]),
                T_per_rotation=int(a['counts'][i]), T_total=int(a['tau'][i]),
                CCZ_total=int(a['batches'][i]), logical_cycle_cycles=float(a['period']),
                logical_cycles=float(a['slots'][i]), cycles=float(a['cycles'][i]),
                run_seconds=float(a['seconds'][i]), copies=int(a['copies'][i]),
                runs_per_copy=int(a['per_copy'][i]), loss_data=float(a['data'][i]),
                loss_magic=float(a['magic'][i]), loss_seed=float(a['seed']),
                loss_engine_idle=float(a['engine_idle'][i]), loss_synthesis=float(a['synthesis'][i]),
                loss_stochastic=float(a['stochastic'][i]), loss_total=float(a['loss'][i]),
                q_no_fault=float(np.exp(-a['stochastic'][i])), quality=float(np.exp(-a['loss'][i])))


def evaluate(w, ensemble, hw, ch, nT, fc=BASELINE):
    """One synthesis candidate in the proposed Frozen/phase-catalysis model."""
    _validate(w, ensemble, ch, fc)
    if support_reason(hw, fc) is not None:
        return None
    if index(nT) < int(NT[0]) or nT > int(NT[-1]):
        raise ValueError('T gates per rotation must be between 5 and 160')
    a = _designs(w, ensemble, hw, ch, [nT], fc)
    return _row(a, 0, w, ensemble, fc) if np.isfinite(a['qubits'][0]) else None


def min_qubits(w, ensemble, hw, ch=SHARED_BASELINE, fc=BASELINE):
    """Fewest qubits over the finite T-count menu; area ties prefer mean runtime."""
    _validate(w, ensemble, ch, fc)
    if support_reason(hw, fc) is not None:
        return None
    a = _designs(w, ensemble, hw, ch, NT, fc)
    i = np.lexsort((a['loss'], a['seconds'], a['qubits']))[0]
    return _row(a, i, w, ensemble, fc) if np.isfinite(a['qubits'][i]) else None


def sweep(hw, ch=SHARED_BASELINE, fc=BASELINE):
    """Passing sampled widths up to N_MAX; no lower-distance Frozen menu is assumed."""
    _validate(2, 'shallow', ch, fc)
    if support_reason(hw, fc) is not None:
        return []
    rows = []
    for ensemble in ENSEMBLES:
        for w in widths():
            row = min_qubits(w, ensemble, hw, ch, fc)
            if row is None or row['qubits'] > N_MAX:
                break
            rows.append(row)
    return rows


def best_at(N, hw, ch, rows, fc=BASELINE):
    """Exact fixed-budget probe over even widths in this finite serial model.

    For each T count, blocks, paired T demand, fixed-slot hazard, retry
    exposure, synthesis loss, runtime, and copies never decrease with w.
    Feasibility is therefore a prefix, including the minimum-three-block
    plateau. This justifies early termination and refinement by bisection.
    """
    _validate(2, 'shallow', ch, fc)
    if not isfinite(N) or N < 0 or N > N_MAX:
        raise ValueError(f'N must lie between zero and {N_MAX:g}')
    if support_reason(hw, fc) is not None:
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
            candidate = min_qubits(mid, ensemble, hw, ch, fc)
            if candidate is not None and candidate['qubits'] <= N:
                lo, row = mid, candidate
            else:
                hi = mid
        if top is None or row['s'] > top['s']:
            top = row
    return top


def capability(N, hw, ch, rows, fc=BASELINE):
    row = best_at(N, hw, ch, rows, fc)
    return row['s'] if row else 0
