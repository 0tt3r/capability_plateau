"""Sequential fast-block PBC with a finite-buffer factory schedule.

All Cliffords are absorbed into measurement axes. A run consumes R*n_T magic
states, then measures w transformed readouts sequentially. Factories start cold,
produce synchronized batches, and transfer them into F*B distance-d buffer tiles.
The next batch is launched just in time for the buffers to empty (or immediately
after the previous production if supply is slower). Thus no unbounded inventory
or uncharged storage at a factory output is needed. Delivery and classical
feed-forward are exposed assumptions; the default charges zero extra cycles.

The shared p_SC law is applied per live tile, including buffers, as an effective
volume model. It is not a calibration of multi-qubit fast-block patches/surgery.
"""
from math import ceil, log1p

import numpy as np

from . import surface
from .config import PBC, WEEK_SECONDS
from .factories import menu
from .quops import ENSEMBLES, layer_counts


def batch_schedule(states, factories, batch_size, batch_cycles, service_cycles, readout_cycles):
    """Completion time in cycles, vectorized over state/factory counts.

    Wave i becomes available at L + i*max(L, F*B*h). It contains F*B states
    consumed sequentially in h cycles each. A partial last wave is retained.
    This is an achievable throttled schedule with F*B output buffer tiles.
    """
    states, factories = np.asarray(states), np.asarray(factories)
    if (np.any(states < 0) or np.any(states != np.floor(states))
            or np.any(factories < 1) or np.any(factories != np.floor(factories))
            or batch_size < 1 or int(batch_size) != batch_size
            or batch_cycles <= 0 or service_cycles <= 0 or readout_cycles < 0):
        raise ValueError('Invalid batch schedule parameters')
    wave = factories*batch_size
    completed_waves = np.maximum(states - 1, 0)//wave
    last = np.maximum(states - 1, 0) % wave + 1
    cycles = (batch_cycles + completed_waves*np.maximum(batch_cycles, wave*service_cycles)
              + last*service_cycles + readout_cycles)
    return np.where(states > 0, cycles, readout_cycles)


def _validate(w, ensemble, ch, architecture):
    if w < 2 or int(w) != w or w % 2 or ensemble not in ENSEMBLES:
        raise ValueError('QUOPS requires even integer width >= 2 and a known ensemble')
    if architecture.schedule != 'pbc' or architecture.layout != 'fast_columns':
        raise ValueError('Sequential PBC requires the fast layout')
    if architecture.supply not in ('buffered', 'throughput_bound'):
        raise ValueError('Unknown PBC supply policy')
    if (not np.isfinite(architecture.delivery_cycles) or architecture.delivery_cycles < 0
            or not np.isfinite(architecture.feedforward_cycles) or architecture.feedforward_cycles < 0
            or ch.t_per_injection_steps < 1 or ch.executions < 1
            or int(ch.executions) != ch.executions or ch.cycle_seconds <= 0):
        raise ValueError('Invalid PBC timing/execution settings')


def _costs(w, ensemble, hw, ch, architecture, d, f, nt, eps, factories):
    rotations = w*layer_counts(w, ensemble)[0]
    states = rotations*nt
    area = surface.data_tiles(w, ch, architecture)
    buffers = factories*f.states_per_batch if architecture.supply == 'buffered' else 0*factories
    tiles = area + buffers
    service = ch.t_per_injection_steps*d + architecture.delivery_cycles + architecture.feedforward_cycles
    if architecture.supply == 'buffered':
        cycles = batch_schedule(states, factories, f.states_per_batch, f.batch_cycles, service, w*d)
    else:
        # The manuscript gives two lower bounds, not a cold-start queue schedule.
        cycles = np.maximum(states*service+w*d, np.ceil(states/factories)*f.cycles_per_output)
    patches = tiles*cycles*(-log1p(-surface.p_sc(d, hw, ch)))
    magic = states*(-log1p(-f.error))
    synthesis = ch.synthesis_scale*rotations*(eps**2 if ch.synthesis_charge == 'squared' else eps)
    loss = patches + magic
    if ch.synthesis_budget is None:
        ok = loss + synthesis <= ch.threshold
    else:
        ok = (loss <= ch.threshold) & (synthesis <= ch.synthesis_budget)
    runs = np.minimum(np.floor(WEEK_SECONDS/(cycles*ch.cycle_seconds)), ch.executions)
    copies = np.ceil(ch.executions/np.maximum(runs, 1))
    per_copy = 2*d*d*tiles + factories*f.qubits
    qubits = copies*per_copy
    return dict(cycles=cycles, loss_patches=patches, loss_magic=magic,
                loss_synthesis=synthesis, runs_per_copy=runs, copies=copies,
                qubits_per_copy=per_copy, qubits=qubits, passes=ok & (runs >= 1),
                data_tiles=area, buffer_tiles=buffers, tiles=tiles, service_cycles=service)


def evaluate_design(w, ensemble, hw, ch, d, f, nt, factories, architecture=PBC):
    """Inspect any design in the synthesis menu, including failing designs."""
    _validate(w, ensemble, ch, architecture)
    if d < 3 or d % 2 == 0 or int(d) != d or factories < 1 or int(factories) != factories:
        raise ValueError('Use odd distance >= 3 and a positive integer factory count')
    index = np.flatnonzero(surface.NT == nt)
    if not len(index):
        raise ValueError('T count is outside the synthesis menu')
    eps = float(surface.EPS[index[0]])
    x = _costs(w, ensemble, hw, ch, architecture, d, f, nt, eps, factories)
    states = w*layer_counts(w, ensemble)[0]*nt
    row = {k:float(v) for k,v in x.items() if k != 'passes'}
    for k in ('runs_per_copy', 'copies', 'data_tiles', 'buffer_tiles', 'tiles'):
        row[k] = int(row[k])
    row.update(passes=bool(x['passes']), distance=d, factory=f.name, factories=factories,
               eps=eps, T_per_rotation=nt, T_total=states,
               run_seconds=row['cycles']*ch.cycle_seconds,
               schedule=architecture.schedule, layout=architecture.layout,
               factory_batch_size=f.states_per_batch, factory_batch_cycles=f.batch_cycles,
               startup_cycles=f.batch_cycles if architecture.supply == 'buffered' else 0,
               final_measurements=w, final_measurement_cycles=w*d,
               cnot_cycles=0, consumption_cycles=states*row['service_cycles']+w*d)
    if architecture.supply == 'throughput_bound':
        row['schedule_status'] = 'optimistic throughput lower bound; not a delivery schedule'
    row['supply_stall_cycles'] = row['cycles'] - row['startup_cycles'] - row['consumption_cycles']
    return row


def min_qubits(w, ensemble, hw, ch, architecture=PBC):
    """Exhaustive distances/T counts/factory types/counts with exact copy sizing.

    F >= ceil(tau/h) eliminates supply stalls. More factories cannot improve
    runtime under this schedule and only increase area and buffer exposure, so
    enumerating through that saturation count suffices. Copy count then follows
    directly from K*floor(W/T) >= M; it is not a separately approximated search.
    """
    _validate(w, ensemble, ch, architecture)
    area = surface.data_tiles(w, ch, architecture)
    nt, eps = surface.NT[:, None], surface.EPS[:, None]
    best = None
    for f in menu(ch.factories):
        for d in range(3, 101, 2):
            if best is not None and area*2*d*d >= best['qubits']:
                break
            service = ch.t_per_injection_steps*d + architecture.delivery_cycles + architecture.feedforward_cycles
            saturation = max(1, ceil(f.cycles_per_output/service))
            counts = np.arange(1, saturation + 1)[None, :]
            x = _costs(w, ensemble, hw, ch, architecture, d, f, nt, eps, counts)
            qubits = np.where(x['passes'], x['qubits'], np.inf)
            index = np.unravel_index(np.argmin(qubits), qubits.shape)
            if np.isfinite(qubits[index]) and (best is None or qubits[index] < best['qubits']):
                i, j = index
                best = evaluate_design(w, ensemble, hw, ch, d, f, int(surface.NT[i]),
                                       int(counts[0, j]), architecture)
    return best
