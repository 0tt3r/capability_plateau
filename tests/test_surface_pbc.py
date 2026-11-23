"""Independent queue simulations, exhaustive searches and supported PBC defaults."""
from dataclasses import replace
from math import ceil, log1p

import numpy as np
import pytest

from capability import quops, surface, surface_pbc
from capability.config import BASELINE, HARDWARE, PBC, PBC_UNITS, SurfaceArchitecture, SURFACE_ARCHITECTURES, WEEK_SECONDS
from capability.factories import FACTORIES, Factory


def simulate_batches(states, factories, batch_size, latency, service, readout):
    """Explicit production/buffer/consumption events; no closed-form recurrence."""
    now, consumed, inventory = 0., 0, 0
    arrival = latency
    while consumed < states:
        if inventory == 0:
            now = max(now, arrival)
            inventory = factories*batch_size
            # Delay production to avoid holding completed outputs in factories.
            arrival = now + max(latency, inventory*service)
        now += service
        consumed += 1
        inventory -= 1
    return now + readout


@pytest.mark.parametrize('batch_size,latency', [(1, 100.), (4, 130.), (4, 157.)])
def test_batch_formula_matches_explicit_finite_buffer_simulation(batch_size, latency):
    for n in (0, 1, 3, 4, 7, 12, 37, 101):
        for f in (1, 2, 5, 12):
            for service in (3., 11., 19.5, 160.):
                got = surface_pbc.batch_schedule(n, f, batch_size, latency, service, 33.)
                assert got == pytest.approx(simulate_batches(n, f, batch_size, latency, service, 33.))


def test_startup_last_partial_batch_and_final_measurements_are_charged():
    # Two factories each deliver four states at cycle 130. Nine states need a
    # second wave at 260, one consumption to 270, and final readout to 290.
    assert surface_pbc.batch_schedule(9, 2, 4, 130, 10, 20) == 290
    assert surface_pbc.batch_schedule(1, 2, 4, 130, 10, 20) == 160
    assert surface_pbc.batch_schedule(9, 4, 4, 130, 10, 20) == 240


def test_factory_saturation_is_a_valid_search_limit():
    for batch, latency in ((1, 100.), (4, 130.), (4, 157.)):
        for service in (3, 11, 37):
            f = ceil(latency/(batch*service))
            for n in (1, 12, 1001):
                expected = latency + n*service + 77
                for extra in (0, 1, 20):
                    assert surface_pbc.batch_schedule(n, f+extra, batch, latency, service, 77) == expected


def test_counted_area_is_exact_on_squares_and_counts_partial_columns():
    for k in range(1, 20):
        assert surface.data_tiles(2*k*k, BASELINE) == (2*k+1)**2
    for w in range(2, 100, 2):
        height = ceil(np.sqrt(w//2))
        columns = ceil((w//2)/height)
        assert surface.data_tiles(w, BASELINE) == 2*w + 2*(height+columns) + 1
    assert surface.data_tiles(6, BASELINE) == 21


def test_pbc_charges_buffers_and_readout_without_cnot_runtime():
    a = surface.min_qubits(6, 'deep', HARDWARE['H1'], BASELINE)
    assert a['cnot_cycles'] == 0
    assert a['final_measurement_cycles'] == 6*a['distance']
    assert a['tiles'] == a['data_tiles'] + a['buffer_tiles']
    assert a['buffer_tiles'] == a['factories']*a['factory_batch_size']
    assert a['loss_patches'] == pytest.approx(a['tiles']*a['cycles']*
                                              -log1p(-surface.p_sc(a['distance'], HARDWARE['H1'], BASELINE)))
    assert a['qubits'] == a['copies']*a['qubits_per_copy']


@pytest.mark.parametrize('hardware', ['H1', 'H3'])
@pytest.mark.parametrize('multiple', [False, True])
def test_optimizer_against_unpruned_independent_small_search(monkeypatch, hardware, multiple):
    nt = np.array([12, 20, 28])
    eps = 2.**(-(nt-surface.FIT['b'])/surface.FIT['a'])
    monkeypatch.setattr(surface, 'NT', nt)
    monkeypatch.setattr(surface, 'EPS', eps)
    factories = [Factory('test', 'single', 2e-9, 1024, 100),
                 Factory('test', 'batch', 1e-10, 2000, 130/4, 4)]
    monkeypatch.setattr(surface_pbc, 'menu', lambda families: factories)
    hw, w, ensemble = HARDWARE[hardware], 6, 'shallow'
    ch = replace(BASELINE, executions=10000, cycle_seconds=1.) if multiple else BASELINE
    arch = replace(PBC, delivery_cycles=2., feedforward_cycles=3.)
    best = float('inf')
    # Includes counts beyond saturation; evaluate no optimizer pruning or helper costs.
    for f in factories:
        for d in range(3, 32, 2):
            for n, epsilon in zip(nt, eps):
                for count in range(1, 40):
                    states = w*w//2*int(n)
                    cycles = simulate_batches(states, count, f.states_per_batch, f.batch_cycles, d+5, w*d)
                    height = ceil(np.sqrt(w//2))
                    tiles = 2*w+2*(height+ceil((w//2)/height))+1+count*f.states_per_batch
                    loss = (tiles*cycles*-log1p(-surface.p_sc(d, hw, ch))
                            + states*-log1p(-f.error) + w*w//2*epsilon**2)
                    runs = min(ch.executions, int(WEEK_SECONDS/(cycles*ch.cycle_seconds)))
                    if loss <= ch.threshold and runs > 0:
                        copies = ceil(ch.executions/runs)
                        best = min(best, copies*(tiles*2*d*d+count*f.qubits))
    got = surface.min_qubits(w, ensemble, hw, ch, arch)
    assert got['qubits'] == best
    assert got['copies']*got['runs_per_copy'] >= ch.executions
    if multiple:
        assert got['copies'] > 1


def test_separate_synthesis_budget_and_failed_deadline():
    hw, f = HARDWARE['H1'], FACTORIES[0]
    ch = replace(BASELINE, synthesis_budget=.01)
    row = surface_pbc.evaluate_design(6, 'shallow', hw, ch, 11, f, 28, 1)
    assert row['passes']
    assert row['loss_synthesis'] <= ch.synthesis_budget
    assert row['loss_magic'] + row['loss_patches'] <= ch.threshold
    failed = surface_pbc.evaluate_design(6, 'shallow', hw, replace(ch, cycle_seconds=1e9), 11, f, 28, 1)
    assert not failed['passes'] and failed['runs_per_copy'] == 0


def test_architecture_is_retained_during_capability_refinement():
    hw = HARDWARE['H3']
    rows = surface.sweep(hw, BASELINE, PBC)
    row = surface.best_at(1e6, hw, BASELINE, rows)
    assert row['surface_architecture']['schedule'] == 'pbc'
    nxt = surface.min_qubits(row['w']+2, row['ensemble'], hw, BASELINE, PBC)
    assert nxt is None or nxt['qubits'] > 1e6
    with pytest.raises(ValueError, match='different surface architecture'):
        surface.best_at(1e6, hw, BASELINE, rows, PBC_UNITS)


def test_width_monotonicity_for_bisection():
    for hw in (HARDWARE['H1'], HARDWARE['H3']):
        for ensemble in quops.ENSEMBLES:
            values = [surface.min_qubits(w, ensemble, hw, BASELINE, PBC)['qubits']
                      for w in range(2, 42, 2)]
            assert np.all(np.diff(values) >= 0)


@pytest.mark.parametrize('change', [dict(layout='paired'), dict(layout='fast'),
                                  dict(schedule='layered'), dict(schedule='typo')])
def test_removed_or_unknown_architectures_fail_at_configuration(change):
    with pytest.raises(ValueError):
        replace(PBC, **change)


def test_invalid_pbc_timing_fails():
    with pytest.raises(ValueError):
        surface.min_qubits(6, 'shallow', HARDWARE['H1'], BASELINE, replace(PBC, feedforward_cycles=-1))


def test_default_is_counted_sequential_pbc():
    hw = HARDWARE['H1']
    default = surface.min_qubits(6, 'deep', hw, BASELINE)
    assert default == surface.min_qubits(6, 'deep', hw, BASELINE, PBC)
    assert SurfaceArchitecture() == PBC
    assert default['data_tiles'] == 21
    assert default['surface_architecture']['schedule'] == 'pbc'
    assert all(a.schedule in ('pbc', 'pbc_units') for a in SURFACE_ARCHITECTURES.values())


def test_retired_rows_cannot_be_silently_refined_as_pbc():
    hw = HARDWARE['H1']
    row = dict(surface.min_qubits(6, 'deep', hw, BASELINE), w=6, s=216, ensemble='deep')
    del row['surface_architecture']
    for architecture in (None, PBC):
        with pytest.raises(ValueError, match='lack architecture metadata'):
            surface.best_at(1e6, hw, BASELINE, [row], architecture)
    row['surface_architecture'] = dict(layout='paired', schedule='layered')
    with pytest.raises(ValueError, match='layered model was removed'):
        surface.best_at(1e6, hw, BASELINE, [row])
