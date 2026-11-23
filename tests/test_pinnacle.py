"""Pinnacle checks against the paper and independently counted resource constraints."""
from dataclasses import FrozenInstanceError, replace
from math import ceil, exp, expm1, floor, log1p

import numpy as np
import pytest

from capability import curves, pinnacle, quops
from capability.config import BASELINE, HARDWARE, WEEK_SECONDS, vary


def code_at(distance):
    return next(code for code in pinnacle.CODES if code.d == distance)


def engine_at(error):
    return next(engine for engine in pinnacle.ENGINES if engine.error == error)


def test_published_processing_blocks():
    # Table 1 of arXiv:2602.11457v2; includes data, checks, gadgets, and bridges.
    expected = [
        (30, 8, 4, 6, 140),
        (62, 10, 6, 8, 244),
        (126, 12, 10, 12, 452),
        (254, 14, 16, 18, 860),
        (510, 16, 24, 26, 1620),
    ]
    assert [(c.n, c.k, c.d, c.code_cycles, c.processing_qubits)
            for c in pinnacle.CODES] == expected


@pytest.mark.parametrize('error,qubits,ancilla', [(1e-9, 4410, 7), (1e-11, 5430, 9)])
def test_published_magic_engines(error, qubits, ancilla):
    engine = engine_at(error)
    assert engine.qubits == qubits
    assert engine.code_distance == 24
    assert engine.ancilla_distance == ancilla
    assert engine.postselection_rounds == 2
    assert engine.rejection == pytest.approx(0.1)
    assert pinnacle.engine_period(engine) == 26


@pytest.mark.parametrize('distance', [4, 6, 10, 16, 24])
def test_processing_block_error_is_the_published_surgery_fit(distance):
    code, hw = code_at(distance), HARDWARE['H2']
    # Eq. 12 and Table 2: this is per block per native logical cycle, not per qubit.
    expected = 6.2 * (hw.p / 0.0158)**(distance / 2 + 0.47)
    assert pinnacle.p_block(code, hw) == pytest.approx(expected, rel=1e-13)


def test_reaction_time_stretches_the_engine_period():
    engine = engine_at(1e-9)
    choices = pinnacle.vary(reaction_cycles=100)
    # Eq. 11 with d_a=7, r=2: max(22, 108, 113), above the native 26-cycle engine.
    assert pinnacle.engine_period(engine, choices) == 113


def test_geometric_retries_have_the_correct_mean_survival():
    code = code_at(4)
    engine = engine_at(1e-9)
    ch = vary(threshold=100.0)
    pc = pinnacle.vary(logical_error_scale=1.0)
    row = pinnacle.evaluate(2, 'shallow', HARDWARE['H2'], ch, code, engine, 5, pc)
    assert row is not None
    period = pinnacle.engine_period(engine, pc)
    hazard = (period / code.code_cycles) * -log1p(-pinnacle.p_block(code, HARDWARE['H2']))
    # Enumerate one successful T consumption after 1, 2, ... attempts, including idle retries.
    survive_one = sum((1 - engine.rejection) * engine.rejection**(attempts - 1)
                      * exp(-attempts * hazard) for attempts in range(1, 200))
    tau = 2 * 5
    fixed = 2 + pc.startup_logical_cycles
    expected = exp(-fixed * hazard) * survive_one**tau
    assert exp(-row['loss_data']) == pytest.approx(expected, rel=1e-12)
    assert row['logical_cycles'] == pytest.approx(fixed + tau / (1 - engine.rejection))
    # Substituting the mean duration into exp(-hazard*time) is not the same expectation.
    assert expected > exp(-hazard * row['logical_cycles'])


@pytest.mark.parametrize('hazard', [0.0, 1e-15, 0.015, 10.0])
@pytest.mark.parametrize('rejection', [0.0, 0.1, 0.25, 0.9])
def test_retry_survival_matches_finite_enumeration(hazard, rejection):
    tau, fixed = 10, 3
    survive_one = sum((1 - rejection) * rejection**(attempts - 1) * exp(-attempts * hazard)
                      for attempts in range(1, 500))
    expected = exp(-fixed * hazard) * survive_one**tau
    actual = exp(-pinnacle.loss_for_retries(tau, fixed, hazard, rejection))
    assert actual == pytest.approx(expected, rel=1e-12, abs=0)


@pytest.mark.parametrize('w,ensemble', [(4, 'shallow'), (12, 'deep'), (100, 'shallow')])
@pytest.mark.parametrize('executions', [1, 10_000])
def test_design_satisfies_independently_counted_resources_and_losses(w, ensemble, executions):
    hw, ch = HARDWARE['H2'], vary(executions=executions)
    row = pinnacle.min_qubits(w, ensemble, hw, ch)
    assert row is not None
    code, engine = code_at(row['distance']), next(e for e in pinnacle.ENGINES if e.name == row['engine'])
    blocks = ceil(w / code.k)
    tau = w * quops.layer_counts(w, ensemble)[0] * row['T_per_rotation']
    fixed = w + pinnacle.BASELINE.startup_logical_cycles
    period = max(code.code_cycles, pinnacle.engine_period(engine))
    expected_cycles = (fixed + tau / (1 - engine.rejection)) * period
    assert row['blocks'] == blocks
    assert row['logical_capacity'] == blocks * code.k
    assert row['T_total'] == tau
    assert row['cycles'] == pytest.approx(expected_cycles)
    assert row['run_seconds'] == pytest.approx(expected_cycles * ch.cycle_seconds)
    runs = min(floor(WEEK_SECONDS / row['run_seconds']), executions)
    copies = ceil(executions / runs)
    assert row['runs_per_copy'] == runs
    assert row['copies'] == copies
    assert row['qubits'] == copies * (blocks * code.processing_qubits + engine.qubits)
    hazard = (blocks * period / code.code_cycles * pinnacle.BASELINE.logical_error_scale
              * -log1p(-pinnacle.p_block(code, hw)))
    expected_data = fixed * hazard + tau * log1p(expm1(hazard) / (1 - engine.rejection))
    expected_magic = tau * -log1p(-engine.error)
    expected_synthesis = w * quops.layer_counts(w, ensemble)[0] * row['eps']**2
    assert row['loss_data'] == pytest.approx(expected_data, rel=1e-6, abs=1e-12)
    assert row['loss_magic'] == pytest.approx(expected_magic)
    assert row['loss_synthesis'] == pytest.approx(expected_synthesis)
    assert row['loss_total'] == pytest.approx(expected_data + expected_magic + expected_synthesis,
                                              rel=1e-6, abs=1e-12)
    assert row['loss_total'] <= ch.threshold + 1e-12
    assert row['q_no_fault'] == pytest.approx(exp(-row['loss_stochastic']))
    assert row['quality'] == pytest.approx(exp(-row['loss_total']))


def test_multiple_copies_obey_the_common_week_budget():
    ch = vary(executions=10_000, cycle_seconds=1e-3)
    row = pinnacle.evaluate(12, 'deep', HARDWARE['H2'], ch, code_at(24), engine_at(1e-11), 40)
    assert row is not None and row['copies'] > 1
    assert row['copies'] * row['runs_per_copy'] >= ch.executions
    assert (row['copies'] - 1) * row['runs_per_copy'] < ch.executions
    assert row['runs_per_copy'] * row['run_seconds'] <= WEEK_SECONDS


def test_large_t_count_does_not_overflow_a_windows_int32():
    row = pinnacle.evaluate(600, 'deep', HARDWARE['H2'], BASELINE, code_at(24), engine_at(1e-11), 100)
    assert row is not None
    expected = 600 * (600**2 // 2) * 100
    assert expected > np.iinfo(np.int32).max
    assert row['T_total'] == expected
    assert row['loss_magic'] == pytest.approx(expected * -log1p(-1e-11))
    assert row['loss_data'] >= 0
    assert row['cycles'] > 0


def test_linear_and_squared_synthesis_charges_have_distinct_pass_criteria():
    args = (4, 'shallow', HARDWARE['H2'])
    code, engine = code_at(16), engine_at(1e-9)
    squared = pinnacle.evaluate(*args, BASELINE, code, engine, 10)
    assert squared is not None
    assert pinnacle.evaluate(*args, vary(synthesis_charge='linear'), code, engine, 10) is None
    linear = pinnacle.evaluate(*args, vary(synthesis_charge='linear', threshold=10), code, engine, 10)
    rotations = 4 * quops.layer_counts(4, 'shallow')[0]
    assert linear['loss_synthesis'] == pytest.approx(rotations * linear['eps'])
    assert squared['loss_synthesis'] == pytest.approx(rotations * squared['eps']**2)


def test_separate_synthesis_budget_preserves_stochastic_no_fault_probability():
    args = (4, 'shallow', HARDWARE['H2'])
    code, engine = code_at(16), engine_at(1e-9)
    shared = vary(threshold=0.001)
    assert pinnacle.evaluate(*args, shared, code, engine, 10) is None
    separate = vary(threshold=0.001, synthesis_budget=0.2)
    row = pinnacle.evaluate(*args, separate, code, engine, 10)
    assert row is not None
    assert row['loss_stochastic'] <= separate.threshold
    assert row['loss_synthesis'] <= separate.synthesis_budget
    assert row['loss_total'] > separate.threshold
    assert row['q_no_fault'] >= exp(-separate.threshold)
    assert row['quality'] < row['q_no_fault']
    assert pinnacle.evaluate(*args, vary(threshold=0.001, synthesis_budget=0.01),
                              code, engine, 10) is None


@pytest.mark.parametrize('key', ['H1', 'H3'])
def test_nearest_neighbor_grid_has_no_free_nonlocal_pinnacle_operations(key):
    hw = HARDWARE[key]
    assert pinnacle.min_qubits(4, 'shallow', hw, BASELINE) is None
    assert pinnacle.sweep(hw, BASELINE) == []


def test_bimodal_hardware_uses_the_same_declared_conservative_fit():
    a = pinnacle.min_qubits(12, 'deep', HARDWARE['H2'], BASELINE)
    b = pinnacle.min_qubits(12, 'deep', HARDWARE['H4'], BASELINE)
    for field in ('qubits', 'distance', 'engine', 'T_per_rotation', 'cycles', 'loss_total'):
        assert a[field] == b[field]


def test_unpublished_physical_error_rate_does_not_reuse_the_engine_menu():
    hw = replace(HARDWARE['H2'], p=1e-4, eps_1q=1e-4, eps_2q=1e-4)
    assert pinnacle.min_qubits(4, 'shallow', hw) is None
    assert pinnacle.sweep(hw) == []


@pytest.mark.parametrize('width', [-2, 0, 1, 3])
def test_invalid_quops_width_is_rejected(width):
    with pytest.raises(ValueError):
        pinnacle.min_qubits(width, 'shallow', HARDWARE['H2'])


@pytest.mark.parametrize('ensemble', ['shallow', 'deep'])
def test_optimizer_matches_exhaustive_discrete_designs(ensemble):
    hw, w = HARDWARE['H2'], 6
    candidates = [pinnacle.evaluate(w, ensemble, hw, BASELINE, code, engine, int(nT))
                  for code in pinnacle.CODES for engine in pinnacle.ENGINES for nT in pinnacle.NT]
    feasible = [row for row in candidates if row is not None]
    expected = min(row['qubits'] for row in feasible)
    actual = pinnacle.min_qubits(w, ensemble, hw, BASELINE)
    assert actual['qubits'] == expected


@pytest.mark.parametrize('ensemble', ['shallow', 'deep'])
def test_minimum_qubits_never_decrease_across_block_boundaries(ensemble):
    rows = [pinnacle.min_qubits(w, ensemble, HARDWARE['H2']) for w in range(2, 82, 2)]
    assert all(row is not None for row in rows)
    assert np.all(np.diff([row['qubits'] for row in rows]) >= 0)


def test_capability_probe_refines_the_budget_boundary_and_works_with_curves():
    hw = HARDWARE['H2']
    rows = pinnacle.sweep(hw, BASELINE)
    row = pinnacle.best_at(10_000, hw, BASELINE, rows)
    assert row is not None and row['qubits'] <= 10_000
    next_design = pinnacle.min_qubits(row['w'] + 2, row['ensemble'], hw, BASELINE)
    assert next_design is None or next_design['qubits'] > 10_000
    assert pinnacle.capability(10_000, hw, BASELINE, rows) == row['s']
    assert row['s'] >= curves.capability(rows, 10_000)
    values = curves.curve(rows)
    assert np.all(np.diff(values) >= 0)


def test_architecture_choices_do_not_mutate_shared_configuration():
    changed = pinnacle.vary(reaction_cycles=20)
    assert changed.reaction_cycles == 20
    assert pinnacle.BASELINE.reaction_cycles == 10
    assert BASELINE.executions == 1
    with pytest.raises(FrozenInstanceError):
        pinnacle.BASELINE.reaction_cycles = 20
