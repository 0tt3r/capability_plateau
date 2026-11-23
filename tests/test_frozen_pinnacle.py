"""Frozen Pinnacle checks of published resources and the declared QUOPS adapter."""
from dataclasses import FrozenInstanceError, replace
from math import ceil, exp, expm1, floor, log1p

import numpy as np
import pytest

from capability import curves, frozen_pinnacle as frozen, quops
from capability.config import BASELINE, HARDWARE, WEEK_SECONDS, vary


def test_published_code_and_fixed_connectivity_resources():
    code = frozen.CODE
    assert (code.n, code.k, code.d, code.processing_k, code.code_cycles) == (510, 16, 24, 8, 30)
    # Section V: the complete block includes its optional memory-sector gadget.
    assert frozen.BASELINE.processing_qubits == 1020 + 9*99 + 101 + 10*51 == 2522
    assert frozen.ENGINE.name == '8T-to-CCZ'
    assert frozen.ENGINE.qubits == 1020 + 4*99 + 3*101 + 12*454 + 8*29 + 3*(101 + 51) == 7855
    assert frozen.ENGINE.error == pytest.approx(2e-10)
    assert frozen.ENGINE.rejection == pytest.approx(0.01)
    assert frozen.ENGINE.toffoli_slots == 5


def test_frozen_logical_error_fit_is_per_full_block_per_complete_cycle():
    # Eq. 1 and Table 2: shift plus bridged measurement; the rescaling to the
    # complete Frozen cycle is already in these fitted component rates.
    expected = 4.01*(0.001/0.0340)**12 + 5.49*(0.001/0.0284)**12
    assert frozen.p_block() == pytest.approx(expected, rel=1e-13, abs=0)
    assert expected == pytest.approx(2.1621585533884903e-17, rel=1e-13, abs=0)
    # Dividing by all 16 logical observables recovers the rounded paper value.
    assert expected/16 == pytest.approx(1.4e-18, rel=0.04, abs=0)


def independent_resources_and_losses(w, ensemble, nT, ch=BASELINE, fc=None):
    """Count resources without invoking the optimizer or its scheduling helpers."""
    if fc is None:
        fc = frozen.BASELINE
    rotations = w * quops.layer_counts(w, ensemble)[0]
    tau = rotations*nT
    assert tau % 2 == 0
    batches = tau//2
    blocks = max(3, ceil((w + 5)/8))
    period = max(30, fc.reaction_cycles)
    fixed = w + fc.seed_transfer_slots + fc.engine_startup_slots + fc.extra_batch_slots*batches
    expected_cycles = period*(fixed + 5*batches/0.99) + fc.seed_attempt_cycles/fc.seed_success
    h = blocks*period/30*fc.logical_error_scale*-log1p(-frozen.p_block())
    # A geometric trial has E[exp(-hG)]=(1-r)exp(-h)/(1-r exp(-h)).
    # This algebraic form is independent of the module's implementation.
    data = (fixed*h + batches*log1p(expm1(5*h)/0.99)
            + log1p(expm1(fc.seed_attempt_cycles/period*h)/fc.seed_success))
    magic = batches*-log1p(-2e-10)
    seed = (-log1p(-fc.seed_error)
            + fc.seed_transfer_slots*period*-log1p(-fc.seed_patch_error_per_cycle))
    engine_idle = (fc.extra_batch_slots*batches*period/30*fc.logical_error_scale
                   * -log1p(-frozen.p_block()))
    return dict(rotations=rotations, tau=tau, batches=batches, blocks=blocks,
                period=period, cycles=expected_cycles, data=data, magic=magic, seed=seed,
                engine_idle=engine_idle)


@pytest.mark.parametrize('w,ensemble', [(2, 'shallow'), (18, 'deep'), (100, 'shallow')])
@pytest.mark.parametrize('executions', [1, 10_000])
def test_design_has_all_ancillas_resources_retry_time_and_error_exposure(w, ensemble, executions):
    ch, hw = vary(executions=executions), HARDWARE['H2']
    row = frozen.min_qubits(w, ensemble, hw, ch)
    assert row is not None
    counted = independent_resources_and_losses(w, ensemble, row['T_per_rotation'], ch)
    blocks = counted['blocks']
    assert row['blocks'] == blocks
    assert row['ancilla_logicals'] == 5
    assert row['logical_capacity'] == 8*blocks
    assert row['workload_capacity'] == 8*blocks - 5 >= w
    assert row['T_total'] == counted['tau']
    assert row['CCZ_total'] == counted['batches']
    assert row['logical_cycle_cycles'] == counted['period']
    assert row['cycles'] == pytest.approx(counted['cycles'], rel=1e-12)
    assert row['run_seconds'] == pytest.approx(counted['cycles']*ch.cycle_seconds)
    runs = min(floor(WEEK_SECONDS/row['run_seconds']), executions)
    copies = ceil(executions/runs)
    assert row['copies'] == copies
    assert row['runs_per_copy'] == runs
    assert row['qubits'] == copies*(blocks*2522 + 7855 + 454 + 29)
    assert row['loss_data'] == pytest.approx(counted['data'], rel=1e-10, abs=1e-22)
    assert row['loss_magic'] == pytest.approx(counted['magic'], rel=1e-12)
    assert row['loss_seed'] == pytest.approx(counted['seed'], rel=1e-12)
    assert row['loss_engine_idle'] == pytest.approx(counted['engine_idle'], rel=1e-12, abs=0)
    synthesis = counted['rotations']*row['eps']**2
    stochastic = counted['data'] + counted['magic'] + counted['seed'] + counted['engine_idle']
    assert row['loss_synthesis'] == pytest.approx(synthesis)
    assert row['loss_stochastic'] == pytest.approx(stochastic, rel=1e-12)
    assert row['loss_total'] == pytest.approx(stochastic + synthesis, rel=1e-12)
    assert row['loss_total'] <= ch.threshold
    assert row['q_no_fault'] == pytest.approx(exp(-stochastic))
    assert row['quality'] == pytest.approx(exp(-stochastic - synthesis))


def test_initial_catalyst_error_is_charged_once_not_once_per_output():
    hw, ch = HARDWARE['H2'], vary(threshold=10)
    a = frozen.evaluate(18, 'deep', hw, ch, 24)
    b = frozen.evaluate(18, 'deep', hw, ch, 48)
    assert a is not None and b is not None
    assert b['CCZ_total'] == 2*a['CCZ_total']
    assert b['loss_magic'] == pytest.approx(2*a['loss_magic'])
    assert b['loss_seed'] == a['loss_seed']


def test_error_exposure_does_not_double_count_the_30_over_26_stretch():
    row = frozen.evaluate(12, 'deep', HARDWARE['H2'], BASELINE, 24)
    assert row is not None
    expected = independent_resources_and_losses(12, 'deep', 24)['data']
    assert row['loss_data'] == pytest.approx(expected, rel=1e-12)
    assert not np.isclose(row['loss_data'], expected*30/26, rtol=1e-3, atol=0)


def test_reaction_limited_slots_stretch_time_and_error_exposure():
    fc = frozen.vary(reaction_cycles=100)
    row = frozen.evaluate(18, 'deep', HARDWARE['H2'], BASELINE, 24, fc)
    assert row is not None
    counted = independent_resources_and_losses(18, 'deep', 24, fc=fc)
    assert row['logical_cycle_cycles'] == 100
    assert row['cycles'] == pytest.approx(counted['cycles'])
    assert row['loss_data'] == pytest.approx(counted['data'], rel=1e-12)


def test_three_distinct_magic_ports_set_the_minimum_footprint():
    row = frozen.min_qubits(2, 'shallow', HARDWARE['H2'])
    assert row is not None
    assert row['blocks'] == 3
    assert row['qubits'] == 3*2522 + 7855 + 454 + 29 == 15904
    rows = frozen.sweep(HARDWARE['H2'])
    assert frozen.best_at(10_000, HARDWARE['H2'], BASELINE, rows) is None
    assert frozen.capability(10_000, HARDWARE['H2'], BASELINE, rows) == 0


def test_integer_counts_exceeding_int32_remain_correct():
    ch = vary(threshold=5.0, cycle_seconds=1e-7)
    row = frozen.evaluate(600, 'deep', HARDWARE['H2'], ch, 100)
    assert row is not None
    tau = 600*(600**2//2)*100
    assert tau > np.iinfo(np.int32).max
    assert tau//2 > np.iinfo(np.int32).max
    assert row['T_total'] == tau
    assert row['CCZ_total'] == tau//2
    assert row['loss_magic'] == pytest.approx(tau//2*-log1p(-2e-10))
    assert row['loss_data'] > 0
    assert row['cycles'] > 0


def test_independent_copies_meet_the_expected_week_throughput():
    ch = vary(executions=10_000, cycle_seconds=1e-3)
    row = frozen.evaluate(18, 'deep', HARDWARE['H2'], ch, 24)
    assert row is not None and row['copies'] > 1
    assert row['copies']*row['runs_per_copy'] >= ch.executions
    assert (row['copies'] - 1)*row['runs_per_copy'] < ch.executions
    assert row['runs_per_copy']*row['run_seconds'] <= WEEK_SECONDS


@pytest.mark.parametrize('ensemble', ['shallow', 'deep'])
def test_optimizer_agrees_with_exhaustive_t_counts(ensemble):
    hw = HARDWARE['H2']
    candidates = [frozen.evaluate(18, ensemble, hw, BASELINE, int(nT)) for nT in frozen.NT]
    feasible = [row for row in candidates if row is not None]
    expected = min(feasible, key=lambda row: (row['qubits'], row['run_seconds'], row['loss_total']))
    row = frozen.min_qubits(18, ensemble, hw)
    assert row['qubits'] == expected['qubits']
    assert row['run_seconds'] == expected['run_seconds']


@pytest.mark.parametrize('ensemble', ['shallow', 'deep'])
def test_footprints_never_decrease_through_reserved_capacity_boundaries(ensemble):
    rows = [frozen.min_qubits(w, ensemble, HARDWARE['H2']) for w in range(2, 162, 2)]
    assert all(row is not None for row in rows)
    assert np.all(np.diff([row['qubits'] for row in rows]) >= 0)


def test_exact_probe_refines_the_sampled_curve_boundary():
    hw = HARDWARE['H2']
    rows = frozen.sweep(hw)
    row = frozen.best_at(50_000, hw, BASELINE, rows)
    assert row is not None and row['qubits'] <= 50_000
    next_row = frozen.min_qubits(row['w'] + 2, row['ensemble'], hw)
    assert next_row is None or next_row['qubits'] > 50_000
    assert frozen.capability(50_000, hw, BASELINE, rows) == row['s']
    assert row['s'] >= curves.capability(rows, 50_000)
    assert np.all(np.diff(curves.curve(rows)) >= 0)


def test_separate_synthesis_budget_preserves_the_stochastic_probability():
    hw = HARDWARE['H2']
    assert frozen.evaluate(4, 'shallow', hw, vary(threshold=1e-5), 10) is None
    ch = vary(threshold=1e-5, synthesis_budget=0.2)
    row = frozen.evaluate(4, 'shallow', hw, ch, 10)
    assert row is not None
    assert row['loss_stochastic'] <= ch.threshold
    assert row['loss_synthesis'] <= ch.synthesis_budget
    assert row['loss_total'] > ch.threshold
    assert row['quality'] < row['q_no_fault']


def test_bimodal_fit_is_the_same_declared_uniform_noise_proxy():
    a = frozen.min_qubits(18, 'deep', HARDWARE['H2'])
    b = frozen.min_qubits(18, 'deep', HARDWARE['H4'])
    assert a == b


@pytest.mark.parametrize('key', ['H1', 'H3'])
def test_degree_eight_does_not_grant_nearest_neighbor_grid_support(key):
    assert frozen.support_reason(HARDWARE[key]) is not None
    assert frozen.min_qubits(4, 'shallow', HARDWARE[key]) is None
    assert frozen.sweep(HARDWARE[key]) == []


def test_unpublished_physical_error_rate_is_not_interpolated():
    hw = replace(HARDWARE['H2'], p=1e-4, eps_1q=1e-4, eps_2q=1e-4)
    assert frozen.support_reason(hw) is not None
    assert frozen.min_qubits(4, 'shallow', hw) is None


@pytest.mark.parametrize('width', [-2, 0, 1, 3])
def test_invalid_width_is_rejected(width):
    with pytest.raises(ValueError):
        frozen.min_qubits(width, 'shallow', HARDWARE['H2'])


def test_architecture_choices_are_frozen_and_do_not_mutate_shared_defaults():
    changed = frozen.vary(extra_batch_slots=3)
    assert changed.extra_batch_slots == 3
    assert frozen.BASELINE.extra_batch_slots == 2
    assert BASELINE.executions == 1
    with pytest.raises(FrozenInstanceError):
        frozen.BASELINE.extra_batch_slots = 3
