"""Integration audit of the seven STAR assumption decisions.

Compare production against the independently saved analysis and exercise
deadline/copy constraints, calibrated domains and angle-dependent supply.
"""
from dataclasses import replace
import json
from math import exp
from pathlib import Path

import numpy as np
import pytest

from capability import curves, star
from capability.config import BASELINE, HARDWARE, STAR_BASELINE, STAR_LEGACY, WEEK_SECONDS
from capability.star_data import tables, transversal_data
from capability.star_states import (counted_quantities, state_profile, uniform_angles,
                                    transversal_correction_resources)

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = json.loads((Path(__file__).parent/'fixtures/star_reference.json').read_text())


def test_strict_lead_intervals_exclude_ties_and_include_last_budget():
    def passing(qubits, s):
        return dict(qubits=qubits, s=s, passes=True)

    rows = [passing(2, 4), passing(6, 8)]
    references = [passing(4, 4), passing(8, 10)]
    assert curves.lead_intervals(rows, references, max_qubits=6) == [[2, 4], [6, 7]]
    assert curves.lead_intervals(rows, references, max_qubits=9) == [[2, 4], [6, 8]]
    assert curves.lead_intervals(references, rows, max_qubits=6) == []


@pytest.mark.parametrize('label', ['H1', 'H3'])
@pytest.mark.parametrize('stress', [False, True])
def test_main_estimator_reproduces_all_saved_optimistic_designs(label, stress):
    name = 'omitted_toshio_sector_stress' if stress else 'adopted'
    reference = REFERENCE['hardware'][label][name]
    settings = replace(STAR_BASELINE, omitted_sector_stress=stress)
    actual = star.sweep(HARDWARE[label], BASELINE, settings)
    expected = reference['rows']
    assert len(actual) == len(expected)
    for row, saved in zip(actual, expected):
        for key in ('ensemble', 'w', 's', 'qubits', 'distance', 'groups', 'stations', 'passes'):
            assert row[key] == saved[key]
        for key in ('cycles', 'run_seconds', 'loss_patches', 'loss_injection', 'loss_total',
                    'preparation_attempts', 'resource_states', 'correction_resources',
                    'correction_preparation_attempts', 'original_state_error', 'original_acceptance'):
            assert row[key] == pytest.approx(saved[key], rel=2e-13, abs=1e-14)
        assert row['loss_preparation_memory'] == row['loss_deadline'] == 0
        assert row['deadline_seconds'] == row['run_seconds']
        assert row['data_routing_tiles'] == 2*row['w']


def test_packaged_input_provenance_matches_counted_sources():
    data = tables()
    folder = ROOT/'analysis/separate_branch_survival_20261113'
    full = json.loads((folder/'results.json').read_text())['rows']
    single = json.loads((folder/'two_fault_terms.json').read_text())['rows']
    original = json.loads((ROOT/'analysis/original_star_channel_20261115/acceptance.json').read_text())['rows']
    for row in full:
        for label in ('H1', 'H3'):
            assert data['hardware'][label]['transversal'][str(row['groups'])]['H_n'] == row['hardware'][label]['H_n']
    for row in single:
        assert data['hardware'][row['hardware']]['transversal'][str(row['distance']//2)]['T1'] == row['exactly_one_active']
    for row in original:
        for label in ('H1', 'H3'):
            for key, value in row['hardware'][label].items():
                assert data['hardware'][label]['original'][str(row['distance'])][key] == value


@pytest.mark.parametrize('label', ['H1', 'H3'])
@pytest.mark.parametrize('k', [3, 5, 7, 9])
def test_retained_single_fault_numerator_keeps_single_qubit_contribution(label, k):
    hw = HARDWARE[label]
    row = transversal_data(hw, k)
    p, q = hw.eps_2q, hw.eps_1q
    # Independently derived zero-signature background times primary-fault odds.
    expected = row['zero_signature_probability']*k*(2*p/(15*(1-p))+q/(3*(1-q)))
    assert row['T1'] == pytest.approx(expected, rel=1e-14)
    assert expected > row['zero_signature_probability']*k*2*p/(15*(1-p))


@pytest.mark.parametrize('mode', ['mean', 'markov'])
def test_deadline_and_multiple_copy_search_match_brute_force(mode):
    hw = HARDWARE['H3']
    ch = replace(BASELINE, executions=100_000, cycle_seconds=1.0)
    sc = replace(STAR_BASELINE, groups=(3,), max_distance=9, angle_points=32,
                 deadline_model=mode, deadline_failure=.01)
    rows = [star.evaluate_design(4, 'shallow', hw, ch, sc, d, 3, f)
            for d in (3, 5, 7, 9) for f in range(1, 5)]
    valid = [r for r in rows if r['passes']]
    if valid:
        expected = min(r['qubits'] for r in valid)
        row = star.min_qubits(4, 'shallow', hw, ch, sc)
        assert row['qubits'] == expected
        assert row['copies']*int(WEEK_SECONDS/row['deadline_seconds']) >= ch.executions
        assert row['success_probability_lower_bound'] == pytest.approx(exp(-row['loss_total']))
        assert row['copies'] > 1
    else:
        assert star.min_qubits(4, 'shallow', hw, ch, sc) is None


def test_mean_runtime_does_not_depend_on_unused_markov_parameter():
    hw = HARDWARE['H1']
    a = star.evaluate_design(8, 'shallow', hw, BASELINE, STAR_BASELINE, 9, 3, 2)
    b = star.evaluate_design(8, 'shallow', hw, BASELINE,
                             replace(STAR_BASELINE, deadline_failure=0.), 9, 3, 2)
    assert a == b


def test_correction_supply_uses_actual_angles_and_clifford_termination():
    sc, hw, k = STAR_BASELINE, HARDWARE['H3'], 3
    phi = np.r_[0., uniform_angles(sc.angle_points)]
    attempts, states, tail = transversal_correction_resources(phi, k, hw, sc)
    assert attempts[0] == states[0] == 0
    assert np.all(attempts >= states)
    assert np.all(states <= 2.)
    assert np.ptp(attempts[1:]) > .1  # not a constant original-STAR cost
    short, _, short_tail = transversal_correction_resources(phi, k, hw, replace(sc, correction_rounds=30))
    assert np.all(np.abs(short-attempts) <= short_tail+tail)
    profile = state_profile(phi[1:], k, 11, hw, BASELINE, sc)
    row = star.evaluate_design(8, 'shallow', hw, BASELINE, sc, 11, k, 2)
    expected = 64*np.mean(profile.correction_probability*attempts[1:])
    assert row['correction_preparation_attempts'] == pytest.approx(expected)
    assert not np.isclose(expected, 64*profile.correction_probability.mean()*attempts[1:].mean(), rtol=1e-4)


@pytest.mark.parametrize('change', [dict(eps_1q=0.), dict(eps_2q=2e-3)])
def test_counted_rates_are_never_silently_extrapolated(change):
    with pytest.raises(ValueError, match='noise rates'):
        star.min_qubits(4, 'shallow', replace(HARDWARE['H1'], **change), BASELINE)


@pytest.mark.parametrize('change', [dict(rejection_scale=2.), dict(max_distance=33),
                                  dict(groups=(11,)), dict(rotation_weight=1)])
def test_counted_circuit_domain_is_explicit(change):
    with pytest.raises(ValueError, match='Counted'):
        star.min_qubits(4, 'shallow', HARDWARE['H1'], BASELINE, replace(STAR_BASELINE, **change))


def test_legacy_markov_and_noisy_corrections_remain_explicit():
    hw = HARDWARE['H1']
    row = star.evaluate_design(8, 'shallow', hw, BASELINE, STAR_LEGACY, 9, 3, 2)
    assert row['deadline_seconds'] == pytest.approx(row['run_seconds']/STAR_LEGACY.deadline_failure)
    assert row['loss_deadline'] > 0 and row['loss_preparation_memory'] > 0
