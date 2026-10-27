"""Independent checks of STAR channel formulas, stochastic timing and design constraints."""
from dataclasses import replace
from math import exp, log

import numpy as np
import pytest

from capability import star
from capability.config import BASELINE, HARDWARE, STAR_BASELINE, STAR_SINGLE_QUBIT, WEEK_SECONDS
from capability.star_states import (fallback_parameters, fold_angle, ideal_tmr,
                                    repository_angle, rus_rotation_loss, state_profile,
                                    uniform_angles)
from capability.star_states import preparation_survival, undetected_probability
from capability.surface import p_sc


def test_rotation_convention_and_clifford_folding():
    theta = np.array([0., .1, -.3, np.pi/2, 2*np.pi])
    residual = repository_angle(theta)
    assert np.all(np.abs(residual) <= np.pi/8)
    multiples = (-theta/2 - residual)/(np.pi/4)
    assert multiples == pytest.approx(np.round(multiples))
    assert residual[:3] == pytest.approx([0., -.05, .15])
    assert residual[3:] == pytest.approx([0., 0.], abs=1e-14)


@pytest.mark.parametrize('k', [3, 5, 7, 9])
def test_tmr_projection_matches_state_amplitudes(k):
    targets = np.array([0., 1e-8, .001, .1, .35])
    alpha, acceptance, erroneous, angle_error = ideal_tmr(targets, k)
    good_amplitudes = np.array([np.cos(alpha)**k, np.sin(alpha)**k])
    assert acceptance == pytest.approx(np.sum(good_amplitudes**2, axis=0))
    assert np.arctan2(good_amplitudes[1], good_amplitudes[0]) == pytest.approx(targets)
    bad_amplitudes = np.array([np.sin(alpha)*np.cos(alpha)**(k-1),
                               np.sin(alpha)**(k-1)*np.cos(alpha)])
    assert erroneous == pytest.approx(np.sum(bad_amplitudes**2, axis=0))
    nonzero = targets > 0
    assert angle_error[nonzero] == pytest.approx(-np.arctan2(bad_amplitudes[1, nonzero],
                                                           bad_amplitudes[0, nonzero]))


def test_corrected_channel_matches_density_matrix_calculation():
    sc = replace(STAR_BASELINE, protocol_switching=False, preparation_floor_scale=0.)
    profile = state_profile(np.array([.001, .1, .3]), 3, 9, HARDWARE['H1'], BASELINE, sc)
    rho = np.array([[.5, .5], [.5, .5]], complex)
    z = np.diag([1., -1.])
    pf = profile.fallback_error
    for b, delta, error in zip(profile.conditional_bad, profile.delta, profile.error):
        u = np.diag(np.exp(1j*delta*np.array([1., -1.])))
        raw = (1-b)*rho + b*u @ rho @ u.conj().T
        corrected = (1-b)*raw + b*u.conj().T @ raw @ u
        # Remove the separately charged bound on imperfections in the correction itself.
        pauli = error - b*2*pf/(1+pf)
        expected = (1-pauli)*rho + pauli*z @ rho @ z
        assert corrected == pytest.approx(expected, abs=1e-14)


def test_paper_virtual_z_coefficient_and_hardware_extension():
    virtual = replace(HARDWARE['H1'], eps_1q=0.)
    sc = replace(STAR_BASELINE, protocol_switching=False)
    phi = np.array([.001])
    p = state_profile(phi, 3, 9, virtual, BASELINE, sc)
    _, ideal, erroneous, _ = ideal_tmr(phi, 3)
    u = 2*3*virtual.eps_2q/15
    assert p.conditional_bad == pytest.approx(u*erroneous/((1-u)*ideal+u*erroneous))
    h1 = state_profile(phi, 3, 9, HARDWARE['H1'], BASELINE, sc)
    assert h1.error[0] > p.error[0]
    assert fallback_parameters(HARDWARE['H1'], 9, sc)[1] == pytest.approx(9e-3/15)


@pytest.mark.parametrize('acceptance', [[1.], [.25], [.25, .75]])
def test_preparation_order_statistics_against_geometric_identities(acceptance):
    a = np.array(acceptance)
    times = star.batch_preparation_rounds(a, 2, tolerance=1e-11)
    # Sum geometric series analytically for b=1,2, including all RUS rounds.
    expected_one = 2*np.mean(1/a)
    cross = a[:, None] + a[None, :] - a[:, None]*a[None, :]
    expected_two = 4*np.mean(1/a) - (4/3)*np.mean(1/cross)
    assert times[1] == pytest.approx(expected_one, abs=1e-10)
    assert times[2] == pytest.approx(expected_two, abs=1e-10)
    assert times[1] >= expected_one - 1e-12


def test_parallel_rus_is_not_two_rounds_per_layer():
    rounds = star.batch_rus_rounds(32)
    assert rounds[1] == pytest.approx(2.)
    assert rounds[2] == pytest.approx(8/3)
    assert 6 < rounds[32] < 7
    assert np.all(np.diff(rounds) > 0)


def test_extreme_rejection_uses_finite_conservative_timing_bound():
    a = 1e-6
    times = star.batch_preparation_rounds([a], 2)
    exact_two = 4/a - (4/3)/(a*(2-a))
    assert times[1] == pytest.approx(2/a)
    assert exact_two <= times[2] <= 4/a


def test_uniform_angles_and_rus_average_agree():
    sc = replace(STAR_BASELINE, preparation_floor_scale=0., angle_points=4096)
    # Irrational offset avoids artificial early Clifford termination on a dyadic grid.
    phi = (np.arange(4096) + np.sqrt(2)/2)*(np.pi/8)/4096
    profile = state_profile(phi, 3, 9, HARDWARE['H1'], BASELINE, sc)
    expected = 2*np.mean(-np.log1p(-profile.error))
    measured = np.mean(rus_rotation_loss(phi, 3, 9, HARDWARE['H1'], BASELINE, sc, rounds=32))
    assert measured == pytest.approx(expected, rel=3e-4)
    assert rus_rotation_loss(0., 3, 9, HARDWARE['H1'], BASELINE, sc) == 0.


@pytest.mark.parametrize('key', ['H1', 'H3'])
def test_design_loss_and_deadline_constraints(key):
    x = star.min_qubits(8, 'shallow', HARDWARE[key], BASELINE)
    assert x is not None
    loss = sum(x[name] for name in ('loss_patches', 'loss_injection',
                                   'loss_preparation_memory', 'loss_deadline'))
    assert loss == pytest.approx(x['loss_total'])
    assert loss <= BASELINE.threshold
    assert exp(-loss) == pytest.approx(x['success_probability_lower_bound'])
    assert x['deadline_seconds'] == pytest.approx(x['run_seconds']/STAR_BASELINE.deadline_failure)
    assert x['copies']*int(WEEK_SECONDS/x['deadline_seconds']) >= BASELINE.executions
    assert x['resource_states'] >= 8*4*2
    assert x['preparation_attempts'] >= x['resource_states']
    area = ((x['data_routing_tiles'] + x['handoff_tiles'])*2*x['distance']**2
            + x['stations']*x['preparation_qubits_per_station'])
    assert x['qubits'] == area*x['copies']


@pytest.mark.parametrize('settings', [STAR_BASELINE, STAR_SINGLE_QUBIT])
def test_search_matches_unpruned_design_enumeration(settings):
    sc = replace(settings, groups=(3, 5), max_distance=9, angle_points=64)
    hw, w = HARDWARE['H1'], 6
    designs = [star.evaluate_design(w, 'shallow', hw, BASELINE, sc, d, k, f)
               for d in range(3, 10, 2) for k in sc.groups for f in range(1, w+1)]
    expected = min(x['qubits'] for x in designs if x['passes'])
    assert star.min_qubits(w, 'shallow', hw, BASELINE, sc)['qubits'] == expected


def test_noise_floor_rejects_large_workload_before_timing():
    hw = HARDWARE['H1']
    assert star.resource_ceiling(hw, BASELINE) < 1200
    assert star.min_qubits(100, 'deep', hw, BASELINE) is None


@pytest.mark.parametrize('settings', [STAR_BASELINE, STAR_SINGLE_QUBIT])
def test_zero_noise_limits(settings):
    hw = replace(HARDWARE['H1'], eps_1q=0., eps_2q=0., p=0.)
    sc = replace(settings, protocol_switching=False)
    phi = uniform_angles(32)
    profile = state_profile(phi, 3, 9, hw, BASELINE, sc)
    assert profile.error == pytest.approx(np.zeros(32))
    assert profile.acceptance == pytest.approx(ideal_tmr(phi, 3)[1])
    assert profile.preparation_floor_loss == 0.


def test_hardware_and_model_domains_are_explicit():
    with pytest.raises(ValueError, match='2D grid'):
        star.min_qubits(4, 'shallow', HARDWARE['H2'], BASELINE)
    with pytest.raises(ValueError, match='odd'):
        ideal_tmr(.1, 4)
    with pytest.raises(ValueError, match='acceptance'):
        star.batch_preparation_rounds([0.], 2)


@pytest.mark.parametrize('k', [3, 5])
def test_single_qubit_projection_in_stabilizer_subspace(k):
    # A small stabilizer code with logical Z = Z^k and code basis |+>^k, |->^k.
    # Apply the physical unitaries and project, independently of the scalar formulas.
    plus, minus = np.array([1., 1.])/np.sqrt(2), np.array([1., -1.])/np.sqrt(2)
    z = np.diag([1., -1.])
    phi = .12
    alpha, a0, ae, phi_error = ideal_tmr(phi, k)
    rotation = np.diag(np.exp(1j*alpha*np.array([1., -1.])))
    basis_plus = basis_minus = state = faulty = np.array([1.], complex)
    for i in range(k):
        basis_plus = np.kron(basis_plus, plus)
        basis_minus = np.kron(basis_minus, minus)
        rotated = rotation @ plus
        state = np.kron(state, rotated)
        faulty = np.kron(faulty, z @ rotated if i == 0 else rotated)
    good = np.array([np.vdot(basis_plus, state), np.vdot(basis_minus, state)])
    bad = np.array([np.vdot(basis_plus, faulty), np.vdot(basis_minus, faulty)])
    phase = 1j**k
    assert np.vdot(good, good).real == pytest.approx(a0)
    assert np.vdot(bad, bad).real == pytest.approx(ae)
    assert (good[1]/good[0]/phase).real == pytest.approx(np.tan(phi))
    assert (bad[1]/bad[0]/phase).real == pytest.approx(np.tan(phi_error))


def test_single_qubit_exposure_model_includes_syndrome_errors():
    # Count marginal Z terms among the 15 nonidentity two-qubit Paulis.
    paulis = [(x, y) for x in 'IXYZ' for y in 'IXYZ' if (x, y) != ('I', 'I')]
    marginal_z = sum(x == 'Z' for x, _ in paulis)/len(paulis)
    sc, k = STAR_SINGLE_QUBIT, 5
    hw = replace(HARDWARE['H3'], eps_1q=0.)
    u = undetected_probability(k, hw, sc)
    assert u == pytest.approx(k*sc.single_qubit_2q_locations*marginal_z*hw.eps_2q)
    assert u > undetected_probability(k, hw, STAR_BASELINE)
    assert undetected_probability(k, hw, replace(sc, single_qubit_undetected_scale=2.)) == pytest.approx(2*u)


@pytest.mark.parametrize('key', ['H1', 'H3'])
def test_single_qubit_preparation_distance_area_and_error_fit(key):
    sc, hw, k, d = STAR_SINGLE_QUBIT, HARDWARE[key], 5, 9
    profile = state_profile(uniform_angles(32), k, d, hw, BASELINE, sc)
    assert profile.preparation_distance == k
    assert profile.preparation_floor_loss == pytest.approx(-4*np.log1p(-p_sc(k, hw, BASELINE)))
    row = star.evaluate_design(4, 'shallow', hw, BASELINE, sc, d, k, 2)
    assert row['rotation_weight'] == 1
    assert row['preparation_qubits_per_station'] == (2*k - 1)**2
    assert row['qubits_per_copy'] == (8 + 2)*2*d*d + 2*(2*k - 1)**2
    assert row['passes']
    assert preparation_survival(9, hw, sc) < preparation_survival(9, hw, replace(sc, postselection='strip'))


@pytest.mark.parametrize('changes, message', [
    ({'rotation_weight': 3}, 'rotation_weight'),
    ({'postselection': 'unknown'}, 'postselection'),
    ({'single_qubit_undetected_scale': -1.}, 'nonnegative'),
    ({'single_qubit_undetected_scale': 100.}, 'leading-order'),
])
def test_single_qubit_model_domains(changes, message):
    sc = replace(STAR_SINGLE_QUBIT, **changes)
    with pytest.raises(ValueError, match=message):
        star.evaluate_design(4, 'shallow', HARDWARE['H1'], BASELINE, sc, 9, 3, 1)
