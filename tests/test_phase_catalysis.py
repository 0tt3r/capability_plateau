"""Check the ideal phase-catalysis circuit independently of resource formulas.

Gidney--Fowler, arXiv:1812.01238v3, Fig. 14, at theta = pi/4.
One AND computation uses one Toffoli, and its uncomputation uses one
scratch X measurement with a conditional CZ. Reusing the scratch also
requires preparation/reset. The remaining gates are Cliffords; treating
them as frame updates is a separate resource-model assumption.

These state-vector checks establish the circuit identity, including for
data entangled with an untouched reference. They do not establish a
fixed-connectivity implementation or the timing of a hardware schedule.
"""
import numpy as np
import pytest


PHASE = np.exp(1j*np.pi/4)


def _controlled_x(state, target, controls=()):
    """Apply X, CNOT, or Toffoli as a computational-basis permutation."""
    indices = np.arange(state.size)
    active = np.ones(state.size, dtype=bool)
    for control in controls:
        active &= (indices & (1 << control)) != 0
    source = indices ^ (active.astype(np.int64) << target)
    return state[source]


def _initial_state(data_reference):
    """Qubit order, low to high: a, b, catalyst, scratch, reference."""
    state = np.zeros(32, dtype=complex)
    for reference in range(2):
        for catalyst in range(2):
            for data in range(4):
                address = data | (catalyst << 2) | (reference << 4)
                state[address] = data_reference[data | (reference << 2)]*PHASE**catalyst/np.sqrt(2)
    return state


def _run_pair(data_reference, outcome):
    """Explicit Fig. 14 gates; return the unnormalized selected branch."""
    state = _initial_state(data_reference)
    state = _controlled_x(state, 2)
    state = _controlled_x(state, 1, (0,))
    state = _controlled_x(state, 2, (0,))
    state = _controlled_x(state, 3, (1, 2))       # compute AND into fresh |0>
    state = _controlled_x(state, 3, (0,))
    state[(np.arange(32) & 8) != 0] *= 1j        # Z^(2 theta) = S
    state = _controlled_x(state, 3, (0,))

    # Project the scratch onto |+> or |->, then remove that qubit.
    # The remaining qubit order is a, b, catalyst, reference.
    addresses = (np.arange(16) & 7) | ((np.arange(16) & 8) << 1)
    branch = (state[addresses] + (-1)**outcome*state[addresses | 8])/np.sqrt(2)
    if outcome:
        indices = np.arange(16)
        branch[((indices & 2) != 0) & ((indices & 4) != 0)] *= -1

    branch = _controlled_x(branch, 2, (1,))
    branch = _controlled_x(branch, 1, (0,))
    branch = _controlled_x(branch, 2, (0,))
    return _controlled_x(branch, 2)


def _ideal_pair(data_reference):
    """T on both data qubits, identity on the reference, original T catalyst."""
    ideal = np.zeros(16, dtype=complex)
    for reference in range(2):
        for catalyst in range(2):
            for data in range(4):
                phases = (data & 1) + ((data >> 1) & 1) + catalyst
                address = data | (catalyst << 2) | (reference << 3)
                ideal[address] = data_reference[data | (reference << 2)]*PHASE**phases/np.sqrt(2)
    return ideal


def _check_pair(data_reference, outcome):
    branch = _run_pair(data_reference, outcome)
    probability = float(np.vdot(branch, branch).real)
    assert probability == pytest.approx(0.5, abs=2e-14)
    actual = branch/np.sqrt(probability)
    expected = _ideal_pair(data_reference)
    # Outcome-dependent global phase cannot affect the quantum channel.
    overlap = np.vdot(expected, actual)
    assert abs(overlap) == pytest.approx(1.0, abs=2e-14)
    np.testing.assert_allclose(actual, expected*overlap/abs(overlap), atol=2e-14, rtol=0)


@pytest.mark.parametrize('outcome', [0, 1])
@pytest.mark.parametrize('basis', range(8))
def test_phase_catalysis_preserves_the_catalyst_for_basis_inputs(basis, outcome):
    data_reference = np.zeros(8, dtype=complex)
    data_reference[basis] = 1
    _check_pair(data_reference, outcome)


@pytest.mark.parametrize('outcome', [0, 1])
@pytest.mark.parametrize('seed', [0, 7, 19, 101, 510])
def test_phase_catalysis_acts_correctly_on_entangled_data_and_reference(seed, outcome):
    rng = np.random.default_rng(seed)
    data_reference = rng.normal(size=8) + 1j*rng.normal(size=8)
    data_reference /= np.linalg.norm(data_reference)
    _check_pair(data_reference, outcome)
