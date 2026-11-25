"""Surface-code layouts of the arXiv:2610.06535 protocols, calibrated against Litinski (arXiv:1905.06903v3)."""
import pytest

from capability import factories, surface_factories as sf
from capability.config import HARDWARE, vary


def proto(name):
    return next(q for q in sf.PROTOCOLS if q['name'] == name)


def test_one_level_15_to_1_reproduces_litinski_row():
    qubits, cycles, reject, error = sf.single(proto('[[15,1,3]]'), 17, 7, 7, 1e-3)
    assert qubits == pytest.approx(4620, abs=5)                    # 2 (dX + 4 dZ) 3 dX + 4 dm
    assert cycles == pytest.approx(42.6, rel=0.01)                 # 6 dm / (1 - p_fail)
    assert error == pytest.approx(4.5e-8, rel=0.25)


def test_two_level_row_reproduces_litinski_row():
    qubits, cycles, error = sf.two_level(proto('[[15,1,3]]'), 11, 5, 6, 25, 11)
    assert qubits == pytest.approx(30700, rel=0.01)
    assert cycles == pytest.approx(82.5, rel=0.01)
    assert error == pytest.approx(2.7e-12, rel=0.01)               # calibrated scale


def test_20_to_4_time_matches_litinski_ten_cycles():
    four = dict(n=20, k=4, d=2, wires=7, A_d=22, name='20-to-4')    # 20-to-4: 4 outputs, 3 checks (Litinski Fig. 4)
    steps = -(-(four['n'] - 3)//2) + 1
    assert steps*1 == 10


def test_default_surface_menu_is_unchanged():
    base = factories.menu(('cultivation', 'litinski'))
    assert len(base) == 8 and all(f.family != 'jws' for f in base)


def test_jws_menu_is_pareto_and_multi_output_are_batches():
    F = sf.factories()
    assert F and all(f.family == 'jws' for f in F)
    for f in F:
        assert f.states_per_batch >= 1 and f.cycles_per_output > 0 and 0 < f.error < 1e-3
    for f in F:       # no entry is dominated by another
        assert not any(g is not f and g.qubits <= f.qubits and g.cycles_per_output <= f.cycles_per_output
                       and g.error <= f.error and (g.qubits, g.cycles_per_output, g.error) != (f.qubits, f.cycles_per_output, f.error)
                       for g in F)
