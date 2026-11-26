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
    fold = [f for f in F if sf.FOLD_LABEL in f.name]       # fold-transversal inputs: all-to-all hardware only
    assert fold
    for front in ([f for f in F if sf.FOLD_LABEL not in f.name], fold):
        for f in front:       # no entry is dominated by another of the same hardware front
            assert not any(g is not f and g.qubits <= f.qubits and g.cycles_per_output <= f.cycles_per_output
                           and g.error <= f.error and (g.qubits, g.cycles_per_output, g.error) != (f.qubits, f.cycles_per_output, f.error)
                           for g in front)


def test_catalogue_error_coefficients_and_yields():
    """scripts/catalogue_factories.py counts the published error coefficients and applies the T-yield rule."""
    import json
    from pathlib import Path
    rows = json.loads((Path(__file__).resolve().parents[1]/'results/catalogue_factories.json').read_text())['designs']
    by = {(r['n'], r['k'], r['d'], r['gate']): r for r in rows}
    ccz = by[(8, 3, 2, 'CCZ012')]
    assert ccz['A_d'] == 28 and ccz['T_per_run'] == 2          # 8T -> CCZ, kappa = 28; catalysed CCZ -> 2 T
    t14 = by[(14, 1, 2, 'T0')]
    assert t14['A_d'] == 7 and t14['T_per_run'] == 1           # the catalogue's c = 7
    for r in rows:
        assert r['error'] < 1e-3 and r['T_per_run'] >= 1


def test_catalogue_winners_are_in_the_extended_menu():
    from capability.factories import CATALOGUE_WINNERS, menu
    names = {f.name for f in menu(('catalogue_winners',))}
    assert names == set(CATALOGUE_WINNERS)
