"""Checks of the two models against brute force or against their own definitions."""
from math import exp, log

import numpy as np
import pytest

from capability import curves, direct, quops, routing, surface
from capability.config import BASELINE, HARDWARE, WEEK_SECONDS, vary


@pytest.mark.parametrize('router', ['snake', 'rcr'])
@pytest.mark.parametrize('w', [2, 4, 6, 8, 12, 20, 36])
def test_routers_make_every_pair_adjacent(router, w):
    rng = np.random.default_rng(w)
    r, c = routing.grid_shape(w)
    placement = np.arange(r*c)
    for _ in range(20):
        pairs = quops.random_pairing(w, rng)
        rounds, placement = routing.ROUTERS[router](placement, pairs, r, c)
        assert routing.pairs_adjacent(placement, pairs, c)
        assert len(rounds) <= (r*c if router == 'snake' else 2*r + c)
        for swaps in rounds:
            tokens = [t for ab in swaps for t in ab]
            assert len(tokens) == len(set(tokens))           # SWAPs in a round are disjoint


def test_grid_shape_fits_and_has_even_columns():
    for w in range(2, 200, 2):
        r, c = routing.grid_shape(w)
        assert r*c >= w and c % 2 == 0


@pytest.mark.parametrize('compile_mode', ['euler', 'native'])
@pytest.mark.parametrize('idle', [True, False])
def test_rotation_layer_matches_brute_force(compile_mode, idle):
    hw = HARDWARE['H1']
    ch = vary(rotation_compile=compile_mode, idle_errors=idle)
    w, rng, eps = 6, np.random.default_rng(0), HARDWARE['H1'].eps_1q
    values = []
    for _ in range(200000):
        steps = [quops.native_steps(a, compile_mode) for a in rng.choice(list('XYZ'), w)]
        locations = w*max(steps) if idle else sum(steps)
        values.append((1 - eps)**locations)
    assert exp(direct.log_rotation_layer(w, hw, ch)) == pytest.approx(np.mean(values), rel=2e-4)


def test_two_qubits_need_no_routing():
    for grid, a2a in (('H1', 'H2'), ('H3', 'H4')):
        for ensemble in quops.ENSEMBLES:
            q_grid, _ = direct.mean_no_fault(2, ensemble, HARDWARE[grid], BASELINE)
            q_a2a, _ = direct.mean_no_fault(2, ensemble, HARDWARE[a2a], BASELINE)
            assert q_grid == pytest.approx(q_a2a, rel=1e-12)


def test_all_to_all_closed_form_by_hand():
    """w = 4 shallow on H2: 4 layers = 2 rotation layers (D = 5 unless no qubit has Y) and 2 CNOT layers of
    2 CNOTs, plus 4 preparations and 4 measurements."""
    p = HARDWARE['H2'].eps_2q
    pD = {1: (1/3)**4, 3: (2/3)**4 - (1/3)**4, 5: 1 - (2/3)**4}
    rot = sum(pr*(1 - p)**(4*D) for D, pr in pD.items())
    expected = (1 - p)**8*(1 - p)**4*rot**2
    assert direct.mean_no_fault(4, 'shallow', HARDWARE['H2'], BASELINE)[0] == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize('key', ['H1', 'H3'])
@pytest.mark.parametrize('w,ensemble', [(4, 'shallow'), (12, 'deep'), (100, 'shallow'), (300, 'deep')])
@pytest.mark.parametrize('executions', [1, 10_000])
def test_surface_design_meets_its_own_constraints(key, w, ensemble, executions):
    hw, ch = HARDWARE[key], vary(executions=executions)
    x = surface.min_qubits(w, ensemble, hw, ch)
    n_rot, n_cnot = quops.layer_counts(w, ensemble)
    d, F = x['distance'], x['factories']
    f = next(f for f in surface.menu(ch.factories) if f.name == x['factory'])
    nT = x['T_per_rotation']
    cycles = n_rot*max(nT*d, np.ceil(w*nT/F)*f.cycles_per_output) + n_cnot*2*d + 2*d
    assert cycles == pytest.approx(x['cycles'])
    loss = (x['tiles']*cycles*-log(1 - surface.p_sc(d, hw, ch)) + w*n_rot*nT*-log(1 - f.error)
            + w*n_rot*x['eps']**2)
    assert loss <= ch.threshold + 1e-12
    assert x['copies']*min(WEEK_SECONDS//(cycles*ch.cycle_seconds), executions) >= executions
    assert x['qubits'] == x['copies']*(x['tiles']*2*d*d + F*f.qubits)


def brute_force_min_qubits(w, ensemble, hw, ch):
    """Every factory type, odd distance, T count, and factory count up to saturation, with the copy count
    set directly by K floor(week / run time) >= executions."""
    n_rot, n_cnot = quops.layer_counts(w, ensemble)
    R, tiles, M = w*n_rot, int(np.ceil(ch.tiles_per_data_patch*w)), ch.executions
    week = WEEK_SECONDS/ch.cycle_seconds
    nT, eps = surface.NT[:, None], surface.EPS[:, None]
    F = np.arange(1, w*surface.NT.max() + 1)[None, :]
    best = np.inf
    for f in surface.menu(ch.factories):
        for d in range(3, 61, 2):
            lam = -log(1 - surface.p_sc(d, hw, ch))
            cycles = n_rot*np.maximum(nT*d, np.ceil(w*nT/F)*f.cycles_per_output) + n_cnot*2*d + 2*d
            loss = tiles*cycles*lam + R*nT*-log(1 - f.error) + R*eps**2
            per_copy = np.minimum(np.floor(week/cycles), M)
            ok = (loss <= ch.threshold) & (per_copy >= 1)
            qubits = np.ceil(M/np.maximum(per_copy, 1))*(tiles*2*d*d + F*f.qubits)
            best = min(best, np.where(ok, qubits, np.inf).min())
    return best


@pytest.mark.parametrize('key,w,ensemble,change', [
    ('H1', 8, 'shallow', {}),
    ('H3', 10, 'deep', {}),
    ('H1', 8, 'deep', dict(executions=100_000, cycle_seconds=1e-3)),
    ('H3', 12, 'deep', dict(executions=100_000, cycle_seconds=1e-3, factories=('litinski',))),
])
def test_design_search_matches_brute_force(key, w, ensemble, change):
    hw, ch = HARDWARE[key], vary(**change)
    x = surface.min_qubits(w, ensemble, hw, ch)
    assert x['qubits'] == brute_force_min_qubits(w, ensemble, hw, ch)
    if change:
        assert x['copies'] > 1                     # the case exercises the copies trade


def test_surface_qubits_never_decrease_with_width():
    for key in ('H1', 'H3'):
        for ensemble in quops.ENSEMBLES:
            q = [surface.min_qubits(w, ensemble, HARDWARE[key], BASELINE)['qubits'] for w in range(2, 161, 2)]
            assert np.all(np.diff(q) >= 0)


def test_capability_probe_is_exact():
    hw = HARDWARE['H3']
    rows = surface.sweep(hw, BASELINE)
    x = surface.best_at(1e6, hw, BASELINE, rows)
    nxt = surface.min_qubits(x['w'] + 2, x['ensemble'], hw, BASELINE)
    assert x['qubits'] <= 1e6 < nxt['qubits']
    assert x['s'] >= max(r['s'] for r in rows if r['qubits'] <= 1e6)


def test_crossover_is_a_step_of_the_surface_curve():
    for key in ('H1', 'H4'):
        hw = HARDWARE[key]
        d, s = direct.sweep(hw, BASELINE), surface.sweep(hw, BASELINE)
        N = curves.crossover(d, s)
        assert curves.curve(s, [N])[0] > curves.curve(d, [N])[0]
        assert curves.curve(s, [N - 1])[0] <= curves.curve(d, [N - 1])[0]


def test_t_count_follows_the_fit():
    eps = np.array([1e-2, 1e-6, 1e-10])
    assert np.all(surface.t_count(eps) == np.ceil(surface.FIT['a']*np.log2(1/eps) + surface.FIT['b']))
    assert surface.FIT['max_residual'] < 1.5


def test_capability_curve_is_nondecreasing():
    rows = [dict(passes=True, qubits=q, s=s) for q, s in ((10, 4), (100, 64), (50, 16), (1000, 8))]
    y = curves.curve(rows)
    assert np.all(np.diff(y) >= 0) and y[-1] == 64
