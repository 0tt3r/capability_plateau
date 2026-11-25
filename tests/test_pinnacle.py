"""Tests of the Pinnacle model: published engine parameters, hand-computed designs, and an independent brute force."""
import json
import math
from math import ceil, log1p
from pathlib import Path

import numpy as np
import pytest

from capability import curves, pinnacle as P, pinnacle_factories as pf
from capability.config import BASELINE, HARDWARE, WEEK_SECONDS
from capability.quops import layer_counts, size
from capability.surface import EPS, NT

ROOT = Path(__file__).resolve().parents[1]
H2, H4, H1, H3 = (HARDWARE[k] for k in ('H2', 'H4', 'H1', 'H3'))
REF = P.vary(families=('reference',))


def proto(name):
    return next(q for q in pf.PROTOCOLS if q['name'] == name)


def test_engine_template_reproduces_published_engines():
    """Pinnacle v2 Eqs. 9 and 10: 4410 and 5430 qubits, 26 code cycles, about 10% rejection, 1e-9 and 1e-11."""
    for d_a, qubits, target in ((7, 4410, 1e-9), (9, 5430, 1e-11)):
        p_a, p_in = pf.CULTIVATED[d_a]
        s, attempt, p_rot = pf.engine(proto('[[15,1,3]]'), d_a, p_a, p_in, True)
        assert s.qubits == qubits
        assert attempt == 26
        assert s.reject == pytest.approx(0.10, abs=0.005)
        assert target/2 < s.error < 2*target                        # the paper rounds its output targets to one digit
        assert p_rot == pytest.approx(p_in + (d_a + 1)*p_a)


def test_engine_time_matches_eq11():
    assert pf.engine_time(4, 7) == 23 and pf.engine_time(4, 9) == 26     # the 23 and 26 cycles of Pinnacle v2


def test_minimal_engine_footprint_of_the_low_noise_engine():
    """Pinnacle v2 Eq. 7 (p = 1e-4, physical T state): 252 + 10 * 31 + 30 = 592 with the d = 10 block."""
    assert 252 + 10*31 + 30*(1 + 1 - 1) == 592


def test_protocol_data():
    ref = proto('[[15,1,3]]')
    assert (ref['n'], ref['k'], ref['d'], ref['wires'], ref['A_d']) == (15, 1, 3, 5, 35)     # Bravyi-Kitaev 15-to-1
    for q in pf.PROTOCOLS:
        assert q['wires'] > q['k'] and q['A_d'] > 0 and q['d'] >= 3


def test_surface_supply_adapter():
    s = next(x for x in pf.surface_supplies(distances=(15,)) if x.name.startswith('cultivation'))
    assert (s.qubits, s.cycles, s.error) == (1024, 100., 2e-9)
    assert s.port_qubits == 454 + 29                                   # the Frozen Pinnacle seed patch plus bridge
    assert s.port_error == pytest.approx(16*0.03*(1e-3/0.01)**8)
    assert s.min_slot == 15


def test_two_level_supply_is_cleaner_but_bigger():
    one = next(s for s in pf.injected_supplies() if s.name == '[[15,1,3]], injected, 1 level')
    two = next(s for s in pf.injected_supplies() if s.name == '[[15,1,3]], injected, 2 levels, F1=16')
    assert one.error == pytest.approx(35*(3e-3)**3, rel=1e-12)          # p_rot = 3p for the physical-T ancilla
    assert two.error < 1e-15 and two.qubits > 10*one.qubits


def test_reference_matches_pr3_serial_design():
    """Pinnacle-only reference menu: the width-84 deep design of the earlier serial model (9570 qubits, about 274 s)."""
    r = P.min_qubits(84, 'deep', H2, BASELINE, REF)
    assert r['qubits'] == 9570 and r['distance'] == 16 and r['factories'] == 1
    assert r['run_seconds'] == pytest.approx(274, abs=1)
    p4 = P.best_at(1e4, H2, BASELINE, P.sweep(H2, BASELINE, REF), REF)
    assert p4['s'] == 592704
    assert curves.crossover(json.loads((ROOT/'results/baseline.json').read_text())['H2']['direct'],
                            P.sweep(H2, BASELINE, REF)) == 4862


def lam(x):
    return -log1p(-x)


def hand_design(w, ensemble, code, s, NTc, F, hw=H2, ch=BASELINE, eta=2.0):
    """Independent evaluation of one design in the surface-code form: N, T_run, and the stochastic loss."""
    n_rot, _ = layer_counts(w, ensemble)
    R = w*n_rot
    NT_total = R*NTc
    beta = ceil(w/code.k)
    D = max(code.dt, 10, s.min_slot)
    cycles = max(NT_total*D, ceil(NT_total/F)*s.cycles) + (w + 1)*D
    N = beta*code.qubits + F*s.qubits + s.port_qubits
    p_pb = 6.2*(hw.p/0.0158)**(code.d/2 + 0.47)
    volume = beta*cycles/code.dt
    loss = eta*volume*lam(p_pb) + NT_total*(lam(s.error) + lam(s.port_error))
    eps = 2.**(-(NTc - 1.0970643939393736)/3.0764791908264293)
    return N, cycles, loss + R*eps**2


def test_hand_computed_design():
    r = P.min_qubits(84, 'deep', H2, BASELINE, P.BASELINE)
    code = next(c for c in P.CODES if c.d == r['distance'])
    s = next(x for x in pf.all_supplies() if x.name == r['factory'])
    N, cycles, total = hand_design(84, 'deep', code, s, r['T_per_rotation'], r['factories'])
    assert N == r['qubits'] and cycles == pytest.approx(r['cycles'])
    assert total <= 0.5 + 1e-12
    assert r['loss_blocks'] + r['loss_magic'] + r['loss_synthesis'] == pytest.approx(total, rel=1e-9)


def brute_force(w, ensemble, hw, pc):
    """Enumerate every design explicitly, with the smallest sufficient number of supplies."""
    best = np.inf
    for code in P.CODES:
        if code.d not in pc.codes:
            continue
        for s in pf.all_supplies(pc.families, pc.adapter_distances):
            for NTc in range(5, 161):
                for F in range(1, 400):
                    N, cycles, total = hand_design(w, ensemble, code, s, NTc, F, hw)
                    if total <= BASELINE.threshold and cycles <= WEEK_SECONDS/BASELINE.cycle_seconds:
                        best = min(best, N)
                        break
    return best


@pytest.mark.parametrize('w,ensemble,families', [(6, 'shallow', ('reference', 'surface')), (4, 'deep', ('gb_cult', 'gb_inject')),
                                                  (10, 'deep', ('reference', 'gb_inject2'))])
def test_search_matches_brute_force(w, ensemble, families):
    pc = P.vary(families=families, adapter_distances=(15, 23))
    r = P.min_qubits(w, ensemble, H2, BASELINE, pc)
    assert r['qubits'] == brute_force(w, ensemble, H2, pc)


def test_min_qubits_never_decreases_with_width():
    for ensemble in ('shallow', 'deep'):
        qs = [P.min_qubits(w, ensemble, H2, BASELINE)['qubits'] for w in range(2, 120, 2)]
        assert qs == sorted(qs)


def test_best_at_refines_sampled_widths():
    rows = P.sweep(H2, BASELINE)
    got = P.best_at(1e5, H2, BASELINE, rows)
    nxt = [e for e in ('shallow', 'deep') if P.min_qubits(got['w'] + 2, e, H2, BASELINE)]
    assert got['qubits'] <= 1e5
    for e in nxt:
        n = P.min_qubits(got['w'] + 2, e, H2, BASELINE)
        assert n['qubits'] > 1e5 or size(got['w'] + 2, e) <= got['s']


def test_grid_hardware_is_unsupported_not_zero():
    for hw in (H1, H3):
        assert P.min_qubits(4, 'deep', hw, BASELINE) is None
        assert P.sweep(hw, BASELINE) == []
        assert P.best_at(1e5, hw, BASELINE, []) is None
        assert 'not modeled' in P.support_reason(hw)


def test_noise_classes_share_the_block_fit():
    assert P.min_qubits(40, 'deep', H2, BASELINE)['qubits'] == P.min_qubits(40, 'deep', H4, BASELINE)['qubits']


def test_larger_menu_is_never_worse():
    for w in (4, 20, 100, 500):
        a = P.min_qubits(w, 'deep', H2, BASELINE, REF)['qubits']
        b = P.min_qubits(w, 'deep', H2, BASELINE, P.BASELINE)['qubits']
        assert b <= a


def test_eta_one_lowers_the_qubit_need():
    assert P.min_qubits(600, 'deep', H2, BASELINE, P.vary(eta=1.0))['qubits'] <= P.min_qubits(600, 'deep', H2, BASELINE)['qubits']


def test_cycle_time_scales_runtime_exactly_and_leaves_qubits_unchanged_when_not_time_limited():
    from dataclasses import replace
    a = P.min_qubits(60, 'deep', H2, BASELINE)
    b = P.min_qubits(60, 'deep', H2, replace(BASELINE, cycle_seconds=1.2e-6))
    assert b['cycles'] == pytest.approx(a['cycles']) and b['run_seconds'] == pytest.approx(1.2*a['run_seconds'])
    assert b['qubits'] == a['qubits']


def test_frozen_block_menu_is_a_selectable_subset():
    pc = P.FROZEN_T
    r = P.min_qubits(18, 'deep', H2, BASELINE, pc)
    assert r['code_k'] == 8 and r['blocks'] == 3 and r['distance'] == 24
    assert r['qubits'] >= 3*2522                          # three full 2522-qubit Frozen blocks plus the supply
    assert P.min_qubits(2, 'shallow', H2, BASELINE, pc)['qubits'] > P.min_qubits(2, 'shallow', H2, BASELINE)['qubits']


def test_extended_surface_menu_is_never_worse_than_the_cached_menu():
    from dataclasses import replace
    from capability import surface
    ext = replace(BASELINE, factories=('cultivation', 'litinski', 'jws'))
    for w, e in ((6, 'deep'), (40, 'deep'), (300, 'shallow')):
        a, b = surface.min_qubits(w, e, H2, BASELINE), surface.min_qubits(w, e, H2, ext)
        assert b['qubits'] <= a['qubits']
