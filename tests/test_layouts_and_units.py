from dataclasses import replace
from math import ceil, sqrt

import numpy as np
import pytest

from capability import layouts, surface, surface_units, surface_pbc, pbc_compile
from capability.config import BASELINE,HARDWARE,PBC,PBC_COLUMNS,PBC_DRAFT_BOUND,PBC_UNITS
from capability.factories import Factory,FACTORIES


def test_explicit_tile_geometry_and_rounding():
    for w in range(2,402,2):
        layout=layouts.fast_layout(w)
        c=layout.certificate()
        assert c['disjoint'] and c['connected_ancilla'] and c['both_pauli_edges_accessible']
        assert c['algorithm_qubits'] == c['data_tiles'] == w
        assert c['tiles'] == 2*w+2*(max(layout.heights)+len(layout.heights))+1
        assert c['tiles']-(2*w+ceil(sqrt(8*w))+1) in (0,1)
    assert layouts.fast_layout(6).area == 21
    assert layouts.fast_layout(10).area == 31
    assert layouts.fast_layout(18).area == 49


def test_old_ladder_does_not_supply_arbitrary_disjoint_parallel_routes():
    assert layouts.ladder_route_depth([(0,1),(2,3)]) == 1
    assert layouts.ladder_route_depth([(0,2),(1,3)]) == 2
    assert layouts.ladder_route_depth([(0,5),(1,4),(2,3)]) == 3


def test_pipeline_matches_independent_wave_events():
    for layers in (1,2,7,30,103):
        for units in (2,3,6):
            for prepare in (2.,100.):
                for reaction in (1.,19.):
                    # Direct event scheduling: a wave cannot launch before the
                    # prior wave's correction slots can be retired at completion.
                    start,readout=0.,0.
                    left=layers
                    while left:
                        ready=start+prepare
                        begin=max(ready,readout)
                        count=min(units-1,left)
                        readout=begin+count*reaction
                        start=begin
                        left-=count
                    assert surface_units.pipeline_cycles(layers,prepare,units,reaction,7) == readout+7


def test_parallel_schedule_has_no_cnot_term_and_charges_all_units():
    f=FACTORIES[0]
    r=surface_units.evaluate_design(8,'shallow',HARDWARE['H1'],BASELINE,13,f,30,2,4)
    assert r['factories'] == 8 and r['parallel_layers'] == 3
    assert r['tiles'] == 4*(r['data_tiles_per_unit']+r['correction_tiles_per_unit']+r['buffer_tiles_per_unit'])
    assert r['qubits_per_copy'] == 2*13**2*r['tiles']+8*f.qubits
    assert r['cnot_cycles'] == 0
    assert r['cycles'] == surface_units.pipeline_cycles(
        r['T_layers'],r['unit_cycles'],r['units'],PBC_UNITS.reaction_cycles,8*13)


@pytest.mark.parametrize('multiple',[False,True])
def test_units_optimizer_matches_unpruned_menu(monkeypatch,multiple):
    nt=np.array([12,20,28])
    monkeypatch.setattr(surface,'NT',nt)
    monkeypatch.setattr(surface,'EPS',2.**(-(nt-surface.FIT['b'])/surface.FIT['a']))
    f=Factory('test','test',2e-9,1024,100)
    monkeypatch.setattr(surface_units,'menu',lambda _: [f])
    a=replace(PBC_UNITS,units=(2,3,4),reaction_cycles=100.)
    ch=replace(BASELINE,executions=1000,cycle_seconds=1.) if multiple else BASELINE
    rows=[surface_units.evaluate_design(4,'shallow',HARDWARE['H3'],ch,d,f,int(n),count,u,a)
          for d in range(3,28,2) for n in nt for count in range(1,38) for u in a.units]
    best=min((r for r in rows if r['passes']),key=lambda r:r['qubits'])
    got=surface.min_qubits(4,'shallow',HARDWARE['H3'],ch,a)
    assert got['qubits'] == best['qubits']
    assert got['copies']*got['runs_per_copy'] >= ch.executions
    if multiple:
        assert got['copies'] > 1


def test_draft_bound_is_separate_and_optimistic_for_a_fixed_design():
    args=(6,'shallow',HARDWARE['H1'],BASELINE,13,FACTORIES[0],28,2)
    old=surface_pbc.evaluate_design(*args,PBC)
    bound=surface_pbc.evaluate_design(*args,PBC_DRAFT_BOUND)
    assert bound['buffer_tiles'] == bound['startup_cycles'] == 0
    assert bound['cycles'] <= old['cycles']
    assert bound['qubits'] <= old['qubits']
    counted=surface_pbc.evaluate_design(*args,PBC_COLUMNS)
    assert counted['data_tiles'] == 21


def test_clifford_elimination_preserves_full_unitary_and_final_readouts():
    layers=[[('H',(0,)),('S',(1,))],[('CNOT',(0,1))],[('T',(0,)),('T',(1,))],
            [('H',(1,))],[('CNOT',(1,0))],[('T',(0,))]]
    result=pbc_compile.compile_layers(layers,2)
    assert result['unitary_error'] < 1e-12
    assert result['max_within_layer_commutator'] < 1e-12
    assert any(sum(c != 'I' for c in label[1:]) > 1 for label in result['axes'])
    rng=np.random.default_rng(123)
    psi=rng.normal(size=4)+1j*rng.normal(size=4)
    psi/=np.linalg.norm(psi)
    expected=abs(result['original'] @ psi)**2
    state=result['rotations'] @ psi
    for b in range(4):
        projector=np.eye(4,dtype=complex)
        for q,op in enumerate(result['final_operators']):
            projector=projector @ (np.eye(4)+(-1)**((b >> (1-q))&1)*op)/2
        assert np.real(np.vdot(state,projector @ state)) == pytest.approx(expected[b],abs=1e-12)
