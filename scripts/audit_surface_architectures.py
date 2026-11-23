"""Tile, circuit and resource audit; writes only architecture_audit_v3 outputs."""
import csv
from copy import deepcopy
from dataclasses import asdict, replace
from datetime import datetime, timezone
import hashlib
import json
from math import ceil, sqrt
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

from capability import curves, direct, layouts, pbc_compile, star, surface, surface_pbc, surface_units
from capability.config import BASELINE,HARDWARE,PBC,PBC_COLUMNS,PBC_DRAFT_BOUND,PBC_UNITS,STAR_BASELINE,STAR_SINGLE_QUBIT
from capability.factories import FACTORIES

OUT=ROOT/'results/architecture_audit_v3'
FIG=ROOT/'figures/architecture_audit_v3'
VARIANTS={
    'draft_bound':('PBC throughput bound',PBC_DRAFT_BOUND,BASELINE,'optimistic bound'),
    'sequential':('Sequential PBC',PBC_COLUMNS,BASELINE,'finite-buffer resource model'),
    'parallel_units':('PBC units (2–8)',PBC_UNITS,BASELINE,'published protocol; modeled factory interface'),
    'sequential_600ns':('Sequential PBC, 600 ns',PBC_COLUMNS,replace(BASELINE,cycle_seconds=600e-9),'clock sensitivity'),
    'parallel_units_600ns':('PBC units, 600 ns',PBC_UNITS,replace(BASELINE,cycle_seconds=600e-9),'clock sensitivity'),
}


def provenance():
    paths=list((ROOT/'capability').glob('*.py'))+[Path(__file__)]
    paths+=list((ROOT/'tests').glob('test*py'))
    return dict(generated_utc=datetime.now(timezone.utc).isoformat(),
                git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                git_status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True).splitlines(),
                source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
                source='https://arxiv.org/html/1808.02892v3',
                source_sections=['2.3 / Figure 13b','5.1–5.2 / Figures 25–30'])


def layout_audit():
    rows=[]
    for w in range(2,1002,2):
        row=layouts.fast_layout(w).certificate()
        row['draft_tiles']=(2*w+ceil(sqrt(8*w))+1)
        row['difference']=row['tiles']-row['draft_tiles']
        assert row['disjoint'] and row['connected_ancilla'] and row['both_pauli_edges_accessible']
        rows.append(row)
    fig,axes=plt.subplots(1,4,figsize=(12,3.4))
    for ax,w in zip(axes,(6,8,10,18)):
        layout=layouts.fast_layout(w)
        for x,y in layout.ancilla:
            ax.add_patch(Rectangle((x,y),1,1,facecolor='#e9dcb8',edgecolor='#aaa',lw=.5))
        for x,y in layout.patches:
            ax.add_patch(Rectangle((x,y),1,2,facecolor='#567ba3',edgecolor='#243b53',lw=1.5))
        ax.autoscale_view()
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title(f'w={w}: {layout.area} tiles\nformula: {(2*w+ceil(sqrt(8*w))+1)}')
    fig.suptitle('Two-qubit patches (blue) and connected ancilla tiles (tan)')
    fig.tight_layout()
    fig.savefig(FIG/'layouts.png',dpi=180)
    fig.savefig(FIG/'layouts.pdf')
    plt.close(fig)
    return dict(layouts=rows,underestimated_widths=[r['width'] for r in rows if r['difference'] > 0],
                scope='Certificate of tile packing and boundary access; no physical syndrome simulation',
                ladder_counterexample=dict(pairing=[[0,2],[1,3]],required_route_rounds=2,
                                            scope='one-row routing realization only'))


def circuit_audit():
    rows=[]
    for w,ensemble in ((2,'shallow'),(2,'deep'),(4,'shallow'),(4,'deep'),(6,'shallow'),(6,'deep')):
        seed=20261109+w+(100 if ensemble == 'deep' else 0)
        print('Circuit audit',w,ensemble,flush=True)
        layers=pbc_compile.synthesized_quops_layers(w,ensemble,seed)
        result=pbc_compile.compile_layers(layers,w)
        assert result['unitary_error'] < 1e-10 and result['max_within_layer_commutator'] < 1e-10
        row={k:v for k,v in result.items() if k not in ('original','rotations','final_operators')}
        row.update(w=w,ensemble=ensemble,seed=seed,synthesis_tolerance=.05,
                   circuit_layers=[[[name,list(q)] for name,q in layer] for layer in layers],
                   max_pauli_weight=max(sum(c != 'I' for c in s[1:]) for s in result['axes']))
        rows.append(row)
    return rows


def resource_audit():
    results={}
    # Surface models depend on noise but not connectivity; copy only after that
    # explicit restriction is stated in each record. Direct results stay separate.
    for hwkey in ('H1','H3'):
        hw=HARDWARE[hwkey]
        record=dict(hardware=asdict(hw),direct=direct.sweep(hw,BASELINE),models={})
        for key,(label,a,ch,status) in VARIANTS.items():
            print(hwkey,key,flush=True)
            rows=surface.sweep(hw,ch,a)
            probes={str(int(n)):surface.best_at(n,hw,ch,rows,a) for n in (1e4,1e5,1e6)}
            record['models'][key]=dict(label=label,architecture=asdict(a),choices=asdict(ch),
                                      status=status,surface=rows,probes=probes,
                                      crossover=curves.crossover(record['direct'],rows))
        record['costed_pbc_envelope']=record['models']['sequential']['surface']+record['models']['parallel_units']['surface']
        record['multi_execution']={}
        for key,a in (('sequential',PBC_COLUMNS),('parallel_units',PBC_UNITS)):
            ch=replace(BASELINE,executions=10000)
            print(hwkey,key,'10000 executions',flush=True)
            rows=surface.sweep(hw,ch,a)
            record['multi_execution'][key]=dict(choices=asdict(ch),best=surface.best_at(1e6,hw,ch,rows,a))
        # Expand the unit menu as a bounded convergence check, without assuming
        # that 8 units are globally optimal.
        a=replace(PBC_UNITS,units=tuple(range(2,17)))
        print(hwkey,'expanded unit menu',flush=True)
        rows=surface.sweep(hw,BASELINE,a)
        record['expanded_unit_menu']=dict(architecture=asdict(a),best=surface.best_at(1e6,hw,BASELINE,rows,a))
        results[hwkey]=record
    for key,source in (('H2','H1'),('H4','H3')):
        record=deepcopy(results[source])
        record['hardware']=asdict(HARDWARE[key])
        record['direct']=direct.sweep(HARDWARE[key],BASELINE)
        record['surface_connectivity_scope']='Same grid-compatible resource model as '+source
        for model in record['models'].values():
            model['crossover']=curves.crossover(record['direct'],model['surface'])
        results[key]=record
    return results


def fixed_design_tradeoff():
    f=next(f for f in FACTORIES if f.qubits == 30700)
    rows=[]
    for count in (1,2,4,8):
        r=surface_pbc.evaluate_design(32,'deep',HARDWARE['H1'],BASELINE,21,f,40,count,PBC_COLUMNS)
        rows.append(dict(model='sequential',**r))
        for u in (2,3,4,6,8,12,16):
            r=surface_units.evaluate_design(32,'deep',HARDWARE['H1'],BASELINE,21,f,40,count,u)
            rows.append(dict(model='parallel_units',**r))
    return dict(w=32,ensemble='deep',distance=21,T_per_rotation=40,rows=rows)


def star_comparison(hardware):
    """Keep STAR's layered architecture explicit, with matched clocks/task."""
    output={}
    for key,record in hardware.items():
        if key not in ('H1','H3'):
            continue
        output[key]={}
        for clock in (1e-6,600e-9):
            ch=replace(BASELINE,cycle_seconds=clock)
            entry={}
            ref=record['models']['sequential' if clock == 1e-6 else 'sequential_600ns']['surface']
            for name,sc in (('m1',STAR_SINGLE_QUBIT),('m2',STAR_BASELINE)):
                print('STAR comparison',key,name,clock,flush=True)
                rows=star.sweep(HARDWARE[key],ch,sc)
                intervals=curves.lead_intervals(rows,ref+record['direct'])
                # The audit historically clips the upper endpoint to its largest probe.
                entry[name]=dict(architecture='gate-based STAR',star_choices=asdict(sc),rows=rows,
                                 best=star.best_at(1e6,rows),above_direct_and_sequential_pbc=[
                                     [lo,min(hi,10**6)] for lo,hi in intervals])
            output[key][str(clock)]=entry
    return output


def plot_results(data):
    fig,axes=plt.subplots(2,2,figsize=(10,8),sharex=True,sharey=True)
    styles={'draft_bound':('#b47919',':'),
            'sequential':('#2678b5','-'),'parallel_units':('#a33943','-')}
    for ax,key in zip(axes.flat,('H1','H2','H3','H4')):
        record=data['hardware'][key]
        direct_y=curves.curve(record['direct']).astype(float)
        direct_y[direct_y == 0]=np.nan
        ax.step(curves.N_GRID,direct_y,where='post',color='#388342',label='Direct')
        for name,(color,style) in styles.items():
            m=record['models'][name]
            y=curves.curve(m['surface']).astype(float)
            y[y == 0]=np.nan
            ax.step(curves.N_GRID,y,where='post',color=color,ls=style,label=m['label'])
        ax.set(xscale='log',yscale='log',xlim=(1e3,1e6),ylim=(1,3e8),
               xlabel='Physical qubits',title=key+': '+HARDWARE[key].name)
        ax.grid(alpha=.15)
    for ax in axes[:,0]:
        ax.set_ylabel('QUOPS capability')
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',ncol=3,frameon=False,fontsize=9)
    fig.tight_layout(rect=(0,.09,1,1))
    for ext in ('png','pdf'):
        fig.savefig(FIG/('comparison.'+ext),dpi=180)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,4))
    for model,color in (('sequential','#2678b5'),('parallel_units','#a33943')):
        rows=[r for r in data['fixed_design_tradeoff']['rows'] if r['model'] == model and r['passes']]
        ax.scatter([r['qubits'] for r in rows],[r['run_seconds'] for r in rows],label=model,color=color)
    ax.set(xscale='log',yscale='log',xlabel='Physical qubits',ylabel='Seconds per run',
           title='Fixed w=32, d=21, 40 T gates per rotation')
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG/'latency_tradeoff.png',dpi=180)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    FIG.mkdir(exist_ok=True)
    data=dict(provenance=provenance(),layout_audit=layout_audit(),circuit_audit=circuit_audit())
    (OUT/'layout_and_circuits.json').write_text(json.dumps(data,indent=1,allow_nan=False)+'\n')
    data.update(hardware=resource_audit(),fixed_design_tradeoff=fixed_design_tradeoff())
    data['star_comparison']=star_comparison(data['hardware'])
    (OUT/'audit.json').write_text(json.dumps(data,indent=1,allow_nan=False)+'\n')
    with (OUT/'summary.csv').open('w',newline='') as f:
        writer=csv.writer(f)
        writer.writerow(['hardware','variant','crossover','s_at_1e6','width','distance','units','factories','seconds'])
        for key,record in data['hardware'].items():
            for name,model in record['models'].items():
                r=model['probes']['1000000']
                writer.writerow([key,name,model['crossover']]+[r.get(k,1) if r else None
                    for k in ('s','w','distance','units','factories','run_seconds')])
    plot_results(data)
    print((OUT/'summary.csv').read_text())


if __name__ == '__main__':
    main()
