"""Compare counted sequential and parallel PBC resource models.

    .venv/bin/python scripts/compare_surface_architectures.py

Writes results/pbc_comparison_v2/{comparison.json,summary.csv} and
figures/pbc_comparison_v2/capability.{png,pdf}. Historical artifacts are read-only.
"""
import csv
from dataclasses import asdict, replace
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from capability import curves, direct, surface
from capability.config import BASELINE, HARDWARE, N_MAX, PBC, PBC_DRAFT_BOUND, PBC_UNITS

VARIANTS = {
    'pbc_1us': ('Sequential PBC, 1 us', PBC, BASELINE),
    'pbc_units': ('Parallel PBC units, 1 us', PBC_UNITS, BASELINE),
    'pbc_bound': ('PBC throughput bound', PBC_DRAFT_BOUND, BASELINE),
    'pbc_600ns': ('Sequential PBC, 600 ns', PBC, replace(BASELINE, cycle_seconds=600e-9)),
}
PROBES = (1e3, 1e4, 1e5, N_MAX)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def provenance():
    sources = sorted((ROOT/'capability').glob('*.py'))
    sources += [Path(__file__), ROOT/'scripts/make_curves.py', ROOT/'results/synthesis_fit.json',
                ROOT/'docs/surface_architectures.md']
    return dict(
        generated_utc=datetime.now(timezone.utc).isoformat(),
        git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        git_status=subprocess.check_output(['git', 'status', '--short'], cwd=ROOT, text=True).splitlines(),
        source_sha256={str(p.relative_to(ROOT)):sha256(p) for p in sources},
    )


def compute():
    data = dict(schema_version=2, provenance=provenance(),
                variants={key:dict(label=label, architecture=asdict(arch), choices=asdict(ch))
                          for key,(label,arch,ch) in VARIANTS.items()}, hardware={})
    for key, hw in HARDWARE.items():
        print(f'{key}: comparing surface architectures', flush=True)
        dr = direct.sweep(hw, BASELINE)
        models = {}
        for name, (_, architecture, ch) in VARIANTS.items():
            rows = surface.sweep(hw, ch, architecture)
            probes = {f'{n:.0e}':surface.best_at(n, hw, ch, rows, architecture) for n in PROBES}
            models[name] = dict(surface=rows, crossover=curves.crossover(dr, rows), probes=probes)
        data['hardware'][key] = dict(hardware=asdict(hw), direct=dr, models=models)
    data['scope'] = (
        'Counted PBC architectures only. Sequential PBC charges deterministic factory startup, '
        'batch supply and output buffers. The throughput bound is not a delivery schedule. '
        'Parallel units charge additional storage and factory interfaces; this is a finite '
        'architecture menu, not an optimum over all PBC implementations. Physical delivery, '
        'factory rejection and decoding remain uncalibrated.')
    return data


def plot(data, destination):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.6), sharex=True, sharey=True)
    colors = {'pbc_units':'#a33943', 'pbc_bound':'#b47919',
              'pbc_1us':'#2878b5', 'pbc_600ns':'#138a72'}
    for ax, (key, record) in zip(axes.flat, data['hardware'].items()):
        entries = [('Direct', record['direct'], '#243b53', '-')]
        entries += [(data['variants'][name]['label'], v['surface'], colors[name],
                     {'pbc_bound':':', 'pbc_600ns':'--'}.get(name, '-'))
                    for name,v in record['models'].items()]
        for label, rows, color, style in entries:
            y = curves.curve(rows).astype(float)
            y[y == 0] = np.nan
            ax.step(curves.N_GRID, y, where='post', label=label, color=color, ls=style, lw=1.7)
        ax.set(xscale='log', yscale='log', xlim=(1, N_MAX), ylim=(1, 4e8),
               title=f"{key}: {record['hardware']['name']}")
        ax.grid(alpha=.15)
        ax.spines[['top', 'right']].set_visible(False)
    for ax in axes[:,0]:
        ax.set_ylabel('QUOPS capability s')
    for ax in axes[1]:
        ax.set_xlabel('Physical qubits N')
    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=3, frameon=False, fontsize=9,
               bbox_to_anchor=(.5,.035))
    fig.suptitle('Surface-code architectural tradeoffs at fixed physical noise', fontsize=13)
    fig.text(.5,.012, 'PBC includes factory startup and output-buffer area/error. '
             'The throughput bound omits startup and buffers; it is not an execution schedule.',
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0,.115,1,.95))
    for ext in ('png','pdf'):
        fig.savefig(destination.with_suffix('.'+ext), dpi=180)
    plt.close(fig)


def main():
    data = compute()
    output = ROOT/'results/pbc_comparison_v2'
    figures = ROOT/'figures/pbc_comparison_v2'
    output.mkdir(exist_ok=True)
    figures.mkdir(exist_ok=True)
    (output/'comparison.json').write_text(json.dumps(data, indent=1, allow_nan=False)+'\n')
    with (output/'summary.csv').open('w', newline='') as handle:
        fields = ('hardware','variant','direct_s','crossover_qubits','s_at_1e4','s_at_1e6',
                  'w_at_1e6','distance_at_1e6','factories_at_1e6','run_seconds_at_1e6')
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for key,v in data['hardware'].items():
            for name,model in v['models'].items():
                lo, hi = model['probes']['1e+04'], model['probes']['1e+06']
                writer.writerow(dict(hardware=key,variant=name,direct_s=curves.capability(v['direct'],N_MAX),
                                     crossover_qubits=model['crossover'],s_at_1e4=lo['s'] if lo else 0,
                                     s_at_1e6=hi['s'] if hi else 0,w_at_1e6=hi['w'] if hi else None,
                                     distance_at_1e6=hi['distance'] if hi else None,
                                     factories_at_1e6=hi['factories'] if hi else None,
                                     run_seconds_at_1e6=hi['run_seconds'] if hi else None))
    plot(data, figures/'capability')
    print((output/'summary.csv').read_text())


if __name__ == '__main__':
    main()
