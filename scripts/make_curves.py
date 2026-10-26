"""Baseline QUOPS capability curves for the direct and surface-code implementations on H1-H4.

    .venv/bin/python scripts/make_curves.py
writes results/baseline.json and figures/quops_direct_surface.{pdf,png}.
"""
import json
import sys
from dataclasses import asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from capability import curves, direct, surface
from capability.config import BASELINE, HARDWARE, N_MAX

PROBES = (1e4, N_MAX)                       # exact capability is reported at these qubit counts
DIRECT, SURFACE = '#2a78d6', '#eb6834'      # categorical slots 1 and 2 of the dataviz reference palette
INK, MUTED = '#2b2b2b', '#8a8a8a'


def compute(ch=BASELINE):
    out = {}
    for key, hw in HARDWARE.items():
        d = direct.sweep(hw, ch)
        s = surface.sweep(hw, ch)
        probes = {f'{N:.0e}': surface.best_at(N, hw, ch, s) for N in PROBES}
        out[key] = dict(hardware=asdict(hw), direct=d, surface=s, crossover=curves.crossover(d, s), probes=probes)
    return out


def best(rows):
    ok = [r for r in rows if r['passes']]
    return max(ok, key=lambda r: r['s']) if ok else None


def plot(data, path):
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.6), sharex=True, sharey=True)
    grid = curves.N_GRID
    for ax, key in zip(axes.flat, ('H1', 'H2', 'H3', 'H4')):
        v = data[key]
        for rows, color, label in ((v['direct'], DIRECT, 'Direct'), (v['surface'], SURFACE, 'Surface code')):
            y = curves.curve(rows, grid).astype(float)
            y[y == 0] = float('nan')
            ax.step(grid, y, where='post', color=color, lw=2, label=label)
        if v['crossover']:
            ax.axvline(v['crossover'], color=MUTED, lw=1, ls=':')
            ax.text(v['crossover']*1.3, 4e8, f"crossover\n{v['crossover']:,.0f} qubits", color=INK, fontsize=8,
                    va='top')
        b = best(v['direct'])
        ax.text(2, b['s']*1.8, f"direct: s = {b['s']:,} ({b['ensemble']}, w = {b['w']})", color=INK,
                fontsize=8, va='bottom')
        ax.set(xscale='log', yscale='log', xlim=(1, N_MAX), ylim=(1, 1e9), title=f"{key}: {v['hardware']['name']}")
        ax.title.set_fontsize(10)
        ax.grid(alpha=.15)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
    for ax in axes[1]:
        ax.set_xlabel('Physical qubits N')
    for ax in axes[:, 0]:
        ax.set_ylabel('QUOPS capability s')
    axes[0, 0].legend(frameon=False, fontsize=9, loc='upper left', bbox_to_anchor=(0, .86))
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(path.with_suffix('.' + ext), dpi=180)
    plt.close(fig)


def main():
    data = compute()
    (ROOT/'results').mkdir(exist_ok=True)
    (ROOT/'results/baseline.json').write_text(json.dumps(dict(choices=asdict(BASELINE), **data), indent=1) + '\n')
    (ROOT/'figures').mkdir(exist_ok=True)
    plot(data, ROOT/'figures/quops_direct_surface')
    for key, v in data.items():
        b = best(v['direct'])
        print(f"{key} direct s={b['s']} ({b['ensemble']}, w={b['w']}); crossover N={v['crossover']:,.0f}; "
              + '; '.join(f"surface s at {N} = {x['s']:,} (w={x['w']}, d={x['distance']}, {x['factories']} x "
                          f"{x['factory']})" for N, x in v['probes'].items()))


if __name__ == '__main__':
    main()
