"""Combined QUOPS capability figure of the draft: direct, surface code, STAR, and Pinnacle on the four hardware classes.

    .venv/bin/python scripts/make_quops_comparison.py
reads the cached results/baseline.json (direct and sequential-PBC surface code), results/quops_comparison_star.json (STAR on the
2D-grid classes H1 and H3, the steps of the draft's figure; see that file for its source), and results/pinnacle.json (Pinnacle on the all-to-all classes H2 and H4), and writes
figures/quops_comparison.{pdf,png}. Nothing is recomputed. The layout, fonts, colors, and line styles reproduce the
comparison figure of the draft (8.1 x 5.0 in, one-row legend above the panels); the Pinnacle curve is the best over the
supply menu of Sec. 3.7 and replaces the single-engine curve that the draft showed.
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'pdf.fonttype': 42,                       # TrueType fonts, as in the draft's figure
                     'axes.linewidth': 0.7, 'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
                     'grid.linewidth': 0.7, 'grid.alpha': 0.18})

DIRECT, SURFACE, STAR, PINNACLE = '#2878b5', '#c95424', '#8143a1', '#087858'
KEYS = ('H1', 'H2', 'H3', 'H4')
TITLES = {'H1': 'H1  2D grid, uniform', 'H2': 'H2  all-to-all, uniform', 'H3': 'H3  2D grid, bimodal', 'H4': 'H4  all-to-all, bimodal'}


def step(rows, n_max=1e6):
    """Vertices of the capability staircase: the best passing size at each distinct qubit count, one vertex pair per
    increase, ending at N = 1e6."""
    best_at = {}
    for r in rows:
        if r['passes']:
            best_at[r['qubits']] = max(best_at.get(r['qubits'], 0), r['s'])
    xs, ys, best = [], [], 0
    for q in sorted(best_at):
        if best_at[q] > best:
            if best:
                xs.append(q); ys.append(best)
            xs.append(q); ys.append(best_at[q]); best = best_at[q]
    xs.append(n_max); ys.append(best)
    return np.array(xs, float), np.array(ys, float)


def with_probes(rows, probes):
    """Sweep rows plus the exact fixed-budget designs, which the plotted curves include."""
    return list(rows) + [p for p in probes if p]


def main(out=ROOT/'figures/quops_comparison'):
    base = json.loads((ROOT/'results/baseline.json').read_text())
    star = json.loads((ROOT/'results/quops_comparison_star.json').read_text())
    pin = json.loads((ROOT/'results/pinnacle.json').read_text())
    fig, axes = plt.subplots(2, 2, figsize=(8.1, 5.0), sharex=True, sharey=True)
    for ax, key in zip(axes.flat, KEYS):
        ax.plot(*step(with_probes(base[key]['direct'], [])), color=DIRECT, lw=1.9, label='Direct')
        ax.plot(*step(with_probes(base[key]['surface'], base[key]['probes'].values())), color=SURFACE, lw=1.9, label='Surface code (sequential PBC)')
        if key in star:
            ax.plot(*step([dict(qubits=q, s=v, passes=True) for q, v in star[key]]), color=STAR, lw=1.9, ls='--', label='STAR')
        if pin[key]['supported']:
            ax.plot(*step(pin[key]['pinnacle']), color=PINNACLE, lw=1.9, ls='-.', label='Pinnacle')
        ax.set(xscale='log', yscale='log', xlim=(1, 10**6.0172), ylim=(3.0, 1.6e9))
        ax.set_title(TITLES[key], loc='left', fontsize=10, pad=4)
        ax.set_xticks([1, 1e2, 1e4, 1e6])
        ax.grid(True, which='major')
        ax.minorticks_off()
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
        ax.tick_params(labelsize=9)
    fig.supxlabel('Physical qubits $N$', fontsize=10, y=0.015)
    fig.supylabel('QUOPS capability $s$', fontsize=10, x=0.025)
    handles, labels = [], []
    for ax in axes.flat:
        for h, l in zip(*ax.get_legend_handles_labels()):
            if l not in labels:
                handles.append(h); labels.append(l)
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.995), ncol=4, frameon=False, fontsize=9, handlelength=1.5)
    fig.subplots_adjust(left=0.115, right=0.980, top=0.870, bottom=0.115, wspace=0.130, hspace=0.270)
    for ext in ('pdf', 'png'):
        fig.savefig(Path(out).with_suffix('.' + ext), dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main()
