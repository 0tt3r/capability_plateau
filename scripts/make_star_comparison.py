"""Compare single-qubit and two-qubit STAR preparation with QUOPS baselines.

    .venv/bin/python scripts/make_star_comparison.py

Writes results/star_comparison.json and figures/star_comparison.{png,pdf},
plus a focused single-qubit comparison and a preparation diagnostic.
"""
import argparse
from dataclasses import asdict, replace
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

from capability import curves, direct, star, surface
from capability.config import BASELINE, HARDWARE, N_MAX, STAR_BASELINE, STAR_SINGLE_QUBIT
from capability.star_states import state_profile

COLORS = {'direct': '#2878b5', 'surface': '#d75c28', 'star': '#7b3fa0',
          'star_single_qubit': '#138a72'}
PROBES = (1e3, 1e4, 1e5, N_MAX)


def compute(angle_points=512, sensitivity=True):
    sc = replace(STAR_BASELINE, angle_points=angle_points)
    single = replace(STAR_SINGLE_QUBIT, angle_points=angle_points)
    try:
        base_commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        base_commit = None
    out = dict(model='STAR m=1 (Choi) and m=2 (Toshio); analytical noise assumptions; no observable PEC',
               base_commit=base_commit, choices=asdict(BASELINE), star_choices=asdict(sc),
               star_single_qubit_choices=asdict(single),
               source='https://doi.org/10.1103/PhysRevX.15.021057', hardware={})
    for key in ('H1', 'H3'):
        print(f'{key}: computing direct, Clifford+T and STAR m=1 / m=2...', flush=True)
        hw = HARDWARE[key]
        dr = direct.sweep(hw, BASELINE)
        sr = surface.sweep(hw, BASELINE)
        tr = star.sweep(hw, BASELINE, sc)
        single_rows = star.sweep(hw, BASELINE, single)
        probes = {}
        for N in PROBES:
            probes[f'{N:.0e}'] = dict(star=star.best_at(N, tr),
                                     star_single_qubit=star.best_at(N, single_rows),
                                     surface=surface.best_at(N, hw, BASELINE, sr))
        variation = {}
        if sensitivity:
            for factor in (.5, 2.):
                print(f'{key}: rejection sensitivity {factor:g}x...', flush=True)
                variation[str(factor)] = star.sweep(hw, BASELINE, replace(sc, rejection_scale=factor))
        # Quantify a second, separate uncertainty without treating it as a confidence band.
        floor_sensitivity = {}
        if sensitivity:
            for factor in (0., 2.):
                alt = replace(sc, preparation_floor_scale=factor)
                rows = star.sweep(hw, BASELINE, alt)
                floor_sensitivity[str(factor)] = dict(
                    resource_ceiling=star.resource_ceiling(hw, BASELINE, alt),
                    crossover=curves.crossover(dr, rows), best=star.best_at(N_MAX, rows))
        single_noise, single_rejection, single_floor = {}, {}, {}
        if sensitivity:
            for factor in (.5, 2.):
                print(f'{key}: m=1 undetected-error sensitivity {factor:g}x...', flush=True)
                single_noise[str(factor)] = star.sweep(
                    hw, BASELINE, replace(single, single_qubit_undetected_scale=factor))
                single_rejection[str(factor)] = star.sweep(
                    hw, BASELINE, replace(single, rejection_scale=factor))
            for factor in (0., 2.):
                alt = replace(single, preparation_floor_scale=factor)
                rows = star.sweep(hw, BASELINE, alt)
                single_floor[str(factor)] = dict(
                    resource_ceiling=star.resource_ceiling(hw, BASELINE, alt),
                    crossover=curves.crossover(dr, rows), best=star.best_at(N_MAX, rows))
        # Isolate postselection policy and original-protocol switching from rotation weight.
        single_strip = star.sweep(hw, BASELINE, replace(single, postselection='strip'))
        single_no_switch = star.sweep(hw, BASELINE, replace(single, protocol_switching=False))
        out['hardware'][key] = dict(
            hardware=asdict(hw), direct=dr, surface=sr, star=tr,
            star_single_qubit=single_rows, probes=probes,
            star_crossover=curves.crossover(dr, tr), surface_crossover=curves.crossover(dr, sr),
            star_single_qubit_crossover=curves.crossover(dr, single_rows),
            star_resource_ceiling=star.resource_ceiling(hw, BASELINE, sc),
            star_single_qubit_resource_ceiling=star.resource_ceiling(hw, BASELINE, single),
            rejection_sensitivity=variation,
            preparation_floor_sensitivity=floor_sensitivity,
            single_qubit_noise_sensitivity=single_noise,
            single_qubit_rejection_sensitivity=single_rejection,
            single_qubit_preparation_floor_sensitivity=single_floor,
            single_qubit_restricted_postselection=single_strip,
            single_qubit_no_switching=single_no_switch,
        )
    return out


def plot(data):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4), sharex=True, sharey=True)
    for ax, (key, v) in zip(axes, data['hardware'].items()):
        for model, label in (('direct', 'Direct'), ('surface', 'Surface code + Clifford+T'),
                             ('star', 'STAR: two-qubit rotations ($m=2$)'),
                             ('star_single_qubit', 'STAR: single-qubit rotations ($m=1$)')):
            y = curves.curve(v[model]).astype(float)
            y[y == 0] = np.nan
            ax.step(curves.N_GRID, y, where='post', color=COLORS[model], lw=2, label=label)
        for rows in v['rejection_sensitivity'].values():
            y = curves.curve(rows).astype(float)
            y[y == 0] = np.nan
            ax.step(curves.N_GRID, y, where='post', color=COLORS['star'], alpha=.45, lw=1, ls=':')
        for rows in v['single_qubit_noise_sensitivity'].values():
            y = curves.curve(rows).astype(float)
            y[y == 0] = np.nan
            ax.step(curves.N_GRID, y, where='post', color=COLORS['star_single_qubit'], alpha=.45, lw=1, ls=':')
        for model, m, offset in (('star', 2, 2.2), ('star_single_qubit', 1, .3)):
            top = star.best_at(N_MAX, v[model])
            if top:
                ax.annotate(f"m = {m}: s = {top['s']:,}", (8e5, top['s']*offset),
                            xytext=(0, 0), textcoords='offset points', ha='right',
                            color=COLORS[model], fontsize=9)
            if v[f'{model}_crossover']:
                ax.axvline(v[f'{model}_crossover'], color=COLORS[model], lw=.8, ls='--', alpha=.6)
        ax.set(xscale='log', yscale='log', xlim=(1, N_MAX), ylim=(1, 4e8),
               xlabel='Physical qubits N', title=f"{key}: {v['hardware']['name']}")
        ax.grid(which='major', alpha=.15)
    axes[0].set_ylabel('QUOPS capability s')
    axes[0].legend(loc='upper left', frameon=False, fontsize=9)
    fig.suptitle('Uniform-angle QUOPS: single-qubit and two-qubit resource-state preparation', fontsize=12)
    fig.text(.5, .015,
             'Both STAR curves include probabilistic coherent error cancellation and original-STAR protocol switching.\n'
             'Dotted: 0.5x / 2x assumptions (m=1 undetected errors; m=2 rejection counts), not confidence intervals.\n'
             'm=1: full postselection; m=2: restricted postselection. Dashed: STAR/direct crossovers. No observable-level PEC.',
             ha='center', fontsize=8, color='#555555')
    fig.tight_layout(rect=(0, .11, 1, .94))
    for extension in ('png', 'pdf'):
        fig.savefig(ROOT/f'figures/star_comparison.{extension}', dpi=200)
    plt.close(fig)


def plot_preparation(data):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.9))
    phi = np.geomspace(1e-5, np.pi/8, 500)
    sc = replace(STAR_BASELINE, angle_points=data['star_choices']['angle_points'])
    single = replace(STAR_SINGLE_QUBIT, angle_points=sc.angle_points)
    for key, color in (('H1', '#2878b5'), ('H3', '#7b3fa0')):
        for settings, style in ((sc, '-'), (single, '--')):
            p = state_profile(phi, 3, 9, HARDWARE[key], BASELINE, settings)
            label = f'{key}, m={settings.rotation_weight}, d_p={p.preparation_distance}'
            axes[0].loglog(phi, p.error, color=color, ls=style, label=label)
            axes[1].semilogx(phi, p.acceptance, color=color, ls=style, label=label)
        axes[0].axhline(p.fallback_error, color=color, lw=.7, ls=':')
    axes[0].set(xlabel=r'Residual paper angle $|\phi|$ (rad)',
                ylabel='Rotation error per teleportation', title='After error cancellation and protocol switching')
    axes[1].set(xlabel=r'Residual paper angle $|\phi|$ (rad)',
                ylabel='Preparation success probability', title='Assumed finite-noise preparation model', ylim=(0, 1))
    for ax in axes:
        ax.grid(alpha=.15)
        ax.legend(frameon=False)
    fig.suptitle('Preparation comparison: k = 3, data distance 9 (preparation distances differ)', fontsize=11)
    fig.text(.5, .01, 'Each curve includes original-STAR switching. Logical errors during preparation are accounted for separately.',
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .04, 1, .96))
    for extension in ('png', 'pdf'):
        fig.savefig(ROOT/f'figures/star_preparation.{extension}', dpi=200)
    plt.close(fig)


def plot_single_qubit_comparison(data):
    """Show how much the m=1 capability depends on original-STAR switching."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
    entries = (
        ('direct', 'Direct', COLORS['direct'], '-'),
        ('star', 'STAR: two-qubit rotations ($m=2$)', COLORS['star'], '-'),
        ('star_single_qubit', 'STAR: single-qubit rotations ($m=1$)', COLORS['star_single_qubit'], '-'),
        ('single_qubit_no_switching', '$m=1$, protocol switching disabled', COLORS['star_single_qubit'], '--'),
    )
    for ax, (key, v) in zip(axes, data['hardware'].items()):
        for model, label, color, style in entries:
            y = curves.curve(v[model]).astype(float)
            y[y == 0] = np.nan
            ax.step(curves.N_GRID, y, where='post', color=color, ls=style, lw=2, label=label)
        ax.set(xscale='log', yscale='log', xlim=(100, N_MAX), ylim=(4, 1e4),
               xlabel='Physical qubits N', title=f"{key}: {v['hardware']['name']}")
        ax.grid(which='major', alpha=.15)
    axes[0].set_ylabel('QUOPS capability s')
    axes[0].legend(loc='upper left', frameon=False, fontsize=8)
    fig.suptitle('Single-qubit transversal rotations: effect of protocol switching', fontsize=12)
    fig.text(.5, .015,
             'Dashed green disables switching for the requested rotations. Additional error-cancellation rotations still use original STAR.\n'
             'Both variants use the same RUS schedule and deadline criterion; m=1 undetected-error coefficients are modeling assumptions.',
             ha='center', fontsize=8, color='#555555')
    fig.tight_layout(rect=(0, .09, 1, .94))
    for extension in ('png', 'pdf'):
        fig.savefig(ROOT/f'figures/star_single_qubit_comparison.{extension}', dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--angle-points', type=int, default=STAR_BASELINE.angle_points)
    parser.add_argument('--no-sensitivity', action='store_true')
    args = parser.parse_args()
    data = compute(args.angle_points, not args.no_sensitivity)
    (ROOT/'results').mkdir(exist_ok=True)
    (ROOT/'figures').mkdir(exist_ok=True)
    (ROOT/'results/star_comparison.json').write_text(json.dumps(data, indent=1, allow_nan=False) + '\n')
    plot(data)
    plot_preparation(data)
    plot_single_qubit_comparison(data)
    for key, v in data['hardware'].items():
        for model in ('star', 'star_single_qubit'):
            top = star.best_at(N_MAX, v[model])
            print(f"{key} {model}: crossover {v[f'{model}_crossover']}; max {top['s'] if top else 0}; "
                  f"resource-state size bound {v[f'{model}_resource_ceiling']:.1f}", flush=True)


if __name__ == '__main__':
    main()
