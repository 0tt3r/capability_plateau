"""Serial Pinnacle QUOPS curves alongside the repository's cached baseline.

    python scripts/make_pinnacle_curves.py
    python scripts/make_pinnacle_curves.py --skip-sensitivity

Writes results/pinnacle.json and figures/quops_direct_surface_pinnacle.{pdf,png}.
The optional sensitivity calculation writes results/pinnacle_sensitivity.json and
figures/pinnacle_sensitivity.{pdf,png}. No legacy model is imported or rerun.
"""
import argparse
import hashlib
import json
import sys
from dataclasses import asdict, replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from capability import curves, pinnacle
from capability.config import BASELINE, HARDWARE, N_MAX

PROBES = (1e4, N_MAX)
KEYS = ('H1', 'H2', 'H3', 'H4')
SUPPORTED_KEYS = ('H2', 'H4')
DIRECT, SURFACE, PINNACLE = '#2a78d6', '#eb6834', '#238b59'
INK, MUTED = '#2b2b2b', '#8a8a8a'
PAPER = 'https://arxiv.org/html/2602.11457v2'
HASH_FORMAT = 'canonical-json-sha256-v1'
MODELING_NOTES = [
    'Serial Pauli-based computation: one processing unit and one magic engine per copy; no memory.',
    'Logical measurement and magic-engine inputs are from arXiv:2602.11457v2, Sections IV-V.',
    'H1 and H3 are not modeled: a 2D nearest-neighbor realization needs placement, routing, and error costs; this is not an impossibility claim.',
    'H4 uses the same p=1e-3 depolarizing-noise fit as H2, without a claimed bimodal-noise gain.',
    'Data survival averages geometric retries; week throughput uses expected runtime, not a deadline guarantee.',
    'The baseline counts both Pauli error sectors using logical_error_scale=2; this is a declared choice.',
    'Synthesis infidelity is charged separately from stochastic no-fault loss within the shared loss budget.',
    'Direct and surface-code curves are cached baseline references and are never recomputed by this script.',
]


def canonical_json_sha256(value):
    """Hash parsed JSON independently of line endings, whitespace, and key order."""
    canonical = json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True, allow_nan=False)
    return hashlib.sha256(canonical.encode('ascii')).hexdigest()


def load_references(ch=BASELINE):
    """Validate baseline choices and hash its parsed, canonically serialized JSON."""
    path = ROOT/'results/baseline.json'
    raw = path.read_bytes()
    data = json.loads(raw)
    choices = json.loads(json.dumps(asdict(ch)))       # JSON normalizes tuples to lists
    if data.get('choices') != choices:
        raise ValueError('Cached results/baseline.json choices do not match the requested shared choices.')
    for key in KEYS:
        if data[key]['hardware'] != asdict(HARDWARE[key]):
            raise ValueError(f'Cached hardware parameters do not match config.py for {key}.')
    return data, canonical_json_sha256(data)


def with_probes(rows, probes):
    """Include exact fixed-budget designs in the plotted staircase."""
    combined = {(r['ensemble'], r['w']): r for r in rows}
    for row in probes.values():
        if row is not None:
            combined[row['ensemble'], row['w']] = row
    return sorted(combined.values(), key=lambda r: (r['ensemble'], r['w']))


def compute_pinnacle(hw, ch=BASELINE, pc=pinnacle.BASELINE):
    reason = pinnacle.support_reason(hw, pc)
    if reason is not None:
        return dict(supported=False, unsupported_reason=reason, pinnacle=[],
                    probes={f'{N:.0e}': None for N in PROBES})
    rows = pinnacle.sweep(hw, ch, pc)
    probes = {f'{N:.0e}': pinnacle.best_at(N, hw, ch, rows, pc) for N in PROBES}
    return dict(supported=True, unsupported_reason=None, pinnacle=with_probes(rows, probes), probes=probes)


def compute(ch=BASELINE, pc=pinnacle.BASELINE):
    reference, digest = load_references(ch)
    out = dict(choices=asdict(ch), pinnacle_choices=asdict(pc),
               source_paper=dict(url=PAPER, version='v2', revised='2026-05-05'),
               provenance_hash_format=HASH_FORMAT,
               baseline_sha256=digest,
               source_of_reference='cached results/baseline.json; legacy models were not rerun',
               modeling_notes=MODELING_NOTES)
    for key in KEYS:
        record = compute_pinnacle(HARDWARE[key], ch, pc)
        record.update(hardware=asdict(HARDWARE[key]), direct=reference[key]['direct'],
                      surface=reference[key]['surface'], reference_probes=reference[key]['probes'])
        record['crossover_vs_direct'] = (curves.crossover(record['direct'], record['pinnacle'])
                                         if record['supported'] else None)
        out[key] = record
    return out


def draw_curve(ax, rows, color, label, **kwargs):
    y = curves.curve(rows).astype(float)
    y[y == 0] = np.nan
    return ax.step(curves.N_GRID, y, where='post', color=color, label=label, **kwargs)


def decorate(ax):
    ax.set(xscale='log', yscale='log', xlim=(1, N_MAX), ylim=(1, 1e9))
    ax.grid(alpha=.15)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=8)


def save_plot(fig, path):
    path.parent.mkdir(exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(path.with_suffix('.' + ext), dpi=180)
    plt.close(fig)


def plot(data, path):
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.6), sharex=True, sharey=True)
    for ax, key in zip(axes.flat, KEYS):
        value = data[key]
        draw_curve(ax, value['direct'], DIRECT, 'Direct', lw=1.8)
        draw_curve(ax, value['surface'], SURFACE, 'Surface code', lw=1.8)
        if value['supported']:
            draw_curve(ax, value['pinnacle'], PINNACLE, 'Pinnacle (serial)', lw=2.1)
            crossing = value['crossover_vs_direct']
            if crossing is not None:
                ax.axvline(crossing, color=PINNACLE, lw=.9, ls=':', alpha=.7)
                ax.text(.045, .89, f'Pinnacle exceeds direct\nat {crossing:,.0f} qubits',
                        transform=ax.transAxes, fontsize=7.4, color=PINNACLE, va='top')
            top = value['probes'][f'{N_MAX:.0e}']
            if top is not None:
                ax.axvspan(top['qubits'], N_MAX, color=PINNACLE, alpha=.045)
                ax.plot([top['qubits'], N_MAX], [top['s'], top['s']], color=PINNACLE, lw=2, ls=':')
                ax.text(.045, .72, f"Pinnacle at $10^6$:\ns = {top['s']:,}", transform=ax.transAxes,
                        fontsize=7.4, va='top', color=PINNACLE)
        else:
            ax.text(.045, .84, '2D grid mapping not modeled\nPlacement and routing costs needed',
                    transform=ax.transAxes, fontsize=8, color=PINNACLE, va='top',
                    bbox=dict(facecolor='white', edgecolor='none', alpha=.9, pad=2))
        direct_best = max((r['s'] for r in value['direct'] if r['passes']), default=0)
        if direct_best:
            ax.text(2, direct_best*1.7, f'direct: s = {direct_best:,}', color=DIRECT, fontsize=7.5)
        ax.set_title(f"{key}: {value['hardware']['name']}", fontsize=10)
        decorate(ax)
    for ax in axes[1]:
        ax.set_xlabel('Physical qubits N', fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel('QUOPS capability s', fontsize=9)
    handles = [Line2D([0], [0], color=color, lw=2, label=label) for color, label in
               ((DIRECT, 'Direct (cached)'), (SURFACE, 'Surface code (cached)'), (PINNACLE, 'Pinnacle (serial)'))]
    fig.suptitle('QUOPS capability: direct, surface code and Pinnacle', fontsize=12, y=.98)
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .945), ncol=3, frameon=False, fontsize=9)
    fig.text(.5, .016, 'One-week budget; Pinnacle uses mean retry time. H4 has no assumed bimodal gain.',
             ha='center', color=MUTED, fontsize=8)
    fig.subplots_adjust(left=.09, right=.985, bottom=.11, top=.84, wspace=.12, hspace=.26)
    save_plot(fig, path)


def sensitivity_variants(pc=pinnacle.BASELINE):
    low_error = tuple(e.name for e in pinnacle.ENGINES if e.error == 1e-11)
    small_engine = tuple(e.name for e in pinnacle.ENGINES if e.error == 1e-9)
    return [
        ('baseline', {}, {}, PINNACLE, '-'),
        ('one Pauli sector', {}, dict(logical_error_scale=1), '#58a65c', '--'),
        ('1 ms cycle', dict(cycle_seconds=1e-3), {}, '#9467bd', '-'),
        ('10,000 executions', dict(executions=10_000), {}, '#c44e52', '-'),
        ('1e-9 engine only', {}, dict(engines=small_engine), '#8c6d31', '--'),
        ('1e-11 engine only', {}, dict(engines=low_error), '#4c78a8', '--'),
        ('measured synthesis charge', dict(synthesis_scale=.14), {}, '#555555', ':'),
    ]


def compute_sensitivity(baseline=None, ch=BASELINE, pc=pinnacle.BASELINE):
    out = dict(source_paper=dict(url=PAPER, version='v2'), choices=asdict(ch), pinnacle_choices=asdict(pc),
               legacy_reference_used=False, modeling_notes=MODELING_NOTES[:-1], variants=[])
    for name, shared, specific, color, style in sensitivity_variants(pc):
        changed_ch, changed_pc = replace(ch, **shared), replace(pc, **specific)
        result = {}
        for key in SUPPORTED_KEYS:
            if name == 'baseline' and baseline is not None:
                value = baseline[key]
                record = {k: value[k] for k in ('supported', 'unsupported_reason', 'pinnacle', 'probes')}
            else:
                record = compute_pinnacle(HARDWARE[key], changed_ch, changed_pc)
            result[key] = dict(hardware=asdict(HARDWARE[key]), **record)
        out['variants'].append(dict(name=name, shared_changes=shared, pinnacle_changes=specific,
                                    choices=asdict(changed_ch), pinnacle_choices=asdict(changed_pc),
                                    color=color, linestyle=style, result=result))
        print('sensitivity:', name, {key: {N: row['s'] if row else 0 for N, row in value['probes'].items()}
                                     for key, value in result.items()}, flush=True)
    return out


def plot_sensitivity(data, path):
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.8), sharex=True, sharey=True)
    handles = []
    for variant in data['variants']:
        for ax, key in zip(axes, SUPPORTED_KEYS):
            draw_curve(ax, variant['result'][key]['pinnacle'], variant['color'], variant['name'],
                       lw=2 if variant['name'] == 'baseline' else 1.5, ls=variant['linestyle'])
        handles.append(Line2D([0], [0], color=variant['color'], lw=2, ls=variant['linestyle'], label=variant['name']))
    for ax, key in zip(axes, SUPPORTED_KEYS):
        decorate(ax)
        ax.set_xlim(1e3, N_MAX)
        ax.set_title(f'{key}: {HARDWARE[key].name}', fontsize=10)
        ax.set_xlabel('Physical qubits N', fontsize=9)
    axes[0].set_ylabel('QUOPS capability s', fontsize=9)
    fig.suptitle('Pinnacle sensitivity: one choice changed at a time', fontsize=12, y=.98)
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .925), ncol=3,
               frameon=False, fontsize=8, columnspacing=1.4)
    fig.text(.5, .022, 'Pinnacle only; H4 uses the H2 noise fit. Runtime includes expected magic-engine retries.',
             ha='center', color=MUTED, fontsize=8)
    fig.subplots_adjust(left=.09, right=.985, bottom=.14, top=.72, wspace=.12)
    save_plot(fig, path)


def write_json(data, path):
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(data, indent=1, allow_nan=False) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-sensitivity', action='store_true', help='Only compute the baseline comparison.')
    args = parser.parse_args()
    data = compute()
    write_json(data, ROOT/'results/pinnacle.json')
    plot(data, ROOT/'figures/quops_direct_surface_pinnacle')
    for key in KEYS:
        value = data[key]
        if not value['supported']:
            print(f"{key}: Pinnacle unsupported ({value['unsupported_reason']})", flush=True)
            continue
        crossing = value['crossover_vs_direct']
        label = f'{crossing:,.0f}' if crossing is not None else 'none'
        probes = '; '.join(f"s at {N} = {row['s']:,} (w={row['w']}, d={row['distance']}, {row['engine']})"
                           if row else f's at {N} = 0' for N, row in value['probes'].items())
        print(f'{key}: Pinnacle/direct crossover N={label}; {probes}', flush=True)
    if not args.skip_sensitivity:
        sensitivity = compute_sensitivity(data)
        write_json(sensitivity, ROOT/'results/pinnacle_sensitivity.json')
        plot_sensitivity(sensitivity, ROOT/'figures/pinnacle_sensitivity')


if __name__ == '__main__':
    main()
