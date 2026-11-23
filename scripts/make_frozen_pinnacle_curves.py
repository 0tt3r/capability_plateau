"""Frozen Pinnacle QUOPS estimates with an explicit phase-catalysis adapter.

    python scripts/make_frozen_pinnacle_curves.py
    python scripts/make_frozen_pinnacle_curves.py --skip-sensitivity

Writes only new results/frozen_pinnacle*.json and figures/quops_frozen_pinnacle
and figures/frozen_pinnacle_sensitivity PNG/PDF files. Direct, surface-code,
and original Pinnacle references are read from cached results, never rerun.

The Frozen paper publishes a CCZ engine, not a QUOPS T-state engine. The new
curve therefore includes a declared phase-catalysis adaptation and is labeled
an estimate, not a published Frozen QUOPS benchmark.
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

from capability import curves, frozen_pinnacle
from capability.config import BASELINE, HARDWARE, N_MAX

KEYS = ('H1', 'H2', 'H3', 'H4')
SUPPORTED_KEYS = ('H2', 'H4')
PROBES = (1e4, N_MAX)
DIRECT, SURFACE, PINNACLE, FROZEN = '#2a78d6', '#eb6834', '#238b59', '#8056b3'
MUTED = '#777777'
FROZEN_LABEL = 'Frozen + catalysis (estimate)'
HASH_FORMAT = 'canonical-json-sha256-v1'
SOURCES = {
    'frozen': dict(url='https://arxiv.org/html/2609.21249v1', version='v1', revised='2026-09-18'),
    'phase_catalysis': dict(url='https://arxiv.org/html/1812.01238v3', version='v3', section='4, Figure 14'),
    'original_pinnacle': dict(url='https://arxiv.org/html/2602.11457v2', version='v2'),
}
MODELING_NOTES = [
    'This is a Frozen Pinnacle plus phase-catalysis QUOPS estimate, not a published Frozen QUOPS result.',
    'Use only the published distance-24 Frozen code, eight processing logicals per block, and CCZ engine.',
    'The CCZ engine has error 2e-10, rejection 1%, and five logical cycles per Toffoli.',
    'Phase catalysis pairs T gates on distinct QUOPS qubits and consumes one CCZ per pair; it is not free CCZ-to-T conversion.',
    'Reserve a catalyst, an AND scratch, and three dedicated magic-port ancillas on three distinct blocks.',
    'The seed preparation and transfer interface is a declared adapter and is included in qubits, timing, and loss.',
    'Charge CCZ error once per batch and seed error once per execution; catalytic output errors are correlated.',
    'The baseline adds two scratch preparation/readout slots per batch and uses serial execution.',
    'Week throughput uses expected seed and CCZ retry time, not a probabilistic deadline guarantee.',
    'H1 and H3 are not modeled: fixed connectivity can include long-range edges, and a 2D grid realization needs placement, routing, and error costs; this is not an impossibility claim.',
    'H4 uses the same p=1e-3 depolarizing-noise fit as H2, without a claimed bimodal-noise gain.',
    'Direct, surface-code, and original Pinnacle curves are cached references and are never recomputed here.',
]


def normalized(value):
    """Normalize dataclass tuples to the list representation stored in JSON."""
    return json.loads(json.dumps(value))


def canonical_json_sha256(value):
    """Hash parsed JSON independently of line endings, whitespace, and key order."""
    canonical = json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True, allow_nan=False)
    return hashlib.sha256(canonical.encode('ascii')).hexdigest()


def load_references(ch=BASELINE):
    """Validate compatible cached inputs and record their unchanged source hashes."""
    references, hashes = {}, {}
    for name, filename in (('baseline', 'baseline.json'), ('pinnacle', 'pinnacle.json')):
        path = ROOT/'results'/filename
        raw = path.read_bytes()
        data = json.loads(raw)
        if data.get('choices') != normalized(asdict(ch)):
            raise ValueError(f'Cached results/{filename} choices do not match the requested shared choices.')
        for key in KEYS:
            if data[key]['hardware'] != normalized(asdict(HARDWARE[key])):
                raise ValueError(f'Cached results/{filename} hardware does not match config.py for {key}.')
        references[name] = data
        hashes[filename] = canonical_json_sha256(data)
    if references['pinnacle'].get('provenance_hash_format') != HASH_FORMAT:
        raise ValueError('Cached results/pinnacle.json has a missing or unsupported provenance hash format. '
                         'Regenerate it with python scripts/make_pinnacle_curves.py before running '
                         'python scripts/make_frozen_pinnacle_curves.py.')
    if references['pinnacle'].get('baseline_sha256') != hashes['baseline.json']:
        raise ValueError('Cached original Pinnacle results used different baseline inputs. '
                         'Regenerate them with python scripts/make_pinnacle_curves.py before running '
                         'python scripts/make_frozen_pinnacle_curves.py.')
    return references['baseline'], references['pinnacle'], hashes


def with_probes(rows, probes):
    """Include exact fixed-budget designs in the sampled plotting staircase."""
    combined = {(r['ensemble'], r['w']): r for r in rows}
    for row in probes.values():
        if row is not None:
            combined[row['ensemble'], row['w']] = row
    return sorted(combined.values(), key=lambda r: (r['ensemble'], r['w']))


def compute_frozen(hw, ch=BASELINE, fc=frozen_pinnacle.BASELINE):
    reason = frozen_pinnacle.support_reason(hw)
    if reason is not None:
        return dict(supported=False, unsupported_reason=reason, frozen_pinnacle=[],
                    probes={f'{N:.0e}': None for N in PROBES})
    rows = frozen_pinnacle.sweep(hw, ch, fc)
    probes = {f'{N:.0e}': frozen_pinnacle.best_at(N, hw, ch, rows, fc) for N in PROBES}
    return dict(supported=True, unsupported_reason=None,
                frozen_pinnacle=with_probes(rows, probes), probes=probes)


def compute(ch=BASELINE, fc=frozen_pinnacle.BASELINE):
    baseline, original, hashes = load_references(ch)
    out = dict(choices=asdict(ch), frozen_choices=asdict(fc), sources=SOURCES,
               model_label=FROZEN_LABEL, provenance_hash_format=HASH_FORMAT, reference_sha256=hashes,
               source_of_reference='cached results/baseline.json and results/pinnacle.json; no reference models rerun',
               original_pinnacle_choices=original['pinnacle_choices'], modeling_notes=MODELING_NOTES)
    for key in KEYS:
        record = compute_frozen(HARDWARE[key], ch, fc)
        record.update(hardware=asdict(HARDWARE[key]), direct=baseline[key]['direct'],
                      surface=baseline[key]['surface'], pinnacle=original[key]['pinnacle'],
                      reference_probes=dict(surface=baseline[key]['probes'], pinnacle=original[key]['probes']))
        record['crossover_vs_direct'] = (curves.crossover(record['direct'], record['frozen_pinnacle'])
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
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 7.2), sharex=True, sharey=True)
    for ax, key in zip(axes.flat, KEYS):
        value = data[key]
        draw_curve(ax, value['direct'], DIRECT, 'Direct', lw=1.6)
        draw_curve(ax, value['surface'], SURFACE, 'Surface code', lw=1.6)
        draw_curve(ax, value['pinnacle'], PINNACLE, 'Original Pinnacle', lw=1.8)
        if value['supported']:
            draw_curve(ax, value['frozen_pinnacle'], FROZEN, FROZEN_LABEL, lw=2.1)
            crossing = value['crossover_vs_direct']
            if crossing is not None:
                ax.axvline(crossing, color=FROZEN, lw=.9, ls=':', alpha=.65)
                ax.text(.045, .89, f'Frozen exceeds direct\nat {crossing:,.0f} qubits',
                        transform=ax.transAxes, fontsize=7.2, color=FROZEN, va='top')
            small = value['probes'][f'{PROBES[0]:.0e}']
            top = value['probes'][f'{N_MAX:.0e}']
            lines = [f"Frozen at $10^4$: {'no design' if small is None else format(small['s'], ',')}"]
            if top is not None:
                lines.append(f"Frozen at $10^6$: s = {top['s']:,}")
                ax.axvspan(top['qubits'], N_MAX, color=FROZEN, alpha=.035)
            ax.text(.045, .70, '\n'.join(lines), transform=ax.transAxes,
                    fontsize=7.2, color=FROZEN, va='top')
        else:
            ax.text(.045, .86, '2D grid mapping not modeled\nPlacement and routing costs needed',
                    transform=ax.transAxes, fontsize=8, color=FROZEN, va='top',
                    bbox=dict(facecolor='white', edgecolor='none', alpha=.9, pad=2))
        direct_best = max((r['s'] for r in value['direct'] if r['passes']), default=0)
        if direct_best:
            ax.text(2, direct_best*1.7, f'direct: s = {direct_best:,}', color=DIRECT, fontsize=7.3)
        ax.set_title(f"{key}: {value['hardware']['name']}", fontsize=10)
        decorate(ax)
    for ax in axes[1]:
        ax.set_xlabel('Physical qubits N', fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel('QUOPS capability s', fontsize=9)
    handles = [Line2D([0], [0], color=color, lw=2, label=label) for color, label in
               ((DIRECT, 'Direct (cached)'), (SURFACE, 'Surface code (cached)'),
                (PINNACLE, 'Original Pinnacle (cached)'), (FROZEN, FROZEN_LABEL))]
    fig.suptitle('QUOPS capability with original and Frozen Pinnacle', fontsize=12, y=.98)
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .945), ncol=2,
               frameon=False, fontsize=8.5, columnspacing=2)
    fig.text(.5, .034, 'Serial modes; one-week mean-time budget. Frozen includes a proposed phase-catalysis adapter.',
             ha='center', color=MUTED, fontsize=7.5)
    fig.text(.5, .014, 'H4 uses the H2 depolarizing-noise fit. Fixed degree eight is not a nearest-neighbor embedding.',
             ha='center', color=MUTED, fontsize=7.5)
    fig.subplots_adjust(left=.09, right=.985, bottom=.12, top=.82, wspace=.12, hspace=.28)
    save_plot(fig, path)


def sensitivity_variants():
    return [
        ('baseline adapter', {}, {}, FROZEN, '-'),
        ('four scratch slots', {}, dict(extra_batch_slots=4), '#d17a22', '--'),
        ('one Pauli sector', {}, dict(logical_error_scale=1), '#58a65c', '--'),
        ('processing-only block (derived)', {}, dict(processing_qubits=2372), '#4c78a8', '--'),
        ('1 ms cycle', dict(cycle_seconds=1e-3), {}, '#9467bd', ':'),
        ('10,000 executions', dict(executions=10_000), {}, '#c44e52', '-'),
        ('measured synthesis charge', dict(synthesis_scale=.14), {}, '#555555', ':'),
    ]


def compute_sensitivity(baseline=None, ch=BASELINE, fc=frozen_pinnacle.BASELINE):
    out = dict(sources=SOURCES, model_label=FROZEN_LABEL, choices=asdict(ch), frozen_choices=asdict(fc),
               cached_reference_used=False, modeling_notes=MODELING_NOTES[:-1], variants=[])
    for name, shared, specific, color, style in sensitivity_variants():
        changed_ch, changed_fc = replace(ch, **shared), replace(fc, **specific)
        result = {}
        for key in SUPPORTED_KEYS:
            if name == 'baseline adapter' and baseline is not None:
                value = baseline[key]
                record = {k: value[k] for k in ('supported', 'unsupported_reason', 'frozen_pinnacle', 'probes')}
            else:
                record = compute_frozen(HARDWARE[key], changed_ch, changed_fc)
            result[key] = dict(hardware=asdict(HARDWARE[key]), **record)
        out['variants'].append(dict(name=name, shared_changes=shared, frozen_changes=specific,
                                    choices=asdict(changed_ch), frozen_choices=asdict(changed_fc),
                                    color=color, linestyle=style, result=result))
        print('sensitivity:', name,
              {key: {N: row['s'] if row else 0 for N, row in value['probes'].items()}
               for key, value in result.items()}, flush=True)
    return out


def plot_sensitivity(data, path):
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 5.0), sharex=True, sharey=True)
    handles = []
    for variant in data['variants']:
        for ax, key in zip(axes, SUPPORTED_KEYS):
            draw_curve(ax, variant['result'][key]['frozen_pinnacle'], variant['color'], variant['name'],
                       lw=2 if variant['name'] == 'baseline adapter' else 1.5, ls=variant['linestyle'])
        handles.append(Line2D([0], [0], color=variant['color'], lw=2,
                              ls=variant['linestyle'], label=variant['name']))
    for ax, key in zip(axes, SUPPORTED_KEYS):
        decorate(ax)
        ax.set_xlim(1e4, N_MAX)
        ax.set_title(f'{key}: {HARDWARE[key].name}', fontsize=10)
        ax.set_xlabel('Physical qubits N', fontsize=9)
    axes[0].set_ylabel('QUOPS capability s', fontsize=9)
    fig.suptitle('Frozen + catalysis: sensitivity of the proposed QUOPS adapter', fontsize=11.5, y=.98)
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .925), ncol=3,
               frameon=False, fontsize=7.7, columnspacing=1.3)
    fig.text(.5, .022, 'Frozen estimates only; expected retry time. The smaller processing-only footprint is derived.',
             ha='center', color=MUTED, fontsize=7.6)
    fig.subplots_adjust(left=.09, right=.985, bottom=.14, top=.71, wspace=.12)
    save_plot(fig, path)


def write_json(data, path):
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(data, indent=1, allow_nan=False) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-sensitivity', action='store_true', help='Only compute the baseline comparison.')
    args = parser.parse_args()
    data = compute()
    write_json(data, ROOT/'results/frozen_pinnacle.json')
    plot(data, ROOT/'figures/quops_frozen_pinnacle')
    for key in KEYS:
        value = data[key]
        if not value['supported']:
            print(f"{key}: Frozen unsupported ({value['unsupported_reason']})", flush=True)
            continue
        crossing = value['crossover_vs_direct']
        label = f'{crossing:,.0f}' if crossing is not None else 'none'
        probes = '; '.join(f"s at {N} = {row['s']:,} (w={row['w']}, qubits={row['qubits']:,.0f})"
                           if row else f's at {N}: no feasible design' for N, row in value['probes'].items())
        print(f'{key}: Frozen/direct crossover N={label}; {probes}', flush=True)
    if not args.skip_sensitivity:
        sensitivity = compute_sensitivity(data)
        write_json(sensitivity, ROOT/'results/frozen_pinnacle_sensitivity.json')
        plot_sensitivity(sensitivity, ROOT/'figures/frozen_pinnacle_sensitivity')


if __name__ == '__main__':
    main()
