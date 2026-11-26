"""QUOPS capability curves for Pinnacle fault tolerance (draft Secs. 3.6, 3.7, 4.10 and the Frozen appendix).

    python scripts/make_pinnacle_curves.py
    python scripts/make_pinnacle_curves.py --skip-sensitivity

Writes results/pinnacle.json, figures/quops_direct_surface_pinnacle.{pdf,png} (main figure), and
figures/quops_pinnacle_appendix.{pdf,png}. The direct and surface-code curves are cached baseline references and are
never recomputed here, except the optimized surface-code curve, which reruns the cached surface-code model with the
surface-code factories of the arXiv:2610.06535 protocols added to its menu (family 'jws').

Pinnacle and Frozen Pinnacle use a code cycle of 12 t_op = 1.2 us (the surface code uses 6 t_op); 1 us (the published
value) and 0.8 us are sensitivities. The shared baseline choices recorded for provenance remain those of the cached
baseline; the Pinnacle cycle is recorded under pinnacle_choices.

The provenance helpers (HASH_FORMAT, canonical_json_sha256, load_references) and the per-hardware record layout
(supported, pinnacle, probes, direct, surface, reference_probes, crossover_vs_direct) are those of the earlier serial
Pinnacle script, so scripts/make_frozen_pinnacle_curves.py and tests/test_pinnacle_provenance.py keep working.
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
import numpy as np

from capability import curves, frozen_pinnacle, pinnacle, surface
from capability.config import BASELINE, HARDWARE, N_MAX

PROBES = (1e4, N_MAX)
KEYS = ('H1', 'H2', 'H3', 'H4')
SUPPORTED_KEYS = ('H2', 'H4')
DIRECT, SURFACE = '#2a78d6', '#eb6834'
PINNACLE, PINNACLE_LIGHT, FROZEN = '#087858', '#6fbfa6', '#8c564b'   # the green of Fig. 7; STAR uses purple, direct blue, surface orange
INK, MUTED = '#2b2b2b', '#8a8a8a'
PAPER = 'https://arxiv.org/html/2602.11457v2'
HASH_FORMAT = 'canonical-json-sha256-v1'
PINNACLE_CYCLE = 12*1e-7                                   # 12 t_op with t_op = 100 ns
PINNACLE_CH = replace(BASELINE, cycle_seconds=PINNACLE_CYCLE)
SURFACE_CYCLE = 6e-7                                        # 6 t_op: the surface-code cycle of Sec. 3.5
SURFACE_CH = replace(BASELINE, cycle_seconds=SURFACE_CYCLE)   # the surface-code comparison curve (cached model, 600 ns)
SURFACE_OPTIMIZED = replace(SURFACE_CH, factories=('cultivation', 'jws_winners', 'catalogue_winners'))   # the options that set a point on the extended curve
FROZEN_QUOPS = frozen_pinnacle.vary(final_measurements=False)  # N_M = 0 for QUOPS, as in the surface-code model
FAMILIES = ('reference', 'surface', 'surface_jws', 'surface_catalogue', 'gb_cat_cult')
PLATEAU_WIDTH = 936                                         # largest deep width that completes within the week at 1.2 us
# One row per supply that can win somewhere on the Pinnacle curve; the filters are substrings of the supply names of
# capability/pinnacle_factories.py.
VERDICT = (('Cultivation behind an adapter', ('cultivation, adapter',), None),
           ('Published 15-to-1 engines', ('15-to-1 (Pinnacle',), None),
           ('[[28,1,3]] engine with cultivated inputs', ('[[28,1,3]] T0, cultivated engine d_a=7',), ('gb_cat_cult',)),
           ('[[32,1,3]] engine with cultivated inputs', ('[[32,1,3]] T0, cultivated engine d_a=7',), ('gb_cat_cult',)),
           ('[[47,5,3]] engine with cultivated inputs', ('[[47,5,3]] T0·T1·T2·T3·T4, cultivated engine d_a=9',), ('gb_cat_cult',)),
           ('[[56,8,3]] engine with cultivated inputs', ('[[56,8,3]] T0·T1·T2·T3·T4·T5·T6·T7, cultivated engine d_a=9',), ('gb_cat_cult',)),
           ('Layout of [[14,1,2]] with faulty T, behind an adapter', ('[[14,1,2]] T0, faulty T, d=(7,5,5)',), ('surface_catalogue',)),
           ('Layout of [[8,3,2]] (CCZ, catalysed to T) with faulty T, behind an adapter', ('[[8,3,2]] CCZ012, faulty T, d=(7,5,5)',), ('surface_catalogue',)),
           ('Two-level layout of [[14,2,2]] fed by 15-to-1 factories, behind an adapter', ('[[14,2,2]] T0,T1, two levels',), ('surface_catalogue',)),
           ('Two-level layout of [[20,3,2]] fed by 15-to-1 factories, behind an adapter', ('[[20,3,2]] T0,T1,T2, two levels',), ('surface_catalogue',)))
MODELING_NOTES = [
    'Serial Pauli-based computation on one processing unit; the supply (magic engine, surface-code factory behind an adapter, or in-place injection) is a searched black box.',
    'Logical measurement and magic-engine inputs are from arXiv:2602.11457v2, Sections IV-V; the code cycle is 12 t_op = 1.2 us.',
    'H1 and H3 are not modeled: a 2D nearest-neighbor realization needs placement, routing, and error costs; this is not an impossibility claim.',
    'H4 uses the same p=1e-3 depolarizing-noise fit as H2, without a claimed bimodal-noise gain.',
    'Rejections stall the unit; the blocks are charged for the mean stall (relative difference from the geometric average below 1e-6).',
    'The baseline counts both Pauli error sectors using eta=2; this is a declared choice.',
    'Synthesis infidelity is charged within the shared loss budget, separately from the stochastic no-fault loss.',
    'Direct and surface-code curves are cached baseline references; the optimized surface-code curve reruns only the cached surface-code model with more factories.',
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


def probe_row(N, hw, rows, pc, ch):
    x = pinnacle.best_at(N, hw, ch, rows, pc)
    keys = ('s', 'w', 'ensemble', 'qubits', 'distance', 'factory', 'family', 'factories', 'T_per_rotation',
            'run_seconds', 'loss_blocks', 'loss_magic', 'loss_synthesis')
    return None if x is None else dict({k: x[k] for k in keys}, passes=True)


def compute_pinnacle(hw, ch=PINNACLE_CH, pc=pinnacle.BASELINE):
    reason = pinnacle.support_reason(hw, pc)
    if reason is not None:
        return dict(supported=False, unsupported_reason=reason, pinnacle=[],
                    probes={f'{N:.0e}': None for N in PROBES})
    rows = pinnacle.sweep(hw, ch, pc)
    probes = {f'{N:.0e}': probe_row(N, hw, rows, pc, ch) for N in PROBES}
    return dict(supported=True, unsupported_reason=None, pinnacle=with_probes(rows, probes), probes=probes)


def frozen_curves(hw, ch, direct):
    """Frozen Pinnacle: its published CCZ engine with catalysis, and T states delivered through an adapter."""
    cc = frozen_pinnacle.sweep(hw, ch, FROZEN_QUOPS)
    tt = pinnacle.sweep(hw, ch, pinnacle.FROZEN_T)
    both = cc + tt
    probes = lambda rows, pc=None: {f'{N:.0e}': (probe_row(N, hw, rows, pc, ch) if pc else frozen_pinnacle.best_at(N, hw, ch, rows, FROZEN_QUOPS))
                                    for N in PROBES}
    sens = {}
    for name, c2, fc in (('four scratch slots per pair', ch, replace(FROZEN_QUOPS, extra_batch_slots=4)),
                         ('code cycle 1 us (published)', replace(ch, cycle_seconds=1e-6), FROZEN_QUOPS),
                         ('code cycle 1 ms', replace(ch, cycle_seconds=1e-3), FROZEN_QUOPS),
                         ('10^4 executions', replace(ch, executions=10**4), FROZEN_QUOPS)):
        r2 = frozen_pinnacle.sweep(hw, c2, fc)
        sens[name] = dict(crossover=curves.crossover(direct, r2),
                          probes={f'{N:.0e}': frozen_pinnacle.best_at(N, hw, c2, r2, fc) for N in PROBES})
    return dict(ccz=dict(rows=cc, crossover=curves.crossover(direct, cc), probes=probes(cc)),
                adapter_t=dict(rows=tt, crossover=curves.crossover(direct, tt), probes=probes(tt, pinnacle.FROZEN_T)),
                optimized=dict(rows=both, crossover=curves.crossover(direct, both)), sensitivity=sens)


def sensitivities(hw, direct):
    cases = {'baseline': (pinnacle.BASELINE, PINNACLE_CH)}
    for fam in FAMILIES:
        cases[f'only {fam}'] = (pinnacle.vary(families=(fam,)), PINNACLE_CH)
    cases['eta = 1'] = (pinnacle.vary(eta=1.0), PINNACLE_CH)
    cases['reaction 20 t_c'] = (pinnacle.vary(reaction=20), PINNACLE_CH)
    cases['adapter d_a = 15 only'] = (pinnacle.vary(adapter_distances=(15,)), PINNACLE_CH)
    for label, seconds in (('code cycle 1 us (published)', 1e-6), ('code cycle 0.8 us', 0.8e-6)):
        cases[label] = (pinnacle.BASELINE, replace(BASELINE, cycle_seconds=seconds))   # t_op fixed: surface code stays at 600 ns
    # A 1 ms Pinnacle cycle means t_op = 1 ms / 12, so the 6-step surface-code cycle of the factories becomes 0.5 ms.
    cases['code cycle 1 ms'] = (pinnacle.vary(surface_cycle_seconds=0.5e-3), replace(BASELINE, cycle_seconds=1e-3))
    cases['linear synthesis charge'] = (pinnacle.BASELINE, replace(PINNACLE_CH, synthesis_charge='linear'))
    cases['10^4 executions'] = (pinnacle.BASELINE, replace(PINNACLE_CH, executions=10**4))
    out = {}
    for name, (pc, ch) in cases.items():
        rows = pinnacle.sweep(hw, ch, pc)
        out[name] = dict(crossover=curves.crossover(direct, rows),
                         probes={f'{N:.0e}': probe_row(N, hw, rows, pc, ch) for N in PROBES})
    return out


def supply_verdict(hw, direct):
    """Each supply used alone: crossover with direct, capability at the probes, qubits for the plateau design."""
    out = {}
    for name, flt, fam in VERDICT:
        pc = pinnacle.vary(supply_filter=flt, families=fam) if fam else pinnacle.vary(supply_filter=flt)
        rows = pinnacle.sweep(hw, PINNACLE_CH, pc)
        plateau = pinnacle.min_qubits(PLATEAU_WIDTH, 'deep', hw, PINNACLE_CH, pc)
        p = {f'{N:.0e}': probe_row(N, hw, rows, pc, PINNACLE_CH) for N in PROBES}
        out[name] = dict(crossover=curves.crossover(direct, rows), s_1e4=p['1e+04'] and p['1e+04']['s'],
                         s_1e6=p['1e+06'] and p['1e+06']['s'], plateau_qubits=plateau and plateau['qubits'],
                         largest_width_deep=max([r['w'] for r in rows if r['ensemble'] == 'deep'], default=None))
    return out


def winners(hw):
    """Supply that wins at every even width of each ensemble, for the best-of-menu curve."""
    counts = {}
    for ensemble in ('shallow', 'deep'):
        for w in range(2, 4000, 2):
            r = pinnacle.min_qubits(w, ensemble, hw, PINNACLE_CH)
            if r is None:
                break
            name = 'cultivation' if r['factory'].startswith('cultivation') else r['factory']
            counts.setdefault(name, []).append((ensemble, w, r['qubits']))
    return {k: dict(widths=len(v), min_qubits=min(x[2] for x in v), max_qubits=max(x[2] for x in v)) for k, v in counts.items()}


def surface_sc(hw, direct):
    """The cached surface-code model at the 600 ns cycle of Sec. 3.5, the curve that Pinnacle is compared with."""
    rows = surface.sweep(hw, SURFACE_CH)
    probes = {f'{N:.0e}': surface.best_at(N, hw, SURFACE_CH, rows) for N in PROBES}
    return dict(rows=with_probes(rows, probes), probes=probes, crossover=curves.crossover(direct, rows))


def surface_optimized(key, hw, base):
    """Cached surface-code model with the surface-code factories of the 2610.06535 protocols added to its menu."""
    rows = surface.sweep(hw, SURFACE_OPTIMIZED)
    probes = {f'{N:.0e}': surface.best_at(N, hw, SURFACE_OPTIMIZED, rows) for N in PROBES}
    return dict(rows=with_probes(rows, probes), probes=probes, crossover=curves.crossover(base[key]['direct'], rows))


def compute(ch=BASELINE, pc=pinnacle.BASELINE, extended=True):
    reference, digest = load_references(ch)
    out = dict(choices=asdict(ch), pinnacle_choices=dict(asdict(pc), cycle_seconds=PINNACLE_CYCLE),
               source_paper=dict(url=PAPER, version='v2', revised='2026-05-05'),
               provenance_hash_format=HASH_FORMAT,
               baseline_sha256=digest,
               source_of_reference='cached results/baseline.json; legacy models were not rerun',
               modeling_notes=MODELING_NOTES)
    for key in KEYS:
        hw = HARDWARE[key]
        record = compute_pinnacle(hw, PINNACLE_CH, pc)
        record.update(hardware=asdict(hw), direct=reference[key]['direct'],
                      surface=reference[key]['surface'], reference_probes=reference[key]['probes'])
        record['surface_sc'] = surface_sc(hw, record['direct'])
        record['crossover_vs_direct'] = (curves.crossover(record['direct'], record['pinnacle'])
                                         if record['supported'] else None)
        if extended:
            ext = dict(surface_optimized=surface_optimized(key, hw, reference))
            if record['supported']:
                rivals = record['direct'] + record['surface_sc']['rows']   # surface code at 600 ns; the extended menu is in the appendix
                record['crossover_vs_all'] = curves.crossover(rivals, record['pinnacle'])
                ext.update(per_family={f: dict(rows=pinnacle.sweep(hw, PINNACLE_CH, pinnacle.vary(families=(f,))))
                                       for f in FAMILIES},
                           frozen=frozen_curves(hw, PINNACLE_CH, record['direct']),
                           sensitivity=sensitivities(hw, record['direct']),
                           supply_verdict=supply_verdict(hw, record['direct']), winners=winners(hw))
                for f, v in ext['per_family'].items():
                    v['crossover'] = curves.crossover(record['direct'], v['rows'])
            record['extended'] = ext
        out[key] = record
    return out


def draw_curve(ax, rows, color, label, **kwargs):
    y = curves.curve(rows).astype(float)
    y[y == 0] = np.nan
    return ax.step(curves.N_GRID, y, where='post', color=color, label=label, **kwargs)


def decorate(ax, title):
    ax.set(xscale='log', yscale='log', xlim=(1, N_MAX), ylim=(1, 1e9), title=title)
    ax.title.set_fontsize(10)
    ax.grid(alpha=.15)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)


def label_direct(ax, rows):
    d = max((r for r in rows if r['passes']), key=lambda r: (r['s'], r['ensemble'] == 'deep'))
    ax.text(2, d['s']*1.8, f"direct: s = {d['s']:,} ({d['ensemble']}, w = {d['w']})", color=INK, fontsize=8, va='bottom')


def crossover_line(ax, x, color, text, y, ha='left'):
    ax.axvline(x, color=color, lw=1, ls=':')
    ax.text(x*(1.2 if ha == 'left' else 1/1.2), y, f'{text}\n{x:,.0f} qubits', color=INK, fontsize=7.5, va='top' if y > 1e5 else 'center',
            ha=ha, bbox=dict(fc='white', ec='none', alpha=.85, pad=1))


def save_plot(fig, path):
    path.parent.mkdir(exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(path.with_suffix('.' + ext), dpi=180)
    plt.close(fig)


def plot_main(data, path):
    fig, axes = plt.subplots(2, 2, figsize=(9, 7.0), sharex=True, sharey=True)
    for ax, key in zip(axes.flat, KEYS):
        v = data[key]
        draw_curve(ax, v['direct'], DIRECT, 'Direct', lw=2)
        draw_curve(ax, v['surface_sc']['rows'], SURFACE, 'Surface code', lw=2)
        if v['supported']:
            draw_curve(ax, v['pinnacle'], PINNACLE, 'Pinnacle, optimized magic engines', lw=2.2)
            crossover_line(ax, v['crossover_vs_all'], PINNACLE, 'Pinnacle crossover', 4e8, ha='right')
        else:
            ax.text(3, 3e5, 'Pinnacle: 2D-grid mapping\nnot modeled', color=PINNACLE, fontsize=8)
        crossover = curves.crossover(v['direct'], v['surface_sc']['rows'])
        if crossover:
            crossover_line(ax, crossover, SURFACE, 'surface-code crossover', 4e8 if not v['supported'] else 30)
        label_direct(ax, v['direct'])
        decorate(ax, f"{key}: {v['hardware']['name']}")
    for ax in axes[1]:
        ax.set_xlabel('Physical qubits N')
    for ax in axes[:, 0]:
        ax.set_ylabel('QUOPS capability s')
    handles, labels = axes[0, 1].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=2, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    save_plot(fig, path)


def plot_appendix(data, path):
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.6), sharex=True, sharey=True)
    for ax, key in zip(axes, SUPPORTED_KEYS):
        v = data[key]; e = v['extended']
        draw_curve(ax, v['direct'], DIRECT, 'Direct', lw=2)
        draw_curve(ax, v['surface_sc']['rows'], SURFACE, 'Surface code', lw=2)
        draw_curve(ax, e['surface_optimized']['rows'], SURFACE, 'Surface code, optimized magic-state factories', lw=1.6, ls='--')
        draw_curve(ax, v['pinnacle'], PINNACLE, 'Pinnacle, optimized magic engines', lw=1.4)
        draw_curve(ax, e['per_family']['reference']['rows'], PINNACLE_LIGHT, 'Pinnacle, published 15-to-1 engines only', lw=1.6, ls='--')
        draw_curve(ax, e['frozen']['ccz']['rows'], FROZEN, 'Frozen Pinnacle, published CCZ engine with catalysis', lw=2.2)
        draw_curve(ax, e['frozen']['optimized']['rows'], FROZEN, 'Frozen Pinnacle, optimized magic engines', lw=1.6, ls='--')
        decorate(ax, f"{key}: {v['hardware']['name']}")
        ax.set_xlabel('Physical qubits N')
        label_direct(ax, v['direct'])
    axes[0].set_ylabel('QUOPS capability s')
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=2, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.17, 1, 1))
    save_plot(fig, path)


def write_json(data, path):
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(data, indent=1, allow_nan=False) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-sensitivity', action='store_true', help='Skip the extended analysis (curves, supply comparison, Frozen).')
    args = parser.parse_args()
    data = compute(extended=not args.skip_sensitivity)
    write_json(data, ROOT/'results/pinnacle.json')
    for key in KEYS:
        v = data[key]
        if not v['supported']:
            print(f"{key}: Pinnacle unsupported ({v['unsupported_reason']})", flush=True)
            continue
        probes = '; '.join(f"s at {N} = {row['s']:,} (w={row['w']}, {row['ensemble']}, d={row['distance']}, {row['factories']} x {row['factory']}, N={row['qubits']:,.0f})"
                           if row else f's at {N} = 0' for N, row in v['probes'].items())
        print(f"{key}: Pinnacle/direct crossover {v['crossover_vs_direct']}; {probes}", flush=True)
    if not args.skip_sensitivity:
        plot_main(data, ROOT/'figures/quops_direct_surface_pinnacle')
        plot_appendix(data, ROOT/'figures/quops_pinnacle_appendix')
        for key in SUPPORTED_KEYS:
            e = data[key]['extended']
            print(key, 'crossover vs direct and cached surface', data[key]['crossover_vs_all'],
                  '| surface (600 ns) crossover', data[key]['surface_sc']['crossover'],
                  '| surface (extended menu) crossover', e['surface_optimized']['crossover'],
                  '| surface extended s at probes', {N: (r['s'] if r else 0) for N, r in e['surface_optimized']['probes'].items()},
                  '| Frozen ccz', e['frozen']['ccz']['crossover'], 'adapter-T', e['frozen']['adapter_t']['crossover'],
                  'optimized', e['frozen']['optimized']['crossover'], flush=True)
            for name, sv in e['supply_verdict'].items():
                print(f"    {name:70s} cross={sv['crossover']} s@1e4={sv['s_1e4']} s@1e6={sv['s_1e6']} plateau qubits={sv['plateau_qubits']}")
            print('    winners:', {k: (x['widths'], x['min_qubits'], x['max_qubits']) for k, x in e['winners'].items()})


if __name__ == '__main__':
    main()
