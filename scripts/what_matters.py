"""Change one choice at a time and report how the curves move.

    .venv/bin/python scripts/what_matters.py
writes results/what_matters.json and docs/what_matters_table.md.
"""
import json
import sys
from dataclasses import asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from capability import curves, direct, surface
from capability.config import BASELINE, HARDWARE, N_MAX, vary

VARIANTS = [
    ('baseline', 'all defaults in config.py', {}),
    ('direct: no idle errors', 'untouched qubits never fault', dict(idle_errors=False)),
    ('direct: native R_P', 'any R_P(theta) in one step instead of up to 5', dict(rotation_compile='native')),
    ('direct: row-column-row router', 'instead of the snake router (grid only)', dict(router='rcr')),
    ('surface: serial CNOTs', 'one long-range CNOT at a time', dict(cnot_layer='serial')),
    ('surface: 2 steps per T', 'consuming a T state takes 2d cycles', dict(t_per_injection_steps=2)),
    ('surface: compact layout', '1.5 tiles per logical qubit', dict(tiles_per_data_patch=1.5)),
    ('surface: roomy layout', '3 tiles per logical qubit', dict(tiles_per_data_patch=3.0)),
    ('surface: synthesis charged linearly', 'eps instead of eps^2 per rotation', dict(synthesis_charge='linear')),
    ('surface: synthesis at measured infidelity', '0.14 eps^2 per rotation (scripts/check_synthesis_accumulation.py)',
     dict(synthesis_scale=0.14)),
    ('surface: synthesis budgeted separately', 'sum of eps^2 <= 1/2 on its own; stochastic failure gets all of 1/2',
     dict(synthesis_budget=0.5)),
    ('surface: cultivation only', 'no distillation factories', dict(factories=('cultivation',))),
    ('surface: distillation only', 'Litinski factories, no cultivation', dict(factories=('litinski',))),
    ('surface: P_SC prefactor 0.1', 'instead of 0.03', dict(sc_prefactor=0.1)),
    ('surface: no bimodal threshold gain', 'p_th = 0.01 for bimodal too', dict(sc_threshold_bimodal=0.01)),
    ('surface: 1 ms cycle', 'instead of 1 us; only the one-week limit sees it', dict(cycle_seconds=1e-3)),
    ('both: 10^4 executions', 'all runs must fit in the week (extra copies)', dict(executions=10_000)),
]
KEYS = ('H1', 'H2', 'H3', 'H4')
PROBES = (1e4, N_MAX)


DESIGN = ('w', 'ensemble', 's', 'qubits', 'distance', 'factory', 'factories', 'copies', 'runs_per_copy',
          'run_seconds', 'eps', 'T_per_rotation', 'loss_patches', 'loss_magic', 'loss_synthesis')


def summarize(ch):
    out = {}
    for key in KEYS:
        hw = HARDWARE[key]
        d, s = direct.sweep(hw, ch), surface.sweep(hw, ch)
        out[key] = dict(direct=curves.capability(d, N_MAX), crossover=curves.crossover(d, s))
        for N in PROBES:
            x = surface.best_at(N, hw, ch, s)
            out[key][f'surface_{N:.0e}'] = x['s'] if x else 0
            out[key][f'design_{N:.0e}'] = {k: x[k] for k in DESIGN} if x else None
    return out


def fmt_n(v):
    return '—' if v is None else f'{v:,.0f}'


def fmt_s(v):
    return f'{v:.1e}'.replace('e+0', 'e').replace('e+', 'e') if v >= 1e4 else f'{v:,}'


def main():
    rows = []
    for name, note, change in VARIANTS:
        rows.append(dict(name=name, note=note, change=change, result=summarize(vary(**change))))
        print(name, {k: v['crossover'] for k, v in rows[-1]['result'].items()}, flush=True)
    (ROOT/'results/what_matters.json').write_text(json.dumps(dict(baseline=asdict(BASELINE), rows=rows), indent=1) + '\n')
    lines = ['| Choice | Direct s (H1 / H2 / H3 / H4) | Crossover N (H1 / H2 / H3 / H4) | '
             'Surface s at 10^4 (uniform / bimodal) | Surface s at 10^6 (uniform / bimodal) |',
             '|---|---|---|---|---|']
    for r in rows:
        x = r['result']
        lines.append(f"| {r['name']} | " + ' / '.join(fmt_s(x[k]['direct']) for k in KEYS) + ' | '
                     + ' / '.join(fmt_n(x[k]['crossover']) for k in KEYS) + ' | '
                     + ' / '.join(fmt_s(x[k]['surface_1e+04']) for k in ('H1', 'H3')) + ' | '
                     + ' / '.join(fmt_s(x[k]['surface_1e+06']) for k in ('H1', 'H3')) + ' |')
    (ROOT/'docs/what_matters_table.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
