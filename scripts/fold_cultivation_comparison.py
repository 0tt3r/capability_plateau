"""All-to-all QUOPS comparison with fold-transversal cultivation (Sahay, Tsai, Chang, arXiv:2509.05212v2) in place of
the color-code cultivation of Gidney, Shutty, and Jones.

Fold-transversal cultivation uses non-local operations, so it is used only on H2 and H4, where it replaces the
'cultivation' factory everywhere:
  surface code   the cached surface-code model at 600 ns with ('fold_cultivation', 'litinski') and, for the optimized curve,
                 with every surface-code layout as well ('jws', 'catalogue'), whose inputs may be faulty T, Pinnacle's
                 cultivated inputs, or fold-transversal f = 3 states ('fold f=3', capability/surface_factories.py)
  Pinnacle       the full menu of magic engines with the fold-transversal factory behind the adapter ('surface_fold')
  Frozen         the adapter path of Frozen Pinnacle with the same factories (the CCZ engine is unchanged)
Writes results/fold_cultivation.json with crossovers, capabilities at 1e4 and 1e6, and the designs on each staircase.

    .venv/bin/python scripts/fold_cultivation_comparison.py
"""
import json
import sys
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'scripts'))

from capability import pinnacle, surface                     # noqa: E402
from capability.config import HARDWARE                        # noqa: E402
import make_pinnacle_curves as m                              # noqa: E402

KEYS = ('H2', 'H4')
SURFACE_BASE = replace(m.SURFACE_CH, factories=('fold_cultivation', 'litinski'))
SURFACE_OPT = replace(m.SURFACE_CH, factories=('fold_cultivation', 'litinski', 'jws', 'catalogue'))
PINNACLE = pinnacle.vary(families=('surface_fold', 'surface_jws', 'surface_catalogue', 'reference', 'gb_cat_cult'),
                         supply_filter=())
FROZEN = replace(pinnacle.FROZEN_T, families=('surface_fold', 'surface_jws', 'surface_catalogue'), supply_filter=())


def staircase(rows):
    best, cur, out = {}, 0, []
    for r in rows:
        if r and r.get('passes'):
            best.setdefault(r['qubits'], []).append(r)
    for q in sorted(best):
        r = max(best[q], key=lambda x: x['s'])
        if r['s'] > cur:
            cur = r['s']
            out.append(r)
    return out


def at(steps, N):
    s = 0
    for r in steps:
        if r['qubits'] <= N:
            s = r['s']
    return s


def crossover(steps, s_direct):
    return next((r['qubits'] for r in steps if r['s'] > s_direct), None)


def brief(r):
    keys = ('w', 's', 'ensemble', 'qubits', 'distance', 'factory', 'factories', 'run_seconds',
            'loss_blocks', 'loss_magic', 'loss_synthesis')
    return {k: r.get(k) for k in keys if k in r} if r else None


def summarize(rows, probes, s_direct):
    steps = staircase(rows + [p for p in probes.values() if p])
    return dict(crossover=crossover(steps, s_direct),
                probes={k: brief(v) for k, v in probes.items()},
                winners=sorted({r['factory'] for r in steps}),
                steps=[(r['qubits'], r['s'], r['factory']) for r in steps])


def main():
    refs, _ = m.load_references()
    out = {}
    for key in KEYS:
        hw = HARDWARE[key]
        s_direct = max(r['s'] for r in refs[key]['direct'] if r and r.get('passes'))
        res = dict(direct=s_direct)
        for name, ch in (('surface', SURFACE_BASE), ('surface_optimized', SURFACE_OPT)):
            rows = surface.sweep(hw, ch)
            probes = {f'{N:.0e}': surface.best_at(N, hw, ch, rows) for N in m.PROBES}
            res[name] = summarize(rows, probes, s_direct)
            print(key, name, res[name]['crossover'], {k: v and v['s'] for k, v in probes.items()}, flush=True)
        rows = pinnacle.sweep(hw, m.PINNACLE_CH, PINNACLE)
        probes = {f'{N:.0e}': m.probe_row(N, hw, rows, PINNACLE, m.PINNACLE_CH) for N in m.PROBES}
        res['pinnacle'] = summarize(rows, probes, s_direct)
        print(key, 'pinnacle', res['pinnacle']['crossover'], {k: v and v['s'] for k, v in probes.items()}, flush=True)
        rows = pinnacle.sweep(hw, m.PINNACLE_CH, FROZEN)
        probes = {f'{N:.0e}': m.probe_row(N, hw, rows, FROZEN, m.PINNACLE_CH) for N in m.PROBES}
        res['frozen_adapter'] = summarize(rows, probes, s_direct)
        print(key, 'frozen', res['frozen_adapter']['crossover'], {k: v and v['s'] for k, v in probes.items()}, flush=True)
        out[key] = res
    (ROOT/'results/fold_cultivation.json').write_text(json.dumps(out, indent=1, default=str) + '\n')
    print('written results/fold_cultivation.json')


if __name__ == '__main__':
    main()
