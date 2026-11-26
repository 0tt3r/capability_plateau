"""Magic engines and factories that set a point on the exact all-to-all curves (H2 and H4), with full menus.

Every even width of both ensembles is solved (min_qubits) with the unfiltered menu, the staircase in N is formed, and the
designs on it are listed. These lists are the winners kept in capability/pinnacle.py (PINNACLE_WINNERS, FROZEN_WINNERS)
and capability/factories.py (FOLD_JWS_WINNERS, FOLD_CATALOGUE_WINNERS). On all-to-all hardware the color-code
cultivation factory is replaced by fold-transversal cultivation ('fold_cultivation', 'surface_fold').

    .venv/bin/python scripts/find_winners.py            # writes results/winners.json
"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from capability import pinnacle, surface                     # noqa: E402
from capability.config import BASELINE, HARDWARE, N_MAX       # noqa: E402
from capability.quops import ENSEMBLES, size                  # noqa: E402

PINNACLE_CH = replace(BASELINE, cycle_seconds=1.2e-6)
SURFACE_CH = replace(BASELINE, cycle_seconds=6e-7, factories=('fold_cultivation', 'litinski', 'jws', 'catalogue'))
PINNACLE_FULL = pinnacle.vary(families=('surface_fold', 'surface_jws', 'surface_catalogue', 'reference', 'gb_cat_cult'),
                              supply_filter=())
FROZEN_FULL = replace(pinnacle.FROZEN_T, families=('surface_fold', 'surface_jws', 'surface_catalogue'), supply_filter=())
MODELS = ('pinnacle', 'frozen', 'surface')


def solve(args):
    model, key, ensemble, w = args
    hw = HARDWARE[key]
    if model == 'surface':
        r = surface.min_qubits(w, ensemble, hw, SURFACE_CH)
    else:
        r = pinnacle.min_qubits(w, ensemble, hw, PINNACLE_CH, PINNACLE_FULL if model == 'pinnacle' else FROZEN_FULL)
    if r is None or r['qubits'] > N_MAX:
        return None
    return dict(model=model, key=key, ensemble=ensemble, w=w, s=size(w, ensemble), qubits=r['qubits'], factory=r['factory'])


def staircase(rows):
    best, cur, out = {}, 0, []
    for r in rows:
        best.setdefault(r['qubits'], []).append(r)
    for q in sorted(best):
        r = max(best[q], key=lambda x: x['s'])
        if r['s'] > cur:
            cur = r['s']
            out.append(r)
    return out


def run(model, key, pool, chunk=64):
    rows = []
    for ensemble in ENSEMBLES:
        w = 2
        while True:
            batch = [(model, key, ensemble, w + 2*i) for i in range(chunk)]
            got = list(pool.map(solve, batch))
            rows += [g for g in got if g]
            if any(g is None for g in got):      # the qubits needed never decrease with w
                break
            w += 2*chunk
    return staircase(rows)


def main():
    out = {}
    with ProcessPoolExecutor() as pool:
        for model in MODELS:
            for key in ('H2', 'H4'):
                steps = run(model, key, pool)
                out[f'{model}_{key}'] = dict(winners=sorted({r['factory'] for r in steps}),
                                             steps=[(r['qubits'], r['s'], r['ensemble'], r['w'], r['factory']) for r in steps])
                print(model, key, len(steps), out[f'{model}_{key}']['winners'], flush=True)
    (ROOT/'results/winners.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
