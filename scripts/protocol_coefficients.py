"""Extract the pure-T factory protocols used by the Pinnacle magic-engine menu.

Source: master_catalog/master_catalog.json of github.com/shubhamj810/magic-state-factory-repo, the database of
Jain, Wills, and Singh, arXiv:2610.06535. Every catalog row is a postselected diagonal CSS protocol that consumes
n noisy T states and outputs T^(x)k, built from N = k + (checks) wires. A_d is the number of weight-d input-error
patterns that go undetected, so the probability that a batch is bad is A_d p_rot^d to leading order.

    python scripts/protocol_coefficients.py /path/to/master_catalog.json

writes results/factory_protocols.json. The catalog file (54 MB) is not committed; download it into its own
directory first. Only the rows listed in SELECTED are kept. For each (n, k, d) the smallest A_d among pure-T rows
is used (pure T: gate family 'T' and a gate that is a sum of single-output terms).
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# (n, k, d): chosen to span the catalogue: the [[15,1,3]] reference, small and large distance-3 multi-output
# protocols, the lowest-A_d distance-3 protocol found in the n <= 130 range, and higher distances.
SELECTED = [(15, 1, 3), (28, 2, 3), (44, 4, 3), (63, 11, 3), (114, 14, 3), (127, 23, 3), (175, 17, 4),
            (429, 83, 4), (125, 3, 5), (208, 16, 5), (213, 11, 6), (222, 2, 8)]


def pure_t(row):
    terms = row['gate'].split(',')[0].split('+')
    return row.get('gate_family') == 'T' and all(t.isdigit() for t in terms) and len(set(terms)) == len(terms)


def main(path):
    rows = json.loads(Path(path).read_text())['factories']
    best = {}
    for r in rows:
        key = (int(r['n']), int(r['k']), int(r['d']))
        if key in SELECTED and pure_t(r) and r['A_d'] is not None:
            if key not in best or r['A_d'] < best[key]['A_d']:
                best[key] = r
    missing = [k for k in SELECTED if k not in best]
    if missing:
        raise SystemExit(f'not found in catalog: {missing}')
    out = [dict(name=f'[[{n},{k},{d}]]', n=n, k=k, d=d, wires=int(best[(n, k, d)]['N']), A_d=int(best[(n, k, d)]['A_d']),
                catalog_label=best[(n, k, d)].get('catalog_label')) for n, k, d in SELECTED]
    (ROOT/'results/factory_protocols.json').write_text(json.dumps(dict(
        source='arXiv:2610.06535 master_catalog.json (github.com/shubhamj810/magic-state-factory-repo)',
        protocols=out), indent=1) + '\n')
    for p in out:
        print(p['name'], 'wires', p['wires'], 'A_d', p['A_d'])


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '/tmp/msf_data/catalog.json')
