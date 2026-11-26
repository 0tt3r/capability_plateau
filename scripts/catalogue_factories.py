"""Surface-code layouts of the distance-2 and distance-3 factory catalogues, costed with the JWS layout rules.

    .venv/bin/python scripts/catalogue_factories.py --master MASTER_CATALOG_JSON --search SEARCH_OUTPUTS_DIR

Sources (read from the authors' repositories; paths are arguments):
  distance 3  master_catalog.json of the arXiv:2610.06535 repository (level-3 rows with n <= 64, which include the
              generalised triorthogonal protocols of arXiv:2609.30860 with n <= 54); A_3 is counted from its columns.
  distance 2  the Magic_state_factory_search repository (Singh, Gidney, Jones): outputs/l3_catalogue_k6.csv (one verified
              circuit per level-3 class with k <= 6) and outputs/large_unified_l3_circuits.txt.zip joined to
              outputs/large_unified_factories_l3.csv (k <= 12); A_2 is the number of input-fault pairs that pass every
              check and act on the outputs, which reproduces the catalogue's error coefficient c (28 for 8T -> CCZ).
Outputs occupy wires 0..k-1 and checks wires k..N-1 in both catalogues.

T-state yield per run [choice]: only outputs with a known conversion into T states are T-state supplies.
  degree-1 outputs (Clifford-equivalent to T states)  yield = exact T count of the class
  products of CCZ on disjoint qubits                   yield = 2 per CCZ (catalysed CCZ -> 3 T, one reused)
  everything else (CS terms, overlapping products)    not a T-state supply; skipped
Each protocol is costed exactly as the 'jws' layouts of capability/surface_factories.py: faulty-T inputs, the
cultivated inputs of the Pinnacle engines (d_a = 7, 9), and two levels fed by faulty-T 15-to-1 factories, over the same
d_X, d_Z grid (d_Z <= d_X), dropping designs with rejection >= 0.99 or error >= 1e-3, with the error and period divided
by the yield instead of k. The designs not dominated in (qubits, cycles per T state, error) are written to
results/catalogue_factories.json.
"""
import argparse
import csv
import io
import itertools
import json
import re
import sys
import zipfile
from math import ceil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from capability import surface_factories as sf
from capability.config import P

N_MAX_D3 = 64


def terms(gate):
    """[(kind, qubits)] of a gate written as T0·CS01·CCZ012, T0,CS01, or \\mathrm{CS}_{12}."""
    g = gate.replace('\\mathrm{', '').replace('}_{', '').replace('}', '')
    return [(kind, tuple(int(c) for c in qs)) for kind, qs in re.findall(r'(CCZ|CS|T)(\d+)', g.replace('·', ' '))]


YIELD_RULE = 'conversion'     # 'tcount' credits every output with its T count: an upper bound used only as a check


def yield_of(gate, t_count, degree):
    """T states delivered per accepted run, or None if the output has no conversion into T states."""
    if YIELD_RULE == 'tcount':
        return int(t_count) if t_count else None
    ts = terms(gate)
    if degree == 1 or (ts and all(kind == 'T' for kind, _ in ts)):
        return t_count if t_count else len(ts)
    if ts and all(kind == 'CCZ' for kind, _ in ts):
        used = [q for _, qs in ts for q in qs]
        if len(used) == len(set(used)):
            return 2*len(ts)
    return None


def undetected(supports, k, N, d):
    """Weight-d input-fault patterns with zero syndrome on the checks that act on the outputs."""
    checks, outs = set(range(k, N)), set(range(k))
    count = 0
    for subset in itertools.combinations(range(len(supports)), d):
        acc = frozenset()
        for i in subset:
            acc ^= supports[i]
        if not acc & checks and acc & outs:
            count += 1
    return count


def distance3(path):
    rows = json.loads(Path(path).read_text())['factories']
    for f in rows:
        if f['d'] != 3 or f.get('level') != 3 or f['n'] > N_MAX_D3:
            continue
        y = yield_of(f['gate_human'], f.get('t_count'), f.get('poly_degree'))
        if not y:
            continue
        supports = [frozenset(c) for c in f['columns']]
        yield dict(source='master', n=f['n'], k=f['k'], d=3, wires=f['N'], gate=f['gate_human'], yield_=y,
                   A_d=undetected(supports, f['k'], f['N'], 3))


def distance2(path):
    path = Path(path)
    for r in csv.DictReader(open(path/'l3_catalogue_k6.csv')):
        k, N = int(r['k']), int(r['N'])
        y = yield_of(r['spec'], int(r['tcount']), None)
        if not y:
            continue
        supports = [frozenset(q) for q, c in json.loads(r['gates']) if c % 2]
        A = undetected(supports, k, N, 2)
        if A:
            yield dict(source='l3_catalogue_k6', n=int(r['n']), k=k, d=2, wires=N, gate=r['spec'], yield_=y, A_d=A)
    table = {(r['k'], r['n'], r['N'], r['search_stage'], r['params']): r
             for r in csv.DictReader(open(path/'large_unified_factories_l3.csv'))}
    text = zipfile.ZipFile(path/'large_unified_l3_circuits.txt.zip').read('large_unified_l3_circuits.txt').decode()
    for block in text.split('# ')[1:]:
        lines = block.strip().split('\n')
        m = re.match(r'\d+ l=3 k=(\d+) n=(\d+) N=(\d+) (\S+) (\S+)', lines[0])
        row = table.get(m.groups()) if m else None
        if not row or row['d'] != '2':
            continue
        k, N = int(row['k']), int(row['N'])
        y = yield_of(row['reduced'], int(row['zcount']) if row.get('zcount') else None, None)
        if not y:
            continue
        supports = []
        for line in lines[1:]:
            qs, c = line.split(':')
            if int(c) % 2:
                supports.append(frozenset(int(x) for x in qs.split(',')))
        A = undetected(supports, k, N, 2)
        if A:
            yield dict(source='large_unified_l3', n=int(row['n']), k=k, d=2, wires=N, gate=row['reduced'],
                       yield_=y, A_d=A)


def designs(proto, p=P):
    """One-level (faulty T, cultivated) and two-level designs of one protocol, per T state delivered."""
    n, k, d, A, y = proto['n'], proto['k'], proto['d'], proto['A_d'], proto['yield_']
    m = proto['wires'] - k
    ref = next(q for q in sf.PROTOCOLS if q['name'] == '[[15,1,3]]')
    inputs = [('faulty T', p, 0., 0)] + [(f'cultivated d_a={da}', pin + (da + 1)*pa, sf.RHO_CULT,
                                         30*(da*da + da - 1) + sf.CULT_POOL) for da, (pa, pin) in sf.CULTIVATED.items()]
    for dX in sf.DX:
        for dZ in sf.DZ:
            if dZ > dX:
                continue
            l = k*dX + m*dZ
            t0 = (ceil((n - m)/2) + (1 if k > 1 else 0))*dZ
            for label, p_rot, rho, extra in inputs + [(sf.FOLD_LABEL, *sf.fold_input(n, t0, p))]:
                reject = min(0.99, n*(rho + p_rot) + t0*m*(dX/(2*dZ))*sf.p_sc(dZ, p))
                t = t0/(1 - reject)
                error = A*p_rot**d/y + t*sf.p_sc(dX, p)*(1 + l/(2*dX))
                if reject < 0.99 and error < 1e-3:
                    yield label, (dX, dZ), 2*l*3*dX + 4*dZ + extra, t/y, error
            for dX1, dZ1 in sf.LEVEL1:
                _, t1, rej1, e1 = sf.single(ref, dX1, dZ1, dZ1, p)
                for nL1 in sf.N_L1:
                    q = (2*l*3*dX + 2*nL1*((dX1 + 4*dZ1)*(3*dX1 + dZ/2) + 2*dZ1)
                         + 2*(20*dZ**2 + 2*dX*dZ)*l/(dX + 4*dZ))
                    t = (n/2)*max(dZ, 6*dZ1/(1 - rej1)/(nL1/2))
                    error = sf.TWO_LEVEL_SCALE*(A*e1**d/y + t*sf.p_sc(dX, p)*(1 + l/(2*dX)))
                    if error < 1e-3:
                        yield f'two levels, (15-to-1)^{nL1}_{{{dX1},{dZ1},{dZ1}}}', (dX, dZ), int(round(q)), t/y, error


def pareto(rows):
    rows = sorted(rows, key=lambda r: (r['qubits'], r['cycles_per_T'], r['error']))
    kept = []
    for r in rows:
        if not any(g['qubits'] <= r['qubits'] and g['cycles_per_T'] <= r['cycles_per_T'] and g['error'] <= r['error']
                   for g in kept):
            kept.append(r)
    return kept


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--master', required=True, help='master_catalog/master_catalog.json of the arXiv:2610.06535 repository')
    ap.add_argument('--search', required=True, help='outputs/ directory of the Magic_state_factory_search repository')
    ap.add_argument('--yield-rule', choices=('conversion', 'tcount'), default='conversion',
                    help="'tcount': credit every output with its T count (an upper bound, for the sensitivity check only)")
    ap.add_argument('--out', default=str(ROOT/'results/catalogue_factories.json'))
    args = ap.parse_args()
    global YIELD_RULE
    YIELD_RULE = args.yield_rule
    protocols = list(distance3(args.master)) + list(distance2(args.search))
    out = []
    for proto in protocols:
        rows = [dict(name=f"[[{proto['n']},{proto['k']},{proto['d']}]] {proto['gate']}, {label}, d=({dX},{dZ},{dZ})",
                     source=proto['source'], n=proto['n'], k=proto['k'], d=proto['d'], wires=proto['wires'],
                     gate=proto['gate'], T_per_run=proto['yield_'], A_d=proto['A_d'], inputs=label, d_X=dX, d_Z=dZ,
                     qubits=q, cycles_per_T=c, error=e)
                for label, (dX, dZ), q, c, e in designs(proto)]
        out += pareto([r for r in rows if r['inputs'] != sf.FOLD_LABEL])
        out += pareto([r for r in rows if r['inputs'] == sf.FOLD_LABEL])
    # fold-transversal inputs need all-to-all hardware, so the 2D-grid designs are kept on their own Pareto front
    kept = (pareto([r for r in out if r['inputs'] != sf.FOLD_LABEL])
            + pareto([r for r in out if r['inputs'] == sf.FOLD_LABEL]))
    summary = dict(protocols=len(protocols), distance2=sum(p['d'] == 2 for p in protocols),
                   distance3=sum(p['d'] == 3 for p in protocols), designs=len(kept))
    Path(args.out).write_text(json.dumps(dict(yield_rule=YIELD_RULE, 
        sources=dict(master='magic-state-factory-search: master_catalog/master_catalog.json (arXiv:2610.06535)',
                     search='Magic_state_factory_search: outputs/l3_catalogue_k6.csv, large_unified_factories_l3.csv, '
                            'large_unified_l3_circuits.txt.zip'),
        summary=summary, designs=kept,
        protocols=[dict(n=q['n'], k=q['k'], d=q['d'], wires=q['wires'], gate=q['gate'], T_per_run=q['yield_'], A_d=q['A_d'],
                        source=q['source']) for q in protocols]), indent=1) + '\n')
    print(summary)


if __name__ == '__main__':
    main()
