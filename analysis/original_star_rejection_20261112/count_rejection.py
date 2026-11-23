"""Original STAR: exact single-Pauli-fault rejection count, not Monte Carlo.

Akahoshi et al., arXiv:2303.13181, Sec. IV B. The implemented nearest-neighbor
schedule is specified in docs/original_star_rejection_count_derivation.md.
All rounds are noisy by default; ideal_final omits noise in the last round.
The finite-distance constants include our explicit SWAP routing of the seed.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
import csv
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location(
    'star_fault_signatures', HERE.parent/'rejection_count_20261111/count_rejection.py')
_engine = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _engine
_spec.loader.exec_module(_engine)
Op, fault_signatures, check_forward = _engine.Op, _engine.fault_signatures, _engine.check_forward


def make_circuit(d, xorder=(0, 1, 2, 3), zorder=(0, 2, 1, 3)):
    assert d >= 3 and d % 2
    # Data on even-even sites, checks at diagonal offsets (+/-1,+/-1).
    data_xy = [(2*x, 2*y) for y in range(d) for x in range(d)]
    checks = {}
    for y in range(-1, d):
        for x in range(-1, d):
            kind = 'X' if (x+y) % 2 == 0 else 'Z'
            bulk = 0 <= x < d-1 and 0 <= y < d-1
            boundary = ((x in (-1, d-1) and 0 <= y < d-1 and kind == 'Z') or
                        (y in (-1, d-1) and 0 <= x < d-1 and kind == 'X'))
            if bulk or boundary:
                checks[2*x+1, 2*y+1] = kind
    assert len(checks) == d*d-1
    # M0 is an extra boundary gauge ancilla, retired after seed verification.
    coords = data_xy + sorted(set(checks) | {(1, -1)})
    idx = {xy: i for i, xy in enumerate(coords)}
    data = [idx[xy] for xy in data_xy]
    anc = [idx[xy] for xy in checks]
    xa = [idx[xy] for xy, kind in checks.items() if kind == 'X']
    seed = [idx[0, 0], idx[0, 2], idx[2, 0], idx[2, 2]]
    m = [idx[1, -1], idx[-1, 1], idx[3, 1], idx[1, 3]]
    fresh = [i for i in data if i not in seed]
    plus = [i for i in fresh if coords[i][1] >= 4]
    offsets = [(-1, -1), (1, -1), (-1, 1), (1, 1)]
    support = {a: [idx[x+dx, y+dy] for dx, dy in offsets
                   if (x+dx, y+dy) in idx and idx[x+dx, y+dy] in data]
               for a in anc for x, y in [coords[a]]}
    assert all(len(s) in (2, 4) for s in support.values())
    # Validate commuting check operators and the seed's placement.
    for a in xa:
        for b in set(anc)-set(xa):
            assert len(set(support[a]) & set(support[b])) % 2 == 0

    raw = {}
    for r in range(2):
        for a in m:
            raw['gauge', r, a] = 1 << len(raw)
    for r in range(2):
        for a in anc:
            raw['growth', r, a] = 1 << len(raw)
    tests = []
    for r in range(2):
        tests.append(raw['gauge', r, m[0]] ^ raw['gauge', r, m[3]])
        tests.append(raw['gauge', r, m[1]] ^ raw['gauge', r, m[2]])

    # A GF(2) basis of known checks immediately before growth. The sign is a
    # parity of recorded gauge outcomes, not an assumption that each is +1.
    bases = {'X': {}, 'Z': {}}

    def insert(kind, wires, ref=0):
        vec = sum(1 << i for i in wires)
        while vec:
            p = vec.bit_length()-1
            if p not in bases[kind]:
                bases[kind][p] = vec, ref
                return
            v, s = bases[kind][p]
            vec ^= v
            ref ^= s
        assert ref == 0

    def known(kind, wires):
        vec, ref = sum(1 << i for i in wires), 0
        while vec:
            p = vec.bit_length()-1
            if p not in bases[kind]:
                return None
            v, s = bases[kind][p]
            vec ^= v
            ref ^= s
        return ref

    insert('X', seed)
    insert('Z', seed[:2], raw['gauge', 1, m[1]])
    insert('Z', seed[2:], raw['gauge', 1, m[2]])
    for i in fresh:
        insert('X' if i in plus else 'Z', [i])
    known_counts = Counter()
    for a in anc:
        kind = checks[coords[a]]
        ref = known(kind, support[a])
        if ref is not None:
            tests.append(raw['growth', 0, a] ^ ref)
            known_counts[kind] += 1
        tests.append(raw['growth', 1, a] ^ raw['growth', 0, a])
    records = defaultdict(int)
    for j, parity in enumerate(tests):
        while parity:
            bit = parity & -parity
            records[bit] ^= 1 << j
            parity ^= bit
    ops, tick = [], 0
    active = set(seed+m)

    def layer(stage, gates):
        nonlocal tick
        used = set()
        for kind, qs in gates:
            assert set(qs) <= active and not used.intersection(qs), (stage, tick, kind, qs)
            used.update(qs)
            rec = 0
            if kind == 'M':
                phase, r = stage.split('_')
                rec = records[raw[phase, int(r), qs[0]]]
            if kind == 'CX':
                x1, y1 = coords[qs[0]]
                x2, y2 = coords[qs[1]]
                assert abs(x1-x2) == abs(y1-y2) == 1
            ops.append(Op(tick, stage, kind, tuple(qs), rec))
        ops.extend(Op(tick, stage, 'I', (i,)) for i in sorted(active-used))
        tick += 1

    def swaps(stage, pairs):
        for rev in (False, True, False):
            layer(stage, [('CX', (b, a) if rev else (a, b)) for a, b in pairs])

    layer('encoding', [('R', (i,)) for i in sorted(active)])
    layer('encoding', [('H', (seed[i],)) for i in (1, 3)])
    pairs = [(seed[0], m[1]), (seed[2], m[2])]
    swaps('encoding', pairs)
    layer('encoding', [('CX', (seed[1], m[1])), ('CX', (seed[3], m[2]))])
    swaps('encoding', pairs)
    swaps('rotation', [(seed[0], m[0])])
    layer('rotation', [('CX', (m[0], seed[2]))])
    layer('rotation', [('RZ', (seed[2],))])
    layer('rotation', [('CX', (m[0], seed[2]))])
    swaps('rotation', [(seed[0], m[0])])

    for r in range(2):
        stage = f'gauge_{r}'
        layer(stage, [('R', (a,)) for a in m])
        layer(stage, [('H', (m[i],)) for i in (0, 3)])
        layer(stage, [('CX', (m[0], seed[0])), ('CX', (m[3], seed[1]))])
        layer(stage, [('CX', (m[0], seed[2])), ('CX', (m[3], seed[3]))])
        layer(stage, [('CX', (seed[0], m[1])), ('CX', (seed[2], m[2]))] +
                     [('H', (m[i],)) for i in (0, 3)])
        layer(stage, [('CX', (seed[1], m[1])), ('CX', (seed[3], m[2]))] +
                     [('M', (m[i],)) for i in (0, 3)])
        layer(stage, [('M', (m[i],)) for i in (1, 2)])

    active = set(data+anc)
    for r in range(2):
        stage = f'growth_{r}'
        layer(stage, [('R', (a,)) for a in anc+(fresh if r == 0 else [])])
        layer(stage, [('H', (a,)) for a in xa+(plus if r == 0 else [])])
        for j in range(4):
            gates = []
            for a in anc:
                x, y = coords[a]
                dx, dy = offsets[(xorder if a in xa else zorder)[j]]
                b = idx.get((x+dx, y+dy))
                if b in support[a]:
                    gates.append(('CX', (a, b) if a in xa else (b, a)))
            layer(stage, gates)
        layer(stage, [('H', (a,)) for a in xa])
        layer(stage, [('M', (a,)) for a in anc])
    return coords, ops, dict(detectors=len(tests), known_checks=dict(known_counts),
                            patch_qubits=2*d*d-1, ticks=tick)


def analyze(d, ideal_final=False, verify=False, dump=False, **schedule):
    coords, ops, meta = make_circuit(d, **schedule)
    locations, branches = fault_signatures(coords, ops)
    # The sole non-Clifford branch is logical Z: rejection is angle-independent.
    assert branches == [0], branches
    if verify:
        meta['forward_checks'] = check_forward(coords, ops, locations)
    groups = defaultdict(lambda: dict(locations=0, risky_locations=0,
                                     rejecting_paulis=0, coefficient=Fraction()))
    details = []
    for op, sigs in locations:
        rate = 'p' if op.kind in ('R', 'M', 'CX') else 'q'
        # Keep ideal-round gates in propagation, but omit their noise locations.
        if ideal_final and op.stage == 'growth_1':
            continue
        reject = sum(s != 0 for s in sigs)
        w = Fraction(reject, len(sigs))
        key = op.stage, op.kind, rate
        row = groups[key]
        row['locations'] += 1
        row['risky_locations'] += bool(reject)
        row['rejecting_paulis'] += reject
        row['coefficient'] += w
        details.append(dict(tick=op.tick, stage=op.stage, kind=op.kind, rate=rate,
                            coordinates=str([coords[i] for i in op.qubits]),
                            rejection_weight=str(w), rejecting_paulis=reject,
                            allowed_paulis=len(sigs)))
    totals = {rate: sum((v['coefficient'] for (*_, r), v in groups.items() if r == rate), Fraction())
              for rate in ('p', 'q')}
    if dump:
        with (HERE/f'fault_locations_d{d}.csv').open('w') as f:
            writer = csv.DictWriter(f, fieldnames=list(details[0]))
            writer.writeheader()
            writer.writerows(details)
    return dict(distance=d, ideal_final=ideal_final, **meta,
                coefficients={r: str(c) for r, c in totals.items()},
                groups=[dict(stage=s, kind=k, rate=r, **{a: str(b) if isinstance(b, Fraction) else b
                             for a, b in v.items()}) for (s, k, r), v in sorted(groups.items())])


def quadratic(rows, field):
    """Exact polynomial interpolation, then require every other distance to agree."""
    mat = [[Fraction(r['distance'])**2, Fraction(r['distance']), Fraction(1), field(r)]
           for r in rows[:3]]
    for j in range(3):
        pivot = mat[j][j]
        mat[j] = [v/pivot for v in mat[j]]
        for i in range(3):
            if i != j:
                a = mat[i][j]
                mat[i] = [x-a*y for x, y in zip(mat[i], mat[j])]
    abc = [mat[i][3] for i in range(3)]
    for r in rows:
        d = r['distance']
        assert abc[0]*d*d+abc[1]*d+abc[2] == field(r), (d, abc, field(r))
    return [str(x) for x in abc]


def main():
    output = dict(method='Deterministic single-fault enumeration; rational Pauli weights.', cases={})
    for ideal in (False, True):
        rows = [analyze(d, ideal, verify=d in (3, 5, 9), dump=d == 11 and not ideal)
                for d in (3, 5, 7, 9, 11, 13, 15, 17, 21, 25, 31)]
        for row in rows:
            d = row['distance']
            assert sum(row['known_checks'].values()) == (d*d-2*d+5)//2
            assert row['detectors'] == (3*d*d-2*d+11)//2
        # At d=3 the seed meets the opposite boundary; keep its direct count.
        polys = {r: quadratic(rows[1:], lambda row: Fraction(row['coefficients'][r])) for r in ('p', 'q')}
        breakdown = {}
        for kind in ('R', 'M', 'CX', 'H', 'RZ', 'I'):
            breakdown[kind] = quadratic(rows[1:], lambda row: sum(
                (Fraction(g['coefficient']) for g in row['groups'] if g['kind'] == kind), Fraction()))
        label = 'ideal_final' if ideal else 'all_noisy'
        output['cases'][label] = dict(polynomials=polys, polynomial_min_distance=5,
                                     by_operation=breakdown, rows=rows)
        print(label, polys, breakdown, flush=True)
    # An independent bulk unit-cell count determines the leading coefficients
    # directly, without polynomial interpolation in distance. This 4x4 data
    # window has 16 data sites and 16 checks, equally split between X and Z.
    coords, ops, _ = make_circuit(17)
    locations, _ = fault_signatures(coords, ops)
    bulk = defaultdict(Counter)
    for op, sigs in locations:
        anchor = next((i for i in op.qubits if coords[i][0] % 2), op.qubits[0])
        x, y = coords[anchor]
        if 8 <= x < 16 and 8 <= y < 16:
            bulk[op.stage, op.kind][sum(s != 0 for s in sigs), len(sigs)] += 1
    output['bulk_unit_cell'] = []
    for (stage, kind), hist in sorted(bulk.items()):
        coef = sum((Fraction(n, den)*count for (n, den), count in hist.items()), Fraction())/16
        output['bulk_unit_cell'].append(dict(stage=stage, kind=kind, coefficient=str(coef),
            fractions=[dict(rejecting=n, allowed=den, locations_per_d_squared=str(Fraction(count, 16)))
                       for (n, den), count in sorted(hist.items())]))
    for label, case in output['cases'].items():
        for rate in ('p', 'q'):
            total = sum((Fraction(row['coefficient']) for row in output['bulk_unit_cell']
                         if (row['kind'] in ('R', 'M', 'CX')) == (rate == 'p') and
                         (label != 'ideal_final' or row['stage'] != 'growth_1')), Fraction())
            assert total == Fraction(case['polynomials'][rate][0])
    output['reflected_hook_schedules'] = []
    for xo, zo in (((0, 1, 2, 3), (0, 2, 1, 3)),
                   ((1, 0, 3, 2), (1, 3, 0, 2)),
                   ((2, 3, 0, 1), (2, 0, 3, 1)),
                   ((3, 2, 1, 0), (3, 1, 2, 0))):
        rows = [analyze(d, xorder=xo, zorder=zo) for d in (5, 7, 9, 11)]
        polys = {r: quadratic(rows, lambda row: Fraction(row['coefficients'][r])) for r in ('p', 'q')}
        assert polys['p'][0] == '179/15' and polys['q'][0] == '22/3'
        output['reflected_hook_schedules'].append(dict(xorder=xo, zorder=zo, polynomials=polys))
    (HERE/'counts.json').write_text(json.dumps(output, indent=2)+'\n')


if __name__ == '__main__':
    main()
