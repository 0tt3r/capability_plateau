"""Single-fault rejection enumeration for a specified Toshio STAR circuit.

Source: Phys. Rev. X 15, 021057 (2025), Algorithm 1, Figs. 5(c), 8.
The source does not specify a syndrome CNOT schedule or idle timing. Here every
native operation takes one tick; all checks use N,E,W,S CNOT order. SWAP is
compiled into three CNOTs by default. Initial data and first-round ancilla
preparations are parallel, as are their initial Hadamards. All unused sites idle.
No growth, delivery, or later error-correction rounds belong to this acceptance
test. Reset/readout flip errors have rate p; CNOT depolarization has rate p;
H/RZ/idle depolarization has rate q, as in the repository hardware definition.

We propagate the signatures of X,Z faults backwards through the Clifford
skeleton, retaining the syndrome signatures of the k physical RZ branches.
A fault outside their span rejects for every rotation branch. A zero signature
leaves acceptance unchanged. A nonzero signature in that span changes which
rotation branches pass. This keeps intrinsic rotation rejection separate.
"""

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import permutations, product
from pathlib import Path
import argparse
import csv
import json


@dataclass(frozen=True)
class Op:
    tick: int
    stage: str
    kind: str
    qubits: tuple
    record: int = 0


def make_circuit(d, native_swap=False, order="NEWS"):
    assert d >= 4 and d % 2 == 0
    side = 2*d-1
    coords = [(x,y) for y in range(side) for x in range(side)]
    index = {xy:i for i,xy in enumerate(coords)}
    data = [i for i,(x,y) in enumerate(coords) if (x+y)%2 == 0]
    anc = [i for i,(x,y) in enumerate(coords) if (x+y)%2 == 1]
    xa = [i for i in anc if coords[i][0]%2]
    za = [i for i in anc if i not in xa]
    band_x = [i for i in xa if coords[i][1] <= 2]
    band_z = [i for i in za if coords[i][1] == 1]
    detector = {}
    for r in range(3):
        for a in band_x:
            detector[r,a] = 1 << len(detector)
    for r in (1,2):
        for a in band_z:
            detector[r,a] = 1 << len(detector)
    records = dict(detector)
    for a in band_z:
        records[0,a] = detector[1,a] ^ detector[2,a]
    ops = []
    tick = 0

    def layer(stage, gates):
        nonlocal tick
        used = set()
        for kind,qubits in gates:
            assert not used.intersection(qubits), (tick,kind,qubits)
            used.update(qubits)
            rec = records.get((int(stage[-1]),qubits[0]),0) if kind == 'M' else 0
            ops.append(Op(tick,stage,kind,tuple(qubits),rec))
        ops.extend(Op(tick,stage,'I',(i,)) for i in range(len(coords)) if i not in used)
        tick += 1

    def syndrome(r):
        stage = f'syndrome_{r}'
        layer(stage, [('R',(i,)) for i in (list(range(len(coords))) if r==0 else anc)])
        layer(stage, [('H',(i,)) for i in (data+xa if r==0 else xa)])
        for direction in order:
            dx,dy = {'N':(0,-1),'E':(1,0),'W':(-1,0),'S':(0,1)}[direction]
            gates = []
            for a in anc:
                x,y = coords[a]
                b = index.get((x+dx,y+dy))
                if b is not None:
                    gates.append(('CX',(a,b) if a in xa else (b,a)))
            layer(stage,gates)
        layer(stage,[('H',(i,)) for i in xa])
        layer(stage,[('M',(i,)) for i in anc])

    syndrome(0)
    triples = [(index[4*j,0],index[4*j+1,0],index[4*j+2,0]) for j in range(d//2)]

    def swap_layer():
        if native_swap:
            layer('rotation',[('SWAP',(a,b)) for a,b,c in triples])
        else:
            for reverse in (False,True,False):
                layer('rotation',[('CX',(b,a) if reverse else (a,b)) for a,b,c in triples])

    swap_layer()
    layer('rotation',[('CX',(b,c)) for a,b,c in triples])
    layer('rotation',[('RZ',(c,)) for a,b,c in triples])
    layer('rotation',[('CX',(b,c)) for a,b,c in triples])
    swap_layer()
    syndrome(1)
    syndrome(2)
    return coords, ops, len(detector)


def fault_signatures(coords, ops):
    xs = [0]*len(coords)
    zs = [0]*len(coords)
    locations = []
    branches = []
    for op in reversed(ops):
        kind = op.kind
        if kind == 'M':
            signatures = [op.record]
        elif kind == 'R':
            signatures = [xs[op.qubits[0]]]
        elif len(op.qubits) == 1:
            i = op.qubits[0]
            signatures = [xs[i],xs[i]^zs[i],zs[i]]
        else:
            i,j = op.qubits
            signatures = [a^b for n,(a,b) in enumerate(product(
                (0,xs[i],xs[i]^zs[i],zs[i]),
                (0,xs[j],xs[j]^zs[j],zs[j]))) if n]
        locations.append((op,signatures))
        if kind == 'RZ':
            branches.append(zs[op.qubits[0]])
        if kind == 'H':
            i = op.qubits[0]
            xs[i],zs[i] = zs[i],xs[i]
        elif kind == 'CX':
            i,j = op.qubits
            xs[i] ^= xs[j]
            zs[j] ^= zs[i]
        elif kind == 'SWAP':
            i,j = op.qubits
            xs[i],xs[j] = xs[j],xs[i]
            zs[i],zs[j] = zs[j],zs[i]
        elif kind == 'M':
            xs[op.qubits[0]] ^= op.record
        elif kind == 'R':
            i = op.qubits[0]
            # Every tested parity must be deterministic at zero rotation and
            # zero noise. Nonzero Z-fault signature at |0> reset disproves it.
            assert zs[i] == 0, ('nondeterministic detector',op,zs[i])
            xs[i] = zs[i] = 0
    assert not any(xs) and not any(zs)
    return list(reversed(locations)), list(reversed(branches))


def check_forward(coords, ops, backwards):
    """Independent forward propagation of all primitive fault frames at once.

    A bit in a wire frame labels an independent fault experiment. Compare the
    resulting detector bits against every signature from the backward pass.
    """
    xf = [0]*len(coords)
    zf = [0]*len(coords)
    detector_frames = defaultdict(int)
    expected = []

    def record(mask, frame):
        while mask:
            bit = mask & -mask
            detector_frames[bit] ^= frame
            mask ^= bit

    for op,(same,sigs) in zip(ops,backwards):
        assert op == same
        kind = op.kind
        i = op.qubits[0]
        if kind == 'H':
            xf[i],zf[i] = zf[i],xf[i]
        elif kind == 'CX':
            j = op.qubits[1]
            xf[j] ^= xf[i]
            zf[i] ^= zf[j]
        elif kind == 'SWAP':
            j = op.qubits[1]
            xf[i],xf[j] = xf[j],xf[i]
            zf[i],zf[j] = zf[j],zf[i]
        elif kind == 'M':
            record(op.record,xf[i])
            zf[i] = 0
        elif kind == 'R':
            xf[i] = zf[i] = 0
        if kind == 'M':
            record(op.record,1 << len(expected))
            expected.append(sigs[0])
        elif kind == 'R':
            xf[i] ^= 1 << len(expected)
            expected.append(sigs[0])
        else:
            for j,q in enumerate(op.qubits):
                if len(op.qubits) == 1:
                    sx,sz = sigs[0],sigs[2]
                elif j == 0:
                    sx,sz = sigs[3],sigs[11]  # XI, ZI, after excluding II
                else:
                    sx,sz = sigs[0],sigs[2]   # IX, IZ
                xf[q] ^= 1 << len(expected)
                expected.append(sx)
                zf[q] ^= 1 << len(expected)
                expected.append(sz)
    for j,want in enumerate(expected):
        got = sum(bit for bit,frame in detector_frames.items() if (frame >> j)&1)
        assert got == want, ('forward/backward mismatch',j,got,want)
    return len(expected)


def span_basis(branches):
    basis = {}
    for j,sig in enumerate(branches):
        bits = 1 << j
        while sig:
            pivot = sig.bit_length()-1
            if pivot in basis:
                s,b = basis[pivot]
                sig ^= s
                bits ^= b
            else:
                basis[pivot] = sig,bits
                break
        if not sig:
            assert bits == (1<<len(branches))-1
    assert len(basis) == len(branches)-1
    return basis


def branch_class(sig, basis, k):
    bits = 0
    while sig:
        pivot = sig.bit_length()-1
        if pivot not in basis:
            return 'reject'
        s,b = basis[pivot]
        sig ^= s
        bits ^= b
    weight = bits.bit_count()
    return str(min(weight,k-weight))


def analyze(d, native_swap=False, order='NEWS', dump=False):
    coords,ops,nd = make_circuit(d,native_swap,order)
    locs,branches = fault_signatures(coords,ops)
    basis = span_basis(branches)
    groups = defaultdict(lambda: dict(locations=0, risky_locations=0,
        certainly_rejecting_locations=0, fault_counts=Counter(), weights=defaultdict(Fraction)))
    details = []
    for op,sigs in locs:
        cls = Counter(branch_class(s,basis,d//2) for s in sigs)
        rate = 'p' if op.kind in ('CX','SWAP','R','M') else 'q'
        category = op.kind
        if op.kind in ('R','H') and op.tick <= 1 and sum(coords[op.qubits[0]])%2==0:
            category += '_data'
        key = op.stage,category,rate
        rec = groups[key]
        rec['locations'] += 1
        rec['risky_locations'] += cls['0'] < len(sigs)
        rec['certainly_rejecting_locations'] += cls['reject'] == len(sigs)
        rec['fault_counts'].update(cls)
        for c,n in cls.items():
            rec['weights'][c] += Fraction(n,len(sigs))
        if dump and cls['0'] < len(sigs):
            labels = (['flip'] if op.kind in ('M','R') else
                      list('XYZ') if len(op.qubits)==1 else
                      [a+b for a,b in product('IXYZ',repeat=2) if a+b!='II'])
            types = {label:branch_class(s,basis,d//2) for label,s in zip(labels,sigs)}
            details.append(dict(tick=op.tick,stage=op.stage,kind=category,rate=rate,
                coords=str([coords[i] for i in op.qubits]),
                total_fault_types=len(sigs),reject_fault_types=cls['reject'],
                unchanged_fault_types=cls['0'],
                modified_fault_types=sum(v for k,v in cls.items() if k not in ('reject','0')),
                reject_paulis=' '.join(p for p,c in types.items() if c=='reject'),
                modified_paulis=' '.join(p for p,c in types.items() if c not in ('reject','0'))))
    rows = []
    total = {r:defaultdict(Fraction) for r in ('p','q')}
    for (stage,kind,rate),rec in sorted(groups.items()):
        for cls,w in rec['weights'].items():
            total[rate][cls] += w
        rows.append(dict(distance=d,stage=stage,kind=kind,rate=rate,
            locations=rec['locations'],risky_locations=rec['risky_locations'],
            certainly_rejecting_locations=rec['certainly_rejecting_locations'],
            reject_weight=str(rec['weights']['reject']),
            modified_weight=str(sum(w for c,w in rec['weights'].items() if c not in ('reject','0'))),
            fault_counts=dict(rec['fault_counts'])))
    return dict(distance=d,native_swap=native_swap,cnot_order=order,
                ticks=ops[-1].tick+1,detectors=nd,
                totals={r:{c:str(w) for c,w in rec.items() if c!='0'} for r,rec in total.items()},
                groups=rows),details


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    checked = 0
    for native in (False,True):
        coords,ops,_ = make_circuit(6,native)
        locs,branches = fault_signatures(coords,ops)
        checked += check_forward(coords,ops,locs)
        # All branch patterns have only the global-complement degeneracy.
        counts = Counter()
        for mask in range(1<<len(branches)):
            sig = 0
            for j,s in enumerate(branches):
                if (mask>>j)&1:
                    sig ^= s
            counts[sig] += 1
        assert set(counts.values()) == {2}
    print('Forward/backward fault checks:',checked,flush=True)
    runs = []
    for native_swap in (False,True):
        for order in ('NEWS','NWES'):
            for d in (6,10,14,18,22):
                result,details = analyze(d,native_swap,order,dump=d==10 and order=='NEWS' and not native_swap)
                assert Fraction(result['totals']['p']['1']) == Fraction(d,15)
                assert Fraction(result['totals']['q']['1']) == Fraction(d,6)
                assert set(result['totals']['p']) == {'reject','1'}
                assert set(result['totals']['q']) == {'reject','1'}
                runs.append(result)
                print(d,native_swap,order,result['totals'],flush=True)
                if details:
                    with (args.out/'fault_locations_d10.csv').open('w') as f:
                        writer=csv.DictWriter(f,fieldnames=list(details[0]))
                        writer.writeheader()
                        writer.writerows(details)
    (args.out/'counts.json').write_text(json.dumps(runs,indent=2)+'\n')
    summary = []
    for native_swap in (False,True):
        rr = [r for r in runs if r['native_swap']==native_swap and r['cnot_order']=='NEWS']
        for rate in ('p','q'):
            values = [Fraction(r['totals'][rate]['reject']) for r in rr]
            slope = (values[1]-values[0])/(rr[1]['distance']-rr[0]['distance'])
            intercept = values[0]-slope*rr[0]['distance']
            assert all(v == slope*r['distance']+intercept for r,v in zip(rr,values))
            summary.append(dict(native_swap=native_swap,rate=rate,
                                reject_slope=str(slope),reject_intercept=str(intercept)))
    sensitivity = []
    for order in permutations('NEWS'):
        order = ''.join(order)
        try:
            a,_ = analyze(6,order=order)
            b,_ = analyze(10,order=order)
        except AssertionError:
            # Reject orders that collide or fail deterministic-check validation.
            continue
        sensitivity.append(dict(order=order,**{
            r+'_slope':str((Fraction(b['totals'][r]['reject'])-
                           Fraction(a['totals'][r]['reject']))/4) for r in ('p','q')}))
    (args.out/'summary.json').write_text(json.dumps(dict(
        forward_backward_primitive_fault_checks=checked,
        linear_coefficients=summary,valid_cnot_order_sensitivity=sensitivity),indent=2)+'\n')


if __name__ == '__main__':
    main()
