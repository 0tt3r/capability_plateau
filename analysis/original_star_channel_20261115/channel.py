"""Original STAR: retained accepted single-fault logical channel.

Extends the existing declared circuit with output stabilizer and logical
observables. No production module or historical result is edited.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from pathlib import Path
from functools import lru_cache
import importlib.util
import json
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('original_count', HERE.parent/'original_star_rejection_20261112/count_rejection.py')
counter = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = counter
spec.loader.exec_module(counter)


def extended(d):
    coords, ops, meta = counter.make_circuit(d)
    nd = meta['detectors']
    index = {xy:i for i,xy in enumerate(coords)}
    data = {i for i,(x,y) in enumerate(coords) if not x%2 and not y%2}
    checks = []
    for i,(x,y) in enumerate(coords):
        if i in data or (x,y)==(1,-1):
            continue
        support = [index[x+dx,y+dy] for dx,dy in ((-1,-1),(1,-1),(-1,1),(1,1))
                   if (x+dx,y+dy) in index and index[x+dx,y+dy] in data]
        kind = 'X' if ((x-1)//2+(y-1)//2)%2==0 else 'Z'
        checks.append((i,kind,support))
    assert len(checks)==d*d-1
    lxbit, lzbit, rotbit = (1<<(nd+len(checks)+j) for j in range(3))
    xs,zs = [0]*len(coords),[0]*len(coords)
    final_bits = {}
    for j,(a,kind,support) in enumerate(checks):
        bit=1<<(nd+j);final_bits[a]=bit
        for i in support:
            (zs if kind=='X' else xs)[i] ^= bit
    for i,(x,y) in enumerate(coords):
        if i in data:
            if x==0: zs[i] ^= lxbit
            if y==0: xs[i] ^= lzbit
    locations=[]
    for op in reversed(ops):
        kind=op.kind;i=op.qubits[0]
        rec=op.record^(final_bits.get(i,0) if kind=='M' and op.stage=='growth_1' else 0)
        if kind=='M': sigs=[rec]
        elif kind=='R': sigs=[xs[i]]
        elif len(op.qubits)==1: sigs=[xs[i],xs[i]^zs[i],zs[i]]
        else:
            j=op.qubits[1]
            sigs=[a^b for n,(a,b) in enumerate(product((0,xs[i],xs[i]^zs[i],zs[i]),(0,xs[j],xs[j]^zs[j],zs[j]))) if n]
        locations.append((op,sigs))
        if kind=='RZ': xs[i] ^= rotbit
        elif kind=='H': xs[i],zs[i]=zs[i],xs[i]
        elif kind=='CX':
            j=op.qubits[1];xs[i]^=xs[j];zs[j]^=zs[i]
        elif kind=='M': xs[i]^=rec
        elif kind=='R': xs[i]=zs[i]=0
    locations.reverse()
    old,_=counter.fault_signatures(coords,ops)
    detector_mask=(1<<nd)-1
    for (_,a),(_,b) in zip(locations,old):
        assert [s&detector_mask for s in a]==b
    outputmask=((1<<len(checks))-1)<<nd
    classified=[]
    for op,sigs in locations:
        kinds=[]
        for s in sigs:
            if s&detector_mask: cls='reject'
            elif s&outputmask: cls='residual_syndrome_idealized'
            else:
                z=bool(s&lxbit); flip=bool(s&lzbit)^bool(s&rotbit)
                cls=('Z' if z else 'I')+('_angle_flip' if flip else '')
            kinds.append(cls)
        classified.append((op,sigs,kinds))
    return coords,ops,meta,classified


def counts(d):
    _,_,meta,locs=extended(d)
    weights=defaultdict(Fraction);groups=defaultdict(Fraction);examples=[]
    for op,sigs,classes in locs:
        rate='p' if op.kind in ('R','M','CX') else 'q'
        labels=['flip'] if op.kind in ('R','M') else list('XYZ') if len(op.qubits)==1 else [a+b for a,b in product('IXYZ',repeat=2) if a+b!='II']
        for label,cls in zip(labels,classes):
            weights[rate,cls]+=Fraction(1,len(classes))
            groups[op.stage,op.kind,rate,cls]+=Fraction(1,len(classes))
            if cls not in ('reject','I','residual_syndrome_idealized'):
                examples.append(dict(tick=op.tick,stage=op.stage,kind=op.kind,qubits=op.qubits,pauli=label,channel=cls))
    return dict(distance=d,**meta,weights={r:{c:str(weights[r,c]) for rr,c in weights if rr==r} for r in ('p','q')},
                groups=[dict(stage=s,kind=k,rate=r,channel=c,weight=str(w)) for (s,k,r,c),w in sorted(groups.items())],examples=examples)


def basis(signatures):
    rows={}
    for s in signatures:
        while s:
            j=s.bit_length()-1
            if j not in rows:
                rows[j]=s;break
            s^=rows[j]
    return tuple(rows[j] for j in sorted(rows))


def acceptance(d, samples=262144, seed=15112026, force_sample=False):
    """Positive Fourier-integrand sampling; sums every physical fault order.

    A=E_x prod_l E_fault[(-1)^(x.signature)]. The dual bits x are uniform.
    This avoids rare-event rejection sampling at large distance. Categorical
    depolarization gives 1 or 1-16p/15 (CX), 1-4q/3 (1q), 1-2p (flip).
    Distance 3 uses exhaustive Fourier summation, not sampling.
    """
    _,_,meta,locs=extended(d)
    nd=meta['detectors'];mask=(1<<nd)-1
    grouped=Counter()
    for op,sigs,_ in locs:
        kind='cx' if op.kind=='CX' else 'flip' if op.kind in ('R','M') else 'one'
        bs=basis([s&mask for s in sigs])
        if bs: grouped[kind,bs]+=1
    exact=d==3 and not force_sample
    if exact:
        samples=1<<nd
        assignments=np.arange(samples,dtype=np.uint64)
        dual=np.array([np.packbits(((assignments>>j)&1).astype(np.uint8),bitorder='little').view(np.uint64) for j in range(nd)])
    else:
        assert samples%64==0
        rng=np.random.default_rng(seed+d)
        dual=rng.integers(0,np.iinfo(np.uint64).max,size=(nd,samples//64),dtype=np.uint64,endpoint=True)
    @lru_cache(maxsize=None)
    def parity(sig):
        out=np.zeros(samples//64,dtype=np.uint64)
        while sig:
            bit=sig&-sig;out^=dual[bit.bit_length()-1];sig^=bit
        return out
    exposures={kind:np.zeros(samples,dtype=np.int32) for kind in ('cx','flip','one')}
    for (kind,bs),multiplicity in grouped.items():
        active=np.zeros(samples//64,dtype=np.uint64)
        for s in bs: active|=parity(s)
        exposures[kind]+=multiplicity*np.unpackbits(active.view(np.uint8),bitorder='little').astype(np.int32)
    result=dict(distance=d,samples=samples,seed=seed,exact=exact,detectors=nd,characteristic_groups=len(grouped),hardware={})
    for hw,p,q in [('H1',.001,.001),('H3',.001,.000001)]:
        log_phi=0.; retained=[]; cp=cq=0.
        for op,sigs,classes in locs:
            rate=p if op.kind in ('CX','R','M') else q
            rejecting=sum(bool(s&mask) for s in sigs)/len(sigs)
            z=1-rate*rejecting
            log_phi+=np.log(z)
            if op.kind in ('CX','R','M'): cp+=rejecting
            else: cq+=rejecting
            w=rate*classes.count('Z')/len(classes)
            assert not any('angle_flip' in c for c in classes)
            if w: retained.append(w/z)
        retained=np.array(retained)
        # Exactly one known logical-Z fault, any zero-signature background,
        # and no other nonzero rejection signature. Additional effects ideal.
        local_one=np.prod(1-retained)*np.sum(retained/(1-retained))
        logs=(exposures['cx']*np.log1p(-16*p/15)+exposures['one']*np.log1p(-4*q/3)
              +exposures['flip']*np.log1p(-2*p))
        ratios=np.exp(logs-log_phi)
        ratio=float(ratios.mean());se=0. if exact else float(ratios.std(ddof=1)/np.sqrt(samples))
        assert ratio+5*se>=1
        A=float(np.exp(log_phi)*ratio)
        result['hardware'][hw]=dict(p=p,q=q,acceptance=A,acceptance_standard_error=float(np.exp(log_phi)*se),
            zero_signature_probability=float(np.exp(log_phi)),acceptance_over_zero_signature=ratio,
            ratio_standard_error=se,retained_local_one=float(local_one),
            accepted_known_single_fault_weight=float(np.exp(log_phi)*local_one),
            conditional_error=float(local_one/ratio),conditional_error_leading=4*p/15+q/3,
            previous_product_acceptance=float(np.exp(cp*np.log1p(-p)+cq*np.log1p(-q))))
    return result


if __name__=='__main__':
    import argparse
    from time import monotonic
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples',type=int,default=262144)
    parser.add_argument('--seed',type=int,default=15112026)
    args=parser.parse_args()
    rows=[]; yields=[]
    for d in range(3,32,2):
        start=monotonic()
        row=counts(d);rows.append(row)
        val=acceptance(d,args.samples,args.seed);yields.append(val)
        print(d,row['weights'],{hw:(v['acceptance'],v['conditional_error'],v['ratio_standard_error']) for hw,v in val['hardware'].items()},'seconds',round(monotonic()-start,2),flush=True)
        (HERE/'single_fault_counts.json').write_text(json.dumps(dict(rows=rows),indent=2)+'\n')
        (HERE/'acceptance.json').write_text(json.dumps(dict(method='Positive Fourier expectation, all orders; exact at d=3, sampled otherwise.',rows=yields),indent=2)+'\n')
