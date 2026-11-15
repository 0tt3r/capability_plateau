"""Exact postselection branch probabilities for the specified Toshio circuit.

Independent categorical Pauli faults are propagated with the existing checked
Clifford signature engine. A spatial Fourier transfer calculation sums ALL
fault multiplicities, without sampling. It determines acceptance-sector
weights, not the complete logical channel of the accepted state.

Production defaults and earlier capability analyses are not modified.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from pathlib import Path
import argparse
import importlib.util
import json
import sys
from time import monotonic

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    'rejection_counter', HERE.parent/'rejection_count_20261111/count_rejection.py')
counter = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = counter
spec.loader.exec_module(counter)


def bits(mask):
    return tuple(i for i in range(mask.bit_length()) if (mask >> i) & 1)


def parity(x):
    x = x.copy()
    for shift in (16, 8, 4, 2, 1):
        x ^= x >> shift
    return x & 1


def local_mask(mask, support):
    return sum(((mask >> bit) & 1) << j for j, bit in enumerate(support))


def build(d):
    coords, ops, nd = counter.make_circuit(d)
    locations, branches = counter.fault_signatures(coords, ops)
    checks = counter.check_forward(coords, ops, locations)
    bitxy = {}
    for op in ops:
        if op.kind == 'M':
            for bit in bits(op.record):
                bitxy[bit] = coords[op.qubits[0]]
    order = sorted(range(nd), key=lambda b: (bitxy[b][0], bitxy[b][1], b))
    remap = {b: i for i, b in enumerate(order)}

    def reorder(sig):
        return sum(1 << remap[b] for b in bits(sig))

    locations = [(op, [reorder(s) for s in sigs]) for op, sigs in locations]
    branches = [reorder(b) for b in branches]
    return locations, branches, nd, checks


def factors_for(locations, branches, p, q, exclude_primary=False):
    """Group identical location supports; retain categorical gate noise exactly."""
    grouped = defaultdict(Counter)
    for op, sigs in locations:
        all_counts = Counter(s for s in sigs if s)
        if not all_counts:
            continue
        union = 0
        for s in all_counts:
            union |= s
        support = bits(union)
        kind = 'p' if op.kind in ('CX', 'R', 'M') else 'q'
        excluded = {b for b in branches} if exclude_primary else set()
        removed = sum(n for s, n in all_counts.items() if s in excluded)
        counts = {s: n for s, n in all_counts.items() if s not in excluded}
        pattern = tuple(sorted((local_mask(s, support), n) for s, n in counts.items()))
        grouped[support][kind, len(sigs), removed, pattern] += 1
    factors = []
    for support, patterns in grouped.items():
        indices = np.arange(1 << len(support), dtype=np.uint32)
        value = np.ones(len(indices))
        for (kind, denom, removed, pattern), multiplicity in patterns.items():
            rate = p if kind == 'p' else q
            odd_weight = sum(n*parity(indices & s) for s, n in pattern)
            # Excluding an outcome removes its probability; it is NOT
            # reassigned to the identity or renormalized.
            value *= (1-rate*removed/denom-2*rate*odd_weight/denom)**multiplicity
        factors.append((support, 'noise', value))
    for b in branches:
        support = bits(b)
        indices = np.arange(1 << len(support), dtype=np.uint32)
        value = 1.-2.*parity(indices)
        factors.append((support, 'branch', value))
    return factors


def transfer(factors, nd):
    """Sum dual detector variables; polynomial degree tracks branch weight.

    P(S=t)=2^-D sum_x (-1)^(x.t) product_l E[(-1)^(x.s_l)].
    Inserting product_j (1+z*(-1)^(x.B_j)) gives coefficient n equal to
    sum_{|b|=n} P(S=xor_{j in b} B_j). Complementary coefficients agree.
    """
    factors = sorted(factors, key=lambda f: (min(f[0]), max(f[0]), f[1]))
    future = Counter(b for support, _, _ in factors for b in support)
    assert set(future) == set(range(nd))
    active = []
    state = np.ones((1, 1))
    max_frontier = 0
    for support, kind, values in factors:
        added = [b for b in support if b not in active]
        if added:
            state = np.tile(state, (1, 1 << len(added)))
            active.extend(added)
        max_frontier = max(max_frontier, len(active))
        assignments = np.arange(1 << len(active), dtype=np.uint32)
        local = np.zeros(len(assignments), dtype=np.uint32)
        for j, bit in enumerate(support):
            local |= ((assignments >> active.index(bit)) & 1) << j
        factor = values[local]
        if kind == 'noise':
            state *= factor[None, :]
        else:
            assert kind == 'branch'
            updated = np.zeros((len(state)+1, state.shape[1]))
            updated[:-1] += state
            updated[1:] += state*factor[None, :]
            state = updated
        for bit in support:
            future[bit] -= 1
        retired = [i for i, bit in enumerate(active) if not future[bit]]
        if retired:
            axes = tuple(len(active)-i for i in retired)
            state = state.reshape((len(state),)+(2,)*len(active)).sum(axis=axes)
            state = state.reshape((len(state), -1))/(1 << len(retired))
            active = [b for i, b in enumerate(active) if i not in retired]
    assert not active
    return state[:, 0], max_frontier


def exact(locations, branches, nd, p, q, exclude_primary=False):
    coefficients, frontier = transfer(factors_for(locations, branches, p, q, exclude_primary), nd)
    np.testing.assert_allclose(coefficients, coefficients[::-1], atol=5e-14, rtol=1e-10)
    assert np.min(coefficients) >= -5e-14
    # Use one representative of each globally complementary pair.
    h_n = coefficients[:(len(branches)+1)//2]
    assert sum(h_n) <= 1+5e-14
    return h_n, frontier


def second_order(locations, branches):
    """Exact rational p^2, pq, q^2 coefficients for P(S=0), sum P(S=B_j)."""
    totals = {kind: defaultdict(Fraction) for kind in ('p', 'q')}
    same = {kind: defaultdict(Fraction) for kind in ('p', 'q')}
    targets = (0, *branches)
    for op, sigs in locations:
        kind = 'p' if op.kind in ('CX', 'R', 'M') else 'q'
        f = {s: Fraction(n, len(sigs)) for s, n in Counter(sigs).items() if s}
        if not f:
            continue
        f[0] = -sum(f.values())
        for s, weight in f.items():
            totals[kind][s] += weight
        for target in targets:
            same[kind][target] += sum(w*f.get(s ^ target, 0) for s, w in f.items())

    def conv(a, b, target):
        return sum(w*b.get(s ^ target, 0) for s, w in a.items())

    answer = {}
    for label, ts in [('ideal', [0]), ('one_group', branches)]:
        pp = sum((conv(totals['p'], totals['p'], t)-same['p'][t])/2 for t in ts)
        qq = sum((conv(totals['q'], totals['q'], t)-same['q'][t])/2 for t in ts)
        pq = sum(conv(totals['p'], totals['q'], t) for t in ts)
        firstp = sum(totals['p'].get(t, 0) for t in ts)
        firstq = sum(totals['q'].get(t, 0) for t in ts)
        answer[label] = {key: str(value) for key, value in
                        zip(('p', 'q', 'pp', 'pq', 'qq'), (firstp, firstq, pp, pq, qq))}
    # Divide out shared extensive first-order depletion. These local terms
    # determine H_1/H_0 = u + J_pp p^2 + J_pq pq + J_qq q^2 + O(eps^3).
    a, e = answer['ideal'], answer['one_group']
    answer['ratio_correction'] = dict(
        pp=str(Fraction(e['pp'])-Fraction(a['p'])*Fraction(e['p'])),
        pq=str(Fraction(e['pq'])-Fraction(a['p'])*Fraction(e['q'])
               -Fraction(a['q'])*Fraction(e['p'])),
        qq=str(Fraction(e['qq'])-Fraction(a['q'])*Fraction(e['q'])))
    return answer


def ratio_formula(d):
    """Bulk and boundary terms in H_1/H_0-u; odd k=d/2 >=3."""
    result = dict(pp=Fraction(661*d+12754, 450),
                  pq=Fraction(635*d+5104, 90),
                  qq=Fraction(97*d+1548, 18))
    if d == 6:
        # Two distinct one-group faults equal the complementary third group.
        for name, extra in zip(('pp', 'pq', 'qq'),
                               (Fraction(4, 75), Fraction(4, 15), Fraction(1, 3))):
            result[name] += extra
    return result


def check_readout_pair(locations, branches):
    """Two readout flips can hide either end of a contiguous branch prefix."""
    prefix = 0
    for n, branch in enumerate(branches[:-1], 1):
        prefix ^= branch
        # y=0, so the row-major qubit index equals x=4n-1.
        ms = [sigs[0] for op, sigs in locations
              if op.kind == 'M' and op.qubits == (4*n-1,)
              and op.stage in ('syndrome_1', 'syndrome_2')]
        assert len(ms) == 2 and ms[0] ^ ms[1] == prefix
        assert all(s not in branches and s != 0 for s in ms)


def validate_transfer():
    """Compare against direct categorical enumeration on a small independent model."""
    rng = np.random.default_rng(78243)
    for nd in (3, 4, 5):
        probabilities = []
        factors = []
        for _ in range(6):
            ss = rng.choice(1 << nd, 3, replace=False)
            ps = rng.dirichlet(np.ones(4))
            probabilities.append([(0, ps[0])]+list(zip(ss, ps[1:])))
            assignments = np.arange(1 << nd, dtype=np.uint32)
            value = np.full(1 << nd, ps[0])
            for s, p in zip(ss, ps[1:]):
                value += p*(1.-2.*parity(assignments & s))
            factors.append((tuple(range(nd)), 'noise', value))
        branches = [1, 2, 3]
        for b in branches:
            support = bits(b)
            factors.append((support, 'branch', 1.-2.*parity(np.arange(1 << len(support), dtype=np.uint32))))
        dist = np.zeros(1 << nd)
        for outcome in product(*probabilities):
            s, p = 0, 1.
            for sig, weight in outcome:
                s ^= int(sig)
                p *= weight
            dist[s] += p
        expected = np.zeros(4)
        for mask in range(8):
            s = 0
            for j, b in enumerate(branches):
                if (mask >> j) & 1:
                    s ^= b
            expected[mask.bit_count()] += dist[s]
        got, _ = transfer(factors, nd)
        np.testing.assert_allclose(got, expected, rtol=2e-14, atol=2e-14)


def run(distances):
    validate_transfer()
    rows = []
    for d in distances:
        start = monotonic()
        locations, branches, nd, checks = build(d)
        k = d//2
        pair = second_order(locations, branches)
        assert all(Fraction(pair['ratio_correction'][name]) == value
                   for name, value in ratio_formula(d).items())
        check_readout_pair(locations, branches)
        cp, cq = Fraction(920*d-664, 15), Fraction(455*d-134, 6)
        assert Fraction(pair['ideal']['p']) == -cp-Fraction(d, 15)
        assert Fraction(pair['ideal']['q']) == -cq-Fraction(d, 6)
        assert Fraction(pair['one_group']['p']) == Fraction(d, 15)
        assert Fraction(pair['one_group']['q']) == Fraction(d, 6)
        zero, _ = exact(locations, branches, nd, 0., 0.)
        np.testing.assert_allclose(zero, [1]+[0]*(len(zero)-1), atol=2e-14)
        row = dict(distance=d, groups=k, detectors=nd, forward_checks=checks,
                   second_order=pair, hardware={})
        for label, p, q in [('H1', .001, .001), ('H3', .001, .000001)]:
            h_n, frontier = exact(locations, branches, nd, p, q)
            no_primary, _ = exact(locations, branches, nd, p, q, exclude_primary=True)
            assert np.all(no_primary <= h_n+5e-14)
            u = d*p/15+d*q/6
            h_prod = (1-p)**float(cp)*(1-q)**float(cq)
            h0 = h_n[0]/(1-u)
            he = h_n[1]/u
            primary_prior = -np.expm1(d*np.log1p(-p/15)+k*np.log1p(-q/3))
            phi = (np.arange(512)+.5)*(np.pi/8)/512
            alpha = np.arctan(np.tan(phi)**(1/k))
            s, c = np.sin(alpha), np.cos(alpha)
            intrinsic = np.array([s**(2*n)*c**(2*(k-n))+s**(2*(k-n))*c**(2*n)
                                  for n in range((k+1)//2)])
            accepted = h_n[:, None]*intrinsic
            baseline_bad = u*intrinsic[1]/((1-u)*intrinsic[0]+u*intrinsic[1])
            two_sector_bad = accepted[1]/(accepted[0]+accepted[1])
            correction = pair['ratio_correction']
            ratio_2 = u+float(Fraction(correction['pp']))*p*p+float(Fraction(correction['pq']))*p*q+float(Fraction(correction['qq']))*q*q
            row['hardware'][label] = dict(p=p, q=q, u=u, h_product=h_prod,
                H_n=h_n.tolist(), h0=h0, he=he, he_over_h0=he/h0,
                H_n_without_primary_faults=no_primary.tolist(),
                he_contribution_with_primary_fault=(h_n[1]-no_primary[1])/u,
                he_contribution_without_primary_fault=no_primary[1]/u,
                primary_fault_prior=primary_prior,
                pass_one_group_given_primary=(h_n[1]-no_primary[1])/primary_prior,
                mean_two_sector_bad_fraction_common=baseline_bad.mean(),
                mean_two_sector_bad_fraction_separate=two_sector_bad.mean(),
                higher_sector_share_of_accepted_nonzero_sectors=accepted[2:].sum()/accepted[1:].sum(),
                higher_sector_share_of_acceptance=accepted[2:].sum()/accepted.sum(),
                error_odds_multiplier_vs_common=(h_n[1]/h_n[0])/(u/(1-u)),
                H1_over_H0_second_order=ratio_2, frontier=frontier)
            print(f'd={d} {label}: h0={h0:.10g} he={he:.10g} he/h0={he/h0:.8g} Hn={h_n}', flush=True)
        # Check the first TWO perturbative coefficients, independently of
        # the Fourier derivation; the residual should scale cubically.
        residuals = []
        for eps in (1e-6, 5e-7):
            values, _ = exact(locations, branches, nd, eps, eps/3)
            errors = []
            for n, name in enumerate(('ideal', 'one_group')):
                z = {key: float(Fraction(v)) for key, v in pair[name].items()}
                approximation = (1 if n == 0 else 0)+(z['p']+z['q']/3)*eps+(z['pp']+z['pq']/3+z['qq']/9)*eps**2
                errors.append(float(values[n]-approximation))
            residuals.append(errors)
        ratio = np.array(residuals[0])/np.array(residuals[1])
        assert np.all((ratio > 7.7) & (ratio < 8.3)), (d, residuals, ratio)
        row['small_noise_remainders'] = residuals
        row['seconds'] = monotonic()-start
        rows.append(row)
        (HERE/'results.json').write_text(json.dumps(dict(
            status='Circuit acceptance-sector calculation; full logical channel not recalibrated; production defaults unchanged.',
            method='Exact Fourier transfer of independent categorical Pauli fault signatures; no Monte Carlo.',
            rows=rows), indent=2)+'\n')
        print(f'd={d} complete in {row["seconds"]:.2f}s', flush=True)
    return rows


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--distances', type=int, nargs='+', default=[6, 10, 14, 18])
    args = parser.parse_args()
    for d in args.distances:
        if d < 6 or d % 4 != 2:
            parser.error('Use d_p=2k with odd k>=3 (6, 10, 14, ...).')
    run(args.distances)
