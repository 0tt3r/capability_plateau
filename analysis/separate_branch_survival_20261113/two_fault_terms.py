"""Quantify the two-fault contribution to the a_e acceptance coefficient.

An active fault means one with nonzero monitored syndrome signature. Faults
with zero signature are marginalized exactly, regardless of their number.
The 'exactly two active faults' probability is distinct from the quadratic
Taylor coefficient and from the all-orders no-primary contribution.
"""
from collections import Counter, defaultdict
from fractions import Fraction
import json
from math import exp, log1p

from branch_survival import HERE, build


def pair_weight(locations, branches, p, q, exclude_primary=True):
    totals_raw, totals_odds = defaultdict(float), defaultdict(float)
    same_raw = same_odds = log_zero = 0.
    for op, sigs in locations:
        rate = p if op.kind in ('CX', 'R', 'M') else q
        full = {s: rate*n/len(sigs) for s, n in Counter(sigs).items() if s}
        active_rate = sum(full.values())
        log_zero += log1p(-active_rate)
        allowed = {s: w for s, w in full.items()
                   if not exclude_primary or s not in branches}
        z = 1-active_rate
        for s, w in allowed.items():
            totals_raw[s] += w
            totals_odds[s] += w/z
        for target in branches:
            local = sum(w*allowed.get(s ^ target, 0.) for s, w in allowed.items())
            same_raw += local
            same_odds += local/z**2
    raw = (sum(w*totals_raw.get(s ^ t, 0.)
               for t in branches for s, w in totals_raw.items())-same_raw)/2
    odds = (sum(w*totals_odds.get(s ^ t, 0.)
                for t in branches for s, w in totals_odds.items())-same_odds)/2
    return raw, exp(log_zero)*odds, exp(log_zero)


def main():
    saved = json.loads((HERE/'results.json').read_text())
    rows = []
    for row in saved['rows']:
        d = row['distance']
        locations, branches, _, _ = build(d)
        for hardware, v in row['hardware'].items():
            p, q = v['p'], v['q']
            raw, exact_two, zero = pair_weight(locations, branches, p, q)
            raw_all, exact_two_all, _ = pair_weight(locations, branches, p, q, False)
            expected = (631*d+12754)/450*p*p+(635*d+5104)/90*p*q+(47*d+774)/9*q*q
            assert abs(raw-expected) < 1e-14
            primary_pair = 3*(2*p/15+q/3)**2 if d == 6 else 0.
            assert abs(raw_all-raw-primary_pair) < 1e-14
            no_primary = v['H_n_without_primary_faults'][1]
            assert 0 < exact_two <= no_primary
            # All primary-gadget Pauli faults have nonzero signatures. The
            # denominators exclude any other active fault at the same gate.
            one = zero*(d*p/(15*(1-p))+d*q/(6*(1-q)))
            up_to_two = one+exact_two_all
            exact_one_group = v['H_n'][1]
            assert 0 < up_to_two <= exact_one_group
            relative_error = 1-up_to_two/exact_one_group
            assert relative_error < (.01 if hardware == 'H1' else .002)
            coefficients = row['second_order']['ratio_correction']
            j = float(Fraction(coefficients['pp']))*p*p+float(Fraction(coefficients['pq']))*p*q+float(Fraction(coefficients['qq']))*q*q
            out = dict(distance=d, hardware=hardware, u=v['u'],
                       ratio_quadratic_correction_J=j,
                       J_over_u=j/v['u'], raw_rejecting_pair_R=raw,
                       no_active_fault_probability=zero,
                       exactly_two_active_no_primary=exact_two,
                       exactly_two_active_all=exact_two_all,
                       exactly_one_active=one,
                       adopted_one_group_up_to_two=up_to_two,
                       adopted_he=up_to_two/v['u'],
                       retained_h0=v['h0'],
                       relative_underestimate_of_one_group=relative_error,
                       exact_one_group_coefficient=exact_one_group,
                       exact_two_share_of_one_group=exact_two_all/v['H_n'][1],
                       all_orders_without_primary=no_primary,
                       exact_pair_share_of_no_primary=exact_two/no_primary,
                       higher_active_no_primary=no_primary-exact_two)
            rows.append(out)
            print(f'{hardware} d={d}: raw pair R={raw:.9g}, exact two active={exact_two_all:.9g}, '
                  f'share of a_e coefficient={100*out["exact_two_share_of_one_group"]:.4f}%, '
                  f'adopted he={out["adopted_he"]:.9g}, '
                  f'underestimate={100*relative_error:.4f}%', flush=True)
    (HERE/'two_fault_terms.json').write_text(json.dumps(dict(
        convention='Exactly two nonzero-signature faults; any number of zero-signature faults marginalized. Raw R and J are quadratic Taylor terms.',
        adopted_model='Retain exact H0 and use H1 from one or two nonzero-signature faults. Preserve the full background zero-signature product. Two-component acceptance model; capability rerun pending.',
        rows=rows), indent=2)+'\n')


if __name__ == '__main__':
    main()
