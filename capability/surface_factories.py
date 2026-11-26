"""Surface-code factories for the pure-T protocols of arXiv:2610.06535, in the layout of Litinski (arXiv:1905.06903v3).

Each factory is returned as a `factories.Factory` row of family 'jws', so the surface-code model of surface.py can use it
unchanged (`config.vary(factories=('cultivation', 'litinski', 'jws'))`) and Pinnacle can use it behind an adapter.

Layout [Litinski Sec. 3, Figs. 11 and 15, generalized]: a protocol with n inputs, k outputs and m = N - k postselected
checks keeps its k outputs in patches of distance d_X and its m checks in patches of distance d_Z, in a region of length
l = k d_X + m d_Z; an ancilla region of the same length mediates the pi/8 rotations, each executed by a faulty T
measurement on a d_m x d_m patch.
  qubits   2 l 3 d_X + 4 d_m                             (15-to-1: 2 (d_X + 4 d_Z) 3 d_X + 4 d_m, Litinski Fig. 11)
  time     (ceil((n - m)/2) + [k > 1]) d_m / (1 - p_fail)   (15-to-1: 6 d_m, 20-to-4: 10 d_m, as in Litinski)
  error    A_d p_rot^d / k + t P_SC(d_X) (1 + l / (2 d_X)) per delivered state: the leading undetected input-error term
           and the storage errors of the output and of the ancilla region; Z errors on the d_Z check patches are
           detected by the postselection and only add to the rejection rate, which is why d_Z < d_X (Litinski Sec. 3.1).
  reject   n (rho + p_rot) + t m (d_X / (2 d_Z)) P_SC(d_Z)
P_SC is the per-cycle surface-code law of the draft (Sec. 3.5), with the uniform-noise threshold.

Inputs:  'faulty T'    p_rot = p, Litinski's faulty T measurement, with no extra qubits and rho = 0;
         'cultivated'  the cultivated inputs of the Pinnacle engines (Pinnacle v2 Eqs. 9 and 10), with the same ancilla
                       pool of 30 (d_a^2 + d_a - 1) + 750 qubits and the same per-input rejection rho;
         'two levels'  n_L1 faulty-T 15-to-1 factories feed a level-2 protocol, costed with Litinski's two-level formula
                       (Sec. 4): space 2 l_2 3 d_X2 + 2 n_L1 ((d_X + 4 d_Z)(3 d_X + d_m2/2) + 2 d_m)
                       + 2 (20 d_m2^2 + 2 d_X2 d_m2) l_2 / (d_X2 + 4 d_Z2), and time (n_2/2) t_L1 with
                       t_L1 = max(d_m2, 6 d_m / (1 - p_fail,1) / (n_L1/2)). The last space term, Litinski's extra ancillas,
                       is scaled with the level-2 region length (a declared assumption for protocols other than 15-to-1).
Calibration (tests/test_surface_factories.py): (15-to-1)_{17,7,7} gives 4,618 qubits and 42.6 cycles (Litinski: 4,620 and
42.6) and output error 3.8e-8 (4.5e-8); (15-to-1)^6_{11,5,5} x (15-to-1)_{25,11,11} gives 30,732 qubits and 82.5 cycles
(30,700 and 82.5); its output error, 5.9e-13 before calibration, is scaled by the factor 4.6 that reproduces
Litinski's 2.7e-12, and the same factor is applied to every two-level design. For each protocol and input only the Pareto-optimal (qubits, period, error) designs are kept.
"""
import json
from functools import lru_cache
from math import ceil
from pathlib import Path

from .config import BASELINE, P
from .factories import Factory

PROTOCOLS = json.loads((Path(__file__).resolve().parents[1]/'results/factory_protocols.json').read_text())['protocols']
CULTIVATED = {7: (2e-5, 1e-4), 9: (2e-6, 5e-5)}        # [paper] Pinnacle v2 Eqs. 9 and 10: d_a -> (p_a, p_in)
CULT_POOL = 750                                           # [paper] 25-qubit cultivation systems of the Pinnacle engine
RHO_CULT = 0.0064                                         # [choice] as in pinnacle_factories: 10% rejection at 15 inputs
DX = tuple(range(7, 33, 2))
DZ = tuple(range(5, 19, 2))                               # Litinski uses d_Z, d_m >= 5 at p = 1e-3
LEVEL1 = ((9, 5), (11, 5), (13, 5), (17, 7))              # (d_X, d_Z = d_m) of the level-1 15-to-1 factories
TWO_LEVEL_SCALE = 2.7e-12/5.894e-13                      # [calibration] Litinski's (15-to-1)^6_{11,5,5} x (15-to-1)_{25,11,11}
N_L1 = (4, 6, 8)


def p_sc(d, p=P):
    return BASELINE.sc_prefactor*(p/BASELINE.sc_threshold_uniform)**((d + 1)/2)


FOLD_LABEL = 'fold f=3'          # inputs from fold-transversal f = 3 cultivation (all-to-all hardware only)


def fold_input(n, t0, p=P):
    """(p_rot, rho, extra qubits) for inputs from a pool of fold-transversal f = 3 cultivation units.

    Each unit is a Rot(13) patch delivering one state per V3/N_u cycles on average, with its discards inside V3; the pool
    has the fewest units whose mean rate meets the n inputs per t0 cycles of one attempt [choice: mean rate, no buffer].
    The state waits at most d_fin + 1 cycles in its distance-13 patch before the rotation."""
    from .factories import FOLD_F3
    units = ceil(n/t0*FOLD_F3['volume']/FOLD_F3['qubits'])
    return FOLD_F3['error'] + 14*p_sc(13, p), 0., units*FOLD_F3['qubits']


def single(proto, dX, dZ, dm, p_rot, rho=0., extra=0, p=P):
    """(qubits, cycles per attempt, rejection, error per delivered state) of a one-level factory."""
    n, k, d, A = proto['n'], proto['k'], proto['d'], proto['A_d']
    m = proto['wires'] - k
    l = k*dX + m*dZ
    qubits = 2*l*3*dX + 4*dm + extra
    steps = ceil((n - m)/2) + (1 if k > 1 else 0)
    t0 = steps*dm
    reject = min(0.99, n*(rho + p_rot) + t0*m*(dX/(2*dZ))*p_sc(dZ, p))
    t = t0/(1 - reject)
    error = A*p_rot**d/k + t*p_sc(dX, p)*(1 + l/(2*dX))
    return qubits, t, reject, error


def two_level(proto, dX1, dZ1, nL1, dX2, dZ2, p=P):
    """(qubits, cycles per attempt, error per delivered state) of n_L1 faulty-T 15-to-1 factories feeding proto."""
    ref = next(q for q in PROTOCOLS if q['name'] == '[[15,1,3]]')
    _, t1, rej1, e1 = single(ref, dX1, dZ1, dZ1, p)
    dm1, dm2 = dZ1, dZ2
    n2, k2, d2, A2 = proto['n'], proto['k'], proto['d'], proto['A_d']
    m2 = proto['wires'] - k2
    l2 = k2*dX2 + m2*dZ2
    qubits = (2*l2*3*dX2 + 2*nL1*((dX1 + 4*dZ1)*(3*dX1 + dm2/2) + 2*dm1)
              + 2*(20*dm2**2 + 2*dX2*dm2)*l2/(dX2 + 4*dZ2))
    tL1 = max(dm2, 6*dm1/(1 - rej1)/(nL1/2))
    t = (n2/2)*tL1
    error = TWO_LEVEL_SCALE*(A2*e1**d2/k2 + t*p_sc(dX2, p)*(1 + l2/(2*dX2)))
    return int(round(qubits)), t, error


def pareto(rows):
    """Keep rows not dominated in (qubits, cycles per output, error)."""
    rows = sorted(rows, key=lambda f: (f.qubits, f.cycles_per_output, f.error))
    kept = []
    for f in rows:
        if not any(g.qubits <= f.qubits and g.cycles_per_output <= f.cycles_per_output and g.error <= f.error
                   for g in kept):
            kept.append(f)
    return kept


@lru_cache(maxsize=None)
def factories(p=P):
    out = []
    for proto in PROTOCOLS:
        k = proto['k']
        groups = {}
        for dX in DX:
            for dZ in DZ:
                if dZ > dX:
                    continue
                cand = [('faulty T', p, 0., 0)] + [(f'cultivated d_a={da}', pin + (da + 1)*pa, RHO_CULT,
                                                     30*(da*da + da - 1) + CULT_POOL)
                                                    for da, (pa, pin) in CULTIVATED.items()]
                m = proto['wires'] - k
                t0 = (ceil((proto['n'] - m)/2) + (1 if k > 1 else 0))*dZ
                cand.append((FOLD_LABEL, *fold_input(proto['n'], t0, p)))
                for label, p_rot, rho, extra in cand:
                    q, t, rej, e = single(proto, dX, dZ, dZ, p_rot, rho, extra, p)
                    if rej >= 0.99 or e >= 1e-3:
                        continue
                    groups.setdefault(label, []).append(
                        Factory('jws', f"{proto['name']}, {label}, d=({dX},{dZ},{dZ})", e, q, t/k, k))
        for dX1, dZ1 in LEVEL1:
            for nL1 in N_L1:
                for dX2 in DX:
                    for dZ2 in DZ:
                        if dZ2 > dX2:
                            continue
                        q, t, e = two_level(proto, dX1, dZ1, nL1, dX2, dZ2, p)
                        if e >= 1e-3:
                            continue
                        groups.setdefault('two levels', []).append(
                            Factory('jws', f"(15-to-1)^{nL1}_{{{dX1},{dZ1},{dZ1}}} x {proto['name']}_{{{dX2},{dZ2},{dZ2}}}",
                                    e, q, t/k, k))
        for rows in groups.values():
            out += pareto(rows)
    # fold-transversal inputs need all-to-all hardware, so the 2D-grid designs are kept on their own Pareto front
    return tuple(pareto([f for f in out if FOLD_LABEL not in f.name]) + pareto([f for f in out if FOLD_LABEL in f.name]))
