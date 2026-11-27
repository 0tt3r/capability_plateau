"""Magic-state supplies for the Pinnacle consumer, each a black box (draft Secs. 3.6 and 3.7).

A supply is a footprint, a mean delivery period, and an error per delivered T state, plus an interface cost shared
by all factories of one copy. Times are in code cycles t_c (1 us), the unit of capability/pinnacle.py. Four families:

  surface    The surface-code factories of capability/factories.py (cultivation and the seven Litinski rows),
             unchanged, delivered through one surface-code patch of distance d_a that is joined to the processing
             unit by generalized surgery. Adapter footprint 2 d_a^2 + 4 plus a bridge of 2 d_a - 1 qubits, which is
             the 454-qubit patch and 29-qubit bridge of the Frozen Pinnacle seed adapter at d_a = 15. Transfer error
             (d_a + 1) P_SC(d_a) per state, from the draft's surface-code law with the uniform-noise threshold,
             since Pinnacle claims no bimodal gain.
  reference  The two published Pinnacle engines (Pinnacle v2, arXiv:2602.11457v2, Eqs. 9 and 10): 15-to-1 on a
             d = 24 GB block, 4410 and 5430 qubits, 26 code cycles, 10% rejection, output targets 1e-9 and 1e-11.
  gb_cult    The Pinnacle engine template of v2 Eqs. 4 to 6 and 11, generalized from 15-to-1 to the other pure-T
             protocols of arXiv:2610.06535 (15-to-1 itself is the 'reference' family), fed by cultivated inputs injected
             through distance-d_a ancilla codes.
  gb_inject  The same template with direct in-place injection: the noisy physical T state is the ancilla
             (n_a = d_a = 1), so p_rot = 3p. One level reaches only ~1e-6 at p = 1e-3; the two-level variant
             injects the logical outputs of F1 level-1 engines into a level-2 engine ('gb_inject2'), with a block-to-block
             transfer error of 2 P_PB(24) per state and no qubits charged for the bridges between engines.

Engine template (n inputs, k outputs, N wires, m = N - k postselected checks, r_ps = 2 repeated rounds):
  qubits   n_me = B n_cb + 10 n_g + 30 (n_a + d_a - 1) + n_alpha, with B = ceil(2N / k_e) GB blocks of the d = 24
           code ([[510,16,24]]: n_cb = 1020, n_g = 99, k_e = 16), n_a = d_a^2, and n_alpha = 25 x 30 qubits of
           cultivation ancillas for the p = 1e-3 cases. This reproduces the published 4410 and 5430 exactly.
  p_rot    p_in + (d_a + 1) p_a                                         (v2 Eq. 5)
  error    A_d p_rot^d / k per delivered state: a batch is bad with probability A_d p_rot^d, and each batch
           supplies k states. The postselection-undetected term of v2 Eq. 5 is below 1e-14 and dropped.
  reject   n rho + n p_rot + m r_ps p^r_ps, with rho the per-input cultivation failure, fixed by the published 10%
           for 15 inputs (v2 Eq. 6 and its footnote). Extrapolated linearly in n for larger protocols.
  period   t_me = max{2 d_a + m r_ps, t_r + m r_ps, d_a + t_r + (m - 1) r_ps}   (v2 Eq. 11, m checks not 4);
           an attempt takes max(26, t_me) code cycles per ceil(n / 16) groups of 16 injections (two batches of at
           most 8 parallel injections per logical cycle). Mean delivery period = attempt / (k (1 - reject)).
"""
import json
from dataclasses import dataclass, replace
from math import ceil
from pathlib import Path

from .config import BASELINE, P
from .factories import menu
from .surface import p_sc

FIT_A, FIT_B, FIT_C = 6.2, 0.0158, 0.47      # [paper] Pinnacle v2 Table 2, logical measurements


def p_block_d(d, p=P):
    """Eq. 12 of Pinnacle v2: failure of a whole GB block of distance d per native logical period."""
    return FIT_A*(p/FIT_B)**(d/2 + FIT_C)


PROTOCOLS = json.loads((Path(__file__).resolve().parents[1]/'results/factory_protocols.json').read_text())['protocols']

# [paper] Pinnacle v2, Table 1 ([[510,16,24]] row) and Sec. V.2: the engine block used at p = 1e-3.
N_CB, N_G, K_E, ENGINE_CYCLE = 1020, 99, 16, 26
PARALLEL_INJECTIONS = 16           # two batches of at most eight injections per logical cycle
CULT_ANCILLA = 25*30               # [paper] 25-qubit cultivation systems, two per injected state, 15 injected
RHO_CULT = 0.0064                  # [choice] per-input cultivation failure; 15 inputs give the published ~10% reject
POSTSELECT_ROUNDS = 2              # [paper] r at p = 1e-3
REACTION = 10                      # [paper] t_r / t_c
# [paper] Pinnacle v2 Eqs. 9 and 10: (d_a, p_a, p_in) of the cultivated ancilla codes at p = 1e-3.
CULTIVATED = {7: (2e-5, 1e-4), 9: (2e-6, 5e-5)}
PUBLISHED_TARGET = {7: 1e-9, 9: 1e-11}
ADAPTER_DISTANCES = (11, 15, 19, 23)


@dataclass(frozen=True)
class Supply:
    family: str
    name: str
    qubits: int          # footprint of one factory
    cycles: float        # mean code cycles between delivered T states, rejection included
    error: float         # logical error charged per delivered T state
    port_qubits: int = 0     # interface shared by all factories of a copy
    port_error: float = 0.   # interface error per delivered state
    min_slot: int = 0        # smallest consumer slot this supply allows, in code cycles
    reject: float = 0.


def engine_time(m, d_a, r=POSTSELECT_ROUNDS):
    """v2 Eq. 11 generalized to m postselected checks."""
    return max(2*d_a + m*r, REACTION + m*r, d_a + REACTION + (m - 1)*r)


def engine(proto, d_a, p_a, p_in, cultivated, p=P, n_alpha=None, published=None, extra_slots=0.):
    """One GB magic engine for a pure-T protocol. published overrides the computed output error."""
    n, k, d, wires, A = proto['n'], proto['k'], proto['d'], proto['wires'], proto['A_d']
    m = wires - k
    blocks = ceil(2*wires/K_E)
    n_a = d_a*d_a
    n_al = (CULT_ANCILLA if cultivated else 0) if n_alpha is None else n_alpha
    qubits = blocks*N_CB + 10*N_G + 30*(n_a + d_a - 1) + n_al
    p_rot = p_in + (d_a + 1)*p_a
    p_m = p**POSTSELECT_ROUNDS
    reject = min(0.99, n*(RHO_CULT if cultivated else 0.) + n*p_rot + m*POSTSELECT_ROUNDS*p_m)
    attempt = max(ENGINE_CYCLE, engine_time(m, d_a))*ceil(n/PARALLEL_INJECTIONS) + extra_slots
    error = published if published is not None else A*p_rot**d/k
    return Supply('gb', proto['name'], qubits, attempt/(k*(1 - reject)), error, 0, 0, ENGINE_CYCLE, reject), attempt, p_rot


def reference_supplies():
    """The two published 15-to-1 engines, with the output error of the engine template (Eq. pinnacle-engine-error), the
    same estimate as every other engine; the published targets PUBLISHED_TARGET (1e-9, 1e-11) are not used."""
    proto = next(q for q in PROTOCOLS if q['name'] == '[[15,1,3]]')
    out = []
    for d_a, (p_a, p_in) in CULTIVATED.items():
        s, _, _ = engine(proto, d_a, p_a, p_in, True)
        out.append(Supply('reference', f'15-to-1 (Pinnacle, d_a={d_a})', s.qubits, s.cycles, s.error, 0, 0,
                          s.min_slot, s.reject))
    return out


def cultivated_supplies():
    out = []
    for proto in PROTOCOLS:
        if proto['name'] == '[[15,1,3]]':
            continue                    # the published engines of reference_supplies() are this template at 15-to-1
        for d_a, (p_a, p_in) in CULTIVATED.items():
            s, _, _ = engine(proto, d_a, p_a, p_in, True)
            out.append(Supply('gb_cult', f"{proto['name']}, cultivated d_a={d_a}", s.qubits, s.cycles, s.error, 0, 0,
                              s.min_slot, s.reject))
    return out


def _with_rejection(s, attempt, k, rej):
    """The same engine with a stated rejection rate per attempt instead of the first-order estimate."""
    return replace(s, cycles=attempt/(k*(1 - rej)), reject=rej)


def injected_supplies(p=P, level1_engines=(2, 4, 8, 16, 32), inplace_reject=None):
    """Direct in-place injection: one level, and two levels with F1 level-1 15-to-1 engines feeding a level-2 engine.
    inplace_reject, if given, replaces the first-order rejection of every in-place attempt (both levels)."""
    ref = next(q for q in PROTOCOLS if q['name'] == '[[15,1,3]]')
    out = []
    one = {}
    for proto in PROTOCOLS:
        s, att, p_rot = engine(proto, 1, p, p, False, p=p)
        if inplace_reject is not None:
            s = _with_rejection(s, att, proto['k'], inplace_reject)
        one[proto['name']] = (s, att)
        out.append(Supply('gb_inject', f"{proto['name']}, injected, 1 level", s.qubits, s.cycles, s.error, 0, 0,
                          s.min_slot, s.reject))
    s1, _ = one[ref['name']]
    port = 2*p_block_d(24, p)                              # engine-block transfer error between levels
    for proto in PROTOCOLS:
        for f1 in level1_engines:
            p_in2 = s1.error + port
            s2, att2, _ = engine(proto, 1, 0., p_in2, False, p=p, n_alpha=0)
            if inplace_reject is not None:
                s2 = _with_rejection(s2, att2, proto['k'], inplace_reject)
            feed = proto['n']*s1.cycles/f1                    # cycles for F1 engines to supply one attempt
            period = max(att2, feed)/(proto['k']*(1 - s2.reject))
            out.append(Supply('gb_inject2', f"{proto['name']}, injected, 2 levels, F1={f1}",
                              f1*s1.qubits + s2.qubits, period, s2.error, 0, 0, ENGINE_CYCLE, s2.reject))
    return out


def surface_supplies(hw=None, distances=ADAPTER_DISTANCES, families=('cultivation', 'litinski'), label='surface'):
    """Surface-code factories behind a distance-d_a adapter patch."""
    out = []
    for f in menu(families):
        for d_a in distances:
            adapter = 2*d_a*d_a + 4 + 2*d_a - 1
            transfer = (d_a + 1)*p_sc(d_a, _UNIFORM, BASELINE)
            out.append(Supply(label, f'{f.name}, adapter d_a={d_a}', f.qubits, f.cycles_per_output, f.error,
                              adapter, transfer, d_a))
    return out


class _Uniform:
    noise, p = 'uniform', P


_UNIFORM = _Uniform()


def _catalogue_protocols():
    """Distance-2 and distance-3 catalogue protocols (scripts/catalogue_factories.py) with their T yield per run."""
    import os
    path = Path(os.environ.get('CATALOGUE_JSON', Path(__file__).resolve().parents[1]/'results/catalogue_factories.json'))
    if not path.exists():
        return []
    return [dict(name=f"[[{q['n']},{q['k']},{q['d']}]] {q['gate']}", n=q['n'], k=q['k'], d=q['d'], wires=q['wires'],
                 A_d=q['A_d'], yield_=q['T_per_run']) for q in json.loads(path.read_text()).get('protocols', [])]


def _per_T(s, attempt, proto, family, name, rej=None):
    """An engine delivering proto['yield_'] T states per accepted run instead of k."""
    y, k = proto['yield_'], proto['k']
    rej = s.reject if rej is None else rej
    return replace(s, family=family, name=name, cycles=attempt/(y*(1 - rej)), error=s.error*k/y, reject=rej)


def catalogue_engines(families, p=P, level1_engines=(2, 4, 8, 16, 32), inplace_reject=None):
    """The catalogue protocols as Pinnacle magic engines, with the template of the published 15-to-1 engine:
      gb_cat_cult     cultivated inputs in ancilla codes, injected by generalized surgery (d_a = 7, 9)
      gb_cat_inject   physical T states injected in place, one level
      gb_cat_inject2  two levels: F1 one-level in-place [[15,1,3]] engines feed the protocol's engine"""
    out = []
    protos = _catalogue_protocols()
    ref = next(q for q in PROTOCOLS if q['name'] == '[[15,1,3]]')
    s1, att1, _ = engine(ref, 1, p, p, False, p=p)
    if inplace_reject is not None:
        s1 = _with_rejection(s1, att1, 1, inplace_reject)
    port = 2*p_block_d(24, p)
    for proto in protos:
        if (proto['n'], proto['k'], proto['d']) == (15, 1, 3):
            continue                          # the published 15-to-1 engines, with their quoted errors
        if 'gb_cat_cult' in families:
            for d_a, (p_a, p_in) in CULTIVATED.items():
                s, att, _ = engine(proto, d_a, p_a, p_in, True, p=p)
                out.append(_per_T(s, att, proto, 'gb_cat_cult', f"{proto['name']}, cultivated engine d_a={d_a}"))
        if 'gb_cat_inject' in families:
            s, att, _ = engine(proto, 1, p, p, False, p=p)
            out.append(_per_T(s, att, proto, 'gb_cat_inject', f"{proto['name']}, in place, 1 level", inplace_reject))
        if 'gb_cat_inject2' in families:
            for f1 in level1_engines:
                s2, att2, _ = engine(proto, 1, 0., s1.error + port, False, p=p, n_alpha=0)
                rej2 = s2.reject if inplace_reject is None else inplace_reject
                feed = proto['n']*s1.cycles/f1
                y = proto['yield_']
                out.append(replace(s2, family='gb_cat_inject2', name=f"{proto['name']}, in place, 2 levels, F1={f1}",
                                   qubits=f1*s1.qubits + s2.qubits, cycles=max(att2, feed)/(y*(1 - rej2)),
                                   error=s2.error*proto['k']/y, reject=rej2))
    return [s for s in out if s.error < 1e-3 and s.reject < 0.99]    # the filter of the surface-code layouts


def all_supplies(families=('surface', 'surface_jws', 'reference', 'gb_cult', 'gb_inject', 'gb_inject2'),
                 distances=ADAPTER_DISTANCES, inplace_reject=None):
    out = []
    if 'surface' in families:
        out += surface_supplies(distances=distances)
    if 'surface_fold' in families:       # Table-1 factories with fold-transversal cultivation (all-to-all), behind the adapter
        out += surface_supplies(distances=distances, families=('fold_cultivation', 'litinski'), label='surface_fold')
    if 'surface_jws' in families:        # surface-code layouts of the arXiv:2610.06535 protocols, behind the same adapter
        out += surface_supplies(distances=distances, families=('jws',), label='surface_jws')
    if 'surface_catalogue' in families:  # surface-code layouts of the distance-2 and distance-3 catalogue factories
        out += surface_supplies(distances=distances, families=('catalogue',), label='surface_catalogue')
    if 'reference' in families:
        out += reference_supplies()
    if 'gb_cult' in families:
        out += cultivated_supplies()
    if any(f in families for f in ('gb_cat_cult', 'gb_cat_inject', 'gb_cat_inject2')):
        out += catalogue_engines(families, inplace_reject=inplace_reject)
    if 'gb_inject' in families or 'gb_inject2' in families:
        out += [s for s in injected_supplies(inplace_reject=inplace_reject) if s.family in families]
    return out
