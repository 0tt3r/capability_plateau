"""Direct implementation (draft Sec. 3.2): QUOPS circuits on physical qubits, no encoding, no mitigation.

Schedule model [choice]: every native operation takes one time step. Preparation and measurement are one
step each. A CNOT layer is one step on all-to-all hardware; on the grid it is preceded by SWAP rounds of
three CNOT steps each. A rotation layer lasts as long as its slowest qubit (R_Y takes 5 steps).
Fault locations: every CNOT (eps_2q), every preparation and measurement (eps_2q), every one-qubit gate
(eps_1q), and, if idle_errors, every step in which a data qubit is untouched (eps_1q). Empty grid sites
carry no information and never fault.

The no-fault probability of one circuit is prod(1 - p_loc). Its mean over the ensemble factorizes over
layers on all-to-all hardware (a closed form); on the grid the routing depends on the placement left by
earlier layers, so we average over sampled circuits.
"""
from math import exp, log

import numpy as np

from .quops import ENSEMBLES, layer_counts, native_steps, random_pairing, size
from .routing import ROUTERS, grid_shape, pairs_adjacent


def log_rotation_layer(w, hw, ch):
    """log of the mean no-fault probability of one rotation layer (exact average over the axes)."""
    l1 = log(1 - hw.eps_1q)
    steps = {a: native_steps(a, ch.rotation_compile) for a in 'XYZ'}
    if not ch.idle_errors:
        per_qubit = np.mean([exp(steps[a]*l1) for a in 'XYZ'])
        return w*log(per_qubit)
    # With idles, every data qubit is a location in every step of the layer: w*D locations.
    depths = sorted(set(steps.values()))
    mean = 0.
    for k, D in enumerate(depths):
        at_most = sum(steps[a] <= D for a in 'XYZ')/3
        below = sum(steps[a] < D for a in 'XYZ')/3
        mean += (at_most**w - below**w)*exp(w*D*l1)      # P(slowest qubit takes D steps)
    return log(mean)


def log_cnot_layer_routed(rounds, w, hw, ch):
    """log no-fault probability of one routed CNOT layer: SWAP rounds (3 CNOT steps each), then the CNOTs."""
    l1, l2 = log(1 - hw.eps_1q), log(1 - hw.eps_2q)
    total = 0.
    for swaps in rounds:
        moved = [(a, b) for a, b in swaps if a < w or b < w]      # swaps of two empty sites cost nothing
        if not moved:
            continue
        busy = {t for ab in moved for t in ab if t < w}
        total += 3*len(moved)*l2 + (3*(w - len(busy))*l1 if ch.idle_errors else 0.)
    return total + (w//2)*l2


def mean_no_fault(w, ensemble, hw, ch, seed=0):
    """(mean no-fault probability, Monte Carlo standard error) for one width and ensemble."""
    n_rot, n_cnot = layer_counts(w, ensemble)
    l2 = log(1 - hw.eps_2q)
    fixed = 2*w*l2 + n_rot*log_rotation_layer(w, hw, ch)       # preparation, measurement, rotations
    if hw.connectivity == 'all_to_all':
        return exp(fixed + n_cnot*(w//2)*l2), 0.
    r, c = grid_shape(w)
    route = ROUTERS[ch.router]
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(ch.samples):
        placement = np.arange(r*c)
        total = fixed
        for _ in range(n_cnot):
            pairs = random_pairing(w, rng)
            rounds, placement = route(placement, pairs, r, c)
            assert pairs_adjacent(placement, pairs, c)
            total += log_cnot_layer_routed(rounds, w, hw, ch)
        values.append(exp(total))
    values = np.array(values)
    return float(values.mean()), float(values.std(ddof=1)/np.sqrt(len(values)))


def required_qubits(w, hw):
    if hw.connectivity == 'all_to_all':
        return w
    r, c = grid_shape(w)
    return r*c


def sweep(hw, ch, w_max=600, give_up=4):
    """Pass/fail for every even width and both ensembles, until `give_up` consecutive widths fail."""
    rows = []
    for ensemble in ENSEMBLES:
        misses = 0
        for w in range(2, w_max + 1, 2):
            q, err = mean_no_fault(w, ensemble, hw, ch, seed=w)
            ok = q >= exp(-ch.threshold)
            rows.append(dict(ensemble=ensemble, w=w, s=size(w, ensemble), q=q, stderr=err, passes=ok,
                             qubits=required_qubits(w, hw)))
            misses = 0 if ok else misses + 1
            if misses >= give_up:
                break
    return rows
