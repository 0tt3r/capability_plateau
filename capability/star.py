"""STAR capability with m=1 or m=2 preparation on a 2D surface-code grid.

The schedule groups rotations over F ancilla patches for preparation. Within a group,
RUS proceeds in synchronized rounds, with geometric preparation retries. We
calculate expected duration, including the slowest preparation and RUS paths.
Jensen's inequality turns expected additive loss into a conservative mean
no-fault bound *within the effective stochastic model*. A separate Markov
deadline bound avoids claiming that expected runtime guarantees completion.

Injection is a leading-order analytic model, with explicit uncalibrated yield
and background-error assumptions. See docs/star_surface_code_estimation.md.
"""
from functools import lru_cache
from math import ceil, exp, log, sqrt

import numpy as np

from .config import N_MAX, STAR_BASELINE, WEEK_SECONDS
from .quops import ENSEMBLES, layer_counts, size
from .star_states import state_profile, uniform_angles
from .surface import p_sc


def _validate(hw, ch, sc):
    if hw.connectivity != 'grid':
        raise ValueError('Toshio STAR currently supports H1/H3 (2D grid) only')
    if not 0 < sc.deadline_failure < min(1 - exp(-ch.threshold), .1):
        raise ValueError('deadline_failure must leave a positive success budget')
    if ch.executions < 1 or ch.cycle_seconds <= 0 or ch.tiles_per_data_patch < 1:
        raise ValueError('Invalid execution, cycle-time, or patch-layout settings')
    if ch.cnot_layer not in ('parallel', 'serial'):
        raise ValueError('Unknown CNOT schedule')
    if sc.max_distance < 3 or not sc.groups:
        raise ValueError('A data distance and at least one preparation design are required')
    if sc.rotation_weight not in (1, 2):
        raise ValueError('rotation_weight must be 1 or 2')
    if sc.postselection not in ('full', 'strip'):
        raise ValueError("postselection must be 'full' or 'strip'")
    if min(sc.preparation_cycles, sc.injection_steps) < 1:
        raise ValueError('Preparation and injection must have positive duration')
    if min(sc.growth_steps, sc.delivery_steps, sc.clifford_steps) < 0:
        raise ValueError('Schedule overheads cannot be negative')
    if not 0 < sc.timing_tail_tolerance < .01:
        raise ValueError('Invalid timing summation tolerance')


def batch_preparation_rounds(acceptance, max_batch, tolerance=1e-9):
    """Expected total preparation trial rounds for batches 1..max_batch.

    At RUS round j, each lane is active with probability 2^-j. Let
    r_g = E_angle[(1-a)^g]. The chance that a batch still needs another
    preparation trial after g trials is 1-(1-2^-j*r_g)^b. Sum over j,g.
    This averages geometric retries before taking the parallel maximum;
    replacing a by its mean would give a different answer.

    Both positive tails are bounded analytically; the returned values include
    the bound, so the numerical truncation does not underestimate runtime.
    Extremely low yields use the analytic serial-preparation upper bound rather
    than allocating an impractically large integration array.
    """
    a = np.asarray(acceptance, dtype=float).reshape(-1)
    if not len(a) or not np.all((a > 0) & (a <= 1)) or max_batch < 1:
        raise ValueError('Need positive acceptance probabilities and batch size')
    amin = float(a.min())
    jmax = max(1, ceil(np.log2(4*max_batch/(tolerance*amin))))
    if amin == 1:
        gmax = 1
    else:
        gmax = max(1, ceil(np.log(tolerance*amin/(4*max_batch))/np.log1p(-amin)))
    if gmax > 100_000:
        return 2*np.arange(max_batch + 1)*np.mean(1/a)
    r = np.ones(gmax)
    # Chunk the angle integration to keep large-yield-tail calculations bounded in memory.
    with np.errstate(divide='ignore'):
        log_reject = np.log1p(-a)
    for start in range(1, gmax, 2048):
        g = np.arange(start, min(start + 2048, gmax))
        r[start:start + len(g)] = np.exp(g[:, None]*log_reject).mean(axis=1)
    q = 2.**(-np.arange(jmax))
    totals = np.zeros(max_batch + 1)
    # Increment 1-(1-x)^b across batch sizes, avoiding cancellation for small x.
    for start in range(0, gmax, 2048):
        x = r[start:start + 2048, None]*q
        increment = x.copy()
        cumulative = 0.
        for b in range(1, max_batch + 1):
            cumulative += increment.sum()
            totals[b] += cumulative
            increment *= 1 - x
    tail_per_lane = 2.*np.mean((1 - a)**gmax/a) + 2.**(1 - jmax)*np.mean(1/a)
    totals += np.arange(max_batch + 1)*tail_per_lane
    return totals


def batch_rus_rounds(max_batch, tolerance=1e-12):
    """E[max of b geometric(1/2) RUS lengths], with a conservative tail."""
    jmax = max(1, ceil(np.log2(2*max_batch/tolerance)))
    q = 2.**(-np.arange(jmax))
    out = np.zeros(max_batch + 1)
    increment, cumulative = q.copy(), 0.
    for b in range(1, max_batch + 1):
        cumulative += increment.sum()
        out[b] = cumulative + b*2.**(1 - jmax)
        increment *= 1 - q
    return out


@lru_cache(maxsize=256)
def _profile(hw, ch, sc, d, k):
    return state_profile(uniform_angles(sc.angle_points), k, d, hw, ch, sc)


@lru_cache(maxsize=256)
def _batch_times(hw, ch, sc, d, k, max_batch):
    profile = _profile(hw, ch, sc, d, k)
    trial = sc.preparation_cycles + sc.growth_steps*d
    consume = (sc.injection_steps + sc.delivery_steps + sc.clifford_steps)*d
    prep = batch_preparation_rounds(profile.acceptance, max_batch, sc.timing_tail_tolerance)
    rus = batch_rus_rounds(max_batch)
    correction_runtime = 2*(trial/profile.fallback_acceptance + consume)
    # Serialize the rare physical coherent-cancellation rotations within each batch.
    corrections = 2*np.arange(max_batch + 1)*profile.correction_probability.mean()*correction_runtime
    return trial*prep + consume*rus + corrections


def rotation_losses(profile):
    """Expected log loss per completed rotation, including correction resources."""
    injection = 2*float(np.mean(-np.log1p(-profile.error)))
    preparation = (2 + 4*float(profile.correction_probability.mean()))*profile.preparation_floor_loss
    return injection, preparation


def _one_rotation_cycles(profile, d, sc):
    """Exact mean service time of one lane; a lower bound for any nonempty batch."""
    trial = sc.preparation_cycles + sc.growth_steps*d
    consume = (sc.injection_steps + sc.delivery_steps + sc.clifford_steps)*d
    correction = 2*(trial/profile.fallback_acceptance + consume)
    return 2*(trial*np.mean(1/profile.acceptance) + consume) + 2*profile.correction_probability.mean()*correction


def evaluate_design(w, ensemble, hw, ch, sc, d, k, stations):
    """Evaluate one layout/schedule; return diagnostics even if it fails."""
    _validate(hw, ch, sc)
    if w < 2 or w % 2 or ensemble not in ENSEMBLES:
        raise ValueError('QUOPS needs an even width >= 2 and a known ensemble')
    if d < 3 or d % 2 == 0 or not 1 <= stations <= w:
        raise ValueError('Need odd data distance >= 3 and 1..w preparation stations')
    profile = _profile(hw, ch, sc, d, k)
    # Power-of-two table sizes share expensive timing integrations across nearby widths.
    max_batch = 1 << (stations - 1).bit_length()
    times = _batch_times(hw, ch, sc, d, k, max_batch)
    batches, remainder = divmod(w, stations)
    rotation_layer = batches*times[stations] + times[remainder]
    n_rot, n_cnot = layer_counts(w, ensemble)
    rotations = w*n_rot
    t_cnot = 2*d if ch.cnot_layer == 'parallel' else w*d
    cycles = float(n_rot*rotation_layer + n_cnot*t_cnot + 2*d)
    data_routing = ceil(ch.tiles_per_data_patch*w)
    live_tiles = data_routing + stations  # one data-distance handoff/storage patch per station
    prep_qubits = (2*profile.preparation_distance - 1)**2  # unrotated preparation patch
    per_copy_qubits = live_tiles*2*d*d + stations*prep_qubits
    injection, preparation = rotation_losses(profile)
    loss_patches = live_tiles*cycles*(-np.log1p(-p_sc(d, hw, ch)))
    loss_injection = rotations*injection
    loss_preparation = rotations*preparation
    expected_loss = float(loss_patches + loss_injection + loss_preparation)
    # Jensen: E[exp(-L_path)] >= exp(-E[L_path]), including time/noise correlations.
    no_fault_bound = exp(-expected_loss)
    # Markov: P(runtime > E[runtime]/delta) <= delta. Union with a logical fault.
    success_bound = max(0., no_fault_bound - sc.deadline_failure)
    total_loss = -log(success_bound) if success_bound else float('inf')
    mean_seconds = cycles*ch.cycle_seconds
    deadline_seconds = mean_seconds/sc.deadline_failure
    per_copy = min(ch.executions, int(WEEK_SECONDS/deadline_seconds))
    copies = ceil(ch.executions/max(1, per_copy))
    correction_probability = float(profile.correction_probability.mean())
    return dict(
        ensemble=ensemble, w=w, s=size(w, ensemble), distance=d, groups=k,
        rotation_weight=sc.rotation_weight, postselection=sc.postselection,
        preparation_distance=profile.preparation_distance, stations=stations,
        qubits=int(copies*per_copy_qubits), qubits_per_copy=int(per_copy_qubits),
        copies=copies, runs_per_copy=per_copy, passes=bool(total_loss <= ch.threshold and per_copy >= 1),
        cycles=cycles, run_seconds=mean_seconds, deadline_seconds=deadline_seconds,
        data_routing_tiles=data_routing, handoff_tiles=stations, preparation_qubits_per_station=prep_qubits,
        loss_patches=float(loss_patches), loss_injection=float(loss_injection),
        loss_preparation_memory=float(loss_preparation),
        loss_deadline=float(total_loss - expected_loss), loss_total=float(total_loss),
        success_probability_lower_bound=success_bound,
        mean_state_acceptance=float(profile.acceptance.mean()),
        min_state_acceptance=float(profile.acceptance.min()),
        fallback_fraction=float(profile.fallback.mean()),
        mean_correction_probability=correction_probability,
        resource_states=float(rotations*(2 + 4*correction_probability)),
        preparation_attempts=float(rotations*(2*np.mean(1/profile.acceptance) +
                                             4*correction_probability/profile.fallback_acceptance)),
        rotation_layer_cycles=float(rotation_layer),
    )


def min_qubits(w, ensemble, hw, ch, sc=STAR_BASELINE):
    """Minimum over declared preparation groups, odd data distances and station counts."""
    _validate(hw, ch, sc)
    if w < 2 or w % 2 or ensemble not in ENSEMBLES:
        raise ValueError('QUOPS needs an even width >= 2 and a known ensemble')
    n_rot, n_cnot = layer_counts(w, ensemble)
    budget = -log(exp(-ch.threshold) + sc.deadline_failure)
    best = None
    for k in sc.groups:
        prep_area = (2*sc.rotation_weight*k - 1)**2
        resource_loss = w*n_rot*sum(rotation_losses(_profile(hw, ch, sc, 3, k)))
        if resource_loss >= budget:
            continue
        for d in range(3, sc.max_distance + 1, 2):
            area_min = (ceil(ch.tiles_per_data_patch*w) + 1)*2*d*d + prep_area
            if best is not None and area_min >= best['qubits']:
                break
            profile = _profile(hw, ch, sc, d, k)
            lane_cycles = _one_rotation_cycles(profile, d, sc)
            fixed = n_cnot*(2*d if ch.cnot_layer == 'parallel' else w*d) + 2*d
            # Even unlimited stations cannot beat one lane per rotation layer.
            if (n_rot*lane_cycles + fixed)*ch.cycle_seconds/sc.deadline_failure > WEEK_SECONDS:
                continue
            for stations in range(1, w + 1):
                area = (ceil(ch.tiles_per_data_patch*w) + stations)*2*d*d + stations*prep_area
                if best is not None and area >= best['qubits']:
                    break
                cycles_min = n_rot*ceil(w/stations)*lane_cycles + fixed
                tiles = ceil(ch.tiles_per_data_patch*w) + stations
                if resource_loss + tiles*cycles_min*(-np.log1p(-p_sc(d, hw, ch))) > budget:
                    continue
                row = evaluate_design(w, ensemble, hw, ch, sc, d, k, stations)
                if row['passes'] and (best is None or row['qubits'] < best['qubits']):
                    best = row
    return best


def resource_ceiling(hw, ch, sc=STAR_BASELINE):
    """Continuous QUOPS-size ceiling from resource noise alone (s=2R)."""
    _validate(hw, ch, sc)
    loss = min(sum(rotation_losses(_profile(hw, ch, sc, 3, k))) for k in sc.groups)
    return 2*(-log(exp(-ch.threshold) + sc.deadline_failure))/loss if loss > 0 else float('inf')


def sweep(hw, ch, sc=STAR_BASELINE):
    """Every even width below the resource-only ceiling; no monotonic pruning assumption."""
    ceiling = resource_ceiling(hw, ch, sc)
    if not np.isfinite(ceiling):
        raise ValueError('An unbounded zero-noise sweep needs an explicit workload limit')
    rows = []
    for ensemble in ENSEMBLES:
        upper = sqrt(ceiling) if ensemble == 'shallow' else np.cbrt(ceiling)
        for w in range(2, int(upper) + 1, 2):
            row = min_qubits(w, ensemble, hw, ch, sc)
            if row is not None and row['qubits'] <= N_MAX:
                rows.append(row)
    return rows


def best_at(N, rows):
    return max((r for r in rows if r['passes'] and r['qubits'] <= N),
               key=lambda r: (r['s'], -r['qubits']), default=None)
