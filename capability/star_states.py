"""STAR resource states in the exp(+i phi P) angle convention.

The adopted m=2 model uses counted acceptance, retained single-fault logical
components, and ideal auxiliary coherent corrections with physical supply costs.
The legacy analytic model remains explicit for historical and m=1 comparisons.
Neither model claims a complete calibrated accepted-state logical channel.
"""
from dataclasses import dataclass

import numpy as np

from .surface import p_sc
from .star_data import original_data, transversal_data, validate_counted


def fold_angle(phi):
    """Remove Clifford multiples of pi/4 in the paper's convention."""
    return (np.asarray(phi) + np.pi/8) % (np.pi/4) - np.pi/8


def repository_angle(theta):
    """Repository exp(-i theta P/2) -> folded paper angle."""
    return fold_angle(-np.asarray(theta)/2)


def uniform_angles(points):
    """Midpoint quadrature of |phi| on [0, pi/8] for uniform QUOPS rotations."""
    if points < 8:
        raise ValueError('At least eight angle points are required')
    return (np.arange(points) + .5)*(np.pi/8)/points


def ideal_tmr(phi, k):
    """Shared m=1/m=2 projection formulas, Eqs. (6), (7), (11).

    tan(phi) = tan(alpha)**k. Only odd k >= 3 is implemented, so no even-k
    Clifford correction is silently omitted. Negative target angles use symmetry.
    """
    if k < 3 or k % 2 == 0:
        raise ValueError('Use an odd number of physical rotations, k >= 3')
    phi = fold_angle(phi)
    magnitude = np.abs(phi)
    alpha = np.arctan(np.tan(magnitude)**(1/k))
    s, c = np.sin(alpha), np.cos(alpha)
    ideal = c**(2*k) + s**(2*k)
    erroneous = s*s*c*c*(s**(2*k - 4) + c**(2*k - 4))
    # atan form is stable at zero, unlike the equivalent arcsin(.../sqrt(perror)).
    error_angle = -np.sign(phi)*np.arctan(np.tan(alpha)**(k - 2))
    return alpha, ideal, erroneous, error_angle


def _survival(p, count):
    return np.exp(count*np.log1p(-p))


def undetected_probability(k, hw, sc):
    """Leading-order prior for the error state, before postselection.

    For m=1, Toshio Sec. III A / Appendix C requires a circuit-dependent estimate.
    We instead expose an effective-exposure approximation: a CNOT's data-qubit
    marginal has Z probability 4 eps_2q/15; one-qubit depolarizing noise has
    Z probability eps_1q/3. These exposures stand for initialization, rotation,
    and syndrome operations. This is not a calibrated or propagated fault count.
    m=2 retains the previous CNOT/Rz coefficient and one-qubit extension.
    """
    if sc.rotation_weight == 1:
        return sc.single_qubit_undetected_scale*k*(
            sc.single_qubit_2q_locations*4*hw.eps_2q/15 +
            sc.single_qubit_1q_locations*hw.eps_1q/3)
    if sc.rotation_weight == 2:
        return k*(2*hw.eps_2q/15 + hw.eps_1q/3)
    raise ValueError('rotation_weight must be 1 or 2')


def preparation_survival(dp, hw, sc):
    """Legacy common-survival surrogate; counted acceptance is angle dependent."""
    if sc.state_model == 'counted':
        raise ValueError('Counted acceptance has no common survival factor; use counted_quantities')
    if sc.postselection == 'full':
        n2 = sc.full_2q_per_distance_squared*dp*dp
        n1 = sc.full_1q_per_distance_squared*dp*dp
    elif sc.postselection == 'strip':
        n2 = sc.strip_2q_per_distance*dp
        n1 = sc.strip_1q_per_distance*dp
    else:
        raise ValueError("postselection must be 'full' or 'strip'")
    if min(n1, n2, sc.rejection_scale) < 0:
        raise ValueError('Rejection counts and scale must be nonnegative')
    return _survival(hw.eps_2q, sc.rejection_scale*n2)*_survival(
        hw.eps_1q, sc.rejection_scale*n1)


def fallback_parameters(hw, d, sc):
    """Original-STAR acceptance and retained resource error for main switching.

    Counted mode uses A_orig and W_1Z/A_orig from the declared circuit. Legacy
    mode keeps the 9 max(p,q)/15 error and product-yield surrogate, including
    their historical use for auxiliary corrections.
    """
    if sc.state_model == 'counted':
        if hw.eps_1q == hw.eps_2q == 0:
            return 1., 0.
        row = original_data(hw, d)
        return row['acceptance'], row['conditional_error']
    error = 9/15*max(hw.eps_1q, hw.eps_2q)
    acceptance = _survival(hw.eps_2q, sc.rejection_scale*
                          sc.fallback_2q_per_distance_squared*d*d)
    acceptance *= _survival(hw.eps_1q, sc.rejection_scale*
                            sc.fallback_1q_per_distance_squared*d*d)
    return float(acceptance), float(error)


def counted_quantities(phi, k, hw):
    """Acceptance, retained bad-state weight, omitted nonzero-sector weight.

T1 includes one nonzero signature in the known one-group sector, with arbitrary
zero-signature background faults. Other accepted logical effects are idealized.
"""
    alpha, a0, ae, _ = ideal_tmr(phi, k)
    if hw.eps_1q == hw.eps_2q == 0:
        return a0, np.zeros_like(a0), np.zeros_like(a0)
    row = transversal_data(hw, k)
    s, c = np.sin(alpha), np.cos(alpha)
    intrinsic = np.array([s**(2*n)*c**(2*(k-n)) + s**(2*(k-n))*c**(2*n)
                          for n in range((k+1)//2)])
    weights = np.asarray(row['H_n'])
    acceptance = np.tensordot(weights, intrinsic, axes=(0, 0))
    numerator = ae*row['T1']
    omitted = (np.tensordot(weights[1:], intrinsic[1:], axes=(0, 0))-numerator)/acceptance
    return acceptance, numerator/acceptance, np.maximum(omitted, 0.)


def transversal_correction_resources(phi, k, hw, sc):
    """Expected preparation attempts and states for ideal auxiliary RUS.

Use the actual -Delta and doubled/folded RUS angles, the same k, no auxiliary
switching and no recursive coherent cancellation. Add a positive tail bound.
"""
    if not 1 <= sc.correction_rounds <= 48:
        raise ValueError('correction_rounds must be between 1 and 48')
    phi = fold_angle(phi)
    _, _, _, bad_angle = ideal_tmr(phi, k)
    angle = fold_angle(-(bad_angle-phi))
    attempts, states = np.zeros_like(phi), np.zeros_like(phi)
    active = np.ones_like(phi, dtype=bool)
    for j in range(sc.correction_rounds):
        active &= np.abs(angle) > 1e-14
        weight = 2.**(-j)*active
        attempts += weight/counted_quantities(angle, k, hw)[0]
        states += weight
        angle = fold_angle(2*angle)
    h0 = 1. if hw.eps_1q == hw.eps_2q == 0 else transversal_data(hw, k)['H_n'][0]
    tail = 2.**(1-sc.correction_rounds)
    attempt_tail = tail/(h0*2.**(1-k))
    return attempts + active*attempt_tail, states + active*tail, attempt_tail


@dataclass(frozen=True)
class StateProfile:
    phi: np.ndarray
    acceptance: np.ndarray
    error: np.ndarray
    correction_probability: np.ndarray
    fallback: np.ndarray
    conditional_bad: np.ndarray
    delta: np.ndarray
    preparation_distance: int
    preparation_floor_loss: float
    fallback_acceptance: float
    fallback_error: float

    @property
    def mean_step_loss(self):
        return float(np.mean(-np.log1p(-self.error)) + self.preparation_floor_loss)


def state_profile(phi, k, d, hw, ch, sc):
    """Conditional delivered-state model and corrected teleportation error.

    Counted m=2 uses A=sum(a_n H_n), b1=a_e T1/A and ideal coherent cancellation
    with dephasing 2 b1(1-b1) sin(delta)^2. Other accepted logical effects are
    idealized, or omitted nonzero sectors are charged in the stress scenario.
    Legacy mode retains common-survival conditioning and noisy original-STAR
    corrections. m=1 remains on its explicit, uncalibrated exposure model.
    """
    if hw.connectivity != 'grid':
        raise ValueError('This STAR implementation is restricted to a 2D grid')
    if not (0 <= hw.eps_1q < 1 and 0 <= hw.eps_2q < 1):
        raise ValueError('Physical error probabilities must be in [0, 1)')
    if sc.rejection_scale < 0 or sc.preparation_floor_scale < 0:
        raise ValueError('Noise scales must be nonnegative')
    if sc.state_model not in ('legacy', 'counted'):
        raise ValueError('Unknown STAR state_model')
    if sc.state_model == 'counted':
        validate_counted(hw, sc)
    phi = fold_angle(phi)
    _, ideal, erroneous, error_angle = ideal_tmr(phi, k)
    dp = sc.rotation_weight*k
    if min(sc.single_qubit_2q_locations, sc.single_qubit_1q_locations,
           sc.single_qubit_undetected_scale) < 0:
        raise ValueError('Single-qubit protocol noise assumptions must be nonnegative')
    undetected = undetected_probability(k, hw, sc)
    if not 0 <= undetected < .1:
        raise ValueError('Physical noise is outside the leading-order model regime')
    delta = error_angle - phi
    af, pf = fallback_parameters(hw, d, sc)
    if sc.state_model == 'counted':
        acceptance, bad, omitted = counted_quantities(phi, k, hw)
        error = 2*bad*(1-bad)*np.sin(delta)**2
        if sc.omitted_sector_stress:
            error = error + omitted
    else:
        keep = preparation_survival(dp, hw, sc)
        norm = (1 - undetected)*ideal + undetected*erroneous
        acceptance = keep*norm
        bad = undetected*erroneous/norm
        # Legacy: noisy original-STAR auxiliary rotations.
        error = 2*bad*(1-bad)*np.sin(delta)**2 + bad*2*pf/(1+pf)
    use_fallback = (error > pf) & sc.protocol_switching
    error = np.where(use_fallback, pf, error)
    acceptance = np.where(use_fallback, af, acceptance)
    correction = np.where(use_fallback, 0., bad)
    # m=1 has odd d_p; only even d_p needs the next lower odd distance.
    fit_distance = dp if dp % 2 else dp - 1
    floor_loss = sc.preparation_floor_scale*sc.preparation_cycles*(-np.log1p(-p_sc(fit_distance, hw, ch)))
    return StateProfile(phi, acceptance, error, correction, use_fallback, bad, delta,
                        dp, float(floor_loss), af, pf)


def rus_rotation_loss(phi, k, d, hw, ch, sc, rounds=40):
    """Expected additive loss over doubled/folded angles, truncated for diagnostics.

    This diagnostic is angle resolved. The QUOPS estimator uses invariance of the
    uniform angle distribution under folding/doubling to integrate the infinite
    sum. The omitted geometric weight is at most 2**(1-rounds).
    """
    if not 1 <= rounds <= 48:
        raise ValueError('RUS diagnostic rounds must be between 1 and 48')
    angle = fold_angle(phi)
    loss = np.zeros_like(angle, dtype=float)
    active = np.ones_like(angle, dtype=bool)
    for j in range(rounds):
        active &= np.abs(angle) > 1e-14
        profile = state_profile(angle, k, d, hw, ch, sc)
        memory = (1 + 2*profile.correction_probability)*profile.preparation_floor_loss
        loss += 2.**(-j)*active*(-np.log1p(-profile.error) + memory)
        angle = fold_angle(2*angle)
    # Exact Clifford angles finish without another resource state.
    return loss
