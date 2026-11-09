"""Every modeling input and choice, in one place.

Inputs fixed by the draft (Overleaf, Secs. 2.1, 2.2, 3.2, 3.5) are marked [draft]. Choices the draft
leaves open are marked [choice]; their defaults are the simplest reasonable option, and
scripts/what_matters.py changes them one at a time to show which ones move the curves.
"""
from dataclasses import dataclass, replace

P = 1e-3                       # [draft] physical error parameter, fixed throughout
WEEK_SECONDS = 7*24*3600       # [draft] time budget
N_MAX = 1e6                    # [choice] largest machine considered; curves and probes stop here


@dataclass(frozen=True)
class Hardware:
    """[draft] Sec. 2.1. eps_1q: H, Z(theta), and idle locations. eps_2q: CNOT, preparation, measurement."""
    label: str
    name: str
    connectivity: str          # 'grid' (2D nearest neighbor) or 'all_to_all'
    noise: str                 # 'uniform' or 'bimodal'
    eps_1q: float
    eps_2q: float
    p: float = P


HARDWARE = {
    'H1': Hardware('H1', '2D grid, uniform', 'grid', 'uniform', P, P),
    'H2': Hardware('H2', 'all-to-all, uniform', 'all_to_all', 'uniform', P, P),
    'H3': Hardware('H3', '2D grid, bimodal', 'grid', 'bimodal', P/1000, P),
    'H4': Hardware('H4', 'all-to-all, bimodal', 'all_to_all', 'bimodal', P/1000, P),
}


@dataclass(frozen=True)
class Choices:
    # ---- shared
    threshold: float = 0.5               # [draft] pass if mean no-fault probability >= exp(-threshold)
    executions: int = 1                  # [choice] runs that must fit in the week; the draft states no shot count
    # ---- direct implementation (Sec. 3.2)
    idle_errors: bool = True             # [choice] each data qubit untouched in a time step is a location at eps_1q
    rotation_compile: str = 'euler'      # [choice] 'euler': R_Z 1, R_X 3 (H Z H), R_Y 5 (S H Z H Sdg) time steps;
                                         #          'native': any R_P(theta) in one step (not in the draft gate set)
    router: str = 'snake'                # [choice] 2D grid only: 'snake' (odd-even transposition along a
                                         #          boustrophedon path) or 'rcr' (row-column-row sorting)
    samples: int = 64                    # Monte Carlo circuits per (width, ensemble) on the grid
    # ---- surface code (Sec. 3.5)
    sc_prefactor: float = 0.03           # [draft] P_SC = A (p/p_th)^((d+1)/2) per live patch-cycle
    sc_threshold_uniform: float = 0.01   # [draft]
    sc_threshold_bimodal: float = 0.015  # [draft]
    qubits_per_patch: str = '2d^2'       # [draft] approximately 2 d^2
    tiles_per_data_patch: float = 2.0    # [choice] one data patch plus one routing patch per logical qubit
    cnot_layer: str = 'parallel'         # [choice] 'parallel': all w/2 long-range CNOTs in 2 surgery steps
                                         #          (no congestion); 'serial': one CNOT at a time
    t_per_injection_steps: int = 1       # [choice] surgery steps (d cycles each) to consume one T state
    synthesis_charge: str = 'squared'    # [choice] each rotation synthesized to diamond-norm tolerance eps counts
                                         #          as eps^2 ('squared', its infidelity) or eps ('linear')
    synthesis_scale: float = 1.0         # [choice] multiplies that charge; scripts/check_synthesis_accumulation.py
                                         #          measures the mean output infidelity of QUOPS circuits
    synthesis_budget: float = None       # [choice] None: synthesis shares the threshold with stochastic failure;
                                         #          a number: stochastic failure gets the whole threshold and the
                                         #          summed synthesis error has this separate budget
    factories: tuple = ('cultivation', 'litinski')   # [choice] which published factory families may be used
    cycle_seconds: float = 1e-6          # [choice] surface-code cycle; only enters the one-week limit


BASELINE = Choices()


@dataclass(frozen=True)
class StarChoices:
    """Toshio STAR on a 2D grid. These are estimator choices, not calibrated device parameters.

    Weight m=1 uses Choi's single-qubit rotations; m=2 uses Toshio's ZZ rotations.
    d_p = m*k; odd k avoids the extra even-k Clifford in the preparation gadget.
    See docs/star_surface_code_estimation.md for the uncalibrated noise assumptions.
    """
    groups: tuple = (3, 5, 7, 9)
    max_distance: int = 31
    angle_points: int = 512
    protocol_switching: bool = True
    # [choice] Effective fault opportunities in the selectively postselected strip.
    rejection_scale: float = 1.0
    strip_2q_per_distance: float = 24.0
    strip_1q_per_distance: float = 8.0
    # [choice] Full-patch rejection surrogate for the original-protocol fallback.
    fallback_2q_per_distance_squared: float = 8.0
    fallback_1q_per_distance_squared: float = 4.0
    # [paper] About four cycles per state preparation attempt (m <= 3).
    preparation_cycles: int = 4
    # [choice] Charge growth on every attempt, including rejected ones (conservative).
    growth_steps: int = 1
    injection_steps: int = 1
    delivery_steps: int = 1
    clifford_steps: int = 1
    # [choice] Background QEC floor during accepted preparation, separately from injection.
    preparation_floor_scale: float = 1.0
    # Markov deadline bound per execution; this probability is charged to success as well.
    deadline_failure: float = 1e-6
    timing_tail_tolerance: float = 1e-9
    rotation_weight: int = 2
    postselection: str = 'strip'         # 'strip' (Toshio) or 'full' (Choi)
    full_2q_per_distance_squared: float = 8.0
    full_1q_per_distance_squared: float = 4.0
    # [choice] Effective physical fault exposures per rotated data qubit for m=1.
    # Eight CNOT exposures (four per round, two rounds), with marginal Z rate 4p/15.
    # Three one-qubit exposures: rotation plus two representing initialization/idling.
    # These are not a propagated syndrome-circuit count; vary the scale to assess them.
    single_qubit_2q_locations: float = 8.0
    single_qubit_1q_locations: float = 3.0
    single_qubit_undetected_scale: float = 1.0


STAR_BASELINE = StarChoices()
STAR_SINGLE_QUBIT = replace(STAR_BASELINE, rotation_weight=1, postselection='full')


def vary(**changes):
    return replace(BASELINE, **changes)
