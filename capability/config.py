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


def vary(**changes):
    return replace(BASELINE, **changes)
