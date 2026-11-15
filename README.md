# capability_plateau

Shared calculations for *The Quantum Capability Plateau* (Overleaf draft, October 2026). This first
version computes two capability curves for the QUOPS problem class on the four hardware classes H1 to
H4: direct implementation and surface-code implementation. The intent is to agree on definitions and
on what matters before adding approaches, problem classes, or sophistication.

Start with the methods note, [docs/quops_direct_and_surface.md](docs/quops_direct_and_surface.md). It
states how each number is obtained, that is, which parts are formulas, which are compiled and counted,
which are simulated, which come from the literature, and which are declared choices. It also shows which
choices move the curves. The notebook [notebooks/quops_walkthrough.ipynb](notebooks/quops_walkthrough.ipynb)
computes one direct and one surface-code design by hand, checks them against the package, and builds the
curves for one hardware class.

## Layout

```
capability/
  config.py      every input and choice, marked [draft] (from the Overleaf) or [choice] (ours, open)
  quops.py       QUOPS ensembles and rotation compilation
  direct.py      direct implementation: closed form (all-to-all), routed Monte Carlo (2D grid)
  routing.py     SWAP-network routers for the 2D grid
  surface.py     surface-code implementation: layout, schedule, loss, design search
  factories.py   published magic-state factories used as black boxes
  curves.py      capability C(N) and the direct/surface-code crossover
scripts/
  fit_synthesis.py   gridsynth T counts -> results/synthesis_fit.json
  make_curves.py     baseline curves -> results/baseline.json, figures/quops_direct_surface.*
  what_matters.py    one-at-a-time choice changes -> results/what_matters.json, docs/what_matters_table.md
  check_synthesis_accumulation.py   state-vector check of how synthesis errors add up in QUOPS circuits
notebooks/
  quops_walkthrough.ipynb   one width through each model by hand, checked against the package
tests/               checks against brute force and against each model's own constraints
```

## Running

```
python -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python scripts/fit_synthesis.py
.venv/bin/python scripts/make_curves.py
.venv/bin/python scripts/what_matters.py
.venv/bin/python scripts/check_synthesis_accumulation.py
.venv/bin/python -m pytest -q
```

To try a different choice, change its default in `capability/config.py` (or pass it through
`config.vary(...)`), rerun `make_curves.py`, and compare. Adding an implementation approach means adding
a module with a `sweep(hardware, choices)` function that returns rows with `s`, `qubits`, and `passes`;
`curves.py` does the rest.

## STAR preparation comparison on a 2D grid

The main comparison adds the adopted optimistic two-qubit STAR model to the
existing direct and Clifford+T implementations on H1 and H3:

```bash
.venv/bin/python -B scripts/make_star_comparison.py
.venv/bin/python -B scripts/package_star_inputs.py --check
.venv/bin/python -B -m pytest -q
```

The comparison is saved in [results/star_comparison.json](results/star_comparison.json)
and [figures/star_comparison.png](figures/star_comparison.png)
([PDF](figures/star_comparison.pdf)). The main figure shows direct, surface-code,
and two-qubit STAR curves. The script also generates preparation and historical
single-qubit diagnostics; those are separate from the main comparison.

`STAR_BASELINE` uses circuit-derived acceptance, retained single-fault logical
errors, ideal auxiliary coherent corrections with supply costs, and the shared
mean-runtime week criterion. `STAR_LEGACY` and `STAR_SINGLE_QUBIT` retain their
historical channel and Markov-deadline choices. Non-STAR implementations and
configuration are unchanged.

| Hardware | Adopted maximum volume | Minimum physical qubits | Stress-test minimum qubits at the same volume |
|---|---:|---:|---:|
| H1 | 1,156 | 31,171 | 73,123 |
| H3 | 3,136 | 30,371 | 119,675 |

The stress test charges omitted accepted nonzero Toshio sectors as failures;
other unmodeled logical effects remain idealized. It is not a complete-channel
bound. Results cover at most one million physical qubits, H1/H3 rates,
$k=3,5,7,9$, and odd computation distances 3–31. Unsupported counted inputs
raise an error rather than being extrapolated.

- [Methods and implementation](docs/star_surface_code_estimation.md).
- [Status and treatment of all seven assumptions](docs/star_two_qubit_assumption_responses.md).
- [Original-STAR channel and coherent corrections](docs/original_star_channel_and_correction.md).
- [Circuit inputs and reproduction](analysis/star_inputs.md).

Assumptions 1–5 and 7 are implemented for the declared optimistic comparison.
Assumption 6 is partial: shared routing, injection, delivery, and frame conventions
are implemented; STAR-specific growth, storage, and adaptive scheduling remain
provisional. Aligning the distance search with the surface-code baseline also
remains open. The main integration is checked against 114 independently saved
designs in [the regression fixture](tests/fixtures/star_reference.json).
