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

## Status

First pass with the simplest reasonable choices. Open questions for the group are collected in
Sec. 7 of the methods note.
