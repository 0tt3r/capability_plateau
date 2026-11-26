# capability_plateau

Shared calculations for *The Quantum Capability Plateau*. The code compares direct execution,
counted sequential Pauli-based computation (PBC), parallel PBC units, and separate STAR,
Pinnacle, and Frozen Pinnacle resource estimates for QUOPS.

The default surface-code model is `PBC` (`PBC_COLUMNS` is the same model). It absorbs all Clifford
gates, including CNOTs, into Pauli axes and charges sequential magic-state consumption, transformed
final measurements, a counted fast-block layout, cold factory startup, and output buffers.

Read [the methods note](docs/quops_direct_and_surface.md),
[the PBC architecture specification](docs/surface_architectures.md), and
[the parallel-unit model](docs/parallel_pbc_audit.md) for assumptions and validation limits.

## Running

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python scripts/make_curves.py
.venv/bin/python scripts/what_matters.py
.venv/bin/python scripts/make_pinnacle_curves.py
.venv/bin/python scripts/make_frozen_pinnacle_curves.py
.venv/bin/python scripts/compare_surface_architectures.py
.venv/bin/python scripts/audit_surface_architectures.py
.venv/bin/python scripts/make_star_comparison.py
.venv/bin/python scripts/make_quops_comparison.py
.venv/bin/python -m pytest -q
```

`make_curves.py` writes the current PBC baseline to `results/baseline.json` and
`figures/quops_direct_surface.*`. `what_matters.py` recomputes PBC sensitivities in
`results/what_matters.json` and `docs/what_matters_table.md`. The comparison and audit scripts write
`pbc_comparison_v2` and `architecture_audit_v3` output directories, respectively.
The Pinnacle script reads the regenerated baseline. Run it before the Frozen Pinnacle script,
which reads both the baseline and the original Pinnacle results.
`make_pinnacle_curves.py` also computes the Frozen Pinnacle curves with an extended supply menu, the
surface-code curve with the same factory menu, the supply comparison, and the sensitivities; see
[the Pinnacle note](docs/pinnacle.md).
`make_quops_comparison.py` draws the combined capability figure of the draft (`figures/quops_comparison.*`)
from the cached direct, surface-code, and Pinnacle results and the STAR steps in
`results/quops_comparison_star.json`; it recomputes nothing.
The notebook [quops_walkthrough.ipynb](notebooks/quops_walkthrough.ipynb) independently calculates
one direct and one sequential PBC design before constructing the curves.

```python
from capability import surface
from capability.config import BASELINE, HARDWARE, PBC_UNITS

rows = surface.sweep(HARDWARE['H1'], BASELINE)  # counted sequential PBC
winner = surface.best_at(1e6, HARDWARE['H1'], BASELINE, rows)
parallel = surface.sweep(HARDWARE['H1'], BASELINE, PBC_UNITS)
```

Every surface row carries architecture metadata. Refinement rejects missing metadata and removed
architectures rather than interpreting old layered results as PBC. `SurfaceArchitecture()` itself
selects counted sequential PBC. `PBC_DRAFT_BOUND` is an explicitly optimistic throughput bound
using the same counted geometry; it is not a finite-buffer execution schedule.

## Code layout

- `capability/surface.py`: common PBC interface, synthesis fit, and capability refinement.
- `capability/surface_pbc.py`: sequential measurements and finite-buffer factory schedule.
- `capability/surface_units.py`: parallel Bell-pair units, correction storage, and factory interfaces.
- `capability/layouts.py`: counted tiles, connected ancilla region, and boundary-access checks.
- `capability/pbc_compile.py`: small-circuit verification of Clifford elimination and final readouts.
- `capability/direct.py` and `routing.py`: direct execution and physical SWAP routing.
- `capability/star.py` and `star_states.py`: STAR preparation, RUS, and capability model.
- `scripts/fit_synthesis.py` and `check_synthesis_accumulation.py`: synthesis calibration and checks.

## STAR scope

`make_star_comparison.py` compares STAR preparation weights 1 and 2 against direct execution and
the current sequential PBC baseline. STAR retains a separate gate-based resource model; its CNOT
schedule and tile ratio are fields of `StarChoices`, not PBC settings. Its routing assumptions have
not been established by the PBC checks. See [the STAR methods note](docs/star_surface_code_estimation.md).

## Retired results

The previous baseline, plots, and documentation are preserved in Git history at `9e5b00a`.
See the [results index](results/README.md) for the current outputs.
