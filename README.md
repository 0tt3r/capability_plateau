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
.venv/bin/python scripts/package_star_inputs.py --check
.venv/bin/python scripts/make_star_comparison.py
.venv/bin/python -m pytest -q
```

`make_curves.py` writes the current PBC baseline to `results/baseline.json` and
`figures/quops_direct_surface.*`. `what_matters.py` recomputes PBC sensitivities in
`results/what_matters.json` and `docs/what_matters_table.md`. The comparison and audit scripts write
`pbc_comparison_v2` and `architecture_audit_v3` output directories, respectively.
The Pinnacle script reads the regenerated baseline. Run it before the Frozen Pinnacle script,
which reads both the baseline and the original Pinnacle results.
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

The main figure compares direct execution, Clifford+T with sequential PBC, and the adopted
optimistic two-qubit STAR model on H1 and H3. STAR uses circuit-derived acceptance and retained
single-fault errors. It charges preparation, retries, and supply time for ideal auxiliary
corrections. Its gate-based layout assumes two data/routing tiles per logical qubit and parallel
CNOT layers. Routing and adaptive scheduling remain provisional. The historical single-qubit model
and the preparation diagnostics have separate figures. [Comparison figure](figures/star_comparison.png) ·
[Saved designs](results/star_comparison.json)

| Hardware | Physical qubits where STAR strictly leads both alternatives | Maximum STAR volume through one million qubits | Fewest qubits reaching that volume |
|---|---:|---:|---:|
| H1 | 1,395–11,299 | 1,156 | 31,171 |
| H3 | 1,787–11,187 | 3,136 | 30,371 |

Both endpoints in the table are inclusive. The estimator uses a 1 µs cycle, H1/H3 noise rates,
$k=3,5,7,9$, and odd STAR computation distances 3–31. Independent reference data cover 114
adopted and omitted-sector stress designs. [STAR methods](docs/star_surface_code_estimation.md) ·
[Assumption decisions](docs/star_two_qubit_assumption_responses.md) ·
[Circuit-input reproduction](analysis/star_inputs.md)

## Retired results

The previous baseline, plots, and documentation are preserved in Git history at `9e5b00a`.
See the [results index](results/README.md) for the current outputs.
