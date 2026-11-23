# Result status

Current entry points produce:

- `baseline.json`: counted sequential PBC and direct execution.
- `what_matters.json`: sensitivities of that PBC baseline.
- `pinnacle.json` and `pinnacle_sensitivity.json`: serial Pinnacle and its comparison with the PBC baseline.
- `frozen_pinnacle.json` and `frozen_pinnacle_sensitivity.json`: Frozen Pinnacle phase catalysis, with the PBC and serial Pinnacle comparisons.
- `star_comparison.json`: separate STAR estimator compared with current PBC and direct curves.
- `pbc_comparison_v2/`: supported PBC architectures and an explicitly optimistic bound.
- `architecture_audit_v3/`: current layout, circuit, and resource checks.
- `synthesis_fit.json` and `synthesis_accumulation.json`: synthesis inputs and finite validation.

Git commit `9e5b00a` preserves the retired layered model, its baseline, and its figures.
