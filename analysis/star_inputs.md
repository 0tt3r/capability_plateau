# Reproducing the counted STAR inputs

The main estimator reads [star_acceptance.json](../capability/data/star_acceptance.json).
This directory retains the circuit calculations and source data needed to audit
and regenerate that table. They compute acceptance and retained logical-error
components, not complete physical logical channels.

| Source | Role |
|---|---|
| [Toshio circuit and propagation](rejection_count_20261111/count_rejection.py) | Native gate schedule, monitored signatures, and independent forward checks |
| [Original-STAR circuit](original_star_rejection_20261112/count_rejection.py) | Four-qubit seed, noisy expansion rounds, and monitored parities |
| [Fourier transfer](separate_branch_survival_20261113/branch_survival.py) | All-orders Toshio sector weights $H_n$ |
| [Fault-sector weights](separate_branch_survival_20261113/two_fault_terms.py) | Zero-signature background and the retained exactly-one-active weight $T_1$ |
| [Original-STAR channel](original_star_channel_20261115/channel.py) | Accepted single-fault classification and all-orders acceptance integration |
| [Distance-21 refinement](original_star_channel_20261115/refine_acceptance.py) | Additional independent Fourier batches near a resource threshold |

The source data are [Toshio sector weights](separate_branch_survival_20261113/results.json),
[single/two-fault weights](separate_branch_survival_20261113/two_fault_terms.json),
and [original-STAR acceptance](original_star_channel_20261115/acceptance.json).
The second file retains historical diagnostics; the current model uses its
`exactly_one_active` and `no_active_fault_probability` fields. It does not adopt
the earlier two-fault logical-channel prescription in that file's metadata.

From the repository root, verify packaging and run the estimator:

```bash
.venv/bin/python -B scripts/package_star_inputs.py --check
.venv/bin/python -B -m pytest -q
.venv/bin/python -B scripts/make_star_comparison.py
```

To regenerate circuit inputs before packaging them:

```bash
.venv/bin/python -B analysis/separate_branch_survival_20261113/branch_survival.py
.venv/bin/python -B analysis/separate_branch_survival_20261113/two_fault_terms.py
.venv/bin/python -B analysis/original_star_channel_20261115/channel.py
.venv/bin/python -B analysis/original_star_channel_20261115/refine_acceptance.py
.venv/bin/python -B scripts/package_star_inputs.py
```

These commands write generated diagnostics beside the source scripts. Only the
three source JSON tables are inputs to packaging. The estimator does not import
analysis scripts at runtime. The counts and Fourier routines use the existing
NumPy environment; no additional simulator is needed.

Original-STAR distance 3 is summed exhaustively. Other distances normally use
262,144 Fourier samples with seed 15112026; distance 21 pools 24 batches, with
seeds incremented by 10,000. Standard errors in the table quantify numerical
integration, not physical-model uncertainty. Supported rates are H1/H3, with
odd computation distances 3–31 and $k=3,5,7,9$.

The [regression fixture](../tests/fixtures/star_reference.json) contains the 114
adopted/stress designs from the independent calculation before main integration.
It records its source checksum and is not regenerated from the production
estimator. See [the methods](../docs/star_surface_code_estimation.md) and
[assumption decisions](../docs/star_two_qubit_assumption_responses.md) for scope.
