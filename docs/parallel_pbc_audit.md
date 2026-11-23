# Parallel PBC units

The active comparison contains counted sequential PBC, a labeled throughput bound, and the
open-chain parallel-unit model of [Litinski, Section 5](https://arxiv.org/html/1808.02892v3).
The unsupported layered Clifford+T estimator and its area-only controls have been removed.
See [surface_architectures.md](surface_architectures.md) for the default sequential model.

## Resource model

Each unit prepares one layer of $w$ commuting Pauli rotations on Bell-pair halves. A chain of
$U$ units supplies $U-1$ layers per wave. The default menu has $U=2,\ldots,8$. Per unit we charge:

- $(2\lceil\sqrt w\rceil+1)^2$ data tiles, including square padding;
- $4w$ correction storage/access tiles;
- $3FB$ factory-interface tiles and $F$ factories;
- factory batch latency $B\tau+2d$ for conversion to post-corrected states;
- serial consumption of $w$ states and $(\lceil\sqrt w\rceil+3)d$ Bell-pair overhead.

Endpoint units are conservatively charged as full interior units. For $D=L_Rn_T$ input T layers,
$q=\lfloor(D-1)/(U-1)\rfloor$, $r=D-q(U-1)$, unit preparation time $t_u$, and reaction time $t_r$,

$$C=t_u+q\max(t_u,(U-1)t_r)+rt_r+wd.$$

All allocated tiles remain live for this duration. The loss and execution constraints are the
same as for sequential PBC. The pipeline is compared against independent wave scheduling in the
tests. The factory interface, correction storage factor, and reaction time remain explicit model
assumptions; no full physical factory/decoder layout has been compiled.

## Interpretation

Parallel units can shorten latency at fixed logical width, but allocate additional storage and
factories. In the tested one-million-qubit cases they achieve lower capability than sequential
PBC. This is evidence about the implemented unit model and finite menu, not a no-go result for
parallel PBC. Mutually commuting Pauli products can overlap on the data block, so their algebraic
commutation alone does not justify parallel surgery there.

## Reproduction

```bash
.venv/bin/python scripts/audit_surface_architectures.py
.venv/bin/python -m pytest -q
```

The audit writes `results/architecture_audit_v3/` and `figures/architecture_audit_v3/`. It checks tile
geometry, synthesized small-circuit Clifford elimination, the architecture menu, an expanded unit
menu, multiple executions, and fixed-design latency/space tradeoffs. STAR comparisons are labeled
as STAR's separate gate-based model; the PBC audit does not validate STAR routing.
Previous v2 outputs are historical provenance, not the active model menu. This code audit does not
write to the manuscript or Overleaf.
