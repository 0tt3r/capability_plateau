# Assumptions in the two-qubit STAR estimate

This note defines the seven uncertainties in the two-qubit transversal STAR
comparison. The [decision record](star_two_qubit_assumption_responses.md) gives
their adopted treatments and implementation status; the
[methods](star_surface_code_estimation.md) give equations, code, and reproduction.
Assumptions 1–5 and 7 are implemented for the optimistic comparison. Assumption 6
remains partial. This does not claim complete physical channel calibration.

## Notation and scope

| Symbol | Meaning |
|---|---|
| $p,q$ | Fault probabilities for CNOT/preparation/readout and H/rotation/idle, respectively |
| H1 / H3 | $p=10^{-3}$ with $q=10^{-3}$ / $10^{-6}$ |
| $k,d_p=2k$ | Physical two-qubit rotation groups and transversal preparation distance |
| $d,F$ | Computation/storage distance and number of preparation stations |
| $\phi$ | Residual preparation angle in $e^{i\phi P}$; repository angle is $\theta=-2\phi$ |
| $A,H_n$ | Total acceptance and all-orders rotation-sector signature weights |
| $T_1,b_1$ | Retained one-group single-fault weight and its accepted fraction $a_eT_1/A$ |
| $\Delta$ | Known single-fault angle offset |
| $A_{\rm orig},p_f$ | Original-STAR acceptance and retained conditional logical error |
| $\ell_p$ | Additional generic logical-error allowance during accepted preparation |

The counted input table covers compiled CNOT–RZ–CNOT preparation on H1/H3,
$k=3,5,7,9$, and odd computation distances 3–31. It does not calibrate the
single-qubit transversal or native-RZZ alternatives. Accurate analog angle control
is an accepted premise.

## 1. Additional logical errors during accepted preparation

Accepted states can contain background logical faults from physical idles,
syndrome-extraction gates, or erroneous syndrome measurements and subsequent
correction. The old generic preparation allowance was not calibrated to the
specified circuit. The optimistic comparison sets it to zero, retaining explicit
main-state errors and memory exposure on live computation/routing/storage patches.

## 2. Separating rejection from conditional state quality

Rejection probability and accepted-state logical quality are different quantities.
A fault may reject, survive with a correctable syndrome, or survive with a logical
error. The implemented model uses counted acceptance and the known single-fault
logical component, and explicitly idealizes other accepted logical effects.
Counting acceptance alone does not classify the full accepted logical channel.

## 3. Circuit-derived restricted-postselection counts

The survival guessed from a count of possible faults was not a calibrated
finite-noise yield. The main model now uses all-orders acceptance for the
specified noisy circuits, retaining rejected-attempt costs. Original-STAR
acceptance has numerical integration uncertainty, recorded with the table.

## 4. Physical fault formula and higher orders

The retained leading mechanism is $k(2p/15+q/3)$, including physical RZ faults.
Its numerator retains zero-signature background survival without expansion.
Other accepted effects are idealized, and omitted nonzero Toshio sectors are
charged as failures in a separate sensitivity. This is not a proof that all
omitted physical mechanisms can only worsen performance.

## 5. Original STAR for protocol switching and coherent corrections

The [circuit-derived retained error](original_star_channel_and_correction.md)
has leading coefficient $4p/15+q/3$, replacing $9\max(p,q)/15$ for main switching.
Auxiliary coherent-cancellation rotations use ideal Toshio output, while their
actual-angle acceptance, RUS, and supply costs remain. Unresolved original-STAR
accepted faults and decoder outcomes are idealized.

## 6. Growth, delivery, layout, and scheduling

STAR and sequential PBC share the patch-memory law, baseline cycle time, one-step
injection cost, free delivery, and software frames. STAR's $2w$ data/routing tiles
and parallel CNOT layers are optimistic gate-based assumptions. Growth time, handoff/storage footprint,
synchronized batches, and serialized auxiliary corrections remain provisional.
Relating the counted circuits to effective factory costs and allowing feasible
adaptive supply overlap remain open.

## 7. A deliberately loose deadline bound

The old Markov rule substantially inflated STAR runtime relative to the baseline.
The main model now compares mean runtime directly with one week, with no
additional deadline-error charge. This is not a deadline-tail guarantee.

## Search limits and validation

Tests cover channel identities, preparation and RUS timing, independent regression
designs, source-data provenance, small exhaustive searches, and deadline/copy
accounting. The finite menu is not an optimum over all STAR implementations.
In particular, alignment with the surface-code computation-distance search
remains open. The main comparison is scoped to the supplied circuit inputs and
explicit optimistic channel assumptions.
