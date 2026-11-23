# Surface-code PBC architectures

Counted sequential PBC is the default for `surface.min_qubits`, `surface.sweep`, and
`scripts/make_curves.py`. The former layered Clifford+T model has been removed: eliminating
entangling Clifford gates does not preserve single-qubit T axes, so its compilation argument did
not justify local parallel consumption at the assumed area and timing.

## Supported choices

| Configuration | Layout and schedule |
|---|---|
| `PBC`, also named `PBC_COLUMNS` | Counted fast-column layout, sequential measurements, buffered factory batches |
| `PBC_UNITS` | Parallel Bell-pair units, explicit correction storage and factory interfaces |
| `PBC_DRAFT_BOUND` | Same counted fast-column layout; optimistic throughput bound without startup or buffers |

`SurfaceArchitecture()` selects `PBC`. Removed schedules/layouts raise an error. Every generated
surface row records its architecture. `best_at` rejects unlabelled historical rows and mixed
architectures, preventing old results from silently entering the current capability refinement.

## Compilation and layout

Clifford+T synthesis produces $N_T=wn_TL_R$ T gates. Track the Clifford frame $C$ and replace each
T gate on qubit $q$ by a rotation about $C^\dagger Z_q C$. The final readout observables are also
conjugated. These axes generally involve multiple qubits. Commutation within a layer does not
establish simultaneous physical access. Measurement-dependent Clifford corrections update the
subsequent axes. See [Litinski, Sections 1 and 2.3](https://arxiv.org/html/1808.02892v3).

The explicit column layout in `layouts.fast_layout` packs $w/2$ two-qubit patches into columns.
For even $w$, let $r=\lceil\sqrt{w/2}\rceil$ and $c=\lceil(w/2)/r\rceil$. Its computational area is

$$A_{\mathrm{col}}(w)=2w+2(r+c)+1.$$

This includes the connected ancilla region and accessible Pauli boundaries. The often quoted
$2w+\lceil\sqrt{8w}\rceil+1$ agrees at square packings but can omit one tile in this construction
at other widths (21 counted tiles versus 20 at $w=6$). It is no longer a selectable layout.
The geometry checks certify tile occupancy, connectivity, and boundary access; they are not a
physical syndrome or decoder simulation.

## Sequential schedule and loss

Each of $F$ identical factories produces $B$ states per batch after declared latency $L=B\tau$.
Allocate $FB$ output-buffer tiles outside the factory footprint. Factories start cold and are
throttled so each new batch arrives when the previous batch has been consumed, or later if supply
is slower. One measurement takes

$$h=d\,n_{\mathrm{injection}}+t_{\mathrm{delivery}}+t_{\mathrm{feedforward}}.$$

With $q=\lfloor(N_T-1)/(FB)\rfloor$ and $r_T=N_T-qFB$, the run time is

$$C=L+q\max(L,FBh)+r_T h+wd.$$

The last term charges $w$ transformed final measurements. There is no CNOT runtime term.
All allocated non-factory tiles remain live throughout the run:

$$A=A_{\mathrm{col}}+FB,\qquad N_{\mathrm{copy}}=2d^2A+FN_{\mathrm{factory}}.$$

The default loss criterion is

$$AC[-\ln(1-p_{\mathrm{SC}})]+N_T[-\ln(1-p_{\mathrm{factory}})]+wL_R\epsilon^2\le\tfrac12.$$

Here $p_{\mathrm{SC}}=0.03(p/p_{\mathrm{th}})^{(d+1)/2}$. The synthesis term is a modeling
interpretation supported by finite small-circuit checks, not a general error-accumulation theorem.
Linear synthesis charging and a separate synthesis budget remain sensitivity choices.
For $M$ requested executions, copies satisfy $K\lfloor W/(C t_{\mathrm{cycle}})\rfloor\ge M$.

The optimizer searches odd distances 3 through 99, T counts 5 through 160, the declared factory
menu, and all factory counts up to consumption saturation. Increasing $F$ beyond
$\lceil\tau/h\rceil$ cannot improve this serial schedule and only adds area. Results are minima
within this finite resource model. The throughput-bound diagnostic instead uses
$\max(N_Th+wd,\lceil N_T/F\rceil\tau)$ and omits buffers/startup; it must not be called an executable
supply schedule.

## Checks and limitations

Tests compare the batch formula with event-by-event consumption, optimization with independent
unpruned searches, copy counts with the exact deadline constraint, and ideal Clifford elimination
with full small-circuit unitaries and final readout probabilities. Entry-point tests protect the
PBC default and rejection of retired results.

Factory periods are treated as deterministic batch latencies. Random rejection stalls, physical
factory-to-block transfer, and decoding/feed-forward costs still need calibration. The default
extra delivery/feed-forward costs are zero and can be varied explicitly. The effective patch-cycle
error law is applied to the entire counted volume; it is not calibrated to every surgery geometry.
Both connectivities use the same grid-compatible surface architecture, so equal surface curves do
not establish that connectivity is irrelevant.

Use `scripts/compare_surface_architectures.py` for the current PBC comparison and
`scripts/audit_surface_architectures.py` for layout, circuit, and parallel-unit checks. Earlier
layered results are archived as described in [the results index](../results/README.md).
