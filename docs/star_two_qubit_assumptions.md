# Assumptions requiring validation in the two-qubit STAR estimate

This note describes the uncertain modeling choices added to the capability estimator
for Toshio et al.'s **transversal multirotation protocol with two-qubit rotations**
on a two-dimensional surface-code grid. It records the implementation as of
27 October 2026. These assumptions are active in the reported estimates.

The single-qubit transversal variant is outside this note's scope. Accurate analog
angle control is taken as an accepted premise here. The main unresolved issues are
the errors and acceptance probabilities of prepared states, their protection during
code deformation, and the conversion of preparation costs into runtime.

The highest priorities for validation are the added preparation-error contribution,
the separation of rejection from conditional state error, and the completeness of
the physical fault model. All three can affect the predicted capability plateau.

## Model and notation

The workload alternates random single-qubit rotation layers with random-pairing
CNOT layers. For even logical width $w$, circuit volume is $s=wL$, with $L=w$
for shallow circuits and $L=w^2$ for deep circuits. There are $R=s/2$ requested
logical rotations. Capability at physical-qubit budget $N$ is the largest $s$
passing the model's error and deadline criteria. A plateau is the largest volume
reached as $N$ increases over the studied range; it is not a throughput measure.

| Symbol | Meaning in this note |
|---|---|
| $p,q$ | Two-qubit and single-qubit physical fault probabilities, respectively. |
| H1 | Grid hardware with $p=q=10^{-3}$. |
| H3 | Grid hardware with $p=10^{-3}$ and $q=10^{-6}$. |
| $k$ | Number of physical two-qubit rotations in one preparation attempt. |
| $d_p=2k$ | Surface-code distance used for transversal resource-state preparation. |
| $d$ | Surface-code distance used for computation and resource-state storage. |
| $F$ | Number of parallel resource-state preparation patches; code parameter `stations`. |
| $\phi$ | Requested residual logical angle in the convention $e^{i\phi P}$, for a Pauli operator $P$. |
| $a_0,a_e$ | Ideal postselection acceptance without a fault, and acceptance conditional on the modeled relevant fault. |
| $u$ | Total probability of that class of physical faults before postselection; already includes the factor $k$. |
| $b$ | Probability of the erroneous resource-state component, conditional on acceptance. |
| $\Delta$ | Angle difference between the erroneous resource state and the target. |
| $h$ | Additional assumed survival probability against rejection. |
| $p_f,a_f$ | Error and preparation acceptance for the original-STAR implementation used for switching and corrections. |

The repository angle convention is $e^{-i\theta P/2}$, so $\phi=-\theta/2$ before
removing Clifford multiples of $\pi/4$. Uniform workload angles give uniform
residual $\phi\in[-\pi/8,\pi/8)$. The estimator uses repeat-until-success (RUS)
teleportation, with success probability $1/2$ and doubled correction angles.

In code, `fallback` means choosing **original STAR preparation** for the requested
angle; the paper's term for choosing between preparation methods is **protocol
switching**. It does not denote an additional preparation failure. The estimator
also uses original STAR for the extra rotations in probabilistic coherent error
cancellation. No observable-level probabilistic error cancellation or associated
sampling overhead is included.

The analytic preparation formulas and ideal coherent-cancellation construction
come from [Toshio et al., Phys. Rev. X 15, 021057, Secs. III–IV](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.15.021057).
For positive $\phi$, the implemented formulas are

$$
\tan\alpha=(\tan\phi)^{1/k},\qquad
a_0=\cos^{2k}\alpha+\sin^{2k}\alpha,
$$

$$
a_e=\sin^2\alpha\cos^{2k-2}\alpha
    +\sin^{2k-2}\alpha\cos^2\alpha,\qquad
\Delta=-\arctan[(\tan\alpha)^{k-2}]-\phi.
$$

$a_e$ does not include the probability of a physical fault. Also, $b$ is a mixture
weight, not an infidelity: its contribution to resource-state infidelity is
$b\sin^2\Delta$ in this model.

**Code:** [angle conversion and projection formulas](../capability/star_states.py#L14),
[workload definitions](../capability/quops.py#L14),
[hardware settings](../capability/config.py#L15), and
[`StarChoices` defaults](../capability/config.py#L70).

## 1. Additional logical errors during accepted preparation

Each accepted resource state is charged an additional negative-log no-fault
contribution

$$
\ell_p=\gamma\,4[-\ln(1-p_{\rm SC}(d_p-1))],\qquad \gamma=1.
$$

The inherited surface-code memory fit is

$$
p_{\rm SC}(d)=0.03(p/p_{\rm th})^{(d+1)/2},
\qquad p_{\rm th}=0.01\text{ (H1) or }0.015\text{ (H3)}.
$$

The two-qubit preparation distances are even; the implementation substitutes the
next lower odd distance into this fit. The four-cycle duration does not establish
that this memory-error fit describes the conditional error after postselection.
Neither the distance substitution nor the overall magnitude is calibrated against
a preparation circuit and decoder.

The term could overlap with faults already represented by $u$, or omit faults not
represented there. There is no derivation partitioning the two contributions.
The code applies the same $d_p$-dependent contribution to accepted original-STAR
states used for switching and corrections, although their preparation differs.

For the uniform-angle workload, the total preparation contribution is
$R(2+4\bar c)\ell_p$, where $\bar c$ is the angle-averaged probability of requesting
an extra coherent-cancellation rotation. The factors count the modeled mean
number of accepted states, including RUS repetitions and corrections.

This assumption directly changes the plateau. Existing sweeps, with all other
settings fixed and designs reoptimized, give:

| Preparation multiplier $\gamma$ | H1 maximum volume | H3 maximum volume |
|---:|---:|---:|
| 0 | 1,156 | 2,744 |
| 1, default | 1,024 | 2,304 |
| 2 | 1,024 | 1,936 |

These are sensitivity scenarios, not uncertainty intervals. They use $N\le10^6$
and the design menu described below.

**Code:** [`state_profile`, `fit_distance` and `floor_loss`](../capability/star_states.py#L165);
[`rotation_losses`](../capability/star.py#L126);
[`p_sc`](../capability/surface.py#L48). The setting is `preparation_floor_scale`.
Saved results are under `hardware.H1/H3.preparation_floor_sensitivity` in
[`results/star_comparison.json`](../results/star_comparison.json).

**Needed validation:** determine the logical error of accepted preparations with
syndrome extraction and decoding, and account separately for subsequent growth
and storage. This would establish which errors remain after postselection and
avoid counting the same mechanism twice.

## 2. Separating rejection from conditional state quality

The estimator assumes

$$
a=h[(1-u)a_0+u a_e],\qquad
b=\frac{u a_e}{(1-u)a_0+u a_e}.
$$

The same $h$ multiplies both ideal and erroneous outcomes. It therefore cancels
when conditioning on acceptance. Changing the rejection coefficients changes
throughput but leaves $b$, $\Delta$, and the modeled error per delivered rotation
unchanged for fixed preparation settings.

That factorization has not been checked for the noisy syndrome circuit.
Rejection can depend on which fault occurred and on correlations between data and
measurement errors. The finite-noise normalization above is an estimator choice;
retaining it exactly does not make it an exact circuit-level noise model.

This is a structural uncertainty, beyond the numerical values of the rejection
coefficients. The weak plateau sensitivity to those coefficients partly follows
from assuming that conditional state quality is unaffected by them.

**Code:** [`state_profile`, `keep`, `norm`, `acceptance`, and `bad`](../capability/star_states.py#L149).

**Needed validation:** calculate acceptance and accepted-state error jointly for
each relevant physical fault class, using the actual postselection rule.

## 3. Assumed rejection counts

For restricted postselection, the default is

$$
h=(1-p)^{24d_p\eta}(1-q)^{8d_p\eta},\qquad \eta=1.
$$

For the original-STAR preparation used in switching and corrections, it is

$$
a_f=(1-p)^{8d^2\eta}(1-q)^{4d^2\eta}.
$$

The constants 24, 8, 8, and 4 are assumed effective counts of fault opportunities.
They were not obtained by counting rejection-causing faults in a compiled circuit
or fitting preparation-yield data. In particular, making $a_f$ depend on the
whole distance-$d$ patch was not derived from an initialization and deformation
schedule. These formulas are not established upper or lower bounds on acceptance.

More rejection increases retries, runtime, and accumulated storage errors. It can
reduce capability at fixed $N$ and can make a design fail the deadline criterion.

**Code:** [`preparation_survival`](../capability/star_states.py#L73),
[`fallback_parameters`](../capability/star_states.py#L89), and
[`StarChoices`](../capability/config.py#L70). Relevant fields are
`strip_2q_per_distance`, `strip_1q_per_distance`,
`fallback_2q_per_distance_squared`, `fallback_1q_per_distance_squared`, and
`rejection_scale`.

### Sensitivity of only the restricted-postselection counts

Scaling only $(24,8)$ by a common factor $f$, keeping $a_f$ and all other settings
fixed, gave the following reoptimized results on 27 October 2026:

| Factor $f$ | H1 maximum volume | H1 qubits to reach it | H3 maximum volume | H3 qubits to reach it |
|---:|---:|---:|---:|---:|
| 0.5 | 1,024 | 22,331 | 2,304 | 29,040 |
| 1 | 1,024 | 22,331 | 2,304 | 31,581 |
| 2 | 1,024 | 22,331 | 2,304 | 31,944 |
| 4 | 1,024 | 22,331 | 2,304 | 32,907 |
| 10 | 1,024 | 40,965 | 2,304 | 32,907 |

Unchanged maxima do not imply unchanged capability curves. At $N=1,000$, H1
volume falls from 64 to 36 when $f$ increases from 1 to 2. At $N=10,000$, it falls
from 1,000 to 512 when $f$ increases from 1 to 10. Holding the baseline designs
fixed, doubling the counts increases mean runtime by about 11% on H1 and 12% on H3.

The existing `rejection_sensitivity` entries in
[`results/star_comparison.json`](../results/star_comparison.json) instead vary
`rejection_scale`, which also changes $a_f$. For example, that broader doubling
requires 44,000 H1 qubits to reach volume 1,024, whereas doubling only the
restricted-postselection counts still requires 22,331. The existing dotted
rejection-sensitivity curves describe the broader change.

The isolated table can be reproduced from the repository root without changing
the saved figures or JSON:

```bash
.venv/bin/python - <<'PY'
from dataclasses import replace
from capability import star
from capability.config import BASELINE, HARDWARE, N_MAX, STAR_BASELINE

for hardware in ('H1', 'H3'):
    for factor in (0.5, 1.0, 2.0, 4.0, 10.0):
        settings = replace(
            STAR_BASELINE,
            strip_2q_per_distance=24.0 * factor,
            strip_1q_per_distance=8.0 * factor,
        )
        rows = star.sweep(HARDWARE[hardware], BASELINE, settings)
        best = star.best_at(N_MAX, rows)
        small = star.best_at(1_000, rows)
        medium = star.best_at(10_000, rows)
        print(hardware, factor, best['s'], best['qubits'],
              small['s'] if small else None,
              medium['s'] if medium else None)
PY
```

**Needed validation:** obtain acceptance versus angle, preparation distance, and
physical error rates from the specified syndrome circuit. Varying the assumed
coefficients measures sensitivity but does not calibrate them.

## 4. Extending the physical fault formula and omitting higher orders

The code assumes independent stochastic depolarizing faults and uses

$$
u=k\left(\frac{2p}{15}+\frac q3\right).
$$

The $2kp/15$ contribution follows the CNOT/Rz implementation with ideal virtual
Z rotations discussed in [Toshio et al., Sec. III B](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.15.021057).
The added $kq/3$ models the Z component of a physical Rz fault propagating to a ZZ
error. The complete circuit with the repository's separate one- and two-qubit
error rates has not been analyzed to establish that all other first-order faults
are detected or corrected.

This extension matters especially on H1, where $q=p$. Multiple physical faults,
their correlations with rejection, and all resulting accepted-state errors are
not systematically enumerated. The check $u<0.1$ is a domain guard, not a bound on
the omitted terms. Nonlinear normalization and the cancellation expression below
do not supply the missing higher-order fault mechanisms.

**Code:** [`undetected_probability`](../capability/star_states.py#L54),
[`state_profile`](../capability/star_states.py#L124), and
[`Hardware`](../capability/config.py#L15).

**Needed validation:** enumerate single faults in the complete preparation
circuit for separate $p$ and $q$, then estimate the leading multiple-fault
corrections over the angles and distances used by the optimizer.

## 5. Original STAR for protocol switching and coherent corrections

The chosen original-STAR error model is

$$
p_f=\frac9{15}\max(p,q).
$$

The coefficient $9/15$ comes from the particular CNOT/Rz implementation in
[Akahoshi et al., Fig. 14, upper circuit](https://arxiv.org/pdf/2303.13181).
It is not a universal error rate for original STAR. Retaining it on H3 by replacing
all rates with their maximum is an intended conservative approximation, without
a separate derivation for unequal rates or a bound on all higher-order errors.

For a transversal preparation, ideal probabilistic coherent error cancellation
gives $p_c=2b(1-b)\sin^2\Delta$. The implementation additionally charges errors
from performing the correction with original STAR:

$$
p_{\rm rot}=2b(1-b)\sin^2\Delta
            +b\frac{2p_f}{1+p_f}.
$$

Here $2p_f/(1+p_f)$ is the probability of at least one faulty resource in a
geometric RUS sequence under the assumed independent constant-$p_f$ model.
The correction is requested with probability $b$. The code switches the requested
rotation to original STAR when $p_{\rm rot}>p_f$.

The ideal cancellation identity is mathematically established for the assumed
mixture channels. The implementation-specific assumptions are the original-STAR
variant, always using it for additional corrections, its yield $a_f$, and the
combination with the separate preparation-error charge in section 1. The code
does not propagate the full imperfect correction channel through the circuit.

**Code:** [`fallback_parameters`](../capability/star_states.py#L89),
[`state_profile`, correction error and switching](../capability/star_states.py#L157),
and [`_batch_times`, correction runtime](../capability/star.py#L114).
`protocol_switching=False` disables switching for requested rotations only;
additional coherent-cancellation rotations still use original STAR.

**Needed validation:** specify a consistent original-STAR circuit, calculate its
error and yield for H1/H3, and evaluate its full preparation, growth, and correction
costs without overlapping error charges.

## 6. Growth, delivery, layout, and scheduling

The estimator charges each preparation attempt $4+d$ code cycles, including the
$d$-cycle growth allowance even when the attempt is rejected. Consumption adds
$3d$ cycles: $d$ each for delivery, injection, and Clifford operations.

The approximately four-cycle attempt comes from
[Toshio et al., Sec. III B 3](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.15.021057).
The added growth and consumption allowances are estimator choices. The paper
allows expansion to begin during the second stabilizer measurement; the present
implementation does not compile that procedure or its overlap with preparation.

The physical-qubit count per copy is

$$
N_{\rm copy}=2d^2(2w+F)+F(2d_p-1)^2.
$$

This adds an unrotated preparation patch and a distance-$d$ storage/transfer
patch for each of the $F$ parallel preparations. All distance-$d$ patches accrue
the inherited logical-error fit throughout runtime. No spatial layout or
deformation circuit verifies that these patches and timing allowances suffice,
or establishes the error rate while the distance changes. The code searches $d$
and $d_p$ separately and does not enforce $d\ge d_p$; an implementation would also
need to specify the transition for any selected design with $d<d_p$.

The chosen schedule synchronizes unfinished rotations in each batch, waits for
their preparations, and starts the next batch after the current one finishes.
Additional coherent-cancellation rotations are serialized. These are explicit
scheduling restrictions that can increase overhead relative to a schedule that
reuses free patches immediately. Their mathematical timing calculation does not
validate the assumed physical layout.

**Code:** [`_batch_times`](../capability/star.py#L114),
[`_one_rotation_cycles`](../capability/star.py#L133),
[`evaluate_design`](../capability/star.py#L141), and
[`StarChoices`](../capability/config.py#L70). Settings are `growth_steps`,
`delivery_steps`, `injection_steps`, `clifford_steps`, and `preparation_cycles`.

**Needed validation:** construct a patch arrangement and syndrome schedule for
preparation, deformation, and teleportation; derive the duration, overlap, and
logical-error exposure from it.

## 7. A deliberately loose deadline bound

The deadline calculation sets $\delta=10^{-6}$ per execution and uses Markov's
inequality:

$$
T_{\rm bound}=\frac{\mathbb E[T]}{\delta},\qquad
\Pr(T>T_{\rm bound})\le\delta.
$$

The inequality is valid for the modeled nonnegative runtime. The concern is its
looseness, not the mathematics. With one execution and a one-week budget, this
criterion requires a mean runtime no greater than $0.6048$ seconds per copy,
even though the available wall time is 604,800 seconds.

The estimator combines its deadline bound with an error-model bound:

$$
q_{\rm lower}=\max(0,e^{-\mu}-\delta),\qquad
\mu=\mathbb E[L_{\rm path}],
$$

where $L_{\rm path}$ is the modeled negative logarithm of no-fault probability
for a realized schedule. The pass threshold is $q_{\rm lower}\ge e^{-1/2}$.
The value $\delta$ limits deadline overruns; it is not a total logical-error
target of $10^{-6}$. Jensen's inequality justifies the error bound within the
assumed conditional stochastic model, not independently of its calibration.

**Code:** [`evaluate_design`](../capability/star.py#L141), with settings
`deadline_failure`, `executions`, `cycle_seconds`, and `WEEK_SECONDS` in
[`config.py`](../capability/config.py).

**Needed improvement:** calculate or tightly bound the runtime distribution for
the selected preparation/RUS schedule and use its $1-\delta$ quantile. This would
test whether designs are excluded mainly by the Markov bound.

## Search limits, inherited assumptions, and existing checks

The default search uses $k\in\{3,5,7,9\}$, hence $d_p\in\{6,10,14,18\}$,
odd $d$ from 3 through 31, and $1\le F\le w$. Only odd $k$ is implemented to avoid
the additional even-$k$ Clifford case. A reported optimum is over this menu, not
all possible STAR implementations. Increasing angle quadrature from 512 to 1,024
points previously left the selected two-qubit designs unchanged; numerical
integration accuracy does not establish physical-model accuracy.

**Code:** [`StarChoices`](../capability/config.py#L70),
[`ideal_tmr`](../capability/star_states.py#L31),
[`min_qubits`](../capability/star.py#L200), and
[`sweep`](../capability/star.py#L244).

The memory-error fit, approximate $2d^2$ qubits per rotated patch, $2w$ data/routing
patches, default $2d$-cycle parallel CNOT layers without congestion, one-microsecond
cycle time, and uniform-angle workload were inherited from the original estimator.
Their reuse does not establish that the new preparation and deformation circuits
have the same error behavior.

The checks in [`tests/test_star.py`](../tests/test_star.py) cover projection formulas,
ideal channel composition, geometric waiting-time identities, and optimization
against a small exhaustive search. They do not constitute a surface-code
preparation/decoder simulation or validate the assumed rejection counts and
postselected logical-error contributions.

The [implementation report](star_surface_code_estimation.md) gives the broader
comparison. The [comparison script](../scripts/make_star_comparison.py#L32)
generates the saved scenarios in
[`results/star_comparison.json`](../results/star_comparison.json). The isolated
rejection-count calculation in section 3 is separate from those broader scenarios.
