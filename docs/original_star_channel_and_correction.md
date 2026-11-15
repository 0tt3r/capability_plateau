# Original STAR and probabilistic coherent cancellation

15 November 2026. This resolves [assumption 5](star_two_qubit_assumptions.md#5-original-star-for-protocol-switching-and-coherent-corrections)
for the deliberately optimistic STAR comparison. The treatment is implemented
in the main estimator; see the [methods and results](star_surface_code_estimation.md)
and [circuit-input reproduction](../analysis/star_inputs.md).

The adopted treatment has two parts:

1. For original STAR, calculate preparation acceptance from the specified noisy
   circuit, retain its identified single-fault logical-Z mechanism, and idealize
   other accepted logical effects. This follows the treatment already adopted
   for Toshio preparation: acceptance and the retained error numerator come
   from the same circuit, without claiming to calibrate its entire output channel.
2. Implement the additional coherent-cancellation rotations through Toshio
   transversal preparation and RUS. Grant those rotations ideal logical output,
   while retaining their angle-dependent preparation yield, rejected attempts,
   RUS repetitions, runtime, and live-patch memory cost.

Original STAR remains available for **main protocol switching**. Setting its
error to zero would make it a perfect implementation of every requested angle;
that is not the adopted treatment.

## Original-STAR acceptance and retained logical error

Use the [existing four-qubit subsystem seed and two noisy expansion rounds](original_star_rejection_count_derivation.md).
The physical noise model assigns probability $p$ to CNOT, initialization, and
readout faults, and $q$ to single-qubit gates and idles. Extend the monitored
syndrome signatures with final ideal stabilizer checks and logical observables.
The final checks are diagnostic: they are not extra physical postselection rounds.

An isolated fault is classified as rejected, accepted with a residual syndrome,
or accepted without a residual syndrome. In the last class, enumeration finds
only an ideal output or a logical-Z error. The nontrivial cases are:

| Location in the compiled seed rotation | Faults producing the retained logical Z | Leading probability |
|---|---|---:|
| First CNOT | IZ, XZ | $2p/15$ |
| Physical RZ | Z | $q/3$ |
| Second CNOT | YY, ZZ | $2p/15$ |

Pauli labels follow the control/target order of the counted circuit. The result is

$$
p_f=\frac{4p}{15}+\frac q3+O((p+q)^2)
$$

at every enumerated odd distance $3\le d\le31$. It recovers $9p/15$ on H1,
where $p=q=10^{-3}$. On H3, where $q=10^{-6}$, the leading error is
$2.67\times10^{-4}$ instead of the previous $6\times10^{-4}$.
These are coefficients for this circuit, not universal original-STAR constants.

Faults leaving a residual output syndrome are **idealized as recovered**; their
actual decoder outcome has not been calibrated. Other accepted logical effects,
including additional mechanisms involving multiple faults, are also idealized.
This is an explicit modeling choice favorable to STAR, not a proof that the next
correction cycle always recovers every accepted physical fault.

### Finite-noise normalization

Let $s_\ell$ be the random monitored signature at fault location $\ell$. The
acceptance calculation includes every physical fault order:

$$
A_{\rm orig}
=\Pr\!\left[\bigoplus_\ell s_\ell=0\right]
=2^{-m}\sum_{x\in\mathbb F_2^m}
       \prod_\ell\mathbb E\left[(-1)^{x\cdot s_\ell}\right].
$$

Here $m$ is the number of monitored detector bits. This Fourier expression is
evaluated exhaustively at $d=3$ and by sampling its positive integrand at larger
distances. Thus it sums all fault orders in the specified signature model, with
a recorded numerical integration error; it is not a product extrapolation of
first-order rejection counts or a full logical-channel simulation.

Define $z_\ell=\Pr(s_\ell=0)$ and $\Phi_0=\prod_\ell z_\ell$.
At the three gadget locations in the table, let $w_\ell$ be the probability of
the listed logical-Z faults and $r_\ell=w_\ell/z_\ell$. The retained numerator is

$$
W_{1Z}=\Phi_0\sum_{\ell\in G}r_\ell
                         \prod_{j\in G\setminus\{\ell\}}(1-r_j),
\qquad p_f=\frac{W_{1Z}}{A_{\rm orig}}.
$$

$G$ contains those three locations. This event contains exactly one listed
logical-Z fault, no nonzero monitored signatures at other locations, and any
number of other zero-signature faults. Their additional logical effects are
idealized. Accepted histories with canceling nonzero signatures contribute to
acceptance but receive ideal output. Consequently this $p_f$ is a **retained
component**, not the complete physical conditional error probability.

At the adopted plateau designs:

| Hardware and distance | $A_{\rm orig}$ | Rejection | Retained $p_f$ |
|---|---:|---:|---:|
| H1, $d=15$ | 0.010980 | 98.902% | $5.9394\times10^{-4}$ |
| H3, $d=11$ | 0.246748 | 75.325% | $2.6662\times10^{-4}$ |

These high rejection probabilities remain charged through preparation retries.
In this idealized model $p_f$ decreases slightly with distance because the
retained numerator is diluted by other accepted histories assigned ideal output.
The optimizer therefore evaluates every distance's own error; it cannot use
the old distance-3 resource-error bound to prune larger distances.

## Ideal correction channel with physical supply costs

[Toshio et al., Sec. IV C](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.15.021057)
construct probabilistic coherent cancellation using a transversal RUS rotation
for the corrective branch. Their correction is itself imperfect; idealizing its
output is the additional favorable assumption adopted here. Sec. IV D discusses
switching the main preparation protocol when original STAR has lower error.

In the repository's $e^{i\phi Z}$ convention, the retained error channel and
adopted correction are

$$
\mathcal E=(1-b_1)\mathcal I+b_1\mathcal U_\Delta,
\qquad
\mathcal C=(1-b_1)\mathcal I+b_1\mathcal U_{-\Delta}.
$$

Their composition is exactly a dephasing channel with

$$
\boxed{p_{\rm tr}=2b_1(1-b_1)\sin^2\Delta.}
$$

Use $\sin^2(\Delta/2)$ when $\Delta$ instead denotes an angle in the manuscript's
$R_Z(\theta)=e^{-i\theta Z/2}$ convention. The old
$b_1\,2p_f/(1+p_f)$ correction-error term is removed. The remaining dephasing
still contributes to the circuit's injection-error budget. An ideal auxiliary
rotation does **not** make the corrected main channel perfect.

For a correction angle $\psi=-\Delta(\phi,k)$, use Toshio acceptance
$A_T(\operatorname{fold}(2^j\psi),k)$ at every RUS round. The selected $k$ is
the same as for the main preparation; correction branches do not switch to
original STAR and do not recursively request coherent cancellation. Until a
Clifford residual terminates the sequence, their expected cost is

$$
T_{\rm corr}(\psi)
=\sum_{j\ge0}2^{-j}
 \left[\frac{4+d}{A_T(\operatorname{fold}(2^j\psi),k)}+d\right].
$$

The current $4+d$ preparation/growth allowance and $d$-cycle injection are
retained. Free delivery and software frame tracking match the surface-code
baseline. Corrections remain serialized under the provisional schedule.
For a main batch of size $B$, their added time is
$2B\langle b_1(\phi)T_{\rm corr}(-\Delta(\phi,k))\rangle$, with $b_1=0$ on
main angles switched to original STAR. Keeping this angle correlation matters;
the cost is not the product of two separately averaged quantities.

## Reoptimized capability

The calculation combines the agreed optimistic treatments and shared
surface-code conventions. It searches both workload families, $N\le10^6$,
odd $d=3,\ldots,31$, $k\in\{3,5,7,9\}$, and $1\le F\le w$.

| Scenario | H1 plateau | H1 minimum qubits | H3 plateau | H3 minimum qubits |
|---|---:|---:|---:|---:|
| Previous 13 November model | 1,156 | 31,171 | 2,744 | 7,139 |
| Ideal transversal corrections only; previous original-STAR inputs | 1,156 | 31,171 | 2,744 | 7,139 |
| Adopted original-STAR inputs and ideal transversal corrections | **1,156** | **31,171** | **3,136** | **30,371** |

H3's plateau rises by 14.3%. This gain requires the revised original-STAR
inputs; ideal corrections alone do not change either plateau. The new H3
maximum is a width-56 shallow circuit, whereas the old maximum was a width-14
deep circuit. Thus the qubit costs in the last two columns reach different
workloads. At $N=10,000$, H3 still reaches 2,744, now with 6,517 qubits.

At the adopted maxima, original STAR supplies 46.9% of consumed resource states
on H1 and 28.9% on H3, all for main protocol switching. Auxiliary corrections
consume only about 0.447 and 0.583 accepted states per execution in expectation;
their preparation attempts and time are still included.

The omitted-Toshio-sector stress test is rerun with these new original-STAR and
correction inputs. It still idealizes unmodeled original-STAR errors and
zero-sector logical effects. Its plateau and resource costs are reported in
the [main results](../results/star_comparison.json).
It is a sensitivity, not a bound on the complete unknown physical channel.
Recomputing error-based protocol switching also changes supply cost, so this
test need not lower capability at every individual qubit budget.

## Scope and checks

The circuit scripts check fault propagation; the main tests check channel
composition, RUS costs, optimizer pruning, and regression against independently
saved designs. Numerical integration uncertainties are retained in the
[source acceptance table](../analysis/original_star_channel_20261115/acceptance.json).
No additional software installation is required for this calculation.

Assumption 5 is settled **for this optimistic comparison**. A calibrated physical
prediction would still require the omitted accepted logical channels and decoder
outcomes. STAR-specific growth, storage, and scheduling in assumption 6 remain
provisional. The finite search menu and unmatched computation-distance cap also
remain scope limitations. This calculation concerns compiled CNOT–RZ–CNOT
preparation; it does not calibrate the separate native-RZZ scenario.

Suggested paper wording:

> For original STAR, we retain the identified single-fault logical errors and
> idealize other accepted logical effects; for probabilistic coherent cancellation,
> we assume exact auxiliary rotations supplied by Toshio's transversal protocol,
> while accounting for rejected preparation attempts, repeat-until-success
> repetitions, and their runtime and memory costs.
