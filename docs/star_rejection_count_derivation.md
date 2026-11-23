# Two-qubit STAR: deriving rejection counts and their capability impact

11 November 2026. This note derives the additional rejection coefficients for
one explicit implementation of two-qubit transversal STAR preparation,
$d_p=2k$ with odd $k\ge3$. The method is deterministic single-fault enumeration
with exact Pauli propagation. It is not Monte Carlo sampling or a fit to yield
data. The coefficients below were calculated here; they are not quoted from
the source paper.

For the circuit and noise model specified below, the result is

$$
C_p(d_p)=\frac{920d_p-664}{15},\qquad
C_q(d_p)=\frac{455d_p-134}{6}.
$$

Thus the leading terms replacing the guessed exponents $24d_p,8d_p$ are
$(184/3)d_p\simeq61.33d_p$ and $(455/6)d_p\simeq75.83d_p$.
These are **Pauli-probability-weighted rejection coefficients**, not raw
numbers of gates. The first-order survival expansion is

$$
h=1-C_p p-C_q q+O\!\left(d_p^2(p+q)^2\right).
$$

These are first-order rejection coefficients. The main estimator uses the
all-orders acceptance-sector weights derived from the same circuit, not a
product extrapolation of these counts. The [methods](star_surface_code_estimation.md)
and [assumption decisions](star_two_qubit_assumption_responses.md) specify the
retained single-fault logical component and idealization of other accepted effects.
Original STAR uses its [own seed/expansion circuit](original_star_rejection_count_derivation.md).

## 1. Circuit, schedule, and physical noise

The construction follows Algorithm 1, Fig. 5(c), and Fig. 8 of
[Toshio et al., Phys. Rev. X 15, 021057 (2025)](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.15.021057).
Those specify the protocol, rotation gadget, and post-selection region,
but do not fix a unique syndrome-CNOT ordering or idle schedule. Those
choices are part of the result, not universal properties of STAR.

Use an unrotated planar patch with coordinates
$0\le x,y\le2d_p-2$. Data qubits have even $x+y$. Among the remaining
sites, odd $x$ labels X-check ancillas and even $x$ labels Z-check ancillas.
The logical Z operator lies along the top data row, $y=0$.
The dark post-selection band contains X checks at $y=0,2$ and Z checks at
$y=1$: $2d_p-2$ X checks and $d_p$ Z checks.

The scheduled attempt is:

1. Initialize data to $|+\rangle$ with computational-basis reset and H.
   Perform one complete syndrome round to prepare the encoded state.
2. Apply $k$ disjoint ZZ rotations along the top row. For the pair
   $(4j,0),(4j+2,0)$ and intermediate ancilla $(4j+1,0)$, use
   SWAP–CNOT–$R_Z$–CNOT–SWAP, with each SWAP compiled into three CNOTs.
   All pairs run in parallel, giving nine native-operation layers.
3. Perform two further complete syndrome rounds and accept only the
   expected monitored outcomes.

Each syndrome round has eight layers: ancilla reset, H on X ancillas,
four CNOT layers in **north, east, west, south** order, H on X ancillas,
and computational-basis measurement. Direction is measured from ancilla
to data. X ancillas control their CNOTs; Z ancillas are targets. Initial
data resets and H gates run alongside the corresponding operations in
the first syndrome round. The complete attempt therefore has 33 layers.

Every native operation lasts one time step. Unused physical qubits idle
in that step. An idle on an ancilla that is reset before it can affect a
monitored measurement has zero rejection contribution. Likewise, faults
after the last relevant interaction may escape this acceptance test;
they are not automatically harmless to the accepted state.

The noise convention matches the repository hardware classes:

| Operations | Failure probability | Fault distribution |
|---|---:|---|
| CNOT | $p$ | Uniform over 15 nonidentity two-qubit Paulis, after the gate |
| Computational-basis preparation | $p$ | Orthogonal-state preparation, represented by X after reset |
| Computational-basis measurement | $p$ | Classical outcome flip |
| H, $R_Z$, idle | $q$ | Uniform over X, Y, Z after the operation |

H1 has $p=q=10^{-3}$; H3 has $p=10^{-3},q=10^{-6}$.
Despite its code name `eps_2q`, the $p$ rate includes preparation and
measurement. No patch growth, delivery, or post-acceptance correction
rounds are included in this rejection test.

## 2. Which measurement outcomes are tested?

Let $m_{r,a}\in\{0,1\}$ encode the measurement result of check $a$ in
round $r=0,1,2$, with round 0 preceding the rotations. For every monitored
X check, acceptance requires $m_{r,a}=0$ in all three rounds. The initial
Z-check results are random; they are references, not outcomes required to
be zero. For each monitored Z check, test

$$
m_{1,a}\oplus m_{0,a}=0,\qquad
m_{2,a}\oplus m_{0,a}=0.
$$

This gives $3(2d_p-2)+2d_p=8d_p-6$ tested parities. Initial Pauli-frame
corrections can change physical rotation signs, but those signs do not
change the syndrome-branch probabilities used in this acceptance count.

## 3. Enumerating faults without confusing them with intrinsic rejection

For each location $\ell$ and each allowed Pauli fault $P$, propagate its
effect to the tested parities. This gives a binary signature
$s_{\ell,P}$. The calculation uses the Clifford skeleton, replacing each
physical $R_Z$ by the identity while retaining its two expansion branches.
For example, CNOT propagation uses

$$
X_c\mapsto X_cX_t,\quad Z_c\mapsto Z_c,\quad
X_t\mapsto X_t,\quad Z_t\mapsto Z_cZ_t,
$$

and H exchanges X and Z. Resets remove previous errors on that wire;
readout faults flip the corresponding recorded outcome. The code uses a
backward propagation of signatures and checks it against an independent
forward propagation of fault frames.

Each physical rotation has expansion
$e^{i\alpha Z_iZ_j}=cI+isZ_iZ_j$, where
$c=\cos\alpha$ and $s=\sin\alpha$. Let $B_j$ be the monitored signature
of the $Z_iZ_j$ branch of rotation $j$. The $k$ vectors $B_j$ have rank
$k-1$: their full product is logical Z and has zero check signature.
Only a branch and its global complement have the same signature.

For a fault signature $s_{\ell,P}$, exactly one of three cases applies:

| Signature | Acceptance effect |
|---|---|
| Outside $\operatorname{span}\{B_j\}$ | No rotation branch cancels it: guaranteed rejection |
| Zero | The accepted branch weights are unchanged |
| Nonzero inside the span | A different pair of complementary rotation branches passes |

Without a physical fault, the two accepted branches have total weight
$a_0=c^{2k}+s^{2k}$. In this enumeration, every single fault in the third
case maps to a single-group branch or its complement, with weight

$$
a_e=s^2c^{2k-2}+s^{2k-2}c^2.
$$

Their summed prior probability is

$$
u=\frac{d_p}{15}p+\frac{d_p}{6}q
 =\frac{2k}{15}p+\frac{k}{3}q.
$$

This reproduces the existing two-CNOT contribution and the physical-$R_Z$
Z-fault contribution. These faults must remain in $u$: they do not always
reject and must not also be included in $C_p$ or $C_q$.

Let $f_\ell$ be the fraction of allowed fault types at location $\ell$
whose signatures lie outside the branch span. Then

$$
C_p=\sum_{\ell\text{ at rate }p}f_\ell,\qquad
C_q=\sum_{\ell\text{ at rate }q}f_\ell.
$$

For a CNOT, $f_\ell$ is an integer divided by 15; for H, $R_Z$, or idle,
it is an integer divided by 3. Preparation and readout have one fault type.
Consequently the complete first-order yield is

$$
P_{\mathrm{accept}}
=a_0(1-C_pp-C_qq-u)+a_eu+O\!\left(d_p^2(p+q)^2\right).
$$

This establishes acceptance coefficients. It does not establish the full
logical channel of all accepted faults or validate a finite-noise common
survival factor multiplying that channel.

## 4. Location counts and boundary terms

Simply counting CNOTs attached to monitored ancillas gives $11d_p-9$
per round. It misses faults propagated from adjacent checks. For the
specified CNOT order, the locations that can change acceptance are:

| Ancilla row and relevant edges | CNOT locations per round |
|---|---:|
| X checks at $y=0$ | $3(d_p-1)$ |
| Z checks at $y=1$ | $4d_p-2$ |
| X checks at $y=2$ | $4(d_p-1)$ |
| Z checks at $y=3$, north/east/west edges | $3d_p-2$ |
| X checks at $y=4$, north edges | $d_p-1$ |
| Total | $15d_p-12$ |

The last two rows lie outside the monitored band. Fault propagation from
them is why the circuit count exceeds the count of gates attached to the
selected checks.

The weighted CNOT contributions also depend on time, because later faults
have fewer remaining opportunities to trigger a test:

| Syndrome round | Contribution to $C_p$ |
|---|---:|
| Initialization | $12d_p-152/15$ |
| First post-rotation | $(196/15)d_p-32/3$ |
| Second post-rotation | $(164/15)d_p-142/15$ |
| Sum | $36d_p-454/15$ |

The full accounting is:

| Location class | Rate | Leading number of potentially rejection-causing locations | Exact weighted coefficient |
|---|---|---:|---:|
| All syndrome CNOTs | $p$ | $45d_p$ | $36d_p-454/15$ |
| Rotation-gadget CNOTs, including compiled SWAPs | $p$ | $4d_p$ | $(10/3)d_p$ |
| Initial data resets | $p$ | $4d_p$ | $4d_p-2$ |
| Ancilla resets and readouts | $p$ | $18d_p$ | $18d_p-12$ |
| Data and ancilla Hadamards | $q$ | $16d_p$ | $(32/3)d_p-28/3$ |
| Physical $R_Z$ gates | $q$ | $d_p/2$ | $d_p/3$ |
| Idles | $q$ | $(147/2)d_p$ | $(389/6)d_p-13$ |

For reproducibility, the idle coefficient is the sum of
$31d_p/3+2$ in initialization, $55d_p/2-15$ during the rotation gadget,
$53d_p/3-2/3$ in the first post-rotation round, and $28d_p/3+2/3$ in the
last round. Summing the table gives $C_p,C_q$ stated above.
The raw numbers of locations where some fault reduces acceptance are
$71d_p-50$ at rate $p$ and $90d_p-31$ at rate $q$.

## 5. Dependence on schedule and conventions

- Omitting idle noise gives $C_q=11d_p-28/3$.
- If preparation and measurement are ideal and $p$ applies only to CNOTs,
  $C_p=(118/3)d_p-454/15$.
- Treating SWAP as one native two-qubit operation of error rate $p$ and
  duration one tick gives $C_p=(896d_p-664)/15$ and
  $C_q=(379d_p-94)/6$. This differs from the repository's native gate set.

Among common directional orders applied to both check types, eight of the
24 permutations have collision-free layers and deterministic intended checks:

| CNOT orders | Leading $C_p/d_p$ | Leading $C_q/d_p$ |
|---|---:|---:|
| NEWS, NWES | $184/3$ | $455/6$ |
| ENSW, ESNW, WNSE, WSNE | $62$ | $463/6$ |
| SEWN, SWEN | $884/15$ | $153/2$ |

This checks schedule dependence within this family; it does not optimize
all possible schedules or establish universal bounds for STAR.

## 6. Reproduction and finite-noise inputs

The [counting script](../analysis/rejection_count_20261111/count_rejection.py)
contains the circuit, exact fault propagation, independent forward checks,
and schedule/count diagnostics. From the repository root:

```bash
.venv/bin/python -B analysis/rejection_count_20261111/count_rejection.py
```

This generates per-distance counts and location diagnostics. The main finite-noise
acceptance is obtained by [Fourier transfer](../analysis/separate_branch_survival_20261113/branch_survival.py),
with [source weights](../analysis/separate_branch_survival_20261113/results.json).
See [input reproduction](../analysis/star_inputs.md) for the complete calculation
and packaging commands. A zero monitored signature does not establish a harmless
accepted output; full logical-channel classification is outside these counts.
