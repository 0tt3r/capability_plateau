# Original STAR: rejection locations, survival, and capability impact

12 November 2026. This calculation concerns **Akahoshi et al.'s original STAR
injection**, using a four-qubit subsystem seed followed by expansion into a
rotated distance-$d$ surface-code patch. It is separate from the
[Toshio two-qubit transversal calculation](star_rejection_count_derivation.md),
whose preparation distance is $d_p=2k$. Here $d$ is the expanded patch distance;
there is no $k$-dependent intrinsic projection penalty.

**15 November follow-up:** the [accepted-state and correction analysis](original_star_channel_and_correction.md)
now evaluates this circuit's acceptance to all fault orders through a Fourier
expectation, with recorded numerical integration uncertainty. It also derives
the retained logical-Z component with leading coefficient $4p/15+q/3$.
That calculation supersedes the product yield below in the adopted optimistic
comparison. The circuit, single-fault counts, and dated results here remain
the record of the 12 November calculation.

For the fully noisy schedule specified below, direct single-fault counting gives

$$
C_p(d)=\frac{179}{15}d^2+O(d),\qquad
C_q(d)=\frac{22}{3}d^2+O(d).
$$

The leading exponents replacing the original-STAR guesses $8d^2,4d^2$ are
therefore **$11.9333d^2$ and $7.3333d^2$**, including preparation, readout, and
idle errors. These are Pauli-probability-weighted rejection coefficients, not
raw numbers of gates. The actual first-order statement is

$$
h_{\rm orig}=1-C_pp-C_qq+O\!\left(d^4(p+q)^2\right).
$$

The historical first-order product surrogate was

$$
\boxed{h_{\rm orig,prod}=(1-p)^{C_p(d)}(1-q)^{C_q(d)}}.
$$

**$h$ is acceptance; the rejection probability is $1-h$.** Unlike transversal
multirotation preparation, no additional angle-dependent $a_0$ factor is needed.
The product matches the first-order expansion but is neither an exact finite-noise
yield nor a proved bound. At the distances below, $C_pp+C_qq$ can exceed one;
the linear expansion must not be used numerically as a probability. The product
values are not the adopted finite-noise inputs; the follow-up now includes multiple-fault acceptance.

## 1. Protocol, noise, and the final measurement convention

The protocol reference is [Akahoshi et al., arXiv:2303.13181, Sec. IV B,
Figs. 9–14](https://arxiv.org/html/2303.13181v1#S4.SS2), published in
PRX Quantum 5, 010337 (2024). It verifies a four-qubit subsystem seed
twice, expands it, and verifies the full patch twice. Random initial gauge or
new-check outcomes serve as references; only prescribed parities cause rejection.
The paper's simulation includes idle noise and uses an ideal final syndrome round
([Sec. VI](https://arxiv.org/html/2303.13181v1#S6)).

Our main count makes **both expansion rounds noisy**, matching the physical
operation convention used in the Toshio audit. For comparison, omitting faults
in the last round, while retaining its error-detection action, gives

$$
C_p^{\rm ideal\ final}=\frac{20}{3}d^2+O(d),\qquad
C_q^{\rm ideal\ final}=4d^2+O(d).
$$

This diagnostic reproduces the ideal-final-round convention, not the paper's
entire simulator or its choice of native gates. It is not our hardware estimate.
Replacing a native ZZ gate by CNOT–$R_Z$–CNOT and routing the seed changes
constant terms only; the expansion determines the $d^2$ terms.

The hardware convention is unchanged:

| Operation | Rate | Error conditional on a fault |
|---|---:|---|
| CNOT | $p$ | Each of the 15 nonidentity two-qubit Paulis equally likely |
| Reset to $|0\rangle$, Z readout | $p$ | Preparation flip, classical readout flip |
| H, $R_Z$, one-tick idle | $q$ | X, Y, Z equally likely |

H1 has $p=q=10^{-3}$; H3 has $p=10^{-3},q=10^{-6}$.
The variable called `eps_2q` in the estimator also covers preparation and readout.

## 2. Direct derivation of the leading coefficients

The patch contains $d^2$ data qubits and $d^2-1$ check ancillas. Each expansion
round uses reset, H on X ancillas, four CNOT layers, H on X ancillas, and readout:
eight ticks. Fresh data resets and H gates run alongside the first round's
corresponding ancilla operations. An inactive live site has one idle location
per tick. Fresh patch qubits are initialized only when expansion begins.

Away from the seed and boundaries, the freshly added data are $|+\rangle$.
Initially X checks are known, while Z checks acquire random reference outcomes.
For every allowed single fault, propagate its Pauli frame to the tested
measurement parities. A nonzero parity signature rejects. A zero signature
means only that this acceptance test does not reject it.

The following table is a direct count in a repeating bulk cell, normalized per
$d^2$. Boundaries and the finite seed contribute $O(d)$ and $O(1)$ terms.

| Locations | First expansion round | Second expansion round | Total coefficient |
|---|---:|---:|---:|
| CNOT, at $p$ | $11/3$ | $49/15$ | $104/15$ |
| Data/ancilla preparation, at $p$ | $2$ | $1$ | $3$ |
| Readout, at $p$ | $1$ | $1$ | $2$ |
| H, at $q$ | $4/3$ | $2/3$ | $2$ |
| Idle, at $q$ | $8/3$ | $8/3$ | $16/3$ |

For clarity, the CNOT entries do not assume every CNOT fault rejects:

- Round 1 has $d^2$ locations with rejecting fraction $12/15$, $2d^2$ with
  fraction $14/15$, and $d^2$ with fraction $15/15$. Their sum is $(11/3)d^2$.
- Round 2 has $d^2$ locations of each fraction $8/15,12/15,14/15,15/15$.
  Their sum is $(49/15)d^2$.
- Round 1 has $d^2$ bulk idle locations with fraction $2/3$ and $2d^2$ with
  fraction one. Round 2 has the same rejecting contributions plus $2d^2$
  idle locations with zero rejection signature.
- Every bulk preparation/readout flip rejects. H faults reject with fraction
  $2/3$, at $2d^2$ H locations in round 1 and $d^2$ in round 2.

Consequently,

$$
C_p/d^2=104/15+3+2=179/15,\qquad
C_q/d^2=2+16/3=22/3.
$$

Removing the second column's noise gives $20/3$ and $4$ for the ideal-final
diagnostic. Excluding idle noise altogether would give $C_q=2d^2+O(d)$ in the
fully noisy circuit. Omitting preparation/readout would give
$C_p=(104/15)d^2+O(d)$; that is not the repository's noise convention.

These coefficients were not fitted to simulated acceptance probabilities. The
enumerator uses exact rational weights, and the bulk-cell count independently
checks the leading terms obtained from finite-distance counts.

## 3. Finite-distance counts and reproducible implementation

For the explicit implementation below, at odd $d\ge5$,

$$
\boxed{C_p(d)=\frac{358d^2-237d+1281}{30}},\qquad
\boxed{C_q(d)=\frac{22d^2+35d+258}{3}}.
$$

The ideal-final diagnostic instead has

$$
C_p^{\rm ideal\ final}(d)=\frac{200d^2-121d+1333}{30},\qquad
C_q^{\rm ideal\ final}(d)=\frac{12d^2+19d+273}{3}.
$$

Distance 3 is counted separately because the seed meets the opposite boundary:
$(C_p,C_q)=(1897/15,560/3)$, or $(462/5,437/3)$ with the last round ideal.
The lower-order terms depend on routing, boundary geometry, and timing. They
are counts for this implementation, not universal constants of the protocol.

The precise construction used by
[`count_rejection.py`](../analysis/original_star_rejection_20261112/count_rejection.py) is:

1. Data coordinates are $(2x,2y)$ for $0\le x,y<d$. A bulk check at
   $(2x+1,2y+1)$ is X for even $x+y$, Z otherwise. Keep Z boundary checks on
   the left/right and X boundary checks on the top/bottom. Each check connects
   to its diagonal neighboring data sites.
2. Seed data $0,1,2,3$ are at $(0,0),(0,2),(2,0),(2,2)$. Gauge ancillas
   $M_0,M_1,M_2,M_3$ are at $(1,-1),(-1,1),(3,1),(1,3)$.
   Initialize these eight sites. H on 1 and 3, then CNOTs $1\to0$ and $3\to2$,
   prepare two Bell pairs. Route each CNOT by swapping the upper data into
   $M_1$ or $M_2$, applying the CNOT, and swapping back. The pair runs in parallel.
3. Implement the logical ZZ rotation on 0 and 2 with
   SWAP$(0,M_0)$–CNOT$(M_0,2)$–$R_Z(2)$–CNOT$(M_0,2)$–SWAP$(0,M_0)$.
   Every SWAP is three CNOTs. This is an explicitly declared nearest-neighbor
   compilation, not a claim that every scheduling detail is fixed by Fig. 13.
4. Run two seven-tick gauge rounds. Reset four ancillas; H on $M_0,M_3$;
   CNOTs $M_0\to0,M_3\to1$; CNOTs $M_0\to2,M_3\to3$;
   CNOTs $0\to M_1,2\to M_2$ alongside H on $M_0,M_3$;
   CNOTs $1\to M_1,3\to M_2$ alongside readout of $M_0,M_3$;
   then read out $M_1,M_2$. Both gauge products must be +1 in each round;
   individual gauge outcomes need not be +1.
5. Retire the extra ancilla $M_0$. For expansion, initialize remaining data in
   the first two rows to $|0\rangle$, and all lower rows to $|+\rangle$.
   Run two eight-tick full-patch rounds. Measured from ancilla to data, X-check
   CNOT order is NW, NE, SW, SE; Z-check order is NW, SW, NE, SE.
   These orders avoid simultaneous gate conflicts and put hook pairs transverse
   to the corresponding logical strings. There are 48 ticks in the whole attempt.

After the second gauge round, the known seed stabilizers are $X_0X_1X_2X_3$
and the two vertical Z pairs, whose signs are the recorded Z-gauge outcomes.
Together with the fresh single-qubit stabilizers, they determine which first
expansion outcomes are predictable. The code finds these by binary elimination
and retains their recorded signs. It tests those first outcomes and the
round-to-round difference of every full-patch check. Random first outcomes
alone never reject. There are $(d^2-2d+5)/2$ initially predictable checks and
$(3d^2-2d+11)/2$ tested parities including the four gauge tests.

The only non-Clifford branch has zero rejection signature: it acts as logical Z.
Thus a single fault's rejection classification is independent of the angle.
We use the Clifford skeleton for parity propagation, not for replacing the
accepted non-Clifford state by a stabilizer state.

Validation includes all odd distances 3, 5, 7, 9, 11, 13, 15, 17, 21, 25, 31;
8,177 independent forward/backward primitive fault comparisons at 3, 5, 9;
deterministic noiseless parity checks; a zero-signature logical rotation branch;
nearest-neighbor connectivity and simultaneous-gate collision checks; and a
separate 16-data-site bulk-cell count. Reflected valid hook schedules give the
same $d^2$ coefficients, with different lower-order terms. Running the counter
writes detailed per-distance counts and location diagnostics beside the source.

## 4. Current finite-noise acceptance and reproduction

The main model uses the [all-orders acceptance and retained logical component](original_star_channel_and_correction.md),
not the historical product surrogate. The acceptance table records exact
summation at distance 3 and numerical integration errors at larger distances.

The [counting script](../analysis/original_star_rejection_20261112/count_rejection.py)
reproduces this circuit's single-fault checks. The
[channel calculation](../analysis/original_star_channel_20261115/channel.py)
extends it with final diagnostic stabilizers/logical observables and evaluates
finite-noise acceptance. Follow [the input reproduction instructions](../analysis/star_inputs.md)
to regenerate and package the production table.

These acceptance calculations do not determine the decoder outcome of every
accepted residual syndrome. That logical-channel uncertainty is explicitly
idealized in the adopted model. Factory timing, footprint, and adaptive
scheduling remain separate unresolved questions.
