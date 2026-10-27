# Plan: STAR resource estimation on a 2D surface-code grid

*Planning only. No estimator implementation has been changed.*

## Scope

Add two implementations—original STAR and Toshio's improved STAR—and compare them with the existing
Clifford+T estimate on **H1 and H3**, the repository's two 2D-grid hardware classes. Both STAR variants
use surface-code patches and lattice surgery. Here, “transversal rotations” refers to Toshio et al.'s
transversal multirotation resource-state preparation, not the separate neutral-atom Transversal STAR
architecture.

The main change is replacing synthesis and T-state supply with **angle-dependent resource-state
preparation, repeat-until-success execution, and residual rotation errors**. The existing surface-code
model provides much of the surrounding framework, but its rotation scheduler and error accounting
need substantial changes.

Sources:

- [Original STAR: Akahoshi et al., arXiv:2303.13181](https://arxiv.org/abs/2303.13181).
- [Improved STAR: Toshio et al., Practical Quantum Advantage on Partially Fault-Tolerant Quantum
  Computer, Phys. Rev. X 15, 021057 (2025)](https://doi.org/10.1103/PhysRevX.15.021057).
  The [arXiv preprint is 2408.14848](https://arxiv.org/abs/2408.14848); equation references below refer
  to the published paper where applicable.

## 1. Replace synthesized rotations with STAR teleportation

Currently, [`surface.py`](../capability/surface.py) assigns each rotation a T count and synthesis
tolerance. STAR would instead request a resource state for the desired angle and consume it through
lattice surgery. An unsuccessful teleportation requires another state at twice the angle; the mean
consumption is two accepted states per generic rotation. See original STAR, Sec. IV.

The STAR branches would remove `T_per_rotation`, the gridsynth fit, and `loss_synthesis`. They would
track requested angles, resource states, retries, and any physical rotation-control error. The
existing Clifford+T implementation remains available as the comparison baseline.

An explicit angle conversion is necessary: the repository uses

$$
R_P(\theta)=e^{-i\theta P/2},
$$

whereas Toshio uses $e^{i\phi P}$, so $\phi=-\theta/2$. This matters for every angle-dependent formula.

## 2. Introduce separate preparation models for the two papers

| Variant | Preparation model needed |
|---|---|
| Original STAR | Prepare in the $[[4,1,1,2]]$ subsystem code, postselect, and expand into a surface-code patch. |
| Toshio STAR | Prepare a logical Clifford state, apply physical multirotations along a logical operator, and use selective postselection plus error correction. |

The original paper's $2p/15+O(p^2)$ resource-state error assumes a native two-qubit rotation. With
CNOTs and single-qubit rotations, its alternative circuits give $9p/15$, or $7p/15$ with an additional
ancilla, to leading order under its uniform noise model. Use a circuit compatible with the
repository's gate set, including its local SWAP overhead. See original STAR, Sec. IV D and
Figs. 13–14.

Toshio introduces choices of rotation weight, number of physical rotations, preparation distance,
and postselection region. Acceptance and output noise depend on those choices and the target angle.
The paper also switches preparation protocols when doubled RUS angles become large. See Toshio
et al., Secs. III–IV.

The replacement for the current fixed `Factory(error, qubits, cycles_per_output)` would return:

- Acceptance probability and conditional output error/channel.
- Footprint and attempt duration, including rejected attempts.
- Expansion, storage, and delivery costs.

H1 and H3 need separate calibration using `eps_1q` and `eps_2q`; the present practice of sharing
factory errors across noise models would be inappropriate here. Published coefficients must be
matched to their physical gate and noise assumptions before use.

## 3. Model preparation retries separately from teleportation retries

These are different costs. A rejected preparation consumes time and workspace without touching the
data. A failed teleportation changes the required correction angle.

For an on-demand preparation policy, a useful accounting identity is

$$
\mathbb E[N_{\mathrm{prep\ attempts}}(\theta)]
=\sum_{j=0}^{\infty}\frac{2^{-j}}{a(\theta_j)},
$$

where $a(\theta_j)$ is preparation acceptance and $\theta_j$ is the doubled angle, reduced by
Clifford equivalences. This assumes generic rotations and ideal equiprobable teleportation outcomes;
exact Clifford termination should be handled explicitly. Speculative preparation would additionally
charge unused states.

**Two teleportations per rotation does not mean two teleportation rounds per layer.** QUOPS rotates
all $w$ qubits in parallel. With independent retries and a layer barrier, the slowest sequence grows
roughly as $\log_2 w$, even with unlimited state supply. Finite preparation capacity adds stalls.

The current deterministic rotation-layer timing therefore needs a stochastic schedule or a validated
approximation, including idle patch errors and a controlled treatment of the retry tail.

## 4. Make angles explicit in the workload

[`quops.py`](../capability/quops.py) defines uniformly random rotations, but the current surface
estimator only needs their count.

For STAR, integrate over—or sample—the actual angle distribution. After extracting Clifford
rotations, the residual angles remain broadly distributed: in Toshio's convention they span
$[-\pi/8,\pi/8]$. Account for the extracted Clifford operations through the chosen implementation
and frame-tracking rules.

Toshio's corrected, RUS-averaged rotation error has small-angle behavior $O(|\phi|p)$; the favorable
single-state infidelity scaling cannot simply be substituted for the error of a completed rotation.
See Toshio et al., Sec. IV C.

The working expectation is that **QUOPS will benefit less from the small-angle improvement than the
Hamiltonian-simulation examples**. This is an inference to test, not an estimated result. Evaluate
the full angle dependence and protocol switching rather than assigning every QUOPS rotation a
representative small angle.

## 5. Replace the error budget and account for preparation space

The current model uses

$$
L=V\lambda_{\mathrm{SC}}+R n_T\lambda_T+R\epsilon^2.
$$

A corresponding effective stochastic STAR model would use

$$
L_{\mathrm{STAR}}
\approx V_{\mathrm{STAR}}\lambda_{\mathrm{SC}}
+\sum_{r=1}^{R}-\log[1-p_{\mathrm{rot}}(\theta_r)],
$$

with $p_{\mathrm{rot}}$ covering the **complete** rotation procedure. When duration and errors are
correlated through retries, average the joint success probability rather than merely inserting mean
duration into this expression.

Toshio's raw channel contains coherent errors; its probabilistic coherent-error cancellation
provides a route to stochastic accounting. Include that procedure and its occasional extra
operations explicitly if using the corrected channel. See Toshio et al., Sec. IV.

Preparation patches also occupy space. Routing patches can be reused only when the schedule makes
them available. A small 2D layout should distinguish data, routing, preparation, and buffered states,
and the error accounting must avoid counting preparation faults twice. A prepared state waiting for
consumption must incur its storage cost.

For the initial comparison, retain the current pass criterion and exclude observable-level
probabilistic error cancellation: recovering expectation values through extra shots would require
a different success metric. This is distinct from the physical coherent-error cancellation used to
obtain Toshio's effective stochastic rotation channel.

## 6. Add a STAR design search and comparison outputs

| Repository component | Proposed change |
|---|---|
| New `capability/star.py` | Shared STAR scheduling, loss calculation, and minimum-qubit search. |
| New preparation module | Original and Toshio preparation models, with documented noise assumptions. |
| [`config.py`](../capability/config.py) | Protocol, preparation capacity/distance, postselection, switching, and retry settings. |
| [`quops.py`](../capability/quops.py) | Angle distribution and Clifford reduction. |
| Scripts and notebook | Compare Clifford+T, original STAR, and Toshio STAR on H1/H3; report area, duration, acceptance, and loss breakdowns. |
| Tests and methods note | Validate preparation limits, retry statistics, error conventions, and optimization against small exhaustive cases. |

The optimizer would search over data distance, preparation configuration, capacity, and buffering.
Its existing monotonicity and pruning assumptions need checking before reuse. The current
surface-code Clifford error law and timing conventions provide the common comparison framework;
any departures should be explicit modeling choices.

## Expected tradeoff

The expected tradeoff is **lower synthesis and factory overhead, but a residual rotation-error
ceiling**. For generic QUOPS angles, increasing code distance suppresses Clifford faults without
necessarily removing that ceiling. The resulting STAR capability curves could improve the
crossover at smaller machines while plateauing below the fully fault-tolerant curve at larger
machines. Numerical conclusions require implementing and validating the preparation and scheduling
models above.
