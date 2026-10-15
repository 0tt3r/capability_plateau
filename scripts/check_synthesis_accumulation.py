"""How do rotation-synthesis errors accumulate in QUOPS circuits? A direct check at small widths.

Every R_P(theta) of a sampled QUOPS circuit is replaced by its gridsynth word at tolerance eps, both circuits
are simulated on a state vector from |0...0>, and the infidelity of the two outputs is averaged over sampled
circuits. Here delta_i is the operator-norm distance of rotation i, minimized over global phase; gridsynth's
tolerance eps bounds the diamond-norm distance, so delta_i <= eps/2. If the errors of different rotations add
incoherently, the mean infidelity is about sum_i delta_i^2; if they add coherently, it is about
(sum_i delta_i)^2, larger by a factor of order the rotation count R. The surface-code model charges R eps^2.
A second, cheaper check measures mean delta^2 / eps^2 for single rotations down to eps = 1e-8, the range the
large designs use.

    .venv/bin/python scripts/check_synthesis_accumulation.py
writes results/synthesis_accumulation.json (about 8 minutes).
"""
import json
import sys
from math import log2, pi
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import mpmath
import numpy as np
from pygridsynth import gridsynth_gates

from capability.quops import layers, random_pairing

CASES = [(8, 'shallow'), (6, 'deep'), (8, 'deep')]
TOLERANCES = [0.05, 0.02, 0.01]
CIRCUITS = 24
SINGLE_TOLERANCES = [1e-2, 1e-4, 1e-6, 1e-8]
SINGLE_ANGLES = 200

H = np.array([[1, 1], [1, -1]])/np.sqrt(2)
S = np.diag([1, 1j])
GATES = dict(H=H, S=S, T=np.diag([1, np.exp(1j*pi/4)]), X=np.array([[0, 1], [1, 0]]), W=np.exp(1j*pi/4)*np.eye(2))
BASIS = dict(Z=np.eye(2), X=H, Y=S @ H)           # R_P(theta) = C R_Z(theta) C^dagger


def rz(theta):
    return np.diag([np.exp(-1j*theta/2), np.exp(1j*theta/2)])


def synthesize(theta, eps):
    mpmath.mp.dps = 30 + int(2*log2(1/eps))
    V = np.eye(2, dtype=complex)
    for g in str(gridsynth_gates(mpmath.mpf(float(theta)), mpmath.mpf(eps))):
        V = V @ GATES[g]                            # the word is an operator product as written
    return V


def distance(U, V):
    """Operator-norm distance, minimized over global phase."""
    t = np.trace(U.conj().T @ V)
    return np.linalg.norm(U - V*np.conj(t)/abs(t), 2)


def apply_1q(psi, U, q):
    return np.moveaxis(np.tensordot(U, psi, axes=([1], [q])), 0, q)


def apply_cnot(psi, c, t):
    psi = psi.copy()
    one = [slice(None)]*psi.ndim
    one[c] = 1
    sub = psi[tuple(one)]
    axis = t if t < c else t - 1
    psi[tuple(one)] = np.flip(sub, axis=axis)
    return psi


def run(w, ensemble, eps, rng):
    ideal = np.zeros((2,)*w, dtype=complex)
    ideal[(0,)*w] = 1
    synth = ideal.copy()
    deltas = []
    for layer in range(layers(w, ensemble)):
        if layer % 2 == 0:                          # rotation layers first, as in capability/quops.py
            for q in range(w):
                C = BASIS[rng.choice(list('XYZ'))]
                theta = rng.uniform(0, 2*pi)
                U, V = rz(theta), synthesize(theta, eps)
                ideal = apply_1q(ideal, C @ U @ C.conj().T, q)
                synth = apply_1q(synth, C @ V @ C.conj().T, q)
                deltas.append(distance(U, V))
        else:
            for a, b in random_pairing(w, rng):
                ideal, synth = apply_cnot(ideal, a, b), apply_cnot(synth, a, b)
    return 1 - abs(np.vdot(ideal, synth))**2, np.array(deltas)


def main():
    rng = np.random.default_rng(20261015)
    rows = []
    for w, ensemble in CASES:
        for eps in TOLERANCES:
            infid, incoh, coh = [], [], []
            for _ in range(CIRCUITS):
                x, delta = run(w, ensemble, eps, rng)
                infid.append(x)
                incoh.append(np.sum(delta**2))
                coh.append(np.sum(delta)**2)
            R = len(delta)
            row = dict(w=w, ensemble=ensemble, eps=eps, rotations=R, circuits=CIRCUITS,
                       infidelity=float(np.mean(infid)), stderr=float(np.std(infid)/np.sqrt(CIRCUITS)),
                       sum_delta_sq=float(np.mean(incoh)), sum_delta_all_sq=float(np.mean(coh)), charge=R*eps**2)
            row['infidelity_over_charge'] = row['infidelity']/row['charge']
            row['infidelity_over_sum_delta_sq'] = row['infidelity']/row['sum_delta_sq']
            rows.append(row)
            print(f"w={w} {ensemble:7s} eps={eps:<5} R={R:3d}  infidelity={row['infidelity']:.4f}"
                  f" +- {row['stderr']:.4f}  sum d^2={row['sum_delta_sq']:.4f}  (sum d)^2={row['sum_delta_all_sq']:.2f}"
                  f"  R eps^2={row['charge']:.4f}", flush=True)
    single = []
    for eps in SINGLE_TOLERANCES:
        theta = rng.uniform(0, 2*pi, SINGLE_ANGLES)
        delta = np.array([distance(rz(t), synthesize(t, eps)) for t in theta])
        single.append(dict(eps=eps, angles=SINGLE_ANGLES, mean_delta_sq_over_eps_sq=float(np.mean(delta**2)/eps**2),
                           max_delta_over_eps=float(delta.max()/eps)))
        print(f"single rotations, eps={eps:.0e}: mean delta^2/eps^2 = {single[-1]['mean_delta_sq_over_eps_sq']:.3f}",
              flush=True)
    (ROOT/'results/synthesis_accumulation.json').write_text(json.dumps(dict(circuits=rows, single=single), indent=1)
                                                            + '\n')


if __name__ == '__main__':
    main()
