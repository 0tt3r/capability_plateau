"""QUOPS problem class (draft Sec. 2.2).

A circuit on an even number w of qubits alternates rotation layers (an independent R_P(theta) on every
qubit, P uniform in {X, Y, Z}, theta uniform) and CNOT layers (an independent uniformly random perfect
pairing). Each layer has size w. The shallow ensemble has w layers (s = w^2); the deep ensemble has
w^2 layers (s = w^3). We start with a rotation layer, so a circuit of L layers has ceil(L/2) rotation
layers and floor(L/2) CNOT layers.
"""
import numpy as np

ENSEMBLES = ('shallow', 'deep')


def layers(w, ensemble):
    return w if ensemble == 'shallow' else w*w


def size(w, ensemble):
    return w*w if ensemble == 'shallow' else w**3


def layer_counts(w, ensemble):
    n = layers(w, ensemble)
    return (n + 1)//2, n//2            # rotation layers, CNOT layers


def native_steps(axis, compile_mode):
    """Native time steps of R_P(theta) with gates {H, Z(theta)}: Z(t); H Z(t) H; Z(-pi/2) H Z(t) H Z(pi/2)."""
    if compile_mode == 'native':
        return 1
    return {'Z': 1, 'X': 3, 'Y': 5}[axis]


def random_pairing(w, rng):
    order = rng.permutation(w)
    return order.reshape(-1, 2)
