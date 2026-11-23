"""QUOPS surface-code interface, using counted sequential PBC by default.

All Clifford gates are absorbed into Pauli axes; no CNOT runtime is charged.
The sequential schedule includes cold factory startup and finite output buffers.
Parallel PBC units remain an explicit alternative with their own storage costs.
See docs/surface_architectures.md for the resource model and its limitations.
"""
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from .config import PBC, N_MAX, SurfaceArchitecture
from .quops import ENSEMBLES, layer_counts, size

FIT = json.loads((Path(__file__).resolve().parents[1]/'results/synthesis_fit.json').read_text())
NT = np.arange(5, 161)                          # T gates per rotation
EPS = 2.**(-(NT - FIT['b'])/FIT['a'])           # loosest tolerance the fit gives exactly NT T gates (0.42 to 3e-16)


def t_count(eps):
    """T gates per rotation from the gridsynth fit (results/synthesis_fit.json), rounded up."""
    return np.ceil(FIT['a']*np.log2(1/eps) + FIT['b'])


def p_sc(d, hw, ch):
    p_th = ch.sc_threshold_uniform if hw.noise == 'uniform' else ch.sc_threshold_bimodal
    return ch.sc_prefactor*(hw.p/p_th)**((d + 1)/2)


def data_tiles(w, ch, architecture=PBC):
    """Count computational tiles; factory output buffers are added by the schedule."""
    if architecture.layout != 'fast_columns':
        raise ValueError('Use the counted fast_columns layout')
    from .layouts import fast_layout
    return fast_layout(w).area


def min_qubits(w, ensemble, hw, ch, architecture=PBC):
    """Fewest qubits within the selected PBC architecture and finite design menu."""
    if architecture.schedule == 'pbc_units':
        from . import surface_units
        result = surface_units.min_qubits(w, ensemble, hw, ch, architecture)
    elif architecture.schedule == 'pbc':
        from . import surface_pbc
        result = surface_pbc.min_qubits(w, ensemble, hw, ch, architecture)
    else:
        raise ValueError('Unknown surface schedule: ' + architecture.schedule)
    if result is not None:
        result['surface_architecture'] = asdict(architecture)
    return result


def widths(w_max=200000):
    """Sampled widths: every even width to 100, then about 6.5% apart."""
    small = list(range(2, 101, 2))
    large = sorted({int(2*round(v/2)) for v in np.geomspace(102, w_max, 120)})
    return small + [v for v in large if v > 100]


def sweep(hw, ch, architecture=PBC):
    """Designs at the sampled widths up to N_MAX qubits. The fewest qubits never decrease with w (see best_at),
    so the first width that fails or needs more than N_MAX ends each ensemble."""
    rows = []
    for ensemble in ENSEMBLES:
        for w in widths():
            design = min_qubits(w, ensemble, hw, ch, architecture)
            if design is None or design['qubits'] > N_MAX:
                break
            rows.append(dict(design, ensemble=ensemble, w=w, s=size(w, ensemble), passes=True))
    return rows


def best_at(N, hw, ch, rows, architecture=None):
    """The largest-s design with at most N qubits, exact over every even width.

    The fewest qubits needed never decreases with w: a design (distance, factory, T count, factory count) that
    passes at w + 2 also passes at w, with no more qubits. So it suffices to bisect the gap between the largest
    sampled width that fits and the next sampled width."""
    if any('surface_architecture' not in row for row in rows):
        raise ValueError('Surface rows lack architecture metadata; recompute them with PBC')
    row_architectures = [SurfaceArchitecture(**row['surface_architecture']) for row in rows]
    if architecture is None:
        architecture = row_architectures[0] if rows else PBC
    if any(a != architecture for a in row_architectures):
        raise ValueError('Cannot refine rows using a different surface architecture')
    top = None
    for ensemble in ENSEMBLES:
        fit = sorted((r['w'], r) for r in rows if r['ensemble'] == ensemble and r['qubits'] <= N)
        if not fit:
            continue
        lo, row = fit[-1]
        hi = next((w for w in widths() if w > lo), lo + 2)     # first sampled width that does not fit
        while hi - lo > 2:
            mid = lo + 2*((hi - lo)//4)
            x = min_qubits(mid, ensemble, hw, ch, architecture)
            if x is not None and x['qubits'] <= N:
                lo, row = mid, dict(x, ensemble=ensemble, w=mid, s=size(mid, ensemble), passes=True)
            else:
                hi = mid
        if top is None or row['s'] > top['s']:
            top = row
    return top


def capability(N, hw, ch, rows, architecture=None):
    row = best_at(N, hw, ch, rows, architecture)
    return row['s'] if row else 0
