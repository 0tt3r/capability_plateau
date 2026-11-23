"""Capability curves C(N): the largest QUOPS size s solved with at most N physical qubits (draft Sec. 2).

A size s is solved if some s' >= s passes in either ensemble, so C(N) is the largest passing size whose
required qubits are at most N.
"""
import numpy as np

from .config import N_MAX

N_GRID = np.unique(np.round(np.geomspace(1, N_MAX, 1201)))


def capability(rows, N):
    fits = [r['s'] for r in rows if r['passes'] and r['qubits'] <= N]
    return max(fits, default=0)


def curve(rows, grid=N_GRID):
    ok = sorted((r['qubits'], r['s']) for r in rows if r['passes'])
    out, best, k = [], 0, 0
    for N in grid:
        while k < len(ok) and ok[k][0] <= N:
            best = max(best, ok[k][1])
            k += 1
        out.append(best)
    return np.array(out)


def crossover(direct_rows, surface_rows):
    """Smallest N at which the surface-code capability exceeds the direct capability. Both curves are step
    functions that change only at the qubit counts of passing rows, so comparing them there is exact."""
    points = sorted({r['qubits'] for r in direct_rows + surface_rows if r['passes']})
    a, b = curve(direct_rows, points), curve(surface_rows, points)
    above = np.nonzero(b > a)[0]
    return float(points[above[0]]) if len(above) else None


def lead_intervals(rows, reference_rows, max_qubits=N_MAX):
    """Half-open integer-qubit intervals where a curve strictly exceeds all pooled references."""
    stop = int(max_qubits) + 1
    points = sorted({int(r['qubits']) for r in rows + reference_rows
                     if r['passes'] and 1 <= r['qubits'] < stop})
    above = curve(rows, points) > curve(reference_rows, points)
    intervals = []
    for i, leads in enumerate(above):
        if leads and (i == 0 or not above[i-1]):
            intervals.append([points[i], stop])
        elif not leads and i > 0 and above[i-1]:
            intervals[-1][1] = points[i]
    return intervals
