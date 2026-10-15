"""Measure gridsynth T counts for Z(theta) and fit n_T(eps) = a log2(1/eps) + b.

gridsynth (Ross and Selinger, arXiv:1403.2975, via pygridsynth) returns an ancilla-free Clifford+T word
whose channel is within eps of Z(theta) in diamond norm: pygridsynth's error is 2 sqrt(1 - Re(z* u)^2), twice
the operator-norm distance minimized over global phase. We synthesize uniformly random angles at each
eps, count T gates, and fit the mean. The surface-code model rounds the fitted law up per rotation.

    .venv/bin/python scripts/fit_synthesis.py
"""
import json
import time
from math import log2, pi
from pathlib import Path

import mpmath
import numpy as np
from pygridsynth import gridsynth_gates

ROOT = Path(__file__).resolve().parents[1]
EPSILONS = [10.**-k for k in range(1, 13)]
ANGLES = 64


def main():
    rng = np.random.default_rng(20261013)
    rows = []
    start = time.time()
    for eps in EPSILONS:
        mpmath.mp.dps = 30 + int(2*log2(1/eps))
        counts = [str(gridsynth_gates(mpmath.mpf(float(theta)), mpmath.mpf(eps))).count('T')
                  for theta in rng.uniform(0, 2*pi, ANGLES)]
        rows.append(dict(eps=eps, mean_T=float(np.mean(counts)), min_T=int(min(counts)), max_T=int(max(counts))))
        print(f'eps={eps:.0e}  mean T={rows[-1]["mean_T"]:.2f}', flush=True)
    x = np.array([log2(1/r['eps']) for r in rows])
    y = np.array([r['mean_T'] for r in rows])
    a, b = np.polyfit(x, y, 1)
    out = dict(law='n_T(eps) = a log2(1/eps) + b, mean over uniformly random angles',
               a=float(a), b=float(b), max_residual=float(np.max(np.abs(a*x + b - y))),
               angles_per_eps=ANGLES, tool='pygridsynth gridsynth_gates', rows=rows,
               seconds=round(time.time() - start, 1))
    (ROOT/'results').mkdir(exist_ok=True)
    (ROOT/'results/synthesis_fit.json').write_text(json.dumps(out, indent=2) + '\n')
    print(f'a = {a:.3f}, b = {b:.3f}')


if __name__ == '__main__':
    main()
