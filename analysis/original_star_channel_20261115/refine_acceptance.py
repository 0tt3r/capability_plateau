"""Refine the distance-21 integration near the H3 stress-test threshold.

Uses independent batches to avoid allocating all Fourier samples at once.
Rerunning is deterministic and replaces, rather than compounds, the refinement.
"""
import json
from pathlib import Path

import numpy as np

from channel import acceptance

HERE = Path(__file__).resolve().parent


def main():
    data = json.loads((HERE / 'acceptance.json').read_text())
    batches = [acceptance(21, samples=262144, seed=15112026 + 10000 * i)
               for i in range(24)]
    result = batches[0].copy()
    result['hardware'] = {}
    result['samples'] = sum(b['samples'] for b in batches)
    result['batch_seeds'] = [b['seed'] for b in batches]
    n = batches[0]['samples']
    for label in ('H1', 'H3'):
        records = [b['hardware'][label] for b in batches]
        values = np.array([r['acceptance_over_zero_signature'] for r in records])
        mean = float(values.mean())
        # Reconstruct the pooled second moment from each batch's sample SE.
        squares = sum((n - 1) * n * r['ratio_standard_error'] ** 2
                      for r in records) + n * float(np.sum((values - mean) ** 2))
        se = float(np.sqrt(squares / (result['samples'] - 1) / result['samples']))
        row = records[0].copy()
        row.update(acceptance_over_zero_signature=mean, ratio_standard_error=se,
                   acceptance=row['zero_signature_probability'] * mean,
                   acceptance_standard_error=row['zero_signature_probability'] * se,
                   conditional_error=row['retained_local_one'] / mean)
        result['hardware'][label] = row
    data['rows'] = [result if r['distance'] == 21 else r for r in data['rows']]
    (HERE / 'acceptance.json').write_text(json.dumps(data, indent=2) + '\n')
    (HERE / 'acceptance_refinement.json').write_text(json.dumps(batches, indent=2) + '\n')
    print({h: (r['acceptance'], r['acceptance_standard_error'])
           for h, r in result['hardware'].items()})


if __name__ == '__main__':
    main()
