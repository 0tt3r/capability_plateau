"""Package the checked circuit inputs used by the main STAR estimator.

Use --check to verify the packaged file without changing it. See
analysis/star_inputs.md for the upstream circuit calculations and their scope.
"""
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCES = (
    'analysis/separate_branch_survival_20261113/results.json',
    'analysis/separate_branch_survival_20261113/two_fault_terms.json',
    'analysis/original_star_channel_20261115/acceptance.json',
)


def package():
    full, one, original = (json.loads((ROOT/path).read_text())['rows'] for path in SOURCES)
    one = {(r['distance'], r['hardware']): r for r in one}
    data = dict(schema_version=1,
                scope='Compiled CNOT-RZ-CNOT m=2; retained single-fault logical components, other accepted logical effects idealized.',
                sources=list(SOURCES), hardware={})
    for label in ('H1', 'H3'):
        transversal = {str(r['groups']): dict(
            H_n=r['hardware'][label]['H_n'],
            T1=one[r['distance'], label]['exactly_one_active'],
            zero_signature_probability=one[r['distance'], label]['no_active_fault_probability'])
            for r in full}
        originals = {str(r['distance']): dict(r['hardware'][label], samples=r['samples'], exact=r['exact'])
                     for r in original}
        data['hardware'][label] = dict(p=.001, q=.001 if label == 'H1' else .000001,
                                      transversal=transversal, original=originals)
    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    data = package()
    output = ROOT/'capability/data/star_acceptance.json'
    if args.check:
        if data != json.loads(output.read_text()):
            parser.exit(1, 'Packaged STAR inputs differ from the circuit source data.\n')
        print('Packaged STAR inputs match the circuit source data.')
    else:
        output.write_text(json.dumps(data, indent=2)+'\n')
        print(f'Wrote {output.relative_to(ROOT)}')
