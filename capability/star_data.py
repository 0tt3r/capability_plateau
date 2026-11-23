"""Circuit-derived inputs for the adopted optimistic two-qubit STAR model.

Tables retain all acceptance orders and only identified logical-error sectors.
Their provenance and original-STAR numerical integration uncertainty are stored
alongside the values. Unsupported noise rates are rejected, not extrapolated.
"""
from functools import lru_cache
import json
from pathlib import Path


@lru_cache(maxsize=1)
def tables():
    return json.loads((Path(__file__).parent / 'data/star_acceptance.json').read_text())


def hardware_data(hw):
    for entry in tables()['hardware'].values():
        if (hw.eps_2q, hw.eps_1q) == (entry['p'], entry['q']):
            return entry
    raise ValueError('Counted STAR inputs cover H1/H3 noise rates only; regenerate '
                     'circuit acceptance for different physical error rates')


def transversal_data(hw, k):
    try:
        return hardware_data(hw)['transversal'][str(k)]
    except KeyError as exc:
        raise ValueError('Counted STAR preparation covers k=3,5,7,9 only') from exc


def original_data(hw, d):
    try:
        return hardware_data(hw)['original'][str(d)]
    except KeyError as exc:
        raise ValueError('Counted original STAR covers odd distances 3–31 only') from exc


def validate_counted(hw, sc):
    if sc.rotation_weight != 2 or sc.postselection != 'strip':
        raise ValueError('Counted STAR requires rotation_weight=2 and strip postselection')
    if sc.rejection_scale != 1:
        raise ValueError('Counted acceptance cannot scale guessed rejection counts; '
                         'use omitted_sector_stress or regenerate circuit inputs')
    if sc.max_distance > 31 or any(k not in (3, 5, 7, 9) for k in sc.groups):
        raise ValueError('Counted STAR inputs cover odd distances 3–31 and k=3,5,7,9')
    if hw.eps_1q != 0 or hw.eps_2q != 0:
        hardware_data(hw)
