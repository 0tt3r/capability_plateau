"""Cached comparison inputs remain portable without accepting stale numerical data."""
import copy
import hashlib
import json
from dataclasses import asdict

import pytest

from capability.config import BASELINE, HARDWARE
from scripts import make_frozen_pinnacle_curves as frozen_runner
from scripts import make_pinnacle_curves as pinnacle_runner


def write_json(path, value, newline=b'\n', indent=2):
    payload = json.dumps(value, indent=indent, ensure_ascii=False).encode('utf-8')
    path.write_bytes(payload.replace(b'\n', newline))


@pytest.fixture
def references(tmp_path, monkeypatch):
    root = tmp_path/'repository'
    results = root/'results'
    results.mkdir(parents=True)
    choices = json.loads(json.dumps(asdict(BASELINE)))
    baseline = dict(choices=choices, description='portable caf\u00e9 reference')
    original = dict(choices=copy.deepcopy(choices), pinnacle_choices={},
                    provenance_hash_format=pinnacle_runner.HASH_FORMAT)
    for key, hw in HARDWARE.items():
        baseline[key] = dict(hardware=asdict(hw), direct=[],
                             surface=[dict(passes=True, qubits=4000, s=64)], probes={})
        original[key] = dict(hardware=asdict(hw), pinnacle=[], probes={})
    original['baseline_sha256'] = pinnacle_runner.canonical_json_sha256(baseline)
    write_json(results/'baseline.json', baseline)
    write_json(results/'pinnacle.json', original)
    monkeypatch.setattr(pinnacle_runner, 'ROOT', root)
    monkeypatch.setattr(frozen_runner, 'ROOT', root)
    return results, baseline, original


@pytest.mark.parametrize('runner', [pinnacle_runner, frozen_runner])
def test_hash_has_an_explicit_canonical_json_format(runner):
    value = {'n': 2, 'label': 'caf\u00e9'}
    expected_bytes = b'{"label":"caf\\u00e9","n":2}'
    assert runner.HASH_FORMAT == 'canonical-json-sha256-v1'
    assert runner.canonical_json_sha256(value) == hashlib.sha256(expected_bytes).hexdigest()


@pytest.mark.parametrize('runner', [pinnacle_runner, frozen_runner])
def test_nonfinite_values_cannot_enter_a_provenance_hash(runner):
    with pytest.raises(ValueError):
        runner.canonical_json_sha256({'invalid': float('nan')})


@pytest.mark.parametrize('newline', [b'\n', b'\r\n'])
@pytest.mark.parametrize('indent', [None, 4])
def test_both_loaders_accept_lf_crlf_and_whitespace_changes(references, newline, indent):
    results, baseline, original = references
    # Reverse top-level key order as well: hashes describe data, not checkout formatting.
    write_json(results/'baseline.json', dict(reversed(list(baseline.items()))), newline, indent)
    write_json(results/'pinnacle.json', dict(reversed(list(original.items()))), newline, indent)
    loaded, digest = pinnacle_runner.load_references()
    loaded_baseline, loaded_original, hashes = frozen_runner.load_references()
    assert loaded == loaded_baseline == baseline
    assert loaded_original == original
    assert digest == hashes['baseline.json'] == original['baseline_sha256']
    assert hashes['pinnacle.json'] == frozen_runner.canonical_json_sha256(original)


def test_changed_numerical_reference_rejects_stale_pinnacle_cache(references):
    results, baseline, _ = references
    baseline['H2']['surface'][0]['s'] += 1
    write_json(results/'baseline.json', baseline)
    with pytest.raises(ValueError, match='different baseline inputs.*make_pinnacle_curves.py'):
        frozen_runner.load_references()


@pytest.mark.parametrize('version', ['missing', None, '', 'raw-byte-sha256-v0'])
def test_missing_or_unknown_hash_format_requires_regeneration(references, version):
    results, _, original = references
    if version == 'missing':
        original.pop('provenance_hash_format')
    else:
        original['provenance_hash_format'] = version
    write_json(results/'pinnacle.json', original)
    with pytest.raises(ValueError, match='provenance hash format.*make_pinnacle_curves.py'):
        frozen_runner.load_references()


@pytest.mark.parametrize('filename', ['baseline.json', 'pinnacle.json'])
def test_cached_choice_mismatch_is_rejected(references, filename):
    results, baseline, original = references
    data = baseline if filename == 'baseline.json' else original
    data['choices']['threshold'] = 0.4
    write_json(results/filename, data)
    with pytest.raises(ValueError, match='choices do not match'):
        frozen_runner.load_references()
    if filename == 'baseline.json':
        with pytest.raises(ValueError, match='choices do not match'):
            pinnacle_runner.load_references()


@pytest.mark.parametrize('filename', ['baseline.json', 'pinnacle.json'])
def test_cached_hardware_mismatch_is_rejected(references, filename):
    results, baseline, original = references
    data = baseline if filename == 'baseline.json' else original
    data['H2']['hardware']['eps_2q'] = 0.002
    write_json(results/filename, data)
    with pytest.raises(ValueError, match='hardware.*does not match'):
        frozen_runner.load_references()
    if filename == 'baseline.json':
        with pytest.raises(ValueError, match='hardware parameters do not match'):
            pinnacle_runner.load_references()


def test_generated_comparisons_declare_their_hash_format(references, monkeypatch):
    _, baseline, original = references

    def no_pinnacle(*args):
        return dict(supported=False, unsupported_reason='test fixture', pinnacle=[], probes={})

    def no_frozen(*args):
        return dict(supported=False, unsupported_reason='test fixture', frozen_pinnacle=[], probes={})

    monkeypatch.setattr(pinnacle_runner, 'compute_pinnacle', no_pinnacle)
    monkeypatch.setattr(frozen_runner, 'compute_frozen', no_frozen)
    p = pinnacle_runner.compute()
    f = frozen_runner.compute()
    assert p['provenance_hash_format'] == f['provenance_hash_format'] == pinnacle_runner.HASH_FORMAT
    assert p['baseline_sha256'] == original['baseline_sha256']
    assert f['reference_sha256']['baseline.json'] == pinnacle_runner.canonical_json_sha256(baseline)
