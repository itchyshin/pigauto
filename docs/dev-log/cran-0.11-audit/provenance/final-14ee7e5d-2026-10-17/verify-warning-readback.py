"""Rehash reviewed actual warning-audit inputs without rerunning fits/tests."""
from pathlib import Path
import hashlib
import json

TASK = Path('/private/tmp/pigauto-final-closure-20261017')
ROOT = TASK / 'final-warning-audit'
result = json.loads((ROOT / 'actual-results-null-error-fix/audit.json').read_text())
assert result['archive_sha256'] == '14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb'
assert result['archive_size_bytes'] == 5137648 and result['direct_run'] == '39082452558'
assert result['passed_cells'] == 4 and len(result['records']) == 4
assert {r['cell'] for r in result['records']} == {'local-macos', 'windows-latest-release', 'windows-latest-devel', 'ubuntu-latest-release'}
for record in result['records']:
    assert record['warning_audit_status'] == 'PASS', record['cell']
    for name, expected in record['input_sha256'].items():
        p = Path(name)
        assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest() == expected, p
    expected_status = 'NOTE' if record['cell'] == 'local-macos' else 'OK'
    assert record['native_verdict']['status'] == expected_status
    assert record['native_counts']['fail'] == 0
    assert sum(record['warning_categories'].values()) == record['native_counts']['warnings']
assert result['native_check_adjudication'] == 'UNADJUDICATED'
print('FOUR_CELL_ACTUAL_WARNING_INPUTS_REVERIFIED_NATIVE_NOTE_RETAINED')
