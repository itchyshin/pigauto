"""Read-only final approval-checkpoint verifier; does not authorize publication."""
import hashlib
import json
import re
import subprocess
from pathlib import Path

REPO = Path('/private/tmp/pigauto-pr236-validator-20261017')
TASK = Path('/private/tmp/pigauto-final-closure-20261017')
SHA = '14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb'
SOURCE = '8e275ca738129c8d608f06ffc407b8141c1a2636'
FINAL = REPO / 'docs/dev-log/cran-0.11-audit/provenance/final-14ee7e5d-2026-10-17'

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def read(p):
    return json.loads(p.read_text())

ledger = read(REPO / 'docs/dev-log/cran-0.11-audit/release-ledger.json')
assert ledger['artifact']['sha256'] == SHA
assert ledger['artifact']['source_commit'] == SOURCE
assert ledger['status_claim'] in ('platform-clean', 'submission-ready')
assert ledger['submission']['state'] == 'not_submitted'
assert ledger['submission']['cran_upload_authorized'] is False
assert ledger['submission']['maintainer_confirmation_authorized'] is False
assert ledger['evidence_pr']['number'] == 236 and ledger['evidence_pr']['merged'] is False

for p in {Path(ledger['artifact']['path']), Path(read(TASK / 'STATE.json')['persistent_archive_path']), Path(read(TASK / 'FINAL.json')['path'])}:
    assert p.is_file() and digest(p) == SHA and p.stat().st_size == 5137648, p
    assert p.stat().st_mode & 0o222 == 0, p

local = Path(read(TASK / 'FINAL.json')['path']).parent / 'local-macos'
native = (local / 'frozen-check/pigauto.Rcheck/00check.log').read_text()
assert re.findall(r'^Status:.*$', native, re.M) == ['Status: 1 NOTE']
assert 'compare_methods 2.532  0.992  18.406' in native
assert 'exit=0' in (local / 'platform-evidence/check-exit.txt').read_text().splitlines()
assert read(local / 'local-receipt.json')['records'][0]['exit_code'] == 1
assert read(local / 'local-note-diagnostic.json')['exit_code'] == 0

audit = read(TASK / 'final-warning-audit/actual-results-null-error-fix/audit.json')
assert audit['passed_cells'] == 4
assert all(c['warning_audit_status'] == 'PASS' for c in audit['records'])

for dirname, wrapped in [('checks-and-integration', True), ('external-platforms', False)]:
    root = FINAL / dirname
    index = read(root / 'publication-files.json')
    entries = index['files'] if wrapped else index
    for name, item in entries.items():
        p = root / name
        assert p.is_file() and not p.is_symlink(), p
        assert digest(p) == item['sha256'] and p.stat().st_size == item['size_bytes'], p

for reviewer in ('grace', 'rose', 'pat'):
    receipt = read(FINAL / 'panel' / reviewer / 'verdict.json')
    assert receipt['artifact_sha256'] == SHA and receipt['source_commit'] == SOURCE
    assert receipt['verdict'] == 'READY', (reviewer, receipt['verdict'])

manifest = FINAL / 'files.sha256'
for line in manifest.read_text().splitlines():
    expected, name = line.split('  ', 1)
    p = FINAL / name
    assert p.is_file() and digest(p) == expected, p

package_paths = ['DESCRIPTION', 'NAMESPACE', 'R', 'man', 'inst', 'data', 'vignettes', 'README.md', 'LICENSE']
subprocess.run(['git', 'diff', '--exit-code', SOURCE, '--', *package_paths], cwd=REPO, check=True)
assert subprocess.check_output(['git', 'branch', '--show-current'], cwd=REPO, text=True).strip() == 'release/cran-0.11-gate-780c084'
print('FINAL_CHECKPOINT_VERIFIED_PENDING_HUMAN_APPROVAL')
