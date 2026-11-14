#!/usr/bin/env python3
"""Build a deterministic, read-only-input package of final local Mac evidence."""
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import stat
import tarfile

from readback import MANIFEST, selftest

TASK = Path('/private/tmp/pigauto-final-closure-20261017')
FINAL = TASK / 'FINAL.json'
FINAL_INFO = json.loads(FINAL.read_text())
ARCHIVE = Path(FINAL_INFO['path']).resolve()
LOCAL = ARCHIVE.parent / 'local-macos'
OUT = Path(__file__).resolve().parent.parent
STAGE = OUT / '.staging-final-local'
TARBALL = OUT / 'pigauto-final-local-macos-evidence.tar.gz'
PUBLICATION = OUT / 'publication'


def sha(data):
    return hashlib.sha256(data).hexdigest()

def file_sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def put(stage, rel, data, mode=0o644, source_path=None):
    rel = Path(rel).as_posix()
    target = stage / rel
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(data)
    target.chmod(mode & 0o777)
    return {'path': rel, 'size_bytes': len(data), 'sha256': sha(data),
            'mode': mode & 0o777, 'source_path': str(source_path) if source_path else None}

def add_original(stage, src, rel, rows):
    src = Path(src)
    st = src.lstat()
    if not stat.S_ISREG(st.st_mode) or src.is_symlink():
        raise ValueError(f'input must be a regular nonsymlink file: {src}')
    rows.append(put(stage, rel, src.read_bytes(), stat.S_IMODE(st.st_mode), src))

def build_tar(stage, dest):
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode='wb', filename='', mtime=0, compresslevel=9) as gz:
        with tarfile.open(fileobj=gz, mode='w', format=tarfile.PAX_FORMAT) as tf:
            files = sorted(p for p in stage.rglob('*') if p.is_file() and not p.is_symlink())
            for path in files:
                rel = path.relative_to(stage).as_posix()
                data = path.read_bytes()
                info = tarfile.TarInfo(rel)
                info.size = len(data)
                info.mode = stat.S_IMODE(path.stat().st_mode)
                info.mtime = 0
                info.uid = info.gid = 0
                info.uname = info.gname = ''
                tf.addfile(info, io.BytesIO(data))
    blob = buf.getvalue()
    dest.write_bytes(blob)
    return blob

def main():
    if STAGE.exists() or TARBALL.exists():
        raise SystemExit('Refusing to overwrite an existing staging directory or evidence archive')
    if sha(ARCHIVE.read_bytes()) != FINAL_INFO['sha256'] or ARCHIVE.stat().st_size != FINAL_INFO['size_bytes']:
        raise SystemExit('Frozen archive identity differs from FINAL.json')
    STAGE.mkdir()
    members = []
    # Complete native platform-evidence folder, byte-for-byte.
    evidence = LOCAL / 'platform-evidence'
    for src in sorted(evidence.rglob('*')):
        if src.is_symlink():
            raise SystemExit(f'Refusing symlink in platform evidence: {src}')
        if src.is_file():
            add_original(STAGE, src, 'local-macos/' + src.relative_to(LOCAL).as_posix(), members)
    # Every regular file directly at the R CMD check root (check outputs, logs, examples, manuals).
    check = LOCAL / 'frozen-check/pigauto.Rcheck'
    for src in sorted(check.iterdir()):
        if src.is_symlink():
            raise SystemExit(f'Refusing symlink at Rcheck root: {src}')
        if src.is_file():
            add_original(STAGE, src, 'local-macos/' + src.relative_to(LOCAL).as_posix(), members)
    # Keep the native test summary, not the full generated source/fixture trees.
    rout = check / 'tests/testthat.Rout'
    add_original(STAGE, rout, 'local-macos/' + rout.relative_to(LOCAL).as_posix(), members)
    for name in ('local-receipt.json', 'native-check.log', 'warning-replay.log', 'local-note-diagnostic.json'):
        add_original(STAGE, LOCAL / name, 'local-macos/' + name, members)

    native_log = (check / '00check.log').read_text(errors='replace')
    status_lines = [line for line in native_log.splitlines() if line.startswith('Status:')]
    check_exit = (evidence / 'check-exit.txt').read_text().splitlines()
    receipt = json.loads((LOCAL / 'local-receipt.json').read_text())
    replay_driver = json.loads((LOCAL / 'local-note-diagnostic.json').read_text())
    replay = json.loads((evidence / 'warning-replay.json').read_text())
    if status_lines != ['Status: 1 NOTE'] or 'exit=0' not in check_exit or receipt['records'][0]['exit_code'] != 1:
        raise SystemExit(f'Unexpected native verdict or wrapper status: {status_lines}, {check_exit}, {receipt}')
    if replay_driver['exit_code'] != 0 or replay_driver['archive_sha256_after'] != FINAL_INFO['sha256']:
        raise SystemExit('Replay did not exit successfully with the exact archive unchanged')
    if replay.get('error') != {}:
        raise SystemExit('Unexpected replay receipt error payload')
    note = '* checking examples ... [13s/63s] NOTE\nExamples with CPU (user + system) or elapsed time > 5s\n                 user system elapsed\ncompare_methods 2.532  0.992  18.406'
    if note not in native_log:
        raise SystemExit('Native timing NOTE differs from the exact retained NOTE')
    native_findings = [line for line in native_log.splitlines()
                       if line.startswith('* ') and line.endswith((' NOTE', ' WARNING', ' ERROR'))]
    if native_findings != ['* checking examples ... [13s/63s] NOTE']:
        raise SystemExit(f'Native log has additional or changed findings: {native_findings}')
    expected_counts = {'fail': 0, 'warnings': 161, 'skips': 86, 'pass': 3195}
    if replay['native_counts'] != expected_counts or replay['replay_counts'] != expected_counts:
        raise SystemExit('Native/replay warning-capture counts differ from retained expected counts')

    # Generated, review-facing readable view. Raw source files above stay untouched.
    summary = (
        'Frozen archive local macOS evidence\n'
        f"Archive SHA-256: {FINAL_INFO['sha256']}\n"
        f"Archive size: {FINAL_INFO['size_bytes']} bytes\n"
        f"Source commit: {FINAL_INFO['source_commit']}\n\n"
        f"Native R CMD check status: {status_lines[0]}\n"
        f"Native R CMD check exit: {check_exit[0]}\n"
        f"Native check wrapper exit: {receipt['records'][0]['exit_code']} (wrapper rejects NOTE findings)\n"
        'Native finding: sole NOTE is compare_methods examples timing, user=2.532s system=0.992s elapsed=18.406s.\n'
        'Native timing NOTE adjudication: UNADJUDICATED. This is not native Status: OK.\n'
        f"Native check driver elapsed: {receipt['records'][0]['elapsed_seconds']} s; R check reporter: {check_exit[1]}\n\n"
        f"Separate warning replay exit: {replay_driver['exit_code']}\n"
        f"Replay driver elapsed: {replay_driver['elapsed_seconds']} s; R-side reporter: {replay['diagnostic_elapsed_seconds']} s\n"
        f"Replay counts: {json.dumps(replay['replay_counts'], sort_keys=True)}\n"
        'Replay scope: diagnostic only; separate from native R CMD check time. Warnings remain unadjudicated.\n'
        'Four-cell warning audit: NOT_DONE; reader guard rejects error={} as non-null. No checker/contract change included.\n'
        'No release READY or CRAN submission claim follows from this package.\n'
    )
    generated = PUBLICATION
    for name in ('README.md', 'build.py', 'readback.py'):
        src = generated / name
        mode = stat.S_IMODE(src.stat().st_mode)
        members.append(put(STAGE, 'publication/' + name, src.read_bytes(), mode, src))
    producer = {
        'format': 'pigauto-final-local-evidence-producer-v1',
        'source_commit': FINAL_INFO['source_commit'],
        'final_archive': {'path': str(ARCHIVE), 'sha256': FINAL_INFO['sha256'], 'size_bytes': FINAL_INFO['size_bytes']},
        'source_hashes': {
            'FINAL.json': file_sha(FINAL),
            'native_check_log': file_sha(LOCAL / 'native-check.log'),
            'native_check_00check': file_sha(check / '00check.log'),
            'local_receipt': file_sha(LOCAL / 'local-receipt.json'),
            'warning_replay_log': file_sha(LOCAL / 'warning-replay.log'),
            'local_note_diagnostic': file_sha(LOCAL / 'local-note-diagnostic.json'),
            'capture_warnings_script': file_sha(Path('/private/tmp/pigauto-pr236-validator-20261017/script/cran-0.11-platform/capture-warnings.R')),
            'direct_verifier_script': file_sha(Path('/private/tmp/pigauto-pr236-validator-20261017/script/cran-0.11-platform/verify_results.py')),
        },
        'producer_sha256': {'build.py': file_sha(generated / 'build.py'), 'readback.py': file_sha(generated / 'readback.py'), 'README.md': file_sha(generated / 'README.md')},
        'native_verdict': {'status': status_lines[0], 'check_exit': 0, 'wrapper_exit': receipt['records'][0]['exit_code'], 'adjudication': 'UNADJUDICATED'},
        'separate_replay': {'exit': replay_driver['exit_code'], 'driver_elapsed_seconds': replay_driver['elapsed_seconds'], 'reporter_elapsed_seconds': replay['diagnostic_elapsed_seconds'], 'counts': replay['replay_counts']},
        'warning_audit': {'status': 'NOT_DONE', 'reason': 'Audit reader requires null error field; replay receipt encodes empty object.'},
    }
    members.append(put(STAGE, 'publication/producer-identity.json', json.dumps(producer, indent=2, sort_keys=True).encode() + b'\n'))
    members.append(put(STAGE, 'publication/readable-summary.txt', summary.encode('utf-8')))
    manifest = {'schema': 'pigauto-final-local-evidence-members-v1', 'source_commit': FINAL_INFO['source_commit'],
                'archive_sha256': FINAL_INFO['sha256'], 'members': sorted(members, key=lambda r: r['path'])}
    manifest_bytes = json.dumps(manifest, indent=2, sort_keys=True).encode() + b'\n'
    put(STAGE, MANIFEST, manifest_bytes)
    blob = build_tar(STAGE, TARBALL)
    result = selftest(blob)
    (OUT / 'roundtrip-validation.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    manifest_sha = sha(manifest_bytes)
    index = {'schema': 'pigauto-final-local-evidence-index-v1', 'archive': TARBALL.name,
             'archive_size_bytes': len(blob), 'archive_sha256': sha(blob), 'member_manifest': MANIFEST,
             'member_manifest_sha256': manifest_sha, 'member_count': len(members) + 1,
             'raw_original_member_count': sum(x['source_path'] is not None and str(x['source_path']).startswith(str(LOCAL)) for x in members),
             'producer_identity_sha256': sha((STAGE / 'publication/producer-identity.json').read_bytes()),
             'roundtrip_validation': str(OUT / 'roundtrip-validation.json'),
             'positive_roundtrip': result['positive_roundtrip']['status'],
             'planted_hash_failure_rejected': result['planted_hash_failure']['rejected'],
             'planted_path_failure_rejected': result['planted_path_failure']['rejected']}
    (OUT / 'publication-index.json').write_text(json.dumps(index, indent=2, sort_keys=True) + '\n')
    print(json.dumps(index, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
