#!/usr/bin/env python3
"""Run only the native strict check against the frozen successor archive."""
from pathlib import Path
import datetime
import hashlib
import json
import os
import signal
import shutil
import subprocess
import sys
import time
import tarfile

TASK = Path('/private/tmp/pigauto-final-closure-20261017/successor-b2e6054')
DRIVER = Path('/private/tmp/pigauto-pr236-validator-20261017')
FINAL = TASK / 'FINAL.json'
EXPECTED_COMMIT = 'b2e6054685e6bdedb073ba46e058af2c0a619a18'
EXPECTED_ARCHIVE_NAME = 'pigauto_0.11.0.tar.gz'
CHECK_SCRIPT = DRIVER / 'script/cran-0.11-platform/check.R'
TIMEOUT_SECONDS = 75 * 60
TERM_GRACE_SECONDS = 15


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def validate_final():
    if not FINAL.is_file():
        raise RuntimeError(f'No frozen successor FINAL.json at {FINAL}; refusing to run')
    receipt = json.loads(FINAL.read_text())
    required = ('path', 'sha256', 'size_bytes', 'source_commit', 'clean_worktree',
                'forbidden_entries', 'rights_unchanged_from_predecessor',
                'extracted_source', 'archive_entries', 'inventory_sha256')
    missing = [key for key in required if key not in receipt]
    if missing:
        raise RuntimeError(f'FINAL.json lacks required freeze fields: {missing}')
    archive = Path(receipt['path'])
    if archive.name != EXPECTED_ARCHIVE_NAME:
        raise RuntimeError(f'Unexpected archive basename: {archive.name}')
    expected_parent = Path('/private/tmp/pigauto-release-candidate') / ('sha256-' + receipt['sha256'])
    if archive.parent != expected_parent:
        raise RuntimeError(f'Archive path is not content-addressed as expected: {archive}')
    sidecar = archive.parent / 'freeze-receipt.json'
    inventory = archive.parent / 'inventory.txt'
    if not sidecar.is_file() or json.loads(sidecar.read_text()) != receipt:
        raise RuntimeError('FINAL.json does not match the immutable-folder freeze receipt')
    if not inventory.is_file():
        raise RuntimeError('Frozen inventory.txt is missing')
    inventory_bytes = inventory.read_bytes()
    if hashlib.sha256(inventory_bytes).hexdigest() != receipt['inventory_sha256']:
        raise RuntimeError('Frozen inventory checksum mismatch')
    with tarfile.open(archive, 'r:gz') as tar:
        actual_names = sorted(m.name.rstrip('/') + ('/' if m.isdir() else '')
                              for m in tar.getmembers())
    if ('\n'.join(actual_names) + '\n').encode() != inventory_bytes:
        raise RuntimeError('Archive member inventory differs from frozen inventory')
    if len(actual_names) != receipt['archive_entries']:
        raise RuntimeError('Archive entry count differs from freeze receipt')
    if receipt['source_commit'] != EXPECTED_COMMIT or receipt['clean_worktree'] is not True:
        raise RuntimeError('Freeze does not bind the expected clean merged source commit')
    if receipt['forbidden_entries'] != []:
        raise RuntimeError('Freeze reports forbidden archive entries')
    if receipt['rights_unchanged_from_predecessor'] is not True:
        raise RuntimeError('Freeze does not affirm unchanged rights/data streams')
    if not archive.is_file():
        raise RuntimeError(f'Frozen archive is missing: {archive}')
    mode = archive.stat().st_mode
    if mode & 0o222:
        raise RuntimeError(f'Frozen archive is writable (mode {mode & 0o777:o})')
    actual_size = archive.stat().st_size
    if actual_size != receipt['size_bytes']:
        raise RuntimeError(f'Archive size mismatch: {actual_size} != {receipt["size_bytes"]}')
    actual_sha = sha256(archive)
    if actual_sha != receipt['sha256']:
        raise RuntimeError(f'Archive SHA-256 mismatch: {actual_sha} != {receipt["sha256"]}')
    extracted_description = Path(receipt['extracted_source']) / 'DESCRIPTION'
    if not extracted_description.is_file():
        raise RuntimeError(f'Frozen extracted DESCRIPTION is missing: {extracted_description}')
    if not CHECK_SCRIPT.is_file():
        raise RuntimeError(f'Native check driver is missing: {CHECK_SCRIPT}')
    return receipt, archive


def stop_process_group(proc):
    """Terminate the complete check process tree, escalating after a short grace."""
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        proc.wait(timeout=TERM_GRACE_SECONDS)
    except subprocess.TimeoutExpired:
        pass
    # Kill the group even if its leader exited: grandchildren may still be alive.
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    proc.wait()


def main():
    receipt, archive = validate_final()
    root = TASK / 'native-macos'
    root.mkdir(exist_ok=False)
    (root / 'frozen-artifact').mkdir()
    (root / 'frozen-artifact' / archive.name).symlink_to(archive)
    (root / 'frozen-deps').mkdir()
    shutil.copyfile(Path(receipt['extracted_source']) / 'DESCRIPTION',
                    root / 'frozen-deps' / 'DESCRIPTION')
    (root / 'platform-evidence').mkdir()

    env = os.environ.copy()
    env.pop('NOT_CRAN', None)
    env.update({
        'OPENBLAS_NUM_THREADS': '1',
        'OMP_NUM_THREADS': '1',
        'MKL_NUM_THREADS': '1',
        'VECLIB_MAXIMUM_THREADS': '1',
        'R_NUM_THREADS': '1',
        'MAKEFLAGS': '-j1',
        '_R_CHECK_FORCE_SUGGESTS_': 'true',
        '_R_CHECK_LIMIT_CORES_': 'true',
        '_R_CHECK_TIMINGS_': '0',
        '_R_CHECK_ALWAYS_LOG_VIGNETTE_OUTPUT_': 'true',
    })
    command = ['Rscript', '--vanilla', str(CHECK_SCRIPT)]
    record = {
        'artifact_sha256': receipt['sha256'],
        'artifact_size_bytes': receipt['size_bytes'],
        'source_commit': receipt['source_commit'],
        'archive_path': str(archive),
        'archive_mode': format(archive.stat().st_mode & 0o777, '04o'),
        'native_check_script': str(CHECK_SCRIPT),
        'native_check_script_sha256': sha256(CHECK_SCRIPT),
        'command': command,
        'cwd': str(root),
        'NOT_CRAN': 'UNSET',
        'force_suggests': True,
        'thread_caps': {key: env[key] for key in (
            'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'R_NUM_THREADS', 'MAKEFLAGS')},
        'timeout_seconds': TIMEOUT_SECONDS,
        'term_grace_seconds': TERM_GRACE_SECONDS,
        'warning_replay': 'not run; requires separate explicit review',
        'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'archive_sha256_before': sha256(archive),
    }
    log_path = root / 'native-check.log'
    started = time.monotonic()
    timed_out = False
    with log_path.open('x') as log:
        proc = subprocess.Popen(command, cwd=root, env=env, stdout=log,
                                stderr=subprocess.STDOUT, start_new_session=True)
        try:
            exit_code = proc.wait(timeout=TIMEOUT_SECONDS)
        except subprocess.TimeoutExpired:
            timed_out = True
            stop_process_group(proc)
            exit_code = 124
    record.update({
        'elapsed_seconds': round(time.monotonic() - started, 3),
        'timed_out': timed_out,
        'exit_code': exit_code,
        'archive_sha256_after': sha256(archive),
        'log': str(log_path),
        'finished_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    })
    (root / 'native-receipt.json').write_text(json.dumps(record, indent=2) + '\n')
    if record['archive_sha256_after'] != receipt['sha256']:
        raise RuntimeError('Frozen archive changed during check; receipt retained')
    if timed_out:
        print(f'NATIVE_CHECK_TIMEOUT after {TIMEOUT_SECONDS}s; process group stopped; logs retained',
              file=sys.stderr)
        return 124
    if exit_code != 0:
        print(f'NATIVE_CHECK_NOT_PASS exit={exit_code}; NOTE/WARNING remains strict failure; logs retained',
              file=sys.stderr)
        return exit_code
    print(f'NATIVE_CHECK_PASS archive_sha256={receipt["sha256"]} elapsed={record["elapsed_seconds"]}s')
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f'PREFLIGHT_REFUSED: {exc}', file=sys.stderr)
        raise SystemExit(2)
