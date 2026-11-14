from pathlib import Path
import datetime, hashlib, json, os, signal, subprocess, sys, time

root = Path(__file__).resolve().parent
meta = json.loads((root / 'input.json').read_text())
archive = Path(meta['artifact']['path'])
started = time.monotonic()
records = []
TOTAL = int(meta['total_timeout_seconds'])
PER = int(meta['per_vignette_timeout_seconds'])

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def kill_group(proc):
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        proc.wait(timeout=10)
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    proc.wait()

assert sha256(archive) == meta['artifact']['sha256']
assert archive.stat().st_size == meta['artifact']['size_bytes']
assert (root / 'source' / 'pigauto').is_dir()
assert (Path(meta['library']) / 'pigauto').is_dir()
env = os.environ.copy()
env.pop('NOT_CRAN', None)
env.update({
    'R_LIBS': meta['library'], 'R_LIBS_USER': meta['library'], 'R_LIBS_SITE': '/nonexistent',
    'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
    'VECLIB_MAXIMUM_THREADS': '1', 'R_NUM_THREADS': '1', 'MAKEFLAGS': '-j1',
    '_R_CHECK_FORCE_SUGGESTS_': 'true',
})
records_path = root / 'execution.json'
for name in sorted(meta['source_rmd_sha256']):
    remaining = TOTAL - (time.monotonic() - started)
    if remaining <= 0:
        records.append({'file': name, 'exit': 'NOT_RUN_TOTAL_CAP'})
        records_path.write_text(json.dumps(records, indent=2) + '\n')
        break
    cmd = ['Rscript', '--vanilla', str(root / 'render-one.R'), str(root), name]
    stamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    begin = time.monotonic()
    record = {'file': name, 'command': cmd,
              'render_script_sha256': sha256(root / 'render-one.R'),
              'source_rmd_sha256': meta['source_rmd_sha256'][name],
              'start_utc': stamp, 'timeout_seconds': min(PER, max(1, int(remaining))),
              'thread_caps': {k: env[k] for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','R_NUM_THREADS','MAKEFLAGS')},
              'NOT_CRAN': 'UNSET', 'cwd': str(root)}
    log_path = root / (name + '.log')
    with log_path.open('x') as out:
        proc = subprocess.Popen(cmd, env=env, cwd=root, stdout=out,
                                stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = proc.wait(timeout=min(PER, max(1, remaining)))
            status = code
        except subprocess.TimeoutExpired:
            kill_group(proc)
            status = 'TIMEOUT'
    record.update({'process_elapsed_seconds': round(time.monotonic() - begin, 3),
                   'exit': status, 'log': str(log_path),
                   'log_sha256': sha256(log_path),
                   'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat()})
    records.append(record)
    records_path.write_text(json.dumps(records, indent=2) + '\n')
    print(json.dumps(record), flush=True)
    if status != 0:
        break
post = sha256(archive)
source_after = {name: sha256(root / 'source' / 'pigauto' / 'vignettes' / name)
                for name in meta['source_rmd_sha256']}
identity = {'archive_sha256_before': meta['artifact']['sha256'], 'archive_sha256_after': post,
            'archive_unchanged': post == meta['artifact']['sha256'],
            'source_rmd_sha256_before': meta['source_rmd_sha256'],
            'source_rmd_sha256_after': source_after,
            'source_rmds_unchanged': source_after == meta['source_rmd_sha256'],
            'total_elapsed_seconds': round(time.monotonic() - started, 3),
            'hard_cap_seconds': TOTAL, 'completed_vignettes': len(records),
            'all_six_success': len(records) == 6 and all(r.get('exit') == 0 for r in records)}
(root / 'identity-after.json').write_text(json.dumps(identity, indent=2) + '\n')
if not identity['archive_unchanged'] or not identity['source_rmds_unchanged']:
    raise SystemExit('IDENTITY_MISMATCH; see identity-after.json')
if not identity['all_six_success']:
    raise SystemExit('VIGNETTE_DIAGNOSTIC_NOT_PASS; preserve logs and report failure')
print('ALL_SIX_SUCCESSOR_VIGNETTE_DIAGNOSTICS_PASS', identity['total_elapsed_seconds'], flush=True)
