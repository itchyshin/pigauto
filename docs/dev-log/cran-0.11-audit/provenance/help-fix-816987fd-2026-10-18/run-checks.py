from pathlib import Path
import datetime
import json
import os
import signal
import subprocess
import time

ROOT = Path('/Users/z3437171/.codex/worktrees/cran-011-final-artifact/pigauto')
OUT = Path(__file__).resolve().parent
env = os.environ.copy()
env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', R_NUM_THREADS='1')
env.pop('NOT_CRAN', None)
env['_R_CHECK_FORCE_SUGGESTS_'] = 'true'
(OUT/'installed-library').mkdir(exist_ok=False)
records = []
for name, command in [
    ('source-AST-and-existing-report-tests', ['Rscript','--vanilla',str(OUT/'check-source.R')]),
    ('provisional-install', ['R','CMD','INSTALL','--library='+str(OUT/'installed-library'),str(ROOT)]),
    ('installed-help-and-novice', ['Rscript','--vanilla',str(OUT/'check-installed.R')]),
]:
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    t = time.monotonic()
    with (OUT/(name+'.log')).open('x') as log:
        proc = subprocess.Popen(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        try:
            code = proc.wait(timeout=300)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
            code = 124
    records.append({'step':name,'command':command,'started_utc':started,'elapsed_seconds':round(time.monotonic()-t,3),'exit_code':code,'hard_cap_seconds':300})
    (OUT/'execution.json').write_text(json.dumps(records,indent=2)+'\n')
    print(name,code,records[-1]['elapsed_seconds'],flush=True)
    if code:
        raise SystemExit(code)
print('HELP_FIX_BOUNDED_SOURCE_AND_INSTALLED_CHECKS_PASS')
