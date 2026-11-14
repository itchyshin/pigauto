"""Bounded source/installed ancillary checks; no source or package writes."""
from pathlib import Path
import datetime, hashlib, json, os, signal, subprocess, time
ROOT=Path(__file__).resolve().parent
IDENTITY=json.loads((ROOT/'identity.json').read_text())
ARCHIVE=Path(IDENTITY['path'])
BASE_LIBRARY=ROOT.parent/'mi-matrix/neither'
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','R_NUM_THREADS','JULIA_NUM_THREADS','TORCH_NUM_THREADS')
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def stop_group(proc):
 for sig in (signal.SIGTERM,signal.SIGKILL):
  try:os.killpg(proc.pid,sig)
  except ProcessLookupError:return
  if sig==signal.SIGTERM:time.sleep(.2)
 if proc.poll() is None:proc.wait(timeout=10)
def run(name):
 assert digest(ARCHIVE)==IDENTITY['sha256'] and ARCHIVE.stat().st_size==IDENTITY['size_bytes']
 script=ROOT/(name+'.R');assert digest(script)==IDENTITY['script_sha256'][name+'.R']
 env=dict(os.environ);env.pop('NOT_CRAN',None);env.update({x:'1' for x in THREADS});env.update(_R_CHECK_FORCE_SUGGESTS_='true',R_LIBS_SITE='/nonexistent')
 lib=ROOT/'absence-library' if name=='absence' else BASE_LIBRARY
 env.update(R_LIBS=str(lib),R_LIBS_USER=str(lib))
 command=['Rscript','--vanilla',str(script),str(ROOT)];start=time.monotonic();proc=None;timed_out=False;exit_code=None
 r={'command':command,'script_sha256':digest(script),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'archive_sha256_before':digest(ARCHIVE),'hard_cap_seconds':180,'environment':{k:env[k] for k in (*THREADS,'R_LIBS','R_LIBS_USER','R_LIBS_SITE','_R_CHECK_FORCE_SUGGESTS_')},'NOT_CRAN':'UNSET'}
 with (ROOT/(name+'.log')).open('x') as log:
  try:
   proc=subprocess.Popen(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True);r['owned_process_group']=proc.pid
   try:exit_code=proc.wait(timeout=180)
   except subprocess.TimeoutExpired:timed_out=True;stop_group(proc);exit_code=proc.wait(timeout=10)
  finally:
   if proc is not None:stop_group(proc)
   r.update(exit_code=exit_code,timed_out=timed_out,elapsed_seconds=round(time.monotonic()-start,3),archive_sha256_after=digest(ARCHIVE),archive_size_bytes_after=ARCHIVE.stat().st_size,script_sha256_after=digest(script))
   r['archive_unchanged']=r['archive_sha256_before']==r['archive_sha256_after']==IDENTITY['sha256'] and r['archive_size_bytes_after']==IDENTITY['size_bytes']
   (ROOT/(name+'-execution.json')).write_text(json.dumps(r,indent=2)+'\n')
 print(name,json.dumps({k:r[k] for k in ['exit_code','timed_out','elapsed_seconds','archive_unchanged']}),flush=True)
 assert r['archive_unchanged']
 return r
if __name__=='__main__':
 records={n:run(n) for n in ['spelling','urls','absence']}
 (ROOT/'execution-summary.json').write_text(json.dumps(records,indent=2)+'\n')
 print('ANCILLARY_CHECKS_COMPLETE_INSPECT_ACTUAL_RESULTS',flush=True)
