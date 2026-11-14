from pathlib import Path
import hashlib, json, os, signal, subprocess, sys, time
root=Path('/private/tmp/pigauto-final-closure-20261017/successor-b2e6054/warning-replay')
script=root/'capture.R'
cmd=['Rscript','--vanilla',str(script),str(root)]
env=os.environ.copy()
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    env[key]='1'
env.pop('NOT_CRAN',None)
env['_R_CHECK_FORCE_SUGGESTS_']='true'
env['R_LIBS_USER']='/Users/z3437171/Library/R/arm64/4.6/library'
env['R_LIBS_SITE']='/nonexistent'
start=time.monotonic()
p=subprocess.Popen(cmd,cwd=root,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,start_new_session=True)
try:
    out,_=p.communicate(timeout=1200)
    code=p.returncode
except subprocess.TimeoutExpired:
    os.killpg(p.pid,signal.SIGTERM)
    try: out,_=p.communicate(timeout=10)
    except subprocess.TimeoutExpired:
        os.killpg(p.pid,signal.SIGKILL); out,_=p.communicate()
    code=124
(root/'capture.log').write_text(out)
record={'command':cmd,'cwd':str(root),'elapsed_seconds':time.monotonic()-start,'timeout_seconds':1200,'exit_code':code,'process_group':p.pid,'killed_own_group_on_timeout':code==124,'capture_script_sha256':hashlib.sha256(script.read_bytes()).hexdigest()}
(root/'supervisor.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
sys.exit(code)
