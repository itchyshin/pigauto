#!/usr/bin/env python3
"""Run the bounded installed-help/novice check against the exact matrix install."""
from __future__ import annotations
import datetime as dt, hashlib, json, os, signal, subprocess, sys, time
from pathlib import Path

ROOT = Path("/private/tmp/pigauto-final-closure-20261017/successor-b2e6054/installed-novice/attempt2")
ARCHIVE = Path("/private/tmp/pigauto-release-candidate/sha256-a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb/pigauto_0.11.0.tar.gz")
ARCHIVE_SHA = "a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb"
SOURCE = "b2e6054685e6bdedb073ba46e058af2c0a619a18"
FINAL = Path("/private/tmp/pigauto-final-closure-20261017/successor-b2e6054/FINAL.json")
LIB = Path("/private/tmp/pigauto-final-closure-20261017/successor-b2e6054/mi-matrix/neither")
SCRIPT = ROOT / "check-installed.R"
CAPS = {"OPENBLAS_NUM_THREADS":"1", "OMP_NUM_THREADS":"1", "MKL_NUM_THREADS":"1", "R_NUM_THREADS":"1", "VECLIB_MAXIMUM_THREADS":"1", "BLIS_NUM_THREADS":"1", "NUMEXPR_NUM_THREADS":"1", "JULIA_NUM_THREADS":"1", "TORCH_NUM_THREADS":"1"}
HARD_CAP = 300

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()

def main():
    if not ARCHIVE.is_file() or sha(ARCHIVE) != ARCHIVE_SHA:
        raise SystemExit("frozen archive absent or changed; refusing installed check")
    final=json.loads(FINAL.read_text())
    if final.get("sha256")!=ARCHIVE_SHA or final.get("source_commit")!=SOURCE or Path(final.get("path","")).resolve()!=ARCHIVE.resolve():
        raise SystemExit("FINAL.json does not bind the exact archive and approved source")
    if LIB.resolve() != Path("/private/tmp/pigauto-final-closure-20261017/successor-b2e6054/mi-matrix/neither/pigauto").parent or not (LIB/"pigauto").is_dir():
        raise SystemExit("expected exact neither-cell installed package is absent")
    if sha(SCRIPT) != "a74c0d4545efcdea5e94f462edb464ba25bd4cf7081de5d1c61a23b28342fe2f":
        raise SystemExit("check script hash differs from supervisor-pinned bytes")
    outputs=[ROOT/"installed-impute-help.txt",ROOT/"installed-pigauto_report-help.txt",ROOT/"installed-workflow-receipt.json",ROOT/"installed-session.txt",ROOT/"novice-run/pigauto_report.html"]
    if any(p.exists() for p in outputs) or (ROOT/"run.log").exists() or (ROOT/"run-receipt.json").exists():
        raise SystemExit("output already exists; refusing overwrite")
    env=os.environ.copy(); env.pop("NOT_CRAN",None); env.update(CAPS)
    env["R_LIBS_USER"]=str(LIB); env["R_LIBS_SITE"]="/nonexistent"
    cmd=["Rscript","--vanilla",str(SCRIPT)]
    before=sha(ARCHIVE); size_before=ARCHIVE.stat().st_size
    try: package_hash=sha(LIB/"pigauto"/"DESCRIPTION")
    except OSError as e: raise SystemExit(f"cannot hash installed package DESCRIPTION: {e}")
    rec={"status":"RUNNING","source_commit":SOURCE,"final_receipt":str(FINAL),"final_receipt_sha256":sha(FINAL),"archive":str(ARCHIVE),"archive_sha256_before":before,"archive_size_bytes_before":size_before,"installed_library":str(LIB),"installed_package_path":str(LIB/"pigauto"),"installed_description_sha256":package_hash,"check_script":str(SCRIPT),"check_script_sha256":sha(SCRIPT),"supervisor_script":str(Path(__file__).resolve()),"supervisor_script_sha256":sha(Path(__file__).resolve()),"command":cmd,"environment_controls":{**CAPS,"NOT_CRAN":"unset","R_LIBS_USER":str(LIB),"R_LIBS_SITE":"/nonexistent","browser_suppression":"R options only, inside check-installed.R"},"hard_cap_seconds":HARD_CAP,"started_utc":dt.datetime.now(dt.timezone.utc).isoformat()}
    receipt=ROOT/"run-receipt.json"; log=ROOT/"run.log"
    receipt.write_text(json.dumps(rec,indent=2)+"\n")
    start=time.monotonic(); code=1; timeout=False; proc=None
    try:
        with log.open("x",encoding="utf-8") as stream:
            proc=subprocess.Popen(cmd,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
            try: code=proc.wait(timeout=HARD_CAP)
            except subprocess.TimeoutExpired:
                timeout=True; os.killpg(proc.pid,signal.SIGTERM)
                try: proc.wait(timeout=10)
                except subprocess.TimeoutExpired: os.killpg(proc.pid,signal.SIGKILL); proc.wait()
                code=124
    except KeyboardInterrupt:
        if proc is not None and proc.poll() is None:
            os.killpg(proc.pid,signal.SIGTERM)
            try: proc.wait(timeout=10)
            except subprocess.TimeoutExpired: os.killpg(proc.pid,signal.SIGKILL); proc.wait()
        code=130
    finally:
        after=sha(ARCHIVE)
        desc_after=sha(LIB/"pigauto"/"DESCRIPTION")
        output_checks={str(p):p.is_file() for p in outputs}
        passed=code==0 and after==ARCHIVE_SHA and desc_after==package_hash and all(output_checks.values())
        rec.update({"status":"TIMEOUT" if timeout else ("PASS" if passed else "FAIL"),"exit_code":code,"timed_out":timeout,"elapsed_seconds":round(time.monotonic()-start,3),"ended_utc":dt.datetime.now(dt.timezone.utc).isoformat(),"archive_sha256_after":after,"archive_size_bytes_after":ARCHIVE.stat().st_size,"archive_unchanged":after==ARCHIVE_SHA,"installed_description_sha256_after":desc_after,"installed_description_unchanged":desc_after==package_hash,"log":str(log),"log_sha256":sha(log) if log.exists() else None,"outputs":{str(p):{"exists":p.is_file(),"sha256":sha(p) if p.is_file() else None,"bytes":p.stat().st_size if p.is_file() else None} for p in outputs},"output_checks":output_checks})
        receipt.write_text(json.dumps(rec,indent=2)+"\n")
    print(f"status={rec['status']}")
    print(f"receipt={receipt}")
    print(f"log={log}")
    return code

if __name__=="__main__": sys.exit(main())
