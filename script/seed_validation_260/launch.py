"""Bounded Totoro launcher: six CPU fits, two torch threads each, 40 min total."""
import json
import csv
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

home = Path.home()
source = home / "pigauto_seed260/source"
runner = home / "pigauto_seed260/runner"
archive = home / "sim-archive/pigauto/seed260-20261108"
out = archive / "production"
out.mkdir(exist_ok=True)
if (source / "SOURCE_SHA").read_text().strip() != "033649b93115c6d8191e791f06b1ddf52eb82d19":
    raise RuntimeError("Wrong source export: no worker launched")
jobs = [dict(fold=f, seed=s, status="queued") for s in [0, *range(51001,51007)] for f in range(1,6)]
active = {}
started = time.monotonic()
environment = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
subprocess.run(["Rscript", str(runner / "summarize.R"),
                str(archive / "input/prepared.rds"), str(out), str(archive / "preflight"),
                "--partial"], env=environment, check=True, timeout=120)
prior = list(csv.DictReader((archive / "preflight/attempts.csv").open()))
completed = {(int(r["fold"]), int(r["seed"])) for r in prior if r["status"] == "ok"}
for j in jobs:
    if (j["fold"], j["seed"]) in completed:
        j.update(status="ok", reused=True)
if "--plan-only" in sys.argv:
    print(json.dumps(dict(reused=len(completed), queued=35-len(completed))))
    sys.exit(0)


def interrupted(signum, frame):
    raise RuntimeError(f"Launcher received signal {signum}")


def save():
    (out / "launcher.json").write_text(json.dumps(jobs, indent=2) + "\n")


signal.signal(signal.SIGTERM, interrupted)
signal.signal(signal.SIGINT, interrupted)
try:
    while active or any(j["status"] == "queued" for j in jobs):
        if time.monotonic() - started > 2400:
            raise RuntimeError("40-minute campaign estimate exceeded; stop and re-estimate")
        for j in jobs:
            if len(active) >= 6:
                break
            if j["status"] != "queued":
                continue
            name = f"fold{j['fold']}-seed{j['seed']}-epochs2000"
            if (out / f"{name}.rds").exists():
                raise RuntimeError(f"Output exists; inspect before resuming: {name}")
            log = (out / f"{name}.log").open("w")
            cmd = ["Rscript", str(runner / "run.R"), str(source),
                   str(archive / "input/prepared.rds"), str(out),
                   str(j["fold"]), str(j["seed"]), "2000"]
            p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT,
                                 env=environment, start_new_session=True)
            j.update(status="running", pid=p.pid)
            active[p.pid] = (p, log, j)
        for pid, (p, log, j) in list(active.items()):
            code = p.poll()
            if code is not None:
                log.close()
                j.update(status="ok" if code == 0 else "failed", exit_code=code)
                del active[pid]
        save()
        time.sleep(2)
    if any(j["status"] != "ok" for j in jobs):
        raise RuntimeError("One or more fits failed; retain all attempts")
    subprocess.run(["Rscript", str(runner / "summarize.R"),
                    str(archive / "input/prepared.rds"), str(out), str(archive / "summary")],
                   env=environment, check=True, timeout=120)
finally:
    for p, log, j in active.values():
        try:
            os.killpg(p.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    for p, log, j in active.values():
        try:
            p.wait(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(p.pid, signal.SIGKILL)
            p.wait()
        log.close()
        j.update(status="stopped", exit_code=p.returncode)
    save()
