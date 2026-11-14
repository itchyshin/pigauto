#!/usr/bin/env python3
"""Fail-closed supervisor for the successor frozen-archive MI matrix.

Preparation is not execution. A run requires both the final freeze receipt and
an explicit --root-authorized flag supplied only after parent authorization.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path


TASK = Path("/private/tmp/pigauto-final-closure-20261017/successor-b2e6054")
FINAL = TASK / "FINAL.json"
OUTPUT = TASK / "mi-matrix"
TEMPLATE = Path(
    "/private/tmp/pigauto-release-candidate/sha256-"
    "14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb/mi-matrix"
)
OLD_SHA = "14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb"
SOURCE_COMMIT = "b2e6054685e6bdedb073ba46e058af2c0a619a18"
DRIVER = Path("/private/tmp/pigauto-pr236-validator-20261017")
MATRIX = DRIVER / "script/cran-0.11-integration/run-exact-artifact-matrix.py"
HARNESS = DRIVER / "script/cran-0.11-integration/check-adapters.R"
GIT = Path("/Library/Developer/CommandLineTools/usr/bin/git")
TIMEOUT_SECONDS = 20 * 60
THREAD_CAPS = {
    "OPENBLAS_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "R_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
    "BLIS_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "JULIA_NUM_THREADS": "1",
    "TORCH_NUM_THREADS": "1",
}
EXECUTED_SOURCES = (
    "script/cran-0.11-integration/check-adapters.R",
    "script/cran-0.11-integration/test-real-backends.R",
    "script/cran-0.11-integration/reload-fit.R",
    "tests/testthat/test-pool-mi-adapter-fixtures.R",
    "tests/testthat/test-mi-provenance.R",
)
EXPECTED_MATRIX_SHA256 = "6d6cc86d7c17f1e5e91c5e4a8e9841eef2dadbc793c2775657874fe79cda2bf7"
EXPECTED_HARNESS_SHA256 = "d2b0b5f01e3470ce1f0f0d2116cadbd1a9a2beffcd864a5644e71c0b7d53fb7d"
EXPECTED_EXECUTED_SOURCE_SHA256 = {
    "script/cran-0.11-integration/check-adapters.R": "d2b0b5f01e3470ce1f0f0d2116cadbd1a9a2beffcd864a5644e71c0b7d53fb7d",
    "script/cran-0.11-integration/test-real-backends.R": "7edb499a09820ecd9c8efc6799f2adf4233c824f064bc88462b9f1720098b6be",
    "script/cran-0.11-integration/reload-fit.R": "8db40649b97f074eefb71cfd591453d9540f8b4518b464ac4a6d13472f68401d",
    "tests/testthat/test-pool-mi-adapter-fixtures.R": "1599e96a85be78afec4fe05fc23cc02f9c0f5ac99e877bb31ff87a358e20a0b1",
    "tests/testthat/test-mi-provenance.R": "585299f93cb4c83013337d1043035bf0dbe7af40389c1dc3075507391591f67a",
}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def git(*args: str) -> str:
    return subprocess.check_output([str(GIT), *args], text=True).strip()


def validate_template() -> dict:
    old_matrix = TEMPLATE
    old_manifest = old_matrix / "manifest.json"
    if not old_manifest.is_file():
        raise RuntimeError(f"missing retained predecessor matrix manifest: {old_manifest}")
    old = json.loads(old_manifest.read_text(encoding="utf-8"))
    if old.get("archive_sha256") != OLD_SHA:
        raise RuntimeError("template manifest is not bound to the withdrawn 14ee predecessor")
    if old.get("status") != "PASS" or len(old.get("cells", [])) != 4:
        raise RuntimeError("retained predecessor template lacks a complete four-cell passing matrix")
    for old_cell in old["cells"]:
        if len(old_cell.get("steps", [])) != 3 or any(
            step.get("exit_code") != 0 or step.get("archive_unchanged") is not True
            for step in old_cell["steps"]
        ):
            raise RuntimeError(f"predecessor template has an incomplete or failed cell: {old_cell.get('cell')}")
    cells = {}
    for name in ("neither", "drm-only", "gllvm-only", "both"):
        cell = old_matrix / name
        old_package = cell / "pigauto"
        if not old_package.is_dir() or old_package.is_symlink():
            raise RuntimeError(f"template does not contain a real predecessor package tree: {old_package}")
        dep_lines = []
        dep_count = 0
        for entry in sorted(cell.iterdir()):
            if entry.name == "pigauto":
                continue
            if not entry.is_symlink():
                raise RuntimeError(f"template has non-symlink non-pigauto entry: {entry}")
            dep_lines.append(f"{entry.name}\t{entry.resolve(strict=True)}")
            dep_count += 1
        cells[name] = {
            "predecessor_package_omitted_from_clone": True,
            "dependency_symlink_count": dep_count,
            "dependency_symlink_manifest_sha256": hashlib.sha256(
                ("\n".join(dep_lines) + "\n").encode("utf-8")
            ).hexdigest(),
        }
    return {"template_root": str(old_matrix), "template_manifest_sha256": digest(old_manifest),
            "template_status": old["status"], "template_cells": len(old["cells"]),
            "template_step_count": sum(len(c["steps"]) for c in old["cells"]),
            "template_archive_sha256": old.get("archive_sha256"), "cells": cells}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root-authorized", action="store_true",
        help="explicit authorization gate; omit during preparation/review",
    )
    args = parser.parse_args()
    if not args.root_authorized:
        parser.error("refusing to launch: wait for FINAL.json and explicit parent authorization, then pass --root-authorized")

    if not GIT.is_file():
        raise SystemExit(f"required Command Line Tools git is absent: {GIT}")
    if not FINAL.is_file():
        raise SystemExit(f"FINAL freeze receipt is absent; no matrix launched: {FINAL}")
    receipt = json.loads(FINAL.read_text(encoding="utf-8"))
    if receipt.get("source_commit") != SOURCE_COMMIT or receipt.get("clean_worktree") is not True:
        raise SystemExit("FINAL receipt does not bind a clean archive to the approved merged source commit")
    archive = Path(receipt.get("path", "")).resolve(strict=True)
    expected_sha = str(receipt.get("sha256", "")).lower()
    if archive.name != "pigauto_0.11.0.tar.gz" or len(expected_sha) != 64:
        raise SystemExit("FINAL receipt has an invalid archive path or SHA-256")
    if expected_sha == OLD_SHA or digest(archive) != expected_sha:
        raise SystemExit("FINAL archive is the withdrawn predecessor or does not match FINAL.json")
    if archive.stat().st_size != receipt.get("size_bytes"):
        raise SystemExit("FINAL archive size does not match FINAL.json")
    if archive.stat().st_mode & 0o222:
        raise SystemExit("FINAL archive is writable; expected the content-addressed frozen copy")
    expected_parent = Path("/private/tmp/pigauto-release-candidate") / f"sha256-{expected_sha}"
    if archive.parent != expected_parent:
        raise SystemExit("FINAL archive is not at its content-addressed freeze path")
    if OUTPUT.exists():
        raise SystemExit(f"successor output already exists; refusing overwrite: {OUTPUT}")
    template_receipt = validate_template()

    source_commit = git("-C", str(DRIVER), "rev-parse", "HEAD")
    dirty = git("-C", str(DRIVER), "status", "--porcelain")
    executed_hashes = {rel: digest(DRIVER / rel) for rel in EXECUTED_SOURCES}
    harness_sha = digest(HARNESS)
    matrix_sha = digest(MATRIX)
    if matrix_sha != EXPECTED_MATRIX_SHA256 or harness_sha != EXPECTED_HARNESS_SHA256:
        raise SystemExit("reviewed matrix runner or harness bytes changed since predecessor review")
    if executed_hashes != EXPECTED_EXECUTED_SOURCE_SHA256:
        raise SystemExit("an executed matrix source differs from the reviewed source bytes")
    wrapper_sha = digest(Path(__file__).resolve())
    env = os.environ.copy()
    env.pop("NOT_CRAN", None)
    env["PIGAUTO_GIT"] = str(GIT)
    env.update(THREAD_CAPS)
    command = [sys.executable, str(MATRIX), str(archive), expected_sha,
               str(OUTPUT), str(TEMPLATE)]
    metadata = {
        "status": "RUNNING",
        "started_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "hard_cap_seconds": TIMEOUT_SECONDS,
        "command": command,
        "environment_controls": {**THREAD_CAPS, "NOT_CRAN": "unset", "PIGAUTO_GIT": str(GIT)},
        "archive": str(archive),
        "archive_sha256_before": digest(archive),
        "archive_size_bytes_before": archive.stat().st_size,
        "final_receipt_sha256": digest(FINAL),
        "archive_source_commit": SOURCE_COMMIT,
        "producer_head": source_commit,
        "producer_dirty_status": dirty,
        "runner_sha256": matrix_sha,
        "harness_sha256": harness_sha,
        "wrapper_sha256": wrapper_sha,
        "executed_sources_sha256": executed_hashes,
        "dependency_template": template_receipt,
        "output_root": str(OUTPUT),
    }
    receipt_path = TASK / "mi-matrix-supervisor.json"
    log_path = TASK / "mi-matrix-supervisor.log"
    if receipt_path.exists() or log_path.exists():
        raise SystemExit("supervisor receipt/log already exists; refusing overwrite")
    receipt_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    started = time.monotonic()
    exit_code = 1
    timed_out = False
    proc = None
    try:
        with log_path.open("x", encoding="utf-8") as log:
            proc = subprocess.Popen(command, cwd=DRIVER, env=env, stdout=log,
                                    stderr=subprocess.STDOUT, start_new_session=True)
            try:
                exit_code = proc.wait(timeout=TIMEOUT_SECONDS)
            except subprocess.TimeoutExpired:
                timed_out = True
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
                exit_code = 124
    except KeyboardInterrupt:
        if proc is not None and proc.poll() is None:
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
        exit_code = 130
        metadata["interrupted"] = True
    finally:
        archive_unchanged = digest(archive) == expected_sha
        code_hashes_after = {rel: digest(DRIVER / rel) for rel in EXECUTED_SOURCES}
        code_unchanged = code_hashes_after == executed_hashes
        producer_head_after = git("-C", str(DRIVER), "rev-parse", "HEAD")
        producer_dirty_after = git("-C", str(DRIVER), "status", "--porcelain")
        producer_identity_unchanged = producer_head_after == source_commit
        producer_status_unchanged = producer_dirty_after == dirty
        matrix_status = "TIMEOUT" if timed_out else ("PASS" if exit_code == 0 else "FAIL")
        if matrix_status == "PASS" and (
            not archive_unchanged or not code_unchanged
            or not producer_identity_unchanged or not producer_status_unchanged
        ):
            matrix_status = "FAIL_PROVENANCE_CHANGED"
        metadata.update({
            "status": matrix_status,
            "exit_code": exit_code,
            "timed_out": timed_out,
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "ended_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "archive_sha256_after": digest(archive),
            "archive_size_bytes_after": archive.stat().st_size,
            "archive_unchanged": archive_unchanged,
            "producer_head_after": producer_head_after,
            "producer_identity_unchanged": producer_identity_unchanged,
            "producer_dirty_status_after": producer_dirty_after,
            "producer_status_unchanged": producer_status_unchanged,
            "executed_sources_sha256_after": code_hashes_after,
            "executed_sources_unchanged": code_unchanged,
            "supervisor_log": str(log_path),
            "supervisor_log_sha256": digest(log_path) if log_path.exists() else None,
            "matrix_manifest": str(OUTPUT / "manifest.json"),
            "matrix_manifest_sha256": digest(OUTPUT / "manifest.json")
            if (OUTPUT / "manifest.json").is_file() else None,
        })
        receipt_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"status={metadata['status']}")
    print(f"supervisor_receipt={receipt_path}")
    print(f"supervisor_log={log_path}")
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
