#!/usr/bin/env python3
"""Install one frozen pigauto tarball into four isolated optional-backend cells."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


CELLS = {
    "neither": "neither",
    "drm-only": "drm-only",
    "gllvm-only": "gllvm-only",
    "both": "real",
}
THREAD_CAPS = {
    "OPENBLAS_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "R_NUM_THREADS": "1",
}
EXECUTED_SOURCES = (
    "script/cran-0.11-integration/check-adapters.R",
    "script/cran-0.11-integration/test-real-backends.R",
    "script/cran-0.11-integration/reload-fit.R",
    "tests/testthat/test-pool-mi-adapter-fixtures.R",
    "tests/testthat/test-mi-provenance.R",
)
EXPECTED_OPTIONAL = {
    "neither": [],
    "drm-only": ["drmTMB"],
    "gllvm-only": ["gllvmTMB"],
    "both": ["drmTMB", "gllvmTMB"],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def clone_cell_dependencies(template: Path, destination: Path) -> None:
    """Clone only dependency symlinks, deliberately excluding template pigauto."""
    if not template.is_dir():
        raise ValueError(f"template cell is not a directory: {template}")
    old_pigauto = template / "pigauto"
    if not old_pigauto.is_dir() or old_pigauto.is_symlink():
        raise ValueError(f"template must contain the prior installed pigauto directory: {old_pigauto}")
    destination.mkdir(parents=True, exist_ok=False)
    for entry in sorted(template.iterdir()):
        if entry.name == "pigauto":
            continue
        if not entry.is_symlink():
            raise ValueError(f"template dependency must be a symlink: {entry}")
        target = entry.resolve(strict=True)
        (destination / entry.name).symlink_to(target, target_is_directory=True)


def run_step(name: str, command: list[str], log_path: Path, env: dict[str, str], archive: Path, expected_sha: str, cwd: Path) -> dict:
    before = sha256(archive)
    if before != expected_sha:
        raise RuntimeError(f"archive hash changed before {name}: {before}")
    started = datetime.now(timezone.utc).isoformat()
    t0 = time.monotonic()
    completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, env=env, check=False, cwd=cwd)
    elapsed = time.monotonic() - t0
    log_path.write_text(completed.stdout, encoding="utf-8")
    after = sha256(archive)
    return {
        "step": name,
        "command": command,
        "exit_code": completed.returncode,
        "started_utc": started,
        "elapsed_seconds": round(elapsed, 3),
        "log": log_path.name,
        "log_sha256": sha256(log_path),
        "archive_sha256_before": before,
        "archive_sha256_after": after,
        "archive_unchanged": before == after == expected_sha,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("expected_sha256")
    parser.add_argument("output_root", type=Path)
    parser.add_argument("template_root", type=Path,
                        help="prior four-cell library tree used only for shared dependency symlinks")
    args = parser.parse_args()
    archive = args.archive.resolve(strict=True)
    output_root = args.output_root.resolve()
    template_root = args.template_root.resolve(strict=True)
    expected_sha = args.expected_sha256.lower()
    if len(expected_sha) != 64 or sha256(archive) != expected_sha:
        raise SystemExit("archive SHA-256 does not match the required frozen-artifact hash")
    if archive.name != "pigauto_0.11.0.tar.gz":
        raise SystemExit(f"unexpected package archive name: {archive.name}")
    if output_root.exists():
        raise SystemExit(f"output root already exists; refusing to overwrite: {output_root}")

    source_root = Path(__file__).resolve().parents[2]
    harness = source_root / "script/cran-0.11-integration/check-adapters.R"
    git = os.environ.get("PIGAUTO_GIT", "git")
    source_commit = subprocess.check_output(
        [git, "-C", str(source_root), "rev-parse", "HEAD"], text=True
    ).strip()
    output_root.mkdir(parents=True)
    manifest = {
        "archive": str(archive),
        "archive_sha256": expected_sha,
        "source_commit": source_commit,
        "working_directory": str(source_root),
        "executed_sources": {
            relative: {"path": relative, "sha256": sha256(source_root / relative)}
            for relative in EXECUTED_SOURCES
        },
        "harness": str(harness),
        "harness_sha256": sha256(harness),
        "runner": str(Path(__file__).resolve()),
        "runner_sha256": sha256(Path(__file__).resolve()),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "r_version": subprocess.check_output(
            ["R", "--version"], text=True
        ).splitlines()[0],
        "cells": [],
    }
    failures = []
    for cell_name, mode in CELLS.items():
        cell = output_root / cell_name
        template = template_root / cell_name
        cell_result = {
            "cell": cell_name,
            "mode": mode,
            "library": str(cell),
            "expected_optional_packages": EXPECTED_OPTIONAL[cell_name],
            "dependency_symlinks": {},
            "steps": [],
        }
        manifest["cells"].append(cell_result)
        try:
            clone_cell_dependencies(template, cell)
            cell_result["dependency_symlinks"] = {
                entry.name: str(entry.resolve(strict=True))
                for entry in sorted(cell.iterdir()) if entry.is_symlink()
            }
        except Exception as error:  # preserve the failure in the manifest
            cell_result["setup_error"] = f"{type(error).__name__}: {error}"
            failures.append(f"{cell_name}: setup failed")
            continue

        env = os.environ.copy()
        env.update(THREAD_CAPS)
        env.update({"R_LIBS_USER": str(cell), "R_LIBS_SITE": "/nonexistent"})
        steps = [
            ("install", ["R", "CMD", "INSTALL", "--no-multiarch", "-l", str(cell), str(archive)]),
            ("fixtures", ["Rscript", "--vanilla", str(harness), "fixtures", "installed", str(cell)]),
            ("real", ["Rscript", "--vanilla", str(harness), mode, "installed", str(cell)]),
        ]
        blocked = False
        for step_name, command in steps:
            log_path = output_root / f"{cell_name}-{step_name}.log"
            if blocked:
                cell_result["steps"].append({"step": step_name, "status": "not_run_due_to_prior_failure"})
                continue
            try:
                result = run_step(step_name, command, log_path, env, archive, expected_sha, source_root)
            except Exception as error:
                result = {"step": step_name, "status": "error", "error": f"{type(error).__name__}: {error}"}
            cell_result["steps"].append(result)
            if result.get("exit_code") != 0 or result.get("archive_unchanged") is False:
                failures.append(f"{cell_name}/{step_name} failed")
                blocked = True

    manifest["status"] = "PASS" if not failures else "FAIL"
    manifest["failures"] = failures
    manifest_path = output_root / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"manifest={manifest_path}")
    print(f"archive_sha256={sha256(archive)}")
    print(f"status={manifest['status']}")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
