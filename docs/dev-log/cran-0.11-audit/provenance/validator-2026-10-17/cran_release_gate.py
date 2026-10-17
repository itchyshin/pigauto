#!/usr/bin/env python3
"""Fail-closed validator for the cross-project CRAN release ledger."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


RUNGS = [
    "source-clean",
    "tarball-clean",
    "platform-clean",
    "submission-ready",
    "submitted",
    "confirmed",
    "incoming-passed",
    "accepted",
    "archived",
    "live-with-check-page",
]

EVIDENCE_BY_RUNG = {
    "source-clean": ["product_contract", "rights_and_consent", "current_cran_policy", "rendered_site"],
    "tarball-clean": ["local_as_cran", "incoming_feasibility", "installed_package", "timing"],
    "platform-clean": ["platform_matrix", "external_logs"],
}

RELEASE_TYPES = {"first_submission", "update", "resubmission"}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def present(value: object) -> bool:
    return isinstance(value, str) and bool(value.strip())


def verify_evidence(name: str, item: object, root: Path | None = None) -> str | None:
    if not isinstance(item, dict):
        return f"evidence.{name} must be a typed evidence object"
    kind, value = item.get("kind"), item.get("value")
    if not present(value):
        return f"evidence.{name} has no value"
    if kind == "file":
        path = Path(value)
        if not path.is_absolute():
            if root is None:
                return f"evidence.{name} relative file path requires a repository root"
            path = root / path
        if not path.is_file():
            return f"evidence.{name} file does not exist: {path}"
        return None
    if kind == "url":
        if not value.startswith("https://") or ".invalid" in value:
            return f"evidence.{name} URL is not an authoritative HTTPS endpoint"
        try:
            request = urllib.request.Request(value, method="HEAD", headers={"User-Agent": "cran-release-gate/1"})
            with urllib.request.urlopen(request, timeout=15) as response:
                if response.status >= 400:
                    return f"evidence.{name} URL returned HTTP {response.status}"
        except (urllib.error.URLError, TimeoutError, ValueError) as exc:
            return f"evidence.{name} URL did not resolve: {exc}"
        return None
    if kind == "rationale":
        if not value.startswith("NOT_APPLICABLE:"):
            return f"evidence.{name} rationale must start NOT_APPLICABLE:"
        return None
    return f"evidence.{name} kind must be file, url, or rationale"


def verify_cran_public_url(name: str, item: object, package: str, root: Path | None = None) -> str | None:
    if not isinstance(item, dict) or item.get("kind") != "url" or not present(item.get("value")):
        return f"{name} must be a canonical CRAN URL evidence object"
    parsed = urllib.parse.urlparse(item["value"])
    if parsed.scheme != "https" or parsed.netloc != "cran.r-project.org":
        return f"{name} must use https://cran.r-project.org"
    allowed = {
        "submission.public_package_evidence": {
            f"/package={package}",
            f"/web/packages/{package}/index.html",
        },
        "submission.public_check_page_evidence": {
            f"/web/checks/check_results_{package}.html",
        },
    }
    if parsed.path not in allowed[name]:
        return f"{name} does not match the canonical CRAN path for {package}"
    return verify_evidence(name, item, root)


def validate(ledger: dict, root: Path | None = None) -> list[str]:
    errors: list[str] = []
    claim = ledger.get("status_claim")
    if claim not in RUNGS:
        return [f"status_claim must be one of: {', '.join(RUNGS)}"]
    rank = RUNGS.index(claim)

    for key in ("package", "version", "release_type"):
        if not present(ledger.get(key)):
            errors.append(f"missing {key}")
    if ledger.get("release_type") not in RELEASE_TYPES:
        errors.append(f"release_type must be one of: {', '.join(sorted(RELEASE_TYPES))}")

    artifact = ledger.get("artifact", {})
    path_text = artifact.get("path", "")
    path = Path(path_text) if present(path_text) else None
    if path is not None and not path.is_absolute() and root is not None:
        path = root / path
    if rank >= RUNGS.index("tarball-clean"):
        if path is None or (not path.is_absolute() and root is None) or not path.is_file():
            errors.append("artifact.path must resolve to an existing file with a root for relative paths")
        else:
            expected_hash = artifact.get("sha256", "")
            actual_hash = digest(path)
            if expected_hash != actual_hash:
                errors.append("artifact SHA-256 does not match the selected file")
            if artifact.get("size_bytes") != path.stat().st_size:
                errors.append("artifact size_bytes does not match the selected file")
        if not present(artifact.get("source_commit")):
            errors.append("missing artifact source_commit")
        if artifact.get("clean_worktree") is not True:
            errors.append("candidate was not built from a recorded clean worktree")
        inventory_error = verify_evidence("artifact.inventory_evidence", artifact.get("inventory_evidence"), root)
        if inventory_error:
            errors.append(inventory_error)
        if artifact.get("forbidden_entries") != []:
            errors.append("tarball inventory contains forbidden entries")

    evidence = ledger.get("evidence", {})
    for rung, keys in EVIDENCE_BY_RUNG.items():
        if rank >= RUNGS.index(rung):
            for key in keys:
                error = verify_evidence(key, evidence.get(key), root)
                if error:
                    errors.append(error)

    profile = ledger.get("profile", {})
    conditional = {
        "reverse_dependencies": "reverse_dependencies",
        "compiled_code": "compiled_diagnostics",
        "system_or_external_services": "external_service_behaviour",
        "large_data_or_vignettes": "size_and_vignette_budget",
    }
    if rank >= RUNGS.index("submission-ready"):
        for flag, key in conditional.items():
            item = evidence.get(key)
            error = verify_evidence(key, item, root)
            if error:
                errors.append(error)
                continue
            if profile.get(flag) is True and item.get("kind") == "rationale":
                errors.append(f"profile.{flag}=true requires evidence.{key}, not NOT_APPLICABLE")
            if profile.get(flag) not in (True, False):
                errors.append(f"profile.{flag} must be true or false")

    panel = ledger.get("panel", {})
    if rank >= RUNGS.index("submission-ready"):
        for reviewer in ("grace", "rose", "pat"):
            if panel.get(reviewer) != "READY":
                errors.append(f"panel.{reviewer} must be READY")

    submission = ledger.get("submission", {})
    if rank >= RUNGS.index("submitted"):
        if submission.get("uploaded_sha256") != artifact.get("sha256"):
            errors.append("uploaded_sha256 does not match frozen artifact")
        if submission.get("uploaded_size_bytes") != artifact.get("size_bytes"):
            errors.append("uploaded_size_bytes does not match frozen artifact")
        for key in ("timestamp", "submission_id"):
            if not present(submission.get(key)):
                errors.append(f"submitted status requires submission.{key}")
    state_evidence = {
        "confirmed": ["confirmation_evidence"],
        "incoming-passed": ["incoming_evidence"],
        "accepted": ["acceptance_evidence"],
        "archived": ["archive_evidence"],
        "live-with-check-page": ["public_package_evidence", "public_check_page_evidence"],
    }
    for rung, keys in state_evidence.items():
        if rank >= RUNGS.index(rung):
            for key in keys:
                full_name = f"submission.{key}"
                if key in ("public_package_evidence", "public_check_page_evidence"):
                    error = verify_cran_public_url(full_name, submission.get(key), ledger.get("package", ""), root)
                else:
                    error = verify_evidence(full_name, submission.get(key), root)
                if error:
                    errors.append(error)
    return errors


def valid_fixture(path: Path, evidence_file: Path) -> dict:
    sha = digest(path)
    return {
        "package": "fixture",
        "version": "0.0.1",
        "release_type": "first_submission",
        "profile": {
            "compiled_code": False,
            "reverse_dependencies": False,
            "system_or_external_services": False,
            "large_data_or_vignettes": False,
        },
        "artifact": {
            "path": str(path),
            "sha256": sha,
            "size_bytes": path.stat().st_size,
            "source_commit": "deadbeef",
            "clean_worktree": True,
            "inventory_evidence": {"kind": "file", "value": str(evidence_file)},
            "forbidden_entries": [],
        },
        "evidence": {
            key: {"kind": "file", "value": str(evidence_file)}
            for keys in EVIDENCE_BY_RUNG.values() for key in keys
        } | {
            "reverse_dependencies": {"kind": "rationale", "value": "NOT_APPLICABLE: fixture"},
            "compiled_diagnostics": {"kind": "rationale", "value": "NOT_APPLICABLE: fixture"},
            "external_service_behaviour": {"kind": "rationale", "value": "NOT_APPLICABLE: fixture"},
            "size_and_vignette_budget": {"kind": "rationale", "value": "NOT_APPLICABLE: fixture"},
        },
        "panel": {"grace": "READY", "rose": "READY", "pat": "READY"},
        "status_claim": "submission-ready",
        "submission": {},
    }


def selftest() -> int:
    with tempfile.TemporaryDirectory() as tmp:
        artifact = Path(tmp) / "fixture_0.0.1.tar.gz"
        artifact.write_bytes(b"frozen candidate")
        evidence_file = Path(tmp) / "evidence.log"
        evidence_file.write_text("verified fixture\n", encoding="utf-8")
        base = valid_fixture(artifact, evidence_file)
        if validate(base):
            print("FAIL: valid fixture rejected")
            return 1
        source_only = copy.deepcopy(base)
        source_only["status_claim"] = "source-clean"
        source_only["artifact"] = {}
        source_only["evidence"] = {
            key: base["evidence"][key] for key in EVIDENCE_BY_RUNG["source-clean"]
        }
        source_only["panel"] = {}
        if validate(source_only):
            print("FAIL: rung-specific source-clean fixture rejected")
            return 1
        print("PASS: source-clean does not require later-rung evidence")
        relative = copy.deepcopy(base)
        relative["artifact"]["path"] = artifact.name
        relative["artifact"]["inventory_evidence"]["value"] = evidence_file.name
        for item in relative["evidence"].values():
            if item.get("kind") == "file":
                item["value"] = evidence_file.name
        if validate(relative, root=Path(tmp)):
            print("FAIL: repository-relative artifact and evidence paths rejected")
            return 1
        print("PASS: repository-relative artifact and evidence paths resolve from root")
        relative_artifact_only = copy.deepcopy(base)
        relative_artifact_only["artifact"]["path"] = artifact.name
        previous_cwd = Path.cwd()
        try:
            os.chdir(tmp)
            if not validate(relative_artifact_only):
                print("FAIL: relative artifact accepted without repository root")
                return 1
        finally:
            os.chdir(previous_cwd)
        print("PASS: relative artifact without repository root fails closed")
        missing_relative = copy.deepcopy(relative)
        missing_relative["evidence"]["rights_and_consent"]["value"] = "missing.log"
        if not validate(missing_relative, root=Path(tmp)):
            print("FAIL: missing repository-relative evidence accepted")
            return 1
        print("PASS: missing repository-relative evidence fails closed")
        if not validate(relative, root=Path(tmp) / "wrong-root"):
            print("FAIL: wrong repository root accepted")
            return 1
        print("PASS: wrong repository root fails closed")
        mutations = {
            "unlicensed-data": lambda x: x["evidence"].update(rights_and_consent={"kind": "file", "value": "/missing"}),
            "predecessor-hash": lambda x: x["artifact"].update(sha256="0" * 64),
            "forbidden-tarball-entry": lambda x: x["artifact"].update(forbidden_entries=["output/x.rds"]),
            "incoming-preflight-missing": lambda x: x["evidence"].update(incoming_feasibility={"kind": "file", "value": "/missing"}),
            "missing-timing": lambda x: x["evidence"].update(timing={"kind": "file", "value": "/missing"}),
            "unconfirmed-upload": lambda x: x.update(status_claim="confirmed"),
            "public-page-404": lambda x: x.update(status_claim="live-with-check-page"),
            "compiled-gate-omitted": lambda x: (
                x["profile"].update(compiled_code=True),
                x["evidence"].update(compiled_diagnostics={"kind": "rationale", "value": "NOT_APPLICABLE: fixture"}),
            ),
            "fictional-live-url": lambda x: (
                x.update(status_claim="live-with-check-page"),
                x["submission"].update(
                    uploaded_sha256=x["artifact"]["sha256"],
                    uploaded_size_bytes=x["artifact"]["size_bytes"],
                    timestamp="now", submission_id="id",
                    confirmation_evidence={"kind": "file", "value": str(evidence_file)},
                    incoming_evidence={"kind": "file", "value": str(evidence_file)},
                    acceptance_evidence={"kind": "file", "value": str(evidence_file)},
                    archive_evidence={"kind": "file", "value": str(evidence_file)},
                    public_package_evidence={"kind": "url", "value": "https://package.invalid/404"},
                    public_check_page_evidence={"kind": "url", "value": "https://checks.invalid/404"},
                ),
            ),
            "reachable-non-cran-live-url": lambda x: (
                x.update(status_claim="live-with-check-page"),
                x["submission"].update(
                    uploaded_sha256=x["artifact"]["sha256"],
                    uploaded_size_bytes=x["artifact"]["size_bytes"],
                    timestamp="now", submission_id="id",
                    confirmation_evidence={"kind": "file", "value": str(evidence_file)},
                    incoming_evidence={"kind": "file", "value": str(evidence_file)},
                    acceptance_evidence={"kind": "file", "value": str(evidence_file)},
                    archive_evidence={"kind": "file", "value": str(evidence_file)},
                    public_package_evidence={"kind": "url", "value": "https://example.com/"},
                    public_check_page_evidence={"kind": "url", "value": "https://example.com/"},
                ),
            ),
            "panel-not-ready": lambda x: x["panel"].update(rose="NOT READY"),
            "invalid-release-type": lambda x: x.update(release_type="preview"),
            "large-vignette-budget-omitted": lambda x: (
                x["profile"].update(large_data_or_vignettes=True),
                x["evidence"].update(size_and_vignette_budget={"kind": "rationale", "value": "NOT_APPLICABLE: fixture"}),
            ),
        }
        for name, mutate in mutations.items():
            candidate = copy.deepcopy(base)
            mutate(candidate)
            if not validate(candidate):
                print(f"FAIL: negative control passed: {name}")
                return 1
            print(f"PASS: negative control failed closed: {name}")
    print("PASS: CRAN release gate selftest")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("ledger", nargs="?")
    parser.add_argument("--root", type=Path, help="repository root for relative artifact and evidence paths")
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    if args.selftest:
        return selftest()
    if not args.ledger:
        parser.error("provide a ledger JSON file or --selftest")
    ledger_path = Path(args.ledger).resolve()
    root = args.root.resolve() if args.root else next(
        (parent for parent in (ledger_path.parent, *ledger_path.parents)
         if (parent / ".git").exists() or (parent / "DESCRIPTION").is_file()),
        ledger_path.parent,
    )
    with ledger_path.open(encoding="utf-8") as stream:
        errors = validate(json.load(stream), root=root)
    if errors:
        print("NOT READY")
        for error in errors:
            print(f"- {error}")
        return 1
    print("READY FOR CLAIMED RUNG")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
