"""Verify the frozen CRAN candidate before and after direct platform checks."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import tarfile
import tempfile
import urllib.request


URL = 'http://builder2.rhub.io/api/-/package/062a7517bd61fe3fd188ab5f0fd12254'
SHA = 'a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb'
SIZE = 5137837
INVENTORY = 'd7b9b0e59ea051b3a22d7f347d8d925e6265b076b238ac893410e79d5b61329d'
ARCHIVE = Path("frozen-artifact/pigauto_0.11.0.tar.gz")
EVIDENCE = Path("platform-evidence")


def identity(path, sha=SHA, size=SIZE, inventory=INVENTORY, entries=250):
    raw = path.read_bytes()
    actual_sha = hashlib.sha256(raw).hexdigest()
    if actual_sha != sha or len(raw) != size:
        raise ValueError("frozen archive SHA-256 or byte size mismatch")
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        names = [m.name.rstrip("/") + "/" if m.isdir() else m.name for m in members]
        inventory_sha = hashlib.sha256(("\n".join(sorted(names)) + "\n").encode()).hexdigest()
        if len(names) != entries or inventory_sha != inventory:
            raise ValueError("frozen archive inventory mismatch")
        desc = archive.extractfile("pigauto/DESCRIPTION").read()
    for field, expected in (("Package", "pigauto"), ("Version", "0.11.0")):
        match = re.search(rf"^{field}:\s*(\S+)\s*$", desc.decode(), re.MULTILINE)
        if match is None or match.group(1) != expected:
            raise ValueError(f"unexpected {field}")
    return {"archive": str(path.resolve()), "sha256": actual_sha, "size_bytes": len(raw),
            "inventory_sha256": inventory_sha, "archive_entries": len(names),
            "package": "pigauto", "version": "0.11.0"}, desc


def record(phase):
    data, desc = identity(ARCHIVE)
    data.update(phase=phase, platform=platform.platform(),
                workflow_sha=os.getenv("GITHUB_SHA"), run_id=os.getenv("GITHUB_RUN_ID"),
                run_attempt=os.getenv("GITHUB_RUN_ATTEMPT"))
    EVIDENCE.mkdir(exist_ok=True)
    (EVIDENCE / f"identity-{phase}.json").write_text(json.dumps(data, indent=2) + "\n")
    print("EXACT_ARTIFACT_" + phase.upper(), json.dumps(data), flush=True)
    return desc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=["download", "pre", "post", "selftest"])
    args = parser.parse_args()
    if args.phase == "selftest":
        # The real candidate is the positive control; mutations never alter it.
        identity(ARCHIVE)
        controls = [{"sha": "0" * 64}, {"size": SIZE + 1},
                    {"inventory": "0" * 64}, {"entries": 249}]
        for wrong in controls:
            try:
                identity(ARCHIVE, **wrong)
            except ValueError:
                continue
            raise AssertionError(f"identity negative control accepted: {wrong}")
        with tempfile.TemporaryDirectory() as tmp:
            damaged = Path(tmp) / "damaged.tar.gz"
            raw = ARCHIVE.read_bytes()
            damaged.write_bytes(raw[:-1] + bytes([raw[-1] ^ 1]))
            try:
                identity(damaged)
            except ValueError:
                pass
            else:
                raise AssertionError("same-size byte mutation accepted")
        print("IDENTITY_CONTROLS_PASS: positive and five negative controls")
        return
    if args.phase == "download":
        ARCHIVE.parent.mkdir(exist_ok=True)
        with urllib.request.urlopen(URL, timeout=120) as response, ARCHIVE.open("xb") as output:
            shutil.copyfileobj(response, output)
        desc = record("download")
        # Only this dependency-discovery copy is used by setup-r-dependencies.
        Path("frozen-deps").mkdir(exist_ok=True)
        Path("frozen-deps/DESCRIPTION").write_bytes(desc)
        ARCHIVE.chmod(0o444)
    else:
        record(args.phase)


if __name__ == "__main__":
    main()
