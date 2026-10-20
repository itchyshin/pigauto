#!/usr/bin/env python3
"""Read-only verifier for the deterministic local macOS evidence archive."""
import argparse
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import tarfile

MANIFEST = "publication/member-manifest.json"

def safe_name(name):
    if not isinstance(name, str) or not name or "\\" in name:
        raise ValueError(f"unsafe member path: {name!r}")
    p = PurePosixPath(name)
    if p.is_absolute() or any(part in ("", ".", "..") for part in p.parts) or ":" in p.parts[0]:
        raise ValueError(f"unsafe member path: {name!r}")
    return name

def verify_bytes(blob):
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as tf:
        files = {}
        for m in tf.getmembers():
            name = safe_name(m.name)
            if not m.isfile() or name in files:
                raise ValueError(f"member is not a unique regular file: {name}")
            f = tf.extractfile(m)
            if f is None:
                raise ValueError(f"cannot read member: {name}")
            files[name] = f.read()
    if MANIFEST not in files:
        raise ValueError("member manifest missing")
    manifest = json.loads(files[MANIFEST])
    expected = {row["path"]: row for row in manifest["members"]}
    if len(expected) != len(manifest["members"]):
        raise ValueError("duplicate paths in member manifest")
    if set(files) != set(expected) | {MANIFEST}:
        raise ValueError("archive path set differs from manifest")
    for name, row in expected.items():
        safe_name(name)
        data = files[name]
        if len(data) != row["size_bytes"] or hashlib.sha256(data).hexdigest() != row["sha256"]:
            raise ValueError(f"member hash/size mismatch: {name}")
    return {"status": "PASS", "member_count": len(files), "manifest_member_count": len(expected),
            "manifest_sha256": hashlib.sha256(files[MANIFEST]).hexdigest(),
            "members": [{"path": n, "size_bytes": len(files[n]), "sha256": hashlib.sha256(files[n]).hexdigest()}
                        for n in sorted(files)]}

def selftest(blob):
    good = verify_bytes(blob)
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as tf:
        data = {m.name: tf.extractfile(m).read() for m in tf.getmembers() if m.isfile()}
    manifest = json.loads(data[MANIFEST])
    bad_hash = json.loads(json.dumps(manifest))
    bad_hash["members"][0]["sha256"] = "0" * 64
    try:
        verify_member_map(data, bad_hash)
    except ValueError as e:
        hash_failure = str(e)
    else:
        raise AssertionError("planted member hash mismatch accepted")
    try:
        safe_name("../escape.txt")
    except ValueError as e:
        path_failure = str(e)
    else:
        raise AssertionError("planted path traversal accepted")
    return {"positive_roundtrip": good, "planted_hash_failure": {"rejected": True, "reason": hash_failure},
            "planted_path_failure": {"rejected": True, "reason": path_failure}}

def verify_member_map(files, manifest):
    expected = {row["path"]: row for row in manifest["members"]}
    if len(expected) != len(manifest["members"]):
        raise ValueError("duplicate paths in member manifest")
    if set(files) != set(expected) | {MANIFEST}:
        raise ValueError("archive path set differs from manifest")
    for name, row in expected.items():
        safe_name(name)
        data = files[name]
        if len(data) != row["size_bytes"] or hashlib.sha256(data).hexdigest() != row["sha256"]:
            raise ValueError(f"member hash/size mismatch: {name}")
    return True

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("archive", type=Path)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    blob = args.archive.read_bytes()
    result = selftest(blob) if args.selftest else verify_bytes(blob)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "members"}, indent=2))

if __name__ == "__main__":
    main()
