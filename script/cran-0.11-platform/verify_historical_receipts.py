"""Check retained records of the 0.11.0 checks, without the original tarball.

This verifies recorded identities and log hashes. It does not reconstruct the
archive, rerun R CMD check, or verify the current release candidate.
"""

import hashlib
import json
from pathlib import Path

ARCHIVE_SHA = "a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb"
INVENTORY_SHA = "d7b9b0e59ea051b3a22d7f347d8d925e6265b076b238ac893410e79d5b61329d"
MATRIX = {("windows-latest", "release"), ("windows-latest", "devel"),
          ("ubuntu-latest", "release")}
SUFFIX = "\nDIRECT_ARCHIVE_MATRIX_VERIFIED_3_OF_3"


def verify(root):
    root = Path(root)
    receipt = json.loads((root / "FINAL.json").read_text())
    if (receipt["sha256"], receipt["size_bytes"], receipt["inventory_sha256"]) != (
        ARCHIVE_SHA, 5137837, INVENTORY_SHA
    ):
        raise ValueError("historical archive receipt does not match the frozen identity")
    text = (root / "native/direct-verification.log").read_text().rstrip()
    if not text.endswith(SUFFIX):
        raise ValueError("historical matrix verification receipt is incomplete")
    records = json.loads(text[:-len(SUFFIX)])
    if len(records) != 3 or {(r["os"], r["r"]) for r in records} != MATRIX:
        raise ValueError("historical platform matrix is incomplete")
    for record in records:
        identity = record["identity"]
        if (identity["sha256"], identity["size_bytes"], identity["inventory_sha256"]) != (
            ARCHIVE_SHA, 5137837, INVENTORY_SHA
        ) or record["check_status"] != "OK":
            raise ValueError("historical platform identity or status does not match")
        path = root / "native" / f"{record['os']}-{record['r']}" / "Rcheck/00check.log"
        content = path.read_bytes()
        if hashlib.sha256(content).hexdigest() != record["check_log_sha256"]:
            raise ValueError(f"historical check log hash mismatch: {path}")
        if "Status: OK" not in content.decode(errors="replace").splitlines():
            raise ValueError(f"historical check log did not report OK: {path}")
    return records


if __name__ == "__main__":
    root = Path("docs/dev-log/cran-0.11-audit/provenance/final-a93fece7-2026-10-18")
    verify(root)
    print("Historical receipts and three platform log hashes verified; archive not rechecked.")
