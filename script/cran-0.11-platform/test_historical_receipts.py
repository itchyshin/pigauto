"""Checks that historical evidence verification fails when records are altered."""
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from verify_historical_receipts import MATRIX, verify

SOURCE = Path("docs/dev-log/cran-0.11-audit/provenance/final-a93fece7-2026-10-18")


class HistoricalReceiptsTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        paths = [Path("FINAL.json"), Path("native/direct-verification.log")]
        paths.extend(Path("native") / f"{os_name}-{version}" / "Rcheck/00check.log"
                     for os_name, version in MATRIX)
        for path in paths:
            dest = self.root / path
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(SOURCE / path, dest)

    def test_valid_records(self):
        self.assertEqual(len(verify(self.root)), 3)

    def test_modified_log_is_rejected(self):
        path = self.root / "native/ubuntu-latest-release/Rcheck/00check.log"
        path.write_bytes(path.read_bytes() + b"altered")
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            verify(self.root)

    def test_wrong_archive_identity_is_rejected(self):
        path = self.root / "FINAL.json"
        receipt = json.loads(path.read_text())
        receipt["sha256"] = "0" * 64
        path.write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError, "frozen identity"):
            verify(self.root)


if __name__ == "__main__":
    unittest.main()
