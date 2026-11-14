"""Offline regression tests for portable CRAN-ledger path resolution."""

import contextlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


TOOL = Path(__file__).with_name("cran_release_gate.py")
SPEC = importlib.util.spec_from_file_location("cran_gate", TOOL)
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)


@contextlib.contextmanager
def cwd(path):
    before = Path.cwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(before)


class PortablePaths(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / "repo"
        self.root.mkdir()
        self.other = Path(self.temp.name) / "elsewhere"
        self.other.mkdir()
        self.artifact = self.root / "candidate.tar.gz"
        self.artifact.write_bytes(b"frozen candidate")
        self.evidence = self.root / "evidence.log"
        self.evidence.write_text("retained evidence\n")
        self.base = gate.valid_fixture(self.artifact, self.evidence)
        self.ledger = self.root / "docs" / "release-ledger.json"
        self.ledger.parent.mkdir()

    def relative(self):
        data = json.loads(json.dumps(self.base))
        data["artifact"]["path"] = self.artifact.name
        data["artifact"]["inventory_evidence"]["value"] = self.evidence.name
        for item in data["evidence"].values():
            if item["kind"] == "file":
                item["value"] = self.evidence.name
        return data

    def cli(self, data, *args, working=None):
        self.ledger.write_text(json.dumps(data))
        return subprocess.run(
            [sys.executable, str(TOOL), str(self.ledger), *map(str, args)],
            cwd=working or self.other, text=True, capture_output=True,
            timeout=10,
        )

    def test_absolute_api_is_unchanged(self):
        for location in (self.root, self.other):
            with cwd(location):
                self.assertEqual(gate.validate(self.base), [])

    def test_relative_artifact_requires_root_even_when_cwd_matches(self):
        self.base["artifact"]["path"] = self.artifact.name
        for location in (self.root, self.other):
            with cwd(location):
                self.assertTrue(gate.validate(self.base))

    def test_relative_evidence_requires_root(self):
        self.base["evidence"]["rights_and_consent"]["value"] = self.evidence.name
        with cwd(self.root):
            self.assertTrue(gate.validate(self.base))

    def test_explicit_root_api_and_cli_ignore_cwd(self):
        for location in (self.root, self.other):
            with cwd(location):
                self.assertEqual(gate.validate(self.relative(), self.root), [])
            result = self.cli(self.relative(), "--root", self.root, working=location)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("READY FOR CLAIMED RUNG", result.stdout)

    def test_cli_discovers_git_directory_git_file_or_description(self):
        for marker_type in ("git-directory", "git-file", "description"):
            marker = self.root / ("DESCRIPTION" if marker_type == "description" else ".git")
            if marker_type == "git-directory":
                marker.mkdir()
            else:
                marker.write_text("fixture marker\n")
            try:
                result = self.cli(self.relative())
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            finally:
                marker.rmdir() if marker.is_dir() else marker.unlink()

    def test_cli_without_marker_uses_ledger_parent(self):
        self.ledger = self.root / "release-ledger.json"
        self.assertEqual(self.cli(self.relative()).returncode, 0)

    def test_wrong_root_missing_paths_and_directories_fail(self):
        self.assertTrue(gate.validate(self.relative(), self.other))
        for target in ("artifact", "inventory", "evidence"):
            for bad in ("missing", "."):
                data = self.relative()
                if target == "artifact":
                    data["artifact"]["path"] = bad
                elif target == "inventory":
                    data["artifact"]["inventory_evidence"]["value"] = bad
                else:
                    data["evidence"]["rights_and_consent"]["value"] = bad
                self.assertTrue(gate.validate(data, self.root), (target, bad))

    def test_bad_hash_and_size_still_fail(self):
        for key, value in (("sha256", "0" * 64), ("size_bytes", 0)):
            data = self.relative()
            data["artifact"][key] = value
            self.assertTrue(gate.validate(data, self.root), key)


if __name__ == "__main__":
    unittest.main()
