"""Audit downloaded GitHub evidence for all three frozen-archive jobs."""

import hashlib
import json
from pathlib import Path
import sys

from artifact import SHA, SIZE, INVENTORY


def verify(root):
    records = []
    for os_name, r_version in [("windows-latest", "release"),
                                ("windows-latest", "devel"),
                                ("ubuntu-latest", "release")]:
        folder = root / f"frozen-14ee7e5d-{os_name}-{r_version}"
        evidence = folder / "platform-evidence"
        identities = []
        for phase in ["download", "pre", "post"]:
            item = json.loads((evidence / f"identity-{phase}.json").read_text())
            assert item["sha256"] == SHA and item["size_bytes"] == SIZE
            assert item["inventory_sha256"] == INVENTORY and item["archive_entries"] == 250
            assert item["package"] == "pigauto" and item["version"] == "0.11.0"
            assert item["phase"] == phase and item["run_id"] and item["workflow_sha"]
            identities.append(item)
        for key in ["archive", "sha256", "size_bytes", "workflow_sha", "run_id", "run_attempt"]:
            assert len({x[key] for x in identities}) == 1, key
        log = folder / "frozen-check/pigauto.Rcheck/00check.log"
        text = log.read_text(errors="replace")
        assert "Status: OK" in text.splitlines(), "Check status not OK"
        assert "exit=0" in (evidence / "check-exit.txt").read_text().splitlines()
        console = (evidence / "check-console.log").read_text(errors="replace")
        assert any("checking for unstated dependencies in" in line and "tests" in line
                   for line in text.splitlines())
        environment = (evidence / "check-environment.txt").read_text().splitlines()
        assert "NOT_CRAN=UNSET" in environment
        assert "_R_CHECK_FORCE_SUGGESTS_=true" in environment
        command = (evidence / "check-command.txt").read_text()
        assert identities[1]["archive"] in command
        for flag in ["--as-cran", "--run-donttest", "--timings"]:
            assert flag in command
        assert "--no-manual" not in command
        records.append({"os": os_name, "r": r_version, "identity": identities[1],
                        "check_log_sha256": hashlib.sha256(log.read_bytes()).hexdigest(),
                        "check_status": "OK", "manual_review_required": True})
    assert len({x["identity"]["run_id"] for x in records}) == 1
    assert len({x["identity"]["workflow_sha"] for x in records}) == 1
    print(json.dumps(records, indent=2))
    print("DIRECT_ARCHIVE_MATRIX_VERIFIED_3_OF_3")
    return records


if __name__ == "__main__":
    verify(Path(sys.argv[1]))
