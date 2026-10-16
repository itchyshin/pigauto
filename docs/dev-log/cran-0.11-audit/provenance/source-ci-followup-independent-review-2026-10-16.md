# Follow-up independent review: GitHub R-devel source-CI receipt

Review date: 2026-10-16

This follow-up reviews the new GitHub R-devel source-CI receipt and its placement in the current gate records. It does not replace the earlier bounded review of the frozen archive, local `--as-cran` result, MI workflow matrix, and deployed website.

## Evidence reviewed

- Frozen archive SHA-256: `4bdf4de71650716b158824f545ad92c4d125fe1601adefab97eeb5f301c5debe`; 5,135,816 bytes.
- Source commit: `8082ef8273177ff664d49ac4bd643f7def671f20`.
- GitHub Actions run: [38843648553](https://github.com/itchyshin/pigauto/actions/runs/38843648553), Ubuntu R-devel job `116789146589`.
- Receipt: `github-r-devel-source-ci-2026-10-16.md`, SHA-256 `30c025b00d7c411a3fc962014b1d9f90f27938ee3e518a6ae2cac3f93e9d0245`.
- Current release ledger and both release and Unlazy gate summaries.

## Independent findings

- **Pat:** The receipt accurately records a source-CI pass on the exact source commit and clearly states that it did not check the frozen tarball. The ledger, source receipt, and Unlazy G3 preserve that distinction. No conflicting current status found.
- **Rose:** The artifact hash and size remain unchanged; the ledger copies agree; G3 remains unmet because the Win-builder R-release log does not echo the archive hash and Windows R-devel is pending. The GitHub source-CI pass cannot substitute for the checksum-bound artifact/platform requirement. G5 and the final evidence PR remain open.
- **Artifact-status reviewer:** Confirmed the archive identity, duplicated ledger fields, source-CI scope, and `NOT_READY` / `tarball-clean` / `platform-clean` status. The final PR remains premature while platform evidence and the strict checksum-binding criterion are unresolved.

## Disposition

The new source-CI receipt is accurate and consistently classified as supplementary source evidence. It does not satisfy G3 or advance the proven rung. Win-builder R-release is associated with the upload by a unique result URL, filename/size match, and local pre/post-upload hashing, but its log does not assert the tarball digest. Win-builder R-devel remains pending. The strict checksum-bound platform gate therefore remains unmet, and no final evidence PR is opened.

These reviews assess record accuracy and gate consistency. They do not claim a new rerun of the exact archive, a Windows result, platform-clean status, or CRAN readiness.
