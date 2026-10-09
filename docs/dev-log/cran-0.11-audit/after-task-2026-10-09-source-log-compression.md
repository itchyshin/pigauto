## 1. Goal

Reduce PR #231's review diff from six large raw GitHub Actions logs while preserving the exact downloaded log bytes, run summaries, and provenance hashes.

## 2. Implemented

Converted the six candidate-source CI logs to deterministic gzip files. Updated `GATES.md` with each compressed file's size and SHA-256 plus the decompressed byte count and original SHA-256. The six logs round-trip byte-for-byte. Their combined size fell from 32,931,259 bytes to 2,603,814 bytes, a 92.1% reduction. Reconciled the `OWNS` manifest after finding six removed raw-log paths and two site-evidence paths that were absent from the current source history and accessible worktrees.

## 3a. Decisions and Rejected Alternatives

Followed the independent source review's finding that the raw logs dominated the PR diff and obscured the source and documentation changes. Compressed and retained all six logs instead of discarding older CI evidence. Kept each source run's URL, test summary, environment limits, uncompressed checksum, and exact raw contents recoverable. No package source, defaults, data, site content, or release-gate status changed.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/provenance/source-ci-37927503719.log` (removed; replaced by `.log.gz`)
- `docs/dev-log/cran-0.11-audit/provenance/source-ci-37927503719.log.gz` (added)
- `docs/dev-log/cran-0.11-audit/provenance/source-ci-37930221718.log` (removed; replaced by `.log.gz`)
- `docs/dev-log/cran-0.11-audit/provenance/source-ci-37930221718.log.gz` (added)
- `docs/dev-log/cran-0.11-audit/provenance/source-ci-37932930957.log` (removed; replaced by `.log.gz`)
- `docs/dev-log/cran-0.11-audit/provenance/source-ci-37932930957.log.gz` (added)
- `docs/dev-log/cran-0.11-audit/provenance/source-ci-37934989398.log` (removed; replaced by `.log.gz`)
- `docs/dev-log/cran-0.11-audit/provenance/source-ci-37934989398.log.gz` (added)
- `docs/dev-log/cran-0.11-audit/provenance/source-ci-37937268183.log` (removed; replaced by `.log.gz`)
- `docs/dev-log/cran-0.11-audit/provenance/source-ci-37937268183.log.gz` (added)
- `docs/dev-log/cran-0.11-audit/provenance/source-ci-37940225845.log` (removed; replaced by `.log.gz`)
- `docs/dev-log/cran-0.11-audit/provenance/source-ci-37940225845.log.gz` (added)
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-source-log-compression.md`

## 5. Checks Run

- Confirmed the six original files matched their recorded raw sizes and SHA-256 values before compression.
- Compressed with fixed gzip metadata (`mtime=0`) and verified each decompressed byte sequence equals its original; recorded the compressed size and checksum in `GATES.md`.
- Combined compressed size: 2,603,814 bytes; original size: 32,931,259 bytes.
- The diff remains documentation and evidence only; no package tests or model runs were required for this representation change.
- Checked that every remaining `OWNS` path exists, that all six compressed CI logs are listed, and that no removed raw-log or unavailable site-evidence path remains listed.
- The after-task structure check passed. The full repository acceptance check exits 1 because five existing leaves in `.unlazy/imputation-sim/gates/` remain unmet; this documentation slice does not own or modify those gates.

## 6. Tests of the Tests

The round-trip check compares every decompressed archive to the original bytes before deleting the uncompressed copy. The raw checksums in the ledger provide an independent identity check after the source files are removed. A separate final verification checks gzip integrity, compressed hashes, decompressed hashes, sizes, and every ledger path.

## 7a. Issue Ledger

Resolved: raw logs previously added 267,919 visible lines and 32.9 MB to PR #231, obscuring the substantive source and reader-surface changes. Preserved: all six exact CI logs and their checksums. Resolved six stale raw-log paths and two unavailable site-evidence paths in `OWNS`. The dated site-review file named by the older after-task report is absent from current source history and accessible worktrees; the current exact-head visual findings remain recorded directly in `GATES.md`. Open: deployment and exact final-artifact gates G7–G9 remain unchanged.

## 8. Consistency Audit

All six GATES entries continue to identify the same Actions runs and summarize the same candidate-source results. Only the storage format and stored-file checksums changed; the uncompressed checksums remain the historical raw-log identity. Every `OWNS` path now resolves in the current source tree, including the six `.log.gz` files and matching after-task report. The current exact-head site review is recorded in `GATES.md`; no gate relies solely on the absent older site-review file. The PR description's statement that raw logs are retained remains accurate because the byte-identical logs are stored as `.log.gz` files.

## 9. What Did Not Go Smoothly

The source PR worktree is outside the default writable roots, so the scoped filesystem approval was used. A ledger recheck caught stale ownership paths after compression and two site paths without current source-history copies; these are now reconciled. The compression and checksum checks completed without errors.

## 10. Known Residuals

This evidence-format change does NOT prove candidate-source CI against a frozen artifact, deployed-site behavior, a final tarball, platform results bound to that final tarball, or CRAN readiness. G7–G9 remain open. The full repository acceptance check also remains unmet on five existing `imputation-sim` leaves outside this slice. The dated 2026-10-07 site-review file listed in an older report is not available in the current source history or checked worktrees; the corresponding current-head visual review is documented in `GATES.md`.

## 11. Team Learning

When large logs are valuable as durable evidence but dominate a code-review diff, deterministic compression can preserve exact bytes and hashes while making the substantive changes easier to inspect.

Memory receipt: `route.py pigauto` loaded the LOAD-FIRST manifest. The separation between source-CI receipts and exact-artifact checks shaped the scope.

Golden Set: no package source or runtime behavior changed; no source known-mistake class was in scope.

## 12. Cross-Product Coverage

Covers six source-CI log artifacts and the ledger references in PR #231.

Does NOT cover package behavior, defaults, MI adapter correctness, rendered or deployed pages, the final frozen tarball, platform checks on that tarball, or CRAN submission.
