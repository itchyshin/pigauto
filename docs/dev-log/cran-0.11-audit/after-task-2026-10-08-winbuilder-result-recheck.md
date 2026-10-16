## 1. Goal

Record and classify the older Win-builder failure against the current exact-artifact release gate.

## 2. Implemented

Added a dated follow-up to the Windows upload receipt. It records the directly inspected R-release error, the paired R-devel error as reported in PR #228, the missing source/checksum binding, and the need for fresh checks on the final tarball.

## 3a. Decisions and Rejected Alternatives

Classified the old run as unbound evidence. It cannot pass or fail the later 10573 artifact. No new upload was made because the artifact is not final while documentation and provenance changes remain open.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/provenance/win-builder-exact-tarball-uploads-2026-10-06.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-winbuilder-result-recheck.md`

## 5. Checks Run

- Opened the R-release log in Chrome. It identifies pigauto 0.11.0, Windows Server 2022, R 4.6.1, and `Status: 1 ERROR`; testthat reports two failures at `test-check-pigauto.R:47-48`.
- Read PR #228's result note for the paired R-devel result and timestamp. Direct R-devel log navigation was blocked by the browser.
- Compared the failed preflight route with `origin/main`; current `check_pigauto()` defaults to `gnn = FALSE` and skips the torch probe unless explicitly requested, as merged in PR #229.
- `git diff --check` on the receipt passed.
- `slop_check.py` on the receipt reported 0 findings across 505 words.
- Lane preflight showed active source and citation lanes. No source or leased page was edited.

## 6. Tests of the Tests

No package code changed, so no R test was run. The negative control was provenance: no checksum or source commit is present in the old result notice, so the receipt explicitly refuses to attribute it to a later tarball.

## 7a. Issue Ledger

- Fixed: the receipt now distinguishes the old Windows error from results that could validate a later artifact.
- Open: exact R-release and R-devel results for the final post-merge tarball.
- Open: documentation/provenance closure and the BirdTree redistribution basis.

## 8. Consistency Audit

Compared the result timestamp with PR #228's upload notes, the 10573 candidate identity, the merged PR #229 default fix, and the current release gate. The evidence remains consistent: the old results predate later uploads and have no cryptographic artifact link. The PR #228 branch is still Draft, undeployed, and reports conflicts with `main`.

## 9. What Did Not Go Smoothly

The first lane-lease attempt lacked permission to write the protected lease registry and was not relied on; the exact-path lease succeeded after escalation. Chrome could not open the R-devel log directly, so its error status is attributed only to the PR conversation. The closeout template initially listed files from the brain repository; I replaced that template content with this pigauto-scoped report; the direct structural validator passed.

## 10. Known Residuals

This review does not identify the source archive used by either old Windows result. It does not establish that the 10573 tarball passes Windows, that all current docs are reconciled, or that redistribution rights are cleared. No new artifact or platform check was produced.

## 11. Team Learning

Record uploaded archive identity and result identity together. A filename, package version, timestamp, and matching byte size do not bind a builder result to a tarball checksum.

Memory receipt: loaded the pigauto route manifest and ran lane preflight; shared-lane and exact-artifact rules shaped this review. Brain files were not changed.

Golden Set: not checked; this was a release-provenance review, not a recurring product defect.

## 12. Cross-Product Coverage

Covers: the older R-release log and the R-devel status note in PR #228, compared with the current `check_pigauto()` behavior.

Does NOT cover: a Windows check of the final tarball, local rendered-site visual review, deployed-site verification after pending documentation changes, or CRAN submission readiness.
