# Exact-head CI refresh on `24d3285`: after-task report

## 1. Goal

Record the completed candidate-source R-CMD-check matrix for PR #231 head `24d3285`.

## 2. Implemented

Added the completed run #37890462685 result to the CRAN audit ledger. The run tested the exact PR source head `24d32858c211c1eed543a9f8ec53f3b3d34567f7` on Ubuntu R release, Ubuntu R-devel, and macOS arm64 R release.

## 3a. Decisions and Rejected Alternatives

The three passing source jobs support candidate-source compatibility only. They do not establish the final tarball's force-Suggests checks, deployed-site behavior, or independent final review. The ledger remains at 8 of 11 gates met, with G7–G9 open.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-ci-refresh-24d3285.md`

## 5. Checks Run

- Chrome Actions summary for run #37890462685: success, 3 of 3 matrix jobs completed, total 17m49s.
- Chrome job detail pages: Ubuntu R release, Ubuntu R-devel, and macOS arm64 R release each reported success.
- GitHub displayed two Ubuntu runner-image migration notices and one macOS arm64 capacity notice.
- The pkgdown pull-request workflow was skipped by repository design.
- PR #231 remains Draft and unmerged; its page showed no deployment. PR #228 remains the unmerged release-evidence PR.
- No implementation or website tests were rerun because this slice only records completed source-CI evidence.

## 6. Tests of the Tests

The run summary and all three job pages agree on success. This verifies the matrix result for the exact PR head. It does not test a post-merge tarball or deployed site.

Golden Set: no package behavior changed in this evidence-recording slice.

## 7a. Issue Ledger

Resolved: PR #231 head `24d3285` now has a complete three-platform source-CI result recorded in the audit ledger.

Open: G7 deployed-site verification; G8 one exact post-merge tarball with its local and platform checks; G9 independent review of the final artifact and site evidence.

## 8. Consistency Audit

The new receipt binds the result to the full commit hash and keeps it separate from G8's exact-artifact evidence. The gate count remains 8 of 11.

## 9. What Did Not Go Smoothly

The first run observations showed the matrix still active. Chrome's run page later confirmed all three jobs completed successfully.

## 10. Known Residuals

This source matrix does not satisfy deployed-site verification, exact post-merge artifact checks, Windows results bound to that artifact, or independent final review. Merge, deployment, and submission remain under Shinichi's control.

## 11. Team Learning

Use the Actions summary together with each matrix job page before recording a complete run; a successful single job does not establish a successful matrix.

## 12. Cross-Product Coverage

This slice covers the exact-head three-platform candidate-source matrix and its local audit record. It does not cover merge, deployment, frozen-tarball checks, or CRAN submission.
