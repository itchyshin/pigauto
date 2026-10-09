# Candidate CI matrix refresh: after-task report

## 1. Goal

Refresh the pigauto CRAN 0.11 audit with the completed source-branch CI matrix, while stating that this is the only active pigauto lane and keeping exact-artifact gates separate.

## 2. Implemented

Added the current one-lane preflight result and the successful Ubuntu release, Ubuntu R-devel, and macOS release CI jobs for candidate commit `7cbbcbc` to `GATES.md`. Corrected the stale PR snapshot by identifying the earlier `621d44b` / 12-commit report as superseded.

## 3a. Decisions and Rejected Alternatives

- Recorded candidate-source CI as source evidence only; it does not validate a frozen tarball.
- Kept G6 and G7 open for local visual review and deployed-site verification, and G8/G9 open for exact-artifact checks and independent review.
- Did not treat the skipped PR pkgdown workflow as a local site-build failure; fresh local build evidence remains recorded under G5/G5b.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-ci-matrix.md`

## 5. Checks Run

- Lane preflight: reported one pigauto lane, the current Codex lane.
- Chrome PR #231 and workflow run #37848354578: Ubuntu R release passed in 9m16s; Ubuntu R-devel passed in 11m49s; macOS R release passed in 19m57s, including MPS prediction (5m14s) and `R CMD check` (11m06s). The pkgdown workflow was skipped on PRs by repository design; one Ubuntu migration notice was present.
- `git diff --check`: passed.

## 6. Tests of the Tests

No package test code changed. The recorded workflow jobs ran package checks on candidate source. Their source commit and run are linked in the ledger; no claim is made that they test an exact release tarball.

## 7a. Issue Ledger

- Fixed: the latest CI matrix and commit count were absent from the current external-state evidence.
- Open: local visual review, deployed-site verification, exact post-merge tarball, and independent exact-artifact review.

## 8. Consistency Audit

Checked the current lane preflight, PR #231, its completed workflow, the existing local G5/G5b site evidence, and the stale external-state paragraph. The evidence now distinguishes candidate-source CI from frozen-artifact and deployed-site checks.

## 9. What Did Not Go Smoothly

The first attempt to invoke a guessed local Unlazy script path failed because that path does not exist. No files were changed by that attempt. `git diff --check` passed after the ledger update.

## 10. Known Residuals

G0, G6, G7, G8, and G9 remain unmet as described in the release ledger. No merge, deployment, or CRAN submission occurred. Candidate-source CI does not establish platform results for the final tarball.

## 11. Team Learning

The sole-lane preflight and the run attached to the candidate source commit shaped this update. The repo's LOAD-FIRST rules require recovery-to-truth and separate source checks from artifact claims. No Golden Set regression class was in scope for this evidence-only update.

Golden Set: not in scope for a CI-ledger refresh.

## 12. Cross-Product Coverage

Covers one-lane ownership status, PR #231 candidate-source CI, and reconciliation with local site-build evidence.

Does NOT cover visual page layout, the deployed site after merge, rights for bundled tree data, a final frozen tarball, platform checks on that artifact, independent exact-artifact review, or CRAN acceptance.
