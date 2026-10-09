# Final candidate-source CI refresh: after-task report

## 1. Goal

Verify the newest PR #231 source-check run and record whether the source candidate has reached the maintainer merge gate.

## 2. Implemented

Recorded completed Ubuntu release, Ubuntu R-devel, and macOS arm64 source checks for PR head `d7ef63a64dd424a768f8ce0dfdb6ff1b2d106af8`, tested through synthetic merge commit `9cb436d` into base `0b0f71f`. Updated PR #231's validation text to report the completed run and retain the force-Suggests and exact-artifact limits.

## 3a. Decisions and Rejected Alternatives

Treat the matrix as candidate-source CI only: it ran with `NOT_CRAN=true` and `_R_CHECK_FORCE_SUGGESTS_=false`. Do not describe it as G8 evidence or as a check of a frozen tarball. Stop at the merge gate and leave the merge to Shinichi.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-final-candidate-ci.md`

## 5. Checks Run

- Run [#37880698572](https://github.com/itchyshin/pigauto/actions/runs/37880698572) completed successfully on PR head `d7ef63a64dd424a768f8ce0dfdb6ff1b2d106af8`; checkout log confirmed the synthetic merge commit.
- Ubuntu R release, Ubuntu R-devel, and macOS arm64 R release each reported `R CMD check Status: OK` and 3,196 passes, 0 failures, 175 warnings, and 83 skips.
- The focused macOS MPS prediction tests reported 201 passes, 0 failures, 50 warnings, and 0 skips.
- The workflow set `NOT_CRAN=true` and `_R_CHECK_FORCE_SUGGESTS_=false`; pkgdown was skipped on the pull request by repository design.
- PR #231 still points to `d7ef63a64dd424a768f8ce0dfdb6ff1b2d106af8`, is Draft, and is unmerged. Local worktree was clean before this evidence update.

## 6. Tests of the Tests

No test code changed. The full three-platform package matrix and focused macOS MPS test exercised the existing suite; no negative control was needed for this evidence-only slice.

## 7a. Issue Ledger

Resolved: the latest source-check matrix passed, supporting the merge-review checkpoint for PR #231.

Open: G7 deployed-site verification follows the maintainer-controlled merge; G8 requires the exact post-merge tarball and checks; G9 requires independent review of the final artifact and site evidence.

## 8. Consistency Audit

The source check binds to PR head `d7ef63a64dd424a768f8ce0dfdb6ff1b2d106af8` through synthetic merge commit `9cb436d` with base `0b0f71f`. Updated the PR description from pending to completed status. Kept the 8-of-11 ledger tally and G7–G9 residuals consistent. The passing checks do not close force-Suggests or frozen-artifact requirements.

## 9. What Did Not Go Smoothly

The macOS job took about 22 minutes overall, including 11 minutes for the package check. GitHub did not expose logs until the job completed; once available, the check summary was verified directly.

## 10. Known Residuals

No merge, deployment, exact final tarball, Windows result, force-Suggests artifact check, independent final verdict, or CRAN submission has occurred. The successful run applies only to candidate source through the stated synthetic merge.

## 11. Team Learning

This completes the candidate-source CI checkpoint. The work is now at the maintainer-controlled source merge gate; later artifact and publication gates remain separate.

## 12. Cross-Product Coverage

Covers candidate-source Linux release, Linux R-devel, and macOS arm64 release checks, including the focused MPS test. Does NOT cover the exact frozen tarball, force-Suggests checks, Windows, deployed site, final independent artifact review, or CRAN acceptance.
