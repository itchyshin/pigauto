# After-task: latest CRAN 0.11 candidate-source CI receipt

## 1. Goal

Record the completed three-platform candidate-source matrix on PR #231 head `d5a5b42` and align both open audit PR descriptions with the verified result.

## 2. Implemented

Added the exact-head result for GitHub Actions run #37894677486 to the CRAN 0.11 ledger. Updated PR #231 and PR #228 descriptions to state that Ubuntu R release, Ubuntu R-devel, and macOS arm64 R release passed on `d5a5b42` in 18m51s. Both descriptions retain the distinction between candidate-source CI and the later force-Suggests exact-artifact gate.

## 3a. Decisions and Rejected Alternatives

Kept the gate tally at 8 of 11. The source matrix does not close the deployed-site, final frozen-artifact, or independent final-review gates. The passing run used `NOT_CRAN=true` and `_R_CHECK_FORCE_SUGGESTS_=false`; it cannot establish the final artifact's required checks. No source behavior, user documentation, website input, scientific default, or campaign result changed in this receipt slice.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-latest-candidate-ci.md`
- External PR descriptions: #231 and #228

## 5. Checks Run

- Chrome Actions run #37894677486: all three matrix jobs succeeded in 18m51s: Ubuntu R release, Ubuntu R-devel, and macOS arm64 R release.
- Chrome PR #231: the updated description names run #37894677486 and head `d5a5b42`; the PR remains Draft and unmerged.
- Chrome PR #228: the updated description names the same run and preserves the 8-of-11 tally; the PR remains Draft and unmerged.
- Candidate workflow settings were `NOT_CRAN=true` and `_R_CHECK_FORCE_SUGGESTS_=false`. The pkgdown pull-request workflow is skipped by repository design.
- The after-task structure checker passed, but its repository-wide acceptance scan exited 1 on five unchecked `.unlazy/imputation-sim/gates/` leaves. These belong to a separate recorded campaign and were not changed here.
- No package or website tests were rerun because this slice changes only release evidence and the two PR descriptions.

## 6. Tests of the Tests

The run summary identifies the tested PR head and reports all three platform jobs as successful. This validates candidate-source CI for that head only. It does not validate a frozen tarball, force-Suggests behavior, deployment, or CRAN acceptance.

Golden set: no implementation behavior changed in this evidence update.

## 7a. Issue Ledger

Resolved: the two PR descriptions no longer call the `d5a5b42` matrix pending, and the local ledger now records its result and source head.

Open: G7 deployed-site verification; G8 one exact post-merge tarball with its required checks; G9 independent review of final artifact and deployed-site evidence.

## 8. Consistency Audit

The ledger and both PR descriptions agree on run #37894677486, source head `d5a5b42`, the three successful platforms, candidate-only scope, and the 8-of-11 tally. Both PRs remain Draft and unmerged. The source branch has not been deployed, and no CRAN submission occurred.

## 9. What Did Not Go Smoothly

The first preflight command was pointed at a script path that does not exist in this repository. Rerunning the documented brain-tool path produced the expected lane report and LOAD-FIRST manifest. No source or unrelated checkout files were changed by the failed lookup.

## 10. Known Residuals

G7, G8, and G9 remain open. Candidate-source CI is not a substitute for the final tarball's exact checks. The repo-wide after-task closeout command also remains red on the five separately tracked imputation-sim gates. Merge, deployment, and CRAN submission remain maintainer-controlled.

## 11. Team Learning

Update the ledger and both open PR descriptions from the same verified run receipt. Keep the tested source head, workflow settings, and artifact limitation together so a green source matrix cannot be mistaken for release evidence.

## 12. Cross-Product Coverage

This slice covers one three-platform candidate-source matrix, its ledger receipt, and synchronized descriptions on PRs #231 and #228. It does not cover deployment, the post-merge tarball, its exact-artifact checks, independent final review, or submission.
