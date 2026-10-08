## 1. Goal

Refresh the release evidence against the current GitHub state of pigauto PRs #228 and #231, using Chrome, while keeping the release ledger bounded to verified facts.

## 2. Implemented

Recorded the fresh PR states and current source-review findings in the release ledger. The release remains not ready. No source changes, merge, deployment, or submission occurred.

## 3a. Decisions and Rejected Alternatives

Kept G0, G4, G6, G7, G8, and G9 open. The latest Chrome view supersedes older PR commit-count snapshots. Kept the BirdTree rights and unsupported ranking claim open rather than treating citation instructions as permission to redistribute data.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-pr231-checkpoint.md`

## 5. Checks Run

- Ran lane preflight and the pigauto route manifest. Preflight found the source/docs lane on `codex/getting-started-birdtree-provenance`; a narrow lease was granted for the two evidence files only.
- Reloaded PR #228 in Chrome: 17 commits, Draft, no conflicts, required checks passed with one skipped check, no review, no deployment.
- Reloaded PR #231 in Chrome: 9 commits, Draft, no conflicts, 3 successful checks and 1 skipped, no review, no deployment.
- Read the current PR #231 description, which identifies its local source assertions, generated-help checks, two-tree Newick regression, and cleaned local pkgdown build. It explicitly leaves deployment verification open.
- Checked the existing independent review at `c6235fa`: it confirms user-directed tree acquisition and citation guidance, and flags the unsupported Robinson-Foulds ranking in `inst/NOTICE` plus the missing upstream byte comparison.
- No test suite or package build was run because this slice changed only release evidence.

## 6. Tests of the Tests

The primary evidence was read directly from the current Chrome-rendered GitHub PR pages. The PR commit counts, check summaries, draft status, conflicts, review state, and deployment state are visible in GitHub's current page. The review finding is recorded as a finding from the named independent review, not as a newly reproduced calculation.

## 7a. Issue Ledger

- Confirmed: PR #228 has 17 commits and PR #231 has 9 at the time of this Chrome check.
- Confirmed: both PRs remain Draft and undeployed; neither has a review on GitHub.
- Open: unsupported Robinson-Foulds ranking claim in `inst/NOTICE` and missing upstream byte comparison in the source/docs lane.
- Open: BirdTree redistribution permission, full site verification, final frozen artifact, platform checks, and exact-artifact independent review.

## 8. Consistency Audit

Compared the new Chrome state with the older commit-count statements already present in GATES.md. The ledger now labels the latest 17- and 9-commit counts as current and treats older snapshots as historical. The PR #231 description's local build and test claims are not described as deployed or as closing release gates.

## 9. What Did Not Go Smoothly

An initial lane-preflight invocation used the repository name instead of its path. The corrected invocation identified the active source/docs lane. The source lane remains active, so the evidence update was limited to explicitly leased evidence files.

## 10. Known Residuals

The source/docs changes are not reviewed on GitHub or deployed. The tree ranking claim has not been substantiated, and no upstream byte comparison is recorded. BirdTree's public citation instructions do not establish redistribution rights. The final source revision, deployed site, frozen artifact, platform checks, exact-hash review, merge, and CRAN submission remain open.

## 11. Team Learning

Refresh PR state directly in Chrome before recording counts. Keep source review findings distinct from independent reproduction, and keep citation permission distinct from redistribution permission.

Memory receipt: used the pigauto CRAN audit memory entry to locate the live release ledger and preserve the lane and artifact boundaries.

Golden Set: not applicable because this slice records PR state and reviewer findings.

## 12. Cross-Product Coverage

This update does NOT cover source-code changes, installed-package behavior, or the local and deployed site renders. It covers only the state shown on two GitHub PR pages and the findings in the completed independent review at `c6235fa`. No shared brain, sibling package, or external publication surface was changed.
