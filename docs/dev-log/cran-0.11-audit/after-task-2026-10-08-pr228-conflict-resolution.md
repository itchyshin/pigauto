## 1. Goal

Resolve the three GitHub conflicts on the draft CRAN 0.11 evidence PR without changing the implementation lane or implying release readiness.

## 2. Implemented

Merged current `origin/main` into the evidence branch in an isolated temporary worktree. Resolved the three overlapping evidence files by retaining the post-merge deployment record, preserving the 2026-10-07 BirdTree source recheck, correcting the exact-candidate status, and explicitly marking stale source/site observations as historical. No source implementation was edited in this slice.

## 3a. Decisions and Rejected Alternatives

Kept the BirdTree redistribution gate open because citation and download access do not establish redistribution permission. Kept the 1bfed5ad tarball as a locally checked candidate from merged main, not as the final artifact. Kept historical observations summarized in the current gate ledger rather than repeating superseded claims as current status.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/provenance/pre-pr-check.md`
- `docs/dev-log/cran-0.11-audit/provenance/rights-and-policy.md`
- `docs/dev-log/cran-0.11-audit/surface-inventory.md`
- This report.

## 5. Checks Run

- `git diff --cached --check` passed for the four owned files after conflict resolution.
- The required naturalness scan found zero hits in this report.
- Confirmed no conflict markers remain in the three resolved files.
- Confirmed the current release ledger identifies the 1bfed5ad candidate and still marks the release `NOT_READY`.
- Chrome review of PR #228 showed Draft status, 14 commits, zero checks, three conflicts, no reviews, and no deployment before this resolution was pushed.

## 6. Tests of the Tests

Not applicable. This documentation-only merge-resolution slice did not alter executable code or checks.

## 7a. Issue Ledger

- Resolved: three GitHub merge conflicts in the release-evidence branch.
- Open: PR description does not mention the separate 1bfed5ad candidate receipt.
- Open: exact-hash Windows results, final post-documentation artifact, rights clearance, deployed-site verification, and exact-artifact independent review.

## 8. Consistency Audit

Compared the resolved claims with `GATES.md`, the candidate receipt, the BirdTree rights evidence, and the post-merge deployment record. No resolved text claims that citation establishes redistribution rights or that the old local candidate validates subsequent source or documentation edits.

## 9. What Did Not Go Smoothly

GitHub reported three conflicts because both the evidence branch and current main had appended dated audit snapshots to the same files. The first attempt to acquire the file lease encountered a busy registry lock; a subsequent authorized lease claim succeeded.

## 10. Known Residuals

The PR description and evidence ledger still need a consistent current-candidate summary. The evidence branch remains a draft, and all open release gates remain open. No source merge, website deployment, or CRAN submission occurred.

## 11. Team Learning

For an evidence branch that trails main, reconcile dated snapshots against the current gate ledger before resolving prose conflicts. Preserve provenance while clearly separating historical observations from the latest verified state.

## 12. Cross-Product Coverage

Checked PR #228's live GitHub state in Chrome and compared the UI's conflicts and branch state with the local merge worktree. The live state will be checked again after the resolution push.
