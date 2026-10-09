## 1. Goal

Verify the published pigauto version before finalizing the CRAN 0.11 version choice.

## 2. Implemented

Recorded current CRAN and GitHub release/tag evidence. The candidate declares 0.11.0; current public CRAN records show 0.10.0, and GitHub's latest release and tag are v0.10.0.

## 3a. Decisions and Rejected Alternatives

Followed the approved rule and retained 0.11.0 because CRAN has not published it. This public check does not establish whether a submission is pending privately. No version field was changed.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/provenance/publication-history-2026-10-08.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-publication-history.md`

## 5. Checks Run

- Opened the official CRAN package page in Chrome and verified version 0.10.0 and publication date 2026-07-30.
- Opened GitHub's releases and tags pages in Chrome; both show v0.10.0 as latest.
- Verified candidate `DESCRIPTION` declares 0.11.0.
- `slop_check.py` reported zero findings on the provenance record.
- `git diff --check` passed.

## 6. Tests of the Tests

Used two independent public release records, the CRAN package page and GitHub's tags/releases, and compared them with the candidate version field. The checks distinguish the published CRAN version from the repository's next candidate.

## 7a. Issue Ledger

- Verified: 0.11.0 is not the currently published CRAN version and is absent from the current GitHub releases/tags listing.
- Open: whether CRAN has a pending submission, bundled-tree redistribution rights, and final source provenance.

## 8. Consistency Audit

Compared the candidate `DESCRIPTION`, the official current CRAN page, GitHub releases, and GitHub tags. They agree that the next candidate is 0.11.0 while the latest public release remains 0.10.0.

## 9. What Did Not Go Smoothly

The read-only `git ls-remote` tag check could not reach GitHub from the shell because SSH/DNS access was blocked. Chrome exposed the public releases and tags pages, so the version check was completed from those primary records.

## 10. Known Residuals

This verifies public version history only. It does not verify an unpublished or pending CRAN submission, rights clearance, the final source commit, a frozen tarball, or platform checks. No merge, deployment, or CRAN submission occurred.

## 11. Team Learning

When shell access to a public source is blocked, verify the same read-only fact in the requested browser and record which source was actually observed.

Assessment: 1/10, high confidence. Facts and references were checked against the live CRAN page, GitHub releases/tags, and candidate `DESCRIPTION`; no scientific claim is made.

## 12. Cross-Product Coverage

Covers the published-version record, public GitHub release/tag state, and candidate version field. Does NOT cover pending submissions, tree-data rights, release checks, final artifact identity, website deployment, or CRAN acceptance.
