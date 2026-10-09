## 1. Goal

Reconcile the current BirdTree rights summary with Shinichi's recorded maintainer confirmation and the CRAN 0.11 gate ledger, while preserving the dated evidence history.

## 2. Implemented

Updated the current component-status row and summary in `provenance/rights-and-policy.md` to record Shinichi's direct confirmation as the maintainer-warranty basis for including the bundled example trees with attribution. Reframed earlier dated conclusions as historical snapshots and pointed to the current 2026-10-09 disposition. The report continues to say that the reviewed BirdTree pages publish no separate redistribution licence and makes no independent legal determination. Exact tree-object and NOTICE inspection remains assigned to G8.

## 3a. Decisions and Rejected Alternatives

Followed the current ledger decision: treat the maintainer's statement as the basis for the CRAN inclusion warranty. Do not claim an independently verified BirdTree licence or rights-holder agreement. Do not repeat the stale open-rights conclusion as current. Preserve it where it documents an earlier dated review.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/provenance/rights-and-policy.md`
- `docs/dev-log/cran-0.11-audit/GATES.md` (OWNS entry only)
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-rights-status-reconcile.md`

## 5. Checks Run

- Re-read the current G0 and G8 ledger evidence and the full rights-summary timeline before editing.
- Confirmed in the GitHub browser that `itchyshin/pigauto` is public; the `gh repo view` API request failed to connect, so no CLI visibility result is claimed.
- The rights-summary style check reported 0 findings; `git diff --check` passed.
- The after-task structure check passed. The repository-wide closeout still reports five unmet `.unlazy/imputation-sim` gates outside this CRAN audit.
- An independent ledger reviewer confirmed the rights-summary wording and found that the closure sentence overstated merge status because PR #229 had merged. The sentence now names the still-unmerged audit PRs and does not imply that no unrelated merge occurred.
- No package tests, push, or new Actions run had occurred when this report was first committed. The package code, bundled objects, installed NOTICE, generated help, and website files were not changed in the slice.

## 6. Tests of the Tests

The slice corrects documentation consistency. The checks cover prose style, patch whitespace, and required report structure. No package behavior changed, so package tests are not applicable. Candidate-source CI at `83117b48b12c864985fc44387b58e47585bb0ce0` predates this audit-record edit and is not claimed as a result for a later source head.

## 7a. Issue Ledger

Resolved for this slice: the active rights summary now agrees with the current G0 maintainer-warranty decision while distinguishing that decision from public BirdTree licensing evidence.

Still open: G7 deployed-site verification; G8 exact post-merge tarball, including tree objects and NOTICE; G9 independent review of the final artifact and deployed site.

## 8. Consistency Audit

Compared the rights summary's current row, its dated historical assessments, the 2026-10-09 warranty record, `inst/NOTICE`, and the current G0/G8 wording. Current statements now agree. Dated earlier open-gate findings remain identified as historical. The independent review confirmed this reconciliation and prompted the corrected merge-status sentence. No gate status changed.

## 9. What Did Not Go Smoothly

The shell GitHub API request could not connect. Repository visibility was verified in the browser instead. The correction does not depend on API data.

## 10. Known Residuals

The public BirdTree pages reviewed do not state a separate redistribution licence. This report records the maintainer's confirmation as the warranty basis and does not make an independent legal determination. Exact shipped tree bytes and NOTICE remain to be inspected in the final frozen tarball. PR #229 had already merged to main; PRs #231 and #228 remain unmerged. This slice did not verify a deployment, freeze the final artifact, or submit to CRAN.

## 11. Team Learning

When a dated evidence note is superseded by a maintainer decision, preserve the dated observation and update the present-status summary so readers do not mistake historical uncertainty for the current gate decision.

## 12. Cross-Product Coverage

Covers the internal CRAN rights summary and its consistency with the current release ledger. It does NOT cover installed documentation, package behavior, website content, shipped data, exact-tarball inspection, deployed-site verification, or CRAN acceptance.
