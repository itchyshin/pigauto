## 1. Goal

Reconcile the pigauto CRAN 0.11 release ledgers with current PR, deployment, and candidate-source evidence while preserving active lane ownership.

## 2. Implemented

- Corrected malformed historical labels in the CRAN audit ledger so its Unlazy parser reads the documented gates.
- Updated both ledgers to record PR #230 merged as `0b0f71f`, Pages run #583, candidate source ancestry, and the distinction between candidate help and deployed pages.
- Rechecked six live documentation URLs with cache-busting query parameters in the Codex in-app browser after Pages run #583 and updated both ledgers with the current discrepancies: Hackett-only/MIT wording in all three tree help pages, missing retrieval instructions in the deployed tree article, stale MCC wording in getting-started, and the old inverse-Wishart boundary in the MI article. A first non-cache-busted `trees300` view differed and was superseded by the direct cache-busted response.
- Left source and help files held by implementation and citation lanes unchanged.

## 3a. Decisions and Rejected Alternatives

Kept the 0.11 recommendation user-controlled: users retrieve BirdTree or `megatrees` trees themselves, cite the source as required, and pass the tree object to pigauto. No downloader or data-rights claim was added. Did not merge, deploy, submit, or freeze a new artifact. Corrected only the owned audit records.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `.unlazy/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-ledger-reconciliation.md`

## 5. Checks Run

- `~/shinichi-brain/tools/lane_preflight.sh "$PWD"`: detected active implementation and citation lanes; no overlapping source edits were made.
- `node ~/.codex/skills/unlazy/scripts/gate-check.mjs --status docs/dev-log/cran-0.11-audit/GATES.md`: parsed 11 gates; 5 met, 6 unmet.
- `node ~/.codex/skills/unlazy/scripts/gate-check.mjs --status .unlazy/cran-0.11-audit/GATES.md`: parsed 11 gates; 10 met, I7 unmet.
- `git diff --check -- docs/dev-log/cran-0.11-audit/GATES.md .unlazy/cran-0.11-audit/GATES.md`: passed.
- Opened the six live routes recorded above with `?audit=20261008-2240` in the Codex in-app browser. Cache-busted `trees300` shows Hackett-only/MIT wording, contrary to an initial cached view; the cache-busted response is recorded as current. All six route observations support the corrected ledger entry.
- The direct structural validator, `Rscript -e 'source("tools/check-after-task.R"); check_after_task(<absolute report path>)'`, exited 0. `closeout.py check <absolute report path>` failed because it also reverified the brain repository's unrelated acceptance ledgers, which contain unmet gates. The pigauto Unlazy ledgers were checked separately in status mode.
- `slop_check.py` reported 1.1 hits per 1,000 words, zero findings, and zero em dashes.
- Current local Git history confirms PR #230 merge commit `0b0f71fee838c6ed51ef832ed819270eeafaf29b`, with candidate HEAD `ea756083e751e37f650f8801a9fc32e81c40a0c5` in its ancestry. Pages run #583 is recorded as successful on that merge.
- The public GitHub page for PR #228 shows Draft, four commits, and predecessor artifact evidence. The PR #230 web page was unavailable in the browser fetch; merge status was checked from local `origin/main` and independently reviewed.
- Independent read-only review confirmed the ledger updates and identified the deployed `trees300` route as unverified after #583.

## 6. Tests of the Tests

Before repair, the Unlazy parser rejected unindented `CHECK` and `EVIDENCE` lines as unattached to a gate. After converting those historical labels to ordinary prose, status mode parsed all 11 gates and reported the open gates without executing their commands. This demonstrates that malformed ledger syntax is detected and that status checks remain non-executing.

## 7a. Issue Ledger

- Fixed: stale statement that PR #230 was still a draft and undeployed.
- Fixed: stale claim that candidate `R/data.R` incorrectly assigned the underlying tree-data licence to `megatrees`.
- Corrected: cache-busted live `trees300` says Hackett-only and still attributes the MIT licence to the `megatrees` package. The initial non-cache-busted browser view differed.
- Open: source/help and live-site mismatches, local visual review, exact final tarball, final artifact review, Windows evidence, and the basis for the maintainer's redistribution warranty.

## 8. Consistency Audit

Compared both ledgers with `origin/main`, candidate source and generated help, the 2026-10-08 site review, current lease registry, public PR #228, and the six cache-busted live URLs. Candidate `tree300` and `tree_full` manuals remain stale; candidate `trees300` help is corrected. Current deployed `trees300` remains Hackett-only and retains the misleading MIT attribution. The tree article lacks the candidate retrieval section, getting-started calls `tree300` an MCC tree, and the MI article lacks the candidate's legacy-only inverse-Wishart boundary. Leased source/help files were not edited.

## 9. What Did Not Go Smoothly

The brain-rooted `closeout.py new` helper resolves relative paths under the brain repository. I initially passed a relative pigauto path, which created a blank template in the wrong repository. I removed that exact file after confirming its contents, then used the pigauto worktree's absolute path for the correct report. The closeout wrapper's validator also rechecked unrelated brain-root gates and failed; the report's direct structural check passed. GitHub API access also failed; PR #228 was checked in the browser, while PR #230's page fetch was unavailable and its merge was verified from local Git history. A non-cache-busted live page view also differed from the cache-busted response, so the latter is the basis for current-page claims.

## 10. Known Residuals

This is a deployed-page evidence refresh, not release completion. The top-level ledger has 6 unmet gates; the implementation ledger has I7 unmet. No final tarball or exact-artifact review exists for the current uncommitted source. Local visual review remains blocked by the browser's file-URL policy. The BirdTree-derived objects remain bundled, and the rights report records the maintainer's confirmation with attribution while the underlying redistribution basis remains undocumented. No merge, deployment, or CRAN submission was performed.

## 11. Team Learning

When a project uses a tool rooted in another repository, inspect how it resolves report paths before asking it to create files. Use an absolute target for external project reports, and remove only the known accidental file if path resolution goes wrong.

Memory receipt: loaded the pigauto `route.py` manifest and ran lane preflight; the prediction, uncertainty, and shared-lane constraints did not require code changes in this slice. Golden Set: not checked because this slice reconciled release records and parser state rather than changing a recurring product defect.

## 12. Cross-Product Coverage

Covers current PR merge identity, Pages run identity, release-ledger syntax, and current deployed text on six named pages. Does NOT cover runtime defaults, optional-backend execution, recovery, local visual rendering, other retained or retired live routes, sitemap contents, rights clearance, final tarball checks, or Windows results.
