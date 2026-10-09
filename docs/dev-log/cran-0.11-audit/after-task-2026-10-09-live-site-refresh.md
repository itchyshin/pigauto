## 1. Goal

Refresh the live-site baseline for G7 using Chrome, compare it with the current candidate source, and preserve the discrepancies on the release-evidence PR.

## 2. Implemented

Recorded a fresh cache-busted review of Getting Started, the `tree300` reference, the multiple-imputation article, the article index, and the retired simulation-study URL. Updated G7 without marking it met, because the deployed site remains partly out of sync with candidate source and still serves the retired page.

## 3a. Decisions and Rejected Alternatives

Followed Shinichi's instruction that the pigauto work is one lane with multiple checkouts. Treated this worktree as the evidence checkout within that lane. Kept G7 open because PR #231 remains unmerged and the live site has unresolved route and text differences. Did not merge, deploy, or infer current sitemap/search state from the local build.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/site-review-chrome-2026-10-09.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-live-site-refresh.md`

## 5. Checks Run

- Opened cache-busted deployed pages in Chrome on 2026-10-09: Getting Started, `tree300`, multiple imputation, articles index, and `simulation-study.html`.
- Confirmed the live article index links to the retired page and the direct URL serves its full article content.
- Confirmed the live `tree300` reference has moved beyond its earlier MCC wording, while Getting Started still uses that label.
- Confirmed the candidate source describes posterior member 69 and limits the inverse-Wishart option to legacy reproduction.
- Sitemap check: Chrome returned `ERR_BLOCKED_BY_CLIENT`; web lookup also could not access it. Search-index membership was not checked.
- Unlazy status remains 8 met and 3 unmet gates (G7-G9). This site review adds evidence but does not close G7.
- The after-task structure check passed. Its wrapper then exited 1 because five `.unlazy/imputation-sim` gates outside this CRAN audit remain unmet.
- `slop_check.py` returned 0 findings on both reports. `git diff --cached --check` passed for the staged text files.

## 6. Tests of the Tests

This slice is a read-only live-site observation, not a test implementation. Each page was opened directly with a cache-busting query. The retired page was independently confirmed through both the article index link and the direct route. The sitemap and search checks are explicitly unverified.

## 7a. Issue Ledger

Confirmed current public mismatches: stale MCC wording in Getting Started; missing member-69 detail on the live `tree300` reference; missing legacy-reproduction boundary in the live MI article; and a retired article still linked and served. Sitemap and search remain unverified.

## 8. Consistency Audit

Compared the live page text with current PR #231 source. The reference-page result partially supersedes the older 2026-10-08 observation. The source change still has not been deployed, so the live site cannot yet demonstrate those edits. The retired-route result matches the earlier pre-deployment concern and confirms it remains present.

## 9. What Did Not Go Smoothly

Chrome blocked the sitemap URL, and the web lookup tool could not retrieve it. The current page set therefore does not establish sitemap or search-index state. The after-task wrapper also reports five unmet `.unlazy/imputation-sim` gates outside this CRAN audit after confirming report structure.

## 10. Known Residuals

This slice does NOT bind live pages to a deployment commit, verify every retained or retired route, confirm search or sitemap contents, or close G7. It does not establish exact-tarball or platform results. Merge and deployment remain under Shinichi's control.

## 11. Team Learning

A retired page can remain in both the live article index and direct-route output even when its own text says it is excluded from the package build. Check both discovery links and direct URLs after deployment.

Memory receipt: loaded the pigauto LOAD-FIRST manifest with `route.py`; its release-boundary and source-versus-deployment distinctions shaped this check.

Golden Set: no package source behavior changed; no source known-mistake class was in scope.

## 12. Cross-Product Coverage

Covers five live browser pages and the comparison with their candidate source counterparts.

Does NOT cover sitemap or search-index contents, every retained or retired route, deployment-to-commit binding, exact-artifact validation, platform checks, or CRAN submission.
