## 1. Goal

Refresh current live public-site evidence for the pigauto CRAN 0.11 audit, focusing on whether the retired simulation-study route remains discoverable and directly served.

## 2. Implemented

Recorded a fresh Chrome check of the live article index and direct retired route in `GATES.md`. Kept G7 open because the deployed site still links to and serves the retired article.

## 3a. Decisions and Rejected Alternatives

Treat the live article-index link and direct page as evidence that deployment retirement is incomplete. Do not infer the live sitemap or search-index contents from the candidate site's local crawl, because Chrome returned `ERR_BLOCKED_BY_CLIENT` for both live endpoints. Do not mark G7 met before the source PR is merged, deployed, and verified against that deployment.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-live-index-refresh.md`

## 5. Checks Run

- Chrome opened `https://itchyshin.github.io/pigauto/articles/index.html?audit=20261009lead2`; the page listed “Four ways to impute a phylogenetic trait matrix” at `/articles/simulation-study.html`.
- Chrome opened `https://itchyshin.github.io/pigauto/articles/simulation-study.html?audit=20261009lead`; the complete historical article was served, including its title, tables, results, reproduction section, and notice that it is excluded from the intended 0.11.0 package build and public site.
- Chrome navigation to `https://itchyshin.github.io/pigauto/sitemap.xml` and `/search.json` returned `net::ERR_BLOCKED_BY_CLIENT`; their live contents remain unverified.
- `node ~/.codex/skills/unlazy/scripts/gate-check.mjs --status /Users/z3437171/.codex/worktrees/cran-011-gate-reconcile/pigauto/docs/dev-log/cran-0.11-audit/GATES.md` reported 8 met and 3 unmet gates: G7, G8, and G9.

## 6. Tests of the Tests

No package code or tests changed. The public article index and direct retired URL independently expose the same stale route. The sitemap and search checks did not reach the server and are not counted as passes.

## 7a. Issue Ledger

Confirmed: the live article index still links to the retired simulation-study article, and the direct URL still serves its full content. Open: sitemap and search-index membership; verification after merge and deployment; exact artifact and independent final review.

## 8. Consistency Audit

Compared the browser results with the current PR #231 state and the release ledger. PR #231 remains Draft and unmerged at `dee8f02888fd2c4c22d426e595ea6451a28a280e`; the candidate source's local crawl is not evidence of current deployed behavior. The public-site result is consistent with G7 remaining open. The release ledger remains 8 of 11 gates met.

## 9. What Did Not Go Smoothly

Chrome blocked the live sitemap and search-index requests with `ERR_BLOCKED_BY_CLIENT`. The closeout helper initially resolved a relative output path against the Shinichi hub; that accidental file was removed, and this report was then generated at its absolute pigauto worktree path.

## 10. Known Residuals

This report does NOT establish live sitemap/search contents, deployment-to-commit identity, retirement across every old URL, the final frozen tarball, platform checks on that tarball, or CRAN readiness. G7, G8, and G9 remain open.

## 11. Team Learning

A local retirement crawl does not establish what users can still reach on the deployed site. Recheck both the article index and direct legacy URL after deployment; treat browser-blocked sitemap or search requests as unverified.

Memory receipt: `route.py pigauto` loaded the LOAD-FIRST manifest, and the repo's after-task protocol was read. The source-versus-deployment boundary and user-facing route checks shaped this work.

Golden Set: no package source behavior changed; no source known-mistake class was in scope.

## 12. Cross-Product Coverage

Covers the live article index, direct retired article, PR state, and release-gate tally.

Does NOT cover live sitemap/search contents, every old or retained route, post-merge deployment, exact artifact checks, or CRAN submission.
