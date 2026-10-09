## 1. Goal

Close local site review for the pigauto 0.11 candidate and refresh the release ledger against the current source commit.

## 2. Implemented

Updated `GATES.md` with the fresh site build and crawl evidence for commit `34e1e3cd5a128d309782bcf43fd92a94298b175e`, Chrome visual review of five reader pages, and the live PR #231 head and CI state. Corrected the earlier G6 note that wrongly treated loopback HTTP review as prohibited. G6 now has machine-readable manual evidence and the ledger reports 7 met and 4 unmet gates.

## 3a. Decisions and Rejected Alternatives

Used Chrome with the existing loopback-only site server after direct `file://` navigation failed. Kept local visual review distinct from live deployment verification. Did not treat candidate-source checks as frozen-artifact evidence. Did not merge, deploy, freeze a post-merge tarball, or submit to CRAN.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-site-visual-review.md`

## 5. Checks Run

- `bash ~/shinichi-brain/tools/lane_preflight.sh /Users/z3437171/Dropbox/Github\ Local/pigauto` reported one active lane and 58 total worktrees.
- Checked the retained exact-head pkgdown build and crawl logs. Build output records 62 HTML pages, 613 search entries, 3,443 local references, 34 retired routes, zero crawl errors, `SITE_CRAWL_OK`, and a clean `pkgdown::check_pkgdown()` result. Both log checksums are recorded in `GATES.md`.
- Reviewed the local rendered home page, Getting Started, multiple-imputation, tree-uncertainty, and `trees300` reference pages in Chrome.
- Chrome showed PR #231 remains Draft at 21 commits and head `34e1e3c`. Its run #37867499185 completed successfully: Ubuntu R release passed in 10m00, Ubuntu R-devel in 11m42, and macOS R release in 17m18, including MPS prediction tests and `R CMD check`.
- `node ~/.codex/skills/unlazy/scripts/gate-check.mjs --status docs/dev-log/cran-0.11-audit/GATES.md` reports 7 met and 4 unmet gates. It parses the ledger; it does not rerun the gate commands.
- `git diff --check` passed.

## 6. Tests of the Tests

The site crawler has a retained planted retirement negative control at `provenance/site-crawler-simulation-retirement-control.txt`; it was not rerun in this slice. The visual review is manual and has no automated layout oracle.

## 7a. Issue Ledger

- G0 remains open: rights for redistributing bundled BirdTree-derived data are not established by the external pages reviewed, and no final post-merge source revision or artifact exists.
- G7 remains open: the corrected pages have not been deployed and checked against the full route, sitemap, search, and retirement inventory.
- G8 remains open: no exact post-merge tarball is frozen and checked across required platforms.
- G9 remains open: independent reviewers have not voted on that final artifact and deployment evidence.

## 8. Consistency Audit

Compared G6 and G5/G5b claims with the fresh build logs, crawler output, and the rendered Chrome pages. The pages show the corrected tree provenance and retrieval guidance, fixed-effect MI boundaries, and current reader-facing navigation. The G0 rights question and G7 deployment requirement remain separate.

## 9. What Did Not Go Smoothly

Chrome refused local `file://` navigation. The temporary loopback-only server allowed the requested Chrome review and resolved that limitation. The closeout wrapper resolved its root to the brain repository rather than this pigauto worktree, so the direct after-task structure validator was used with its unrelated ledger recheck disabled.

## 10. Known Residuals

This slice does NOT establish BirdTree redistribution rights, deployed-page correctness, frozen-tarball validity, checksum-bound platform results, or CRAN submission readiness. The site crawler does not check CSS URLs, `srcset`, or JavaScript-generated URLs.

## 11. Team Learning

A browser's refusal of `file://` does not by itself prevent local visual review; a loopback-only HTTP server can provide a normal rendered-page review when the browser permits it. The routed `pigauto` LOAD-FIRST manifest and repository instructions were consulted. Their prediction-path and uncertainty guards did not affect this documentation-only slice. Golden Set regression was not in scope.

## 12. Cross-Product Coverage

Covers: exact-commit local site build evidence, route and search retirement checks, and Chrome visual review of five rendered pages.

Does NOT cover: deployed website state, package runtime behavior, a post-merge tarball, checksum-bound external artifact checks, full CSS/srcset/JavaScript link behavior, redistribution rights, or CRAN submission.
