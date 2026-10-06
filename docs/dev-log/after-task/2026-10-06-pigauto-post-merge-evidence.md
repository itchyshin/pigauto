# After-task: pigauto 0.11.0 post-merge evidence

## 1. Goal

Verify the merged 0.11.0 deployment and exact artifact evidence, refresh the audit checkpoint, and identify the next unmet release gate.

## 2. Implemented

Confirmed PR #226 merged as `bb5835d1b214d783da7b0b414df99aa6ba926bc7`. Verified Pages run #570 succeeded on that commit and inspected the live site. All 34 retired benchmark URLs and four retired walkthrough URLs returned 404. The sitemap and search index contained no retired route.

Confirmed the exact merged-source tarball identity and clean path scan. Its full macOS/R 4.6.0 `R CMD check --as-cran --no-manual` completed with `Status: OK`. Actions run #37526171411 passed Ubuntu R release, Ubuntu R-devel, and macOS R release on the merged source commit.

Updated the checkpoint and pre-merge records to distinguish historical pending claims from current proof. G6 and the bounded G8 audit criterion pass. G7 remains partial while exact-tarball Windows R-release and R-devel checks are pending. The release authorization remains NOT_READY because BirdTree redistribution rights remain open.

Reconciled the recovery child gates against the retained per-seed receipts and summary tables: both registered ten-seed adapter-source campaigns meet their bounded mechanism and bias-margin criteria. Corrected a stale compute estimate to include the two-worker feasibility result that reduced the projected remaining gllvm campaign below the three-hour line. Reconciled the surface child gates with the completed local-build record and post-merge live receipts. A later redundant local rebuild stalled during `mixed-types.Rmd`; its partial output failed the crawler and is recorded as incomplete, without replacing the earlier completed build evidence.

## 3a. Decisions and Rejected Alternatives

Kept the release ledger `NOT_READY` and its release artifact/rights fields blank, as the rights record directs while BirdTree redistribution rights remain unresolved. Recorded the exact merged-source audit tarball and checks separately; this does not declare CRAN readiness.

Used a fresh clone of public `main` at the verified merge commit for the evidence branch. The prior temporary clone had no remote history and was unsuitable as a PR base. Did not change the dirty shared checkout, merge anything, contact an upstream data holder, or submit to CRAN.

## 4. Files Touched

Evidence payload prepared for a separate PR:

- `LOOP/lanes/cran-0.11-audit/checkpoint.md`
- `docs/dev-log/cran-0.11-audit/post-merge-verification-2026-10-06.md`
- `docs/dev-log/cran-0.11-audit/provenance/post-merge-artifact-identity.json`
- `docs/dev-log/cran-0.11-audit/provenance/post-merge-R-CMD-check-macos.log`
- `docs/dev-log/cran-0.11-audit/provenance/post-merge-testthat-macos.Rout`
- `docs/dev-log/cran-0.11-audit/provenance/post-review-live-page-recheck-2026-10-06.md`
- `docs/dev-log/cran-0.11-audit/provenance/release-gate-selftest-2026-10-06.log`
- `docs/dev-log/cran-0.11-audit/provenance/review-panel-2026-10-06.md`
- `docs/dev-log/cran-0.11-audit/provenance/win-builder-exact-tarball-uploads-2026-10-06.md`
- `docs/dev-log/cran-0.11-audit/provenance/public-deployment-live-check-2026-10-06.tsv`
- `docs/dev-log/cran-0.11-audit/provenance/pre-pr-check.md`
- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `docs/dev-log/cran-0.11-audit/provenance/rights-and-policy.md`
- `docs/dev-log/cran-0.11-audit/retirement-manifest.md`
- `docs/dev-log/cran-0.11-audit/surface-inventory.md`
- `docs/dev-log/cran-0.11-audit/provenance/compute-estimates.md`
- `docs/dev-log/after-task/2026-10-06-pigauto-post-merge-evidence.md`

The generated reports and retained logs below are ignored by the broad `docs/` rule and need explicit staging if approved for the PR. The local `.unlazy` gate receipts stay outside the PR:

- `docs/dev-log/after-task/2026-10-06-pigauto-post-merge-evidence.md`
- `docs/dev-log/cran-0.11-audit/post-merge-verification-2026-10-06.md`
- `docs/dev-log/cran-0.11-audit/provenance/post-merge-artifact-identity.json`
- `docs/dev-log/cran-0.11-audit/provenance/post-merge-R-CMD-check-macos.log`
- `docs/dev-log/cran-0.11-audit/provenance/post-merge-testthat-macos.Rout`
- `docs/dev-log/cran-0.11-audit/provenance/post-review-live-page-recheck-2026-10-06.md`
- `docs/dev-log/cran-0.11-audit/provenance/release-gate-selftest-2026-10-06.log`
- `docs/dev-log/cran-0.11-audit/provenance/review-panel-2026-10-06.md`
- `docs/dev-log/cran-0.11-audit/provenance/win-builder-exact-tarball-uploads-2026-10-06.md`
- `docs/dev-log/cran-0.11-audit/provenance/public-deployment-live-check-2026-10-06.tsv`
- `.unlazy/cran-0.11-audit/GATES.md` (G7 partial; G8 complete)
- `.unlazy/cran-0.11-audit/adapters/GATES.md`

## 5. Checks Run

- Lane preflight and `route.py pigauto`: completed in the isolated clone; narrow lane lease granted.
- GitHub Actions: Pages run #570 succeeded on the merge commit. R-CMD-check run #37526171411 passed all three jobs: Ubuntu R release, Ubuntu R-devel, and macOS R release.
- Public deployment: home, getting-started, multiple-imputation article, and `multi_impute()` reference returned HTTP 200. The expected navbar labels appeared and no per-trait benchmark links were present.
- Fresh reader-surface recheck: cache-busted primary-browser loads of the homepage, getting-started article, multiple-imputation article, and `multi_impute()` reference showed current 0.11.0 content. The browser does not expose response-body hashes; the distinction and observation markers are recorded in `post-review-live-page-recheck-2026-10-06.md`. An independent web-index result was stale; direct browser rechecks resolved it.
- Public discovery and retirement: sitemap returned 200 with 64 locations; search index returned 200 with 574 paths. No retired route appeared in either. All 34 retired `/dev/` HTML URLs and four old walkthrough URLs returned HTTP 404.
- Artifact: `pigauto_0.11.0.tar.gz`, source commit `bb5835d1b214d783da7b0b414df99aa6ba926bc7`, SHA-256 `d34d469981386546b4ad03b9c7aab81674556708eeaf524c1c7365b6c68dec26`, 5,122,013 bytes, 247 members. The forbidden development-path scan found zero hits.
- Exact artifact check: macOS Tahoe 26.7, R 4.6.0, `R CMD check --as-cran --no-manual`, `Status: OK`. Testthat summary: 3,066 PASS, 0 FAIL, 161 WARN, 86 SKIP. Full log SHA-256: `657a19dc6584f1475f78f35d96824304be45b1279459f72bfa5f1c72113d41a0`.
- CRAN validator selftest: passed, including 13 planted negative controls; exact output is retained in `provenance/release-gate-selftest-2026-10-06.log` (SHA-256 `21dff08b96b98af6bb6b6c17729898e9381f872820c060285916851bfbd49b13`).
- Fresh read-only panel: Grace, Rose, and Pat all returned READY for the bounded audit evidence. Their evidence scopes and the initial stale-index discrepancy are in `provenance/review-panel-2026-10-06.md`.
- Exact-artifact Windows checks: win-builder R-release and R-devel forms accepted the same filename and size; results are pending. The upload receipt and artifact SHA are in `provenance/win-builder-exact-tarball-uploads-2026-10-06.md`.
- After-task structure validator: passed. The full closeout compiler then stopped on unmet gates in this lane's G7, the audit child lanes (adapters, defaults, recovery, surfaces), and the separate imputation-simulation lane. It reported a well-formed report with unfinished work. Those other lane receipts were left unchanged.
- Adapter matrix receipt from the adapter lane: all four isolated library modes passed 20 adapter and 46 provenance expectations each; real-mode counts were 0 with two expected skips, 8, 8, and 16.

## 6. Tests of the Tests

The earlier site crawler's deliberate broken-link and retired-search-path controls passed under G5; they were not rerun in this slice. This slice directly compared the live sitemap and indexed paths against the 34 retired routes, then checked each route's HTTP status. The CRAN release-gate selftest passed the planted negative controls for unlicensed data, bad predecessor hash, forbidden tarball entry, omitted incoming/timing/upload evidence, broken or fictional live URLs, omitted compiled gate, not-ready panel, invalid release type, and omitted large-vignette budget.

In the continued evidence pass, a separate package-site rebuild stalled while rendering `mixed-types.Rmd`. Its partial output was intentionally sent through the site crawler; the crawler rejected the incomplete output for missing articles, sitemap/search files, and internal pages. This failed repeat is not treated as a pass; the earlier complete local build remains the evidence recorded under G5.

## 7a. Issue Ledger

- G6, merged deployment and live verification: complete.
- G7, exact artifact identity/check and platform checks: partial. The exact tarball passed on macOS and merged-commit source CI passed on Ubuntu release, Ubuntu devel, and macOS release. Exact-tarball Windows R-release and R-devel results are pending.
- G8, fresh independent panel, negative-validator controls, and re-verification: complete at the bounded audit-evidence rung. All three reviewers returned READY; direct cache-busted checks resolved the stale-index discrepancy.
- BirdTree data redistribution rights: open; no permission claim or upstream contact.
- CRAN submission: not performed.

## 8. Consistency Audit

Qualified the pre-merge deployment observations in the surface inventory and retirement manifest. Updated the checkpoint, pre-PR check note, and rights note so older pending statements cannot be mistaken for current status. The current release ledger remains `NOT_READY` by design. The package source and API were not changed.

This continuation also removed stale pending statements from the surface child ledger and aligned its S3/S4 statuses with the completed G5/G6 receipts. The recovery child ledger now records the ten-seed results and their source-workflow limits. The compute estimate now records both the initial over-three-hour serial projection and the measured parallel projection used before the continuation.

## 9. What Did Not Go Smoothly

The browser blocked direct XML/JSON reads, so the live sitemap, search index, and retired URLs required a read-only HTTP check. The first isolated clone was a root snapshot with no usable PR history; a fresh clone of public `main` fixed that. An earlier status message described the exact R CMD check as partial because it followed the network diagnostic log; inspecting the retained `00check.log` and `testthat.Rout` showed the full check had completed with `Status: OK`. The audit record now reflects the completed check.

## 10. Known Residuals

The exact tarball's full check was run on macOS; the three Actions platform jobs ran on the exact merged source commit. The exact tarball was then uploaded to win-builder for Windows R-release and R-devel; results remain pending. Until both logs are inspected, G7 remains partial. The latest separate local site rebuild did not complete; the earlier completed site build and the deployed pages were verified. The fresh panel resolved the bounded audit criterion. BirdTree redistribution rights are unresolved, so the release ledger remains `NOT_READY`. The tarball is retained locally under `/private/tmp/pigauto-cran-011-gate-artifacts/`; win-builder received it for prechecks, and CRAN has received no submission. This evidence branch has not yet been pushed or opened as a PR.

## 11. Team Learning

Memory receipt: loaded the repository manifest with `route.py pigauto` and the project instructions; they shaped the scope, lane ownership, and evidence boundaries. Searched the brain for prior pigauto lane guidance. Golden Set: not run; no model implementation or statistical-method change was in scope.

## 12. Cross-Product Coverage

Covers the merged source's Ubuntu R release, Ubuntu R-devel, and macOS R release checks; the exact merged-source tarball's full macOS/R 4.6.0 check and pending win-builder Windows R-release/R-devel checks; the fresh independent release-evidence panel; validator negative controls; and the published pkgdown reader, discovery and retirement surfaces.

Does NOT cover the results of the two pending Windows exact-tarball checks, a byte-identical tarball check on every CI platform, BirdTree redistribution permission, CRAN acceptance, or CRAN submission.
