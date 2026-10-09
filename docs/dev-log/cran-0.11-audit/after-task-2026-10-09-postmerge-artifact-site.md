## 1. Goal

Record the current deployed-site retirement result and freeze one post-merge 0.11.0 artifact with exact local-check evidence. Keep PR #228 as an unmerged evidence PR.

## 2. Implemented

Updated the gate ledger, retirement manifest, and release-ledger JSON to bind the deployment to `d76804e` and the frozen tarball to its SHA-256. Added the build receipt. Added the corrected live-site verifier and its full route capture after an independent review found that the first sitemap filter omitted three named coordination routes. The corrected recheck closes G7; G8 and G9 remain open.

## 3a. Decisions and Rejected Alternatives

Kept the accepted artifact as the single candidate and discarded the earlier archive containing `.-00build.log`. Did not rerun the site verifier after shell DNS failed; retained the prior successful HTTP evidence and recorded the rerun failure. No scientific defaults changed. Assumptions made without asking: none.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/retirement-manifest.md`
- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `docs/dev-log/cran-0.11-audit/provenance/post-merge-build-receipt-2026-10-09.md`
- `docs/dev-log/cran-0.11-audit/provenance/post-merge-artifact-identity.json`
- `docs/dev-log/cran-0.11-audit/provenance/post-merge-tarball-inventory.txt`
- `docs/dev-log/cran-0.11-audit/provenance/post-merge-R-CMD-check.log`
- `docs/dev-log/cran-0.11-audit/provenance/post-merge-testthat.Rout`
- `docs/dev-log/cran-0.11-audit/provenance/live-retired-routes-2026-10-09.tsv`
- `docs/dev-log/cran-0.11-audit/provenance/live-site-discovery-2026-10-09.json` (initial filter, superseded)
- `docs/dev-log/cran-0.11-audit/provenance/live-retired-routes-2026-10-09-v2.tsv`
- `docs/dev-log/cran-0.11-audit/provenance/live-retained-routes-2026-10-09-v2.tsv`
- `docs/dev-log/cran-0.11-audit/provenance/live-site-discovery-2026-10-09-v2.json`
- `docs/dev-log/cran-0.11-audit/provenance/live-site-verification-2026-10-09-v2.log`
- `docs/dev-log/cran-0.11-audit/provenance/live-site-verifier-tests-2026-10-09.log`
- `docs/dev-log/cran-0.11-audit/provenance/live-site-recheck-2026-10-09-v2.md`
- `script/cran-0.11-site/verify_live_site.py`
- `script/cran-0.11-site/test_verify_live_site.py`
- This after-task report.

## 5. Checks Run

- Exact tarball: `env -u NOT_CRAN _R_CHECK_FORCE_SUGGESTS_=true R CMD check --as-cran --run-donttest pigauto_0.11.0.tar.gz`; macOS arm64/R 4.6.0; `Status: OK`; 3,190 passes, 0 failures, 161 warnings, 86 skips.
- Artifact inspection: 250 archive entries; exact SHA-256 and inventory recorded; hidden-path, forbidden-path, optional-backend dependency, and shipped NOTICE/tree-object scans passed.
- Live site: the corrected verifier checked 44 retired routes (all 404), all 62 sitemap URLs (60 content pages, the intentionally retained noindex `validation_suite.html` tombstone, and the expected `/404.html` utility route), and 613 search entries, with no retired sitemap or search targets. The first summary had an incomplete sitemap filter; the corrected filter and fresh successful capture are retained in `provenance/live-site-recheck-2026-10-09-v2.md`.
- Deployment: Pages workflow run #37975894397 succeeded on merged commit `d76804e`.
- `git diff --cached --check` passed before this update.

## 6. Tests of the Tests

The exact-artifact R CMD check exercised package installation, examples, documentation, code checks, and the full testthat suite. The site verifier requires HTTP 404 for each retired route and rejects retired paths in sitemap/search targets; the saved successful result has no failures. A fresh shell rerun could not exercise those assertions because DNS failed before the first response.

## 7a. Issue Ledger

- Resolved: the first post-merge archive included the build log as a hidden file and produced a NOTE; discarded and rebuilt cleanly.
- Open: checksum-bound Win-builder R-release and R-devel result logs.
- Open: independent exact-artifact and deployed-site verdicts from Grace, Rose, and Pat.
- Open: any separate rights evidence beyond the maintainer warranty; no separate published redistribution grant was found.

## 8. Consistency Audit

Cross-checked GATES, release-ledger, retirement-manifest, archive identity, check log, and corrected deployed-site receipt. The first empty sitemap-target list was superseded after its filter omission was found; the corrected capture checks all 62 live sitemap URLs. The source version is 0.11.0 at merged commit `d76804e`. Both optional model packages remain absent from DESCRIPTION dependency fields. The artifact includes `inst/NOTICE`; its hashes identify the shipped tree objects. Historical page and search text remains in Git/changelog history while retired route targets are absent from live sitemap/search.

## 9. What Did Not Go Smoothly

The initial build embedded `.-00build.log` and had to be discarded. The sandboxed shell could not resolve DNS, but the approved read-only network check succeeded and produced the corrected live-site capture. The managed worktree is outside the default writable root, so file updates require the approved sandbox escalation path. The closeout compiler also reported unmet checks in separate imputation-simulation ledgers in the brain root; those belong to another project lane and were left untouched.

## 10. Known Residuals

G7 is met for the checked deployment, with 44 retired routes and all 62 sitemap routes verified. G8 is incomplete until the exact artifact's Windows results are available and bound to its hash. G9 has not passed. The release-evidence PR remains Draft and unmerged. Nothing has been submitted to CRAN.

## 11. Team Learning

Keep the archive-cleanliness scan before freezing: build logging can add a hidden file that creates a CRAN NOTE even when package tests pass. Preserve one exact hash through local and platform checks.

Memory receipt: loaded the pigauto LOAD-FIRST manifest via `route.py pigauto`; its focus on defaults, tree provenance, and exact release checks shaped this work.

Golden Set: Not run; no package behavior changed and no known-mistake code class was in scope.

## 12. Cross-Product Coverage

This update covers one merged macOS/R 4.6.0 artifact check and one deployed pkgdown site, including all 62 sitemap routes and 44 retired routes. It does NOT cover checksum-bound Windows output, future site deployments, CRAN acceptance, or statistical validity outside the previously bounded MI recovery regimes.


## Current reconciliation update (2026-10-09)

Recorded the current R-release and R-devel Win-builder logs, preserving their raw bytes in deterministic gzip files with both raw and compressed SHA-256 values. Both report `Status: 1 ERROR` and five failed expectations; neither identifies the uploaded archive hash. Updated the release ledger to keep G8 and G9 open and to record Grace/Rose as NOT READY and Pat as NOT ASSESSED for durable site evidence. The live-site v2 receipts were ignored by the repository ignore rule, so this update force-adds those specific receipts to the evidence branch. The homepage warning rendering remains an unresolved reader-facing discrepancy. No source change, merge, or CRAN submission was made.


## Pat re-review update (2026-10-09)

Pat independently reviewed the now-committed site receipts and passes G7 for the checked deployment. This supersedes the prior NOT ASSESSED verdict on receipt durability. The review confirms the receipt hashes, Pages run `37975894397` at deployment commit `d76804e768bf77f430f43dcf591633f9cd900dab`, the 44 retired routes, 62 sitemap URLs, 613 search entries, and three verifier tests. Chrome still shows a literal `[!WARNING]` token on the homepage; Pat could not directly fetch the live sitemap in this re-review, so that check is supported by the committed verifier capture. G8/G9 remain open because Windows logs fail without binding to the frozen artifact hash and Grace/Rose remain NOT READY. PR #228 stays Draft and unmerged; no CRAN submission occurred.
