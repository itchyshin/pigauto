## 1. Goal

Record the current deployed-site retirement result and freeze one post-merge 0.11.0 artifact with exact local-check evidence. Keep PR #228 as an unmerged evidence PR.

## 2. Implemented

Updated the gate ledger, retirement manifest, and release-ledger JSON to bind the deployment to `d76804e` and the frozen tarball to its SHA-256. Added the build receipt. G7 was reopened after an independent review found that the first sitemap filter omitted three named coordination routes; G8 and G9 also remain open.

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
- `docs/dev-log/cran-0.11-audit/provenance/live-site-discovery-2026-10-09.json`
- `script/cran-0.11-site/verify_live_site.py`
- This after-task report.

## 5. Checks Run

- Exact tarball: `env -u NOT_CRAN _R_CHECK_FORCE_SUGGESTS_=true R CMD check --as-cran --run-donttest pigauto_0.11.0.tar.gz`; macOS arm64/R 4.6.0; `Status: OK`; 3,190 passes, 0 failures, 161 warnings, 86 skips.
- Artifact inspection: 250 archive entries; exact SHA-256 and inventory recorded; hidden-path, forbidden-path, optional-backend dependency, and shipped NOTICE/tree-object scans passed.
- Live site: 44 direct retired-route checks returned 404, and search reported 613 entries with no retired route targets. The saved discovery summary reports 62 sitemap entries, but its filter omitted `/AGENTS.html`, `/CLAUDE.html`, and `/goodagents.html`; its empty sitemap-target result is not sufficient. A corrected offline filter test passes, but live sitemap verification remains pending because browser access is blocked and shell DNS failed before requests.
- Deployment: Pages workflow run #37975894397 succeeded on merged commit `d76804e`.
- `git diff --cached --check` passed before this update.

## 6. Tests of the Tests

The exact-artifact R CMD check exercised package installation, examples, documentation, code checks, and the full testthat suite. The site verifier requires HTTP 404 for each retired route and rejects retired paths in sitemap/search targets; the saved successful result has no failures. A fresh shell rerun could not exercise those assertions because DNS failed before the first response.

## 7a. Issue Ledger

- Resolved: the first post-merge archive included the build log as a hidden file and produced a NOTE; discarded and rebuilt cleanly.
- Open: corrected live sitemap capture for the three coordination routes.
- Open: checksum-bound Win-builder R-release and R-devel result logs.
- Open: independent exact-artifact and deployed-site verdicts from Grace, Rose, and Pat.
- Open: any separate rights evidence beyond the maintainer warranty; no separate published redistribution grant was found.

## 8. Consistency Audit

Cross-checked GATES, release-ledger, retirement-manifest, archive identity, check log, and the partial deployed-site receipt. The original empty sitemap-target list was downgraded because its filter omitted three named routes. The source version is 0.11.0 at merged commit `d76804e`. Both optional model packages remain absent from DESCRIPTION dependency fields. The artifact includes `inst/NOTICE`; its hashes identify the shipped tree objects. Historical page and search text remains in Git/changelog history while retired route targets are absent from live sitemap/search.

## 9. What Did Not Go Smoothly

The initial build embedded `.-00build.log` and had to be discarded. Shell DNS then prevented rerunning the live verifier. The managed worktree is outside the default writable root, so file updates require the approved sandbox escalation path. The closeout compiler also reported unmet checks in separate imputation-simulation ledgers in the brain root; those belong to another project lane and were left untouched.

## 10. Known Residuals

G7 remains open until the three coordination routes are checked against a captured live sitemap with the corrected filter. G8 is incomplete until the exact artifact's Windows results are available and bound to its hash. G9 has not passed. The release-evidence PR remains Draft and unmerged. Nothing has been submitted to CRAN.

## 11. Team Learning

Keep the archive-cleanliness scan before freezing: build logging can add a hidden file that creates a CRAN NOTE even when package tests pass. Preserve one exact hash through local and platform checks.

Memory receipt: loaded the pigauto LOAD-FIRST manifest via `route.py pigauto`; its focus on defaults, tree provenance, and exact release checks shaped this work.

Golden Set: Not run; no package behavior changed and no known-mistake code class was in scope.

## 12. Cross-Product Coverage

This update covers one merged macOS/R 4.6.0 artifact check and one deployed pkgdown site. It does NOT cover checksum-bound Windows output, future site deployments, CRAN acceptance, or statistical validity outside the previously bounded MI recovery regimes.
