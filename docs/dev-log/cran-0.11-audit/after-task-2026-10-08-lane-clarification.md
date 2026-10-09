## 1. Goal

Reconcile the single-lane pigauto CRAN audit state, refresh reader-site evidence, and verify current defaults checks.

## 2. Implemented

Updated the public BirdTree guidance and removed unsupported example-data performance claims. Clarified that BirdTree retrieval is user-directed, added attribution instructions, and marked the inverse-Wishart option as historical reproduction only. Added a NEWS entry, updated the site retirement controls, and refreshed the audit ledger to record the current sole pigauto lane and the latest local results.

The fresh offline pkgdown output and route crawler now pass on this candidate. Rendered-content assertions cover the three edited reader pages and verify that the retired simulation-study routes are absent from output, sitemap, and search.

## 3a. Decisions and Rejected Alternatives

Kept the bundled trees as cited examples while describing how users can obtain their own BirdTree sample. Did not claim a general redistribution licence: the official BirdTree pages reviewed require citation and describe research use, but do not establish the CRAN redistribution warranty for bundled objects. Preserved the retired study source and evidence in Git while removing its generated routes and discovery entries.

## 4. Files Touched

- `README.md`
- `NEWS.md`
- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-chrome-refresh.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-ledger-reconciliation.md`
- `docs/dev-log/cran-0.11-audit/defaults-inventory.md`
- `docs/dev-log/cran-0.11-audit/provenance/rights-and-policy.md`
- `docs/dev-log/cran-0.11-audit/site-review-2026-10-08.md`
- `docs/dev-log/cran-0.11-audit/provenance/site-crawler-simulation-retirement-control.txt`
- `inst/NOTICE`
- `pkgdown/clean-internal-pages.R`
- `script/cran-0.11-integration/README.md`
- `script/cran-0.11-integration/check-adapters.R`
- `script/cran-0.11-site/build-offline.R`
- `script/cran-0.11-site/crawl.py`
- `script/cran-0.11-site/seed_pkgdown_cache.py`
- `tests/testthat/test-cran-audit-defaults.R`
- `vignettes/getting-started.Rmd`
- `vignettes/multiple-imputation.Rmd`
- `vignettes/tree-uncertainty.Rmd`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-lane-clarification.md`

## 5. Checks Run

- `UNLAZY_APPROVAL_DIR=/private/tmp/pigauto-cran011-unlazy-approvals node .../gate-check.mjs --reverify --approve --timeout 600 docs/dev-log/cran-0.11-audit/GATES.md`: G5 and G5b passed. The crawler reported 62 HTML pages, 3,441 local references, 34 retired routes, and zero errors. The search index has 613 entries and 57 unique non-empty paths.
- Rendered HTML assertions: `RENDERED_CONTENT_AND_RETIREMENT_OK`.
- `Rscript --vanilla -e 'pkgdown::check_pkgdown()'`: no problems found.
- `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript --vanilla -e 'devtools::test(filter = "cran-audit-defaults", stop_on_failure = TRUE)'`: 198 passed, 0 failures, 0 warnings, 0 skips.
- `git diff --check` and `git diff --cached --check`: passed.
- Unlazy status: 11 gates, 5 met and 6 unmet (G0, G4, G6, G7, G8, G9).

## 6. Tests of the Tests

The site-crawler retirement control previously planted a retired Markdown route and confirmed the crawler failed as expected; its retained receipt is `docs/dev-log/cran-0.11-audit/provenance/site-crawler-simulation-retirement-control.txt`. The current rendered-content checks inspect generated HTML and discovery files rather than source text alone.

## 7a. Issue Ledger

- Fixed: Getting Started described all four example traits as lognormal and highly heritable without support.
- Fixed: BirdTree retrieval and attribution instructions were absent from the main reader paths.
- Fixed: The inverse-Wishart option was not clearly limited to reproducing an earlier analysis.
- Fixed: the retired simulation-study page could reappear in generated output and discovery files.
- Open: CRAN redistribution-warranty basis for bundled BirdTree-derived objects.
- Open: local visual inspection, deployment verification, final exact tarball checks, and independent artifact review.

## 8. Consistency Audit

Reviewed README, NOTICE, Getting Started, tree uncertainty, multiple imputation, generated site pages, retirement manifest, crawler behavior, sitemap, and search. The updated statements agree across the checked local reader pages. Current site evidence is local and structural; it does not certify the deployed site. The live site must be checked after the source is merged and deployed.

## 9. What Did Not Go Smoothly

The offline build needed SRI-verified cached JavaScript and temporary font and network overrides. The first Unlazy reverify used its 120-second default and stopped during sitemap generation; the subsequent approved run used 600 seconds and passed. The report generator resolved a relative output path under the brain vault, so that erroneous draft was discarded and this repo-specific report was written directly. Independent reviewers caught one overstated NEWS attribution and one stale saved-site crawler result; both were corrected, and the latest site build and route check were rerun.

## 10. Known Residuals

G0, G4, G6, G7, G8, and G9 remain unmet. The exact release artifact is still the earlier candidate and does not include this follow-up. No merge, deployment, or CRAN submission occurred. The local site build disables the homepage sidebar and CRAN date annotations and uses system fonts, so it is structural evidence rather than a visual review of production styling.

## 11. Team Learning

Memory receipt: loaded the pigauto LOAD-FIRST manifest through `route.py`; trust recovery-to-truth, compare against main, audit defaults and prediction paths, and preserve `r_cal = 0` shaped the work. Retrieved the existing CRAN audit record from memory. No Golden Set regression class was in scope; this was a documentation, route-retirement, and audit-ledger reconciliation.

Golden Set: not in scope for this reader-surface and ledger slice.

## 12. Cross-Product Coverage

Covers the README and three rendered reader articles, NOTICE, local pkgdown output, route retirement, sitemap and search checks, and the focused defaults suite.

Does NOT cover live visual rendering, deployed-site behavior after merge, redistribution rights for the bundled tree data, every optional-backend model class, the final frozen tarball, platform checks on that artifact, or CRAN acceptance.
