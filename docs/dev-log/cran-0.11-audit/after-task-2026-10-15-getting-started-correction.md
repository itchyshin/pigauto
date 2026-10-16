# After-task report: getting-started dimensions and species matching

## 1. Goal

Correct reader-facing claims in Getting Started so the example output, object dimensions, and species-matching description agree with pigauto 0.11.0 behavior.

## 2. Implemented

- Clarified that AVONET includes seven original traits, while categorical expansion yields 14 latent columns.
- Removed hard-coded fit and prediction examples that described only four continuous traits.
- Corrected prediction-object descriptions and explained that discrete `se` values are uncertainty scores.
- Corrected species matching: input-only species remain in `completed` but are not modeled or filled; tree-only tips may contribute to fitting.

## 3a. Decisions and Rejected Alternatives

- Replaced stale numeric fit output with a structural description rather than presenting seed-dependent values as fixed output.
- Kept the change limited to the Getting Started source and its tracked purl output.
- Assumptions made without asking: none.

## 4. Files Touched

- `vignettes/getting-started.Rmd`
- `vignettes/getting-started.R`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-15-getting-started-correction.md`

## 5. Checks Run

- `git diff --check`: passed.
- `knitr::knit("vignettes/getting-started.Rmd", ...)`: passed; all 55 chunks processed.
- `pkgdown::build_article("getting-started", pkg = ..., new_process = TRUE)`: passed and wrote the article HTML into a fresh site snapshot.
- Rendered-page assertions by text search: found seven original traits, 300 × 14 latent dimensions, the separate prediction dimensions, and corrected unmatched-input wording; no old 300 × 4 or dropped-species wording remained.
- Full `pkgdown::build_site()`: stopped before article generation when pkgdown tried to request CRAN metadata and write to a protected shared cache. This is not a full site-build pass.
- Chrome local-file open: rejected by browser URL policy. No workaround was attempted; visual review remains open.

## 6. Tests of the Tests

The rendered article includes computed preprocessing output showing seven traits and 14 latent columns. The source and rendered HTML search checks would detect the prior four-trait matrix claim and dropped-species claim. No negative-control edit was introduced because this is a documentation-only correction.

## 7a. Issue Ledger

- Fixed: Getting Started described four traits and a 300 × 4 matrix although the example contains seven traits and 14 latent columns.
- Fixed: stale prediction output omitted categorical and ordinal traits.
- Fixed: species-matching prose incorrectly said unmatched species are dropped.
- Open: full fresh site build and visual review.
- Open: merge, deployment, and post-deployment verification.
- Open: exact-artifact checks must be rerun after final source is merged.

## 8. Consistency Audit

Checked the source article, tracked purl output, and rendered Getting Started article. The baseline is described as using 14 latent columns; the prediction object is described separately as seven original-trait columns and a 14-column latent prediction. The unmatched-input description now matches the implementation behavior. Broader website routes and deployment were not rechecked in this slice.

## 9. What Did Not Go Smoothly

The fresh site build attempted a network request for CRAN metadata and encountered protected cache writes. Chrome refused the local file URL. The individual article renderer completed successfully, but it does not replace the full site or visual checks.

## 10. Known Residuals

- This branch is not merged and the correction is not deployed.
- The previously frozen tarball predates this documentation change and is not the final candidate.
- Windows, network-capable `--as-cran`, exact-artifact independent review, and redistribution-rights gates remain separate open release gates.
- No claim is made that CRAN submission is ready.

## 11. Team Learning

The rendered example itself is a useful check on prose dimensions: compare original trait count, latent expansion, and returned prediction shapes separately. An audit of the live page can expose inconsistencies that source-level package tests do not cover.

## 12. Cross-Product Coverage

- Covers: Getting Started prose, tracked purl output, and direct rendered article text for mixed-type trait dimensions and species matching.
- Does NOT cover: all package defaults, every public page, full site navigation/search/sitemap, deployed rendering, platform checks, exact-artifact release checks, or statistical validity of downstream analyses.

### Memory receipt

Loaded the pigauto LOAD-FIRST manifest through `route.py pigauto`, the repository AGENTS instructions, and the after-task protocol. The relevant guard was to verify neighboring reader surfaces and distinguish original trait columns from latent columns. No cross-project scouting was done. The Golden Set was not run because no code path changed.
