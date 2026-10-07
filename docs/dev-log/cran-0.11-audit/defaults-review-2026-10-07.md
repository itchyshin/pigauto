# Effective-default audit: 2026-10-07

Read-only review at `bb5835d1b214d783da7b0b414df99aa6ba926bc7`. No tests, builds, benchmarks, or source edits were part of the audit. The existing `defaults-inventory.md` lists all 33 `NAMESPACE` exports and their signatures; an independent parse of source formals found the same 33. The registered `predict.pigauto_fit` S3 method is not in that export list and was checked separately.

## Confirmed defaults and routes

The `0.11.0` defaults remain as documented for trait encodings, baseline selection, `gnn = FALSE`, `safety_floor = FALSE`, `phylo_signal_gate = FALSE`, posterior MI eligibility, and `r_cal = 0`. Intentional exceptions include benchmark wrappers that turn the GNN on. `multi_impute(draws_method = "auto")` uses posterior draws only for all-continuous, single-observation data without covariates or compositions; other inputs receive explicitly diagnostic conformal draws that the analysis-aware pooling route refuses. Saved-fit reconstruction retains legacy fallbacks. No scientific default change is indicated.

## Mismatches requiring correction

1. `fit_pigauto()`'s formal `k_eigen` default is `"auto"`, but its roxygen argument text says integer/default 8. Correct the help text to describe `"auto"` and its tree-size resolution.
2. `impute(..., k_eigen = 3L)` is accepted through `...` and documented as reaching `fit_pigauto()`, but `impute()` has already built a graph with `k_eigen = "auto"`; `fit_pigauto()` then replaces the supplied value with the graph coordinate count. Add an explicit `k_eigen = "auto"` formal to `impute()`, pass it to `build_phylo_graph()`, document it, and test that the selected dimension reaches the fitted configuration. This does not change the default.
3. The `fit_pigauto()` description calls the supported trait set “all five trait types”; pigauto supports eight. Correct this adjacent factual error in the same help block.
4. `defaults-inventory.md` already labels the 2026-10-02 audit as a dated historical contradiction. Preserve that dated evidence and ensure it is never linked as a current defaults statement.

## Proposed verification

- `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript --vanilla script/cran-0.11-defaults/run.R`; the prior recorded runtime was 3.9 seconds, so allow 1 minute. It checks 126 existing assertions but does not cover the `impute()` `k_eigen` override or roxygen text.
- `Rscript --vanilla -e 'devtools::test(filter = "cran-audit-defaults", stop_on_failure = TRUE)'`; allow 1 minute after adding behavior and help assertions.
- `Rscript --vanilla -e 'devtools::document()'`; allow 1 minute to regenerate help after roxygen changes.
- If routing changes extend beyond the graph argument, run `Rscript --vanilla -e 'devtools::test(filter = "cran-audit-defaults|route-choice|lambda-default|discrete-lambda|multi-impute|multi-impute-trees|fit-predict", stop_on_failure = TRUE)'`; allow 5 minutes for bounded unit tests.

All are proposals, not evidence of passing runs. The full audit does not justify repeating historical comparison campaigns.

## Implementation addendum: 2026-10-07

The `impute()` wrapper now declares `k_eigen = "auto"` and passes it to `build_phylo_graph()`. `fit_pigauto()` help now reports the auto default and eight supported trait types. The targeted defaults batch passed 127 assertions, including explicit `k_eigen = 2L`; the defaults runner passed 127 assertions and printed `CRAN_DEFAULTS_AUDIT_OK`. Roxygen regenerated both affected help files.

The adjacent route/fit/predict suite then completed in 53.7 seconds: 860 passed, 0 failed, 27 expected small-validation warnings, and one smcfcs integration skip because that package is installed. The single explicit `k_eigen` test passed.

The complete source test suite passed 3,464 assertions with zero failures, 182 warnings, and 8 environment skips in 335.7 seconds. The warnings are predominantly intentional small-validation and conformal-ceiling diagnostics from compact fixtures. A local `devtools::check()` built and tested the package archive in 7m17.5s: 0 errors, 0 warnings, 1 NOTE (the host could not verify remote system time); the embedded test run passed 3,318, warned 182 times, and skipped 32.

## Independent coverage review of candidate source: 2026-10-07

Read-only review at candidate head `3788c64`. The reviewer confirmed that the inventory records all 33 exports and `predict.pigauto_fit()` separately, and that the bounded tests cover the eight trait types, mixed observed-value preservation, named overrides, selected MI routes, saved-current-fit round-trip, and wrapper forwarding. The review classified the overall defaults gate as **partial**: the inventory is a manual source audit and the executable suite checks selected release-critical defaults rather than every formal and every route; automatic analysis-aware posterior MI and the complete downstream tree-analysis workflow are not run; the conformal-MI fallback message and legacy saved-GNN reconstruction are not asserted; installed-library adapter evidence does not establish the full defaults matrix; and release-artifact validation is separate.

The reviewer also identified that the mixed default test checked `gnn = FALSE` in stored configuration but lacked an explicit assertion on the effective calibrated GNN weight. Added assertions that `result$fit$r_cal_gnn` is nonempty and exactly zero. Added `tests/testthat/fixtures/cran-audit-exported-formals.tsv`, a reviewed snapshot of all 286 formal names/default expressions for the 33 exports and registered `predict.pigauto_fit()` method, with an assertion that compares every row to the loaded namespace. The first run found a TSV quoting defect in the fixture generator; after correcting it, the comparison passed. A bounded real two-draw automatic conformal fallback now checks the user-facing warning, diagnostic object, and refusal by `with_imputations()`.

The independent reviewer confirmed the fixture meaningfully closes the manual-inventory gap for formal names/defaults and withdrew an initial concern about the `...` sentinel after observing the runnable test. The final focused command `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript --vanilla -e 'devtools::test(filter = "cran-audit-defaults", stop_on_failure = TRUE)'` passed 169 assertions, 0 failures, 0 warnings, and 0 skips in 4.9 seconds. The same test passed against pigauto installed from clean commit `d2c2736` in an isolated library without `drmTMB` or `gllvmTMB`: 169 passed, 0 failed, 0 warnings, 0 skips. Installation output is `provenance/installed-defaults-install-2026-10-07.log` (SHA-256 `5860f89a589bb6fc84e3be923188690550c0f67bb5f51c18228c3255925e3bef`); the installed test output is `provenance/installed-defaults-run-2026-10-07.log` (SHA-256 `0e32efebb93d4eb0666dad1c276cc807ce0ebb4e7501caa144da53e48a8c2dfa`). The overall G1 gate remains partial because the full automatic posterior and tree-analysis routes, legacy saved-GNN reconstruction, and exact frozen artifact are not covered by these checks.

## Independent evidence follow-up

The reviewer rechecked the installed evidence and found no factual mismatch. Both log hashes match; the test log identifies source commit `d2c2736`, R 4.6.0, the isolated library, and the absence of both optional backends. The install log itself does not include its source path, command, or commit, so those provenance details rely on this report and the test log. This evidence is not frozen-artifact validation.
