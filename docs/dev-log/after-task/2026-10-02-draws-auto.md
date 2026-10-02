# After-task: `multi_impute(draws_method = "auto")` becomes the default

2026-10-02. Branch `feat/draws-auto` (worktree `pigauto-draws-auto`), based on origin/main b565cad. Platform: Claude Code. Third PR of the defaults work (after #199 audit and #200 `gnn = FALSE`), on Shinichi's instruction "now do draws_method = "auto" in another PR".

## 1. Goal

Make the default `multi_impute()` route end in a supported analysis when the data allow it: choose posterior draws for continuous, one-row-per-species data without covariates, and otherwise fall back to conformal draws with a one-line message saying why and that they cannot be pooled.

## 2. Implemented

- `draws_method = c("auto", "conformal", "mc_dropout", "posterior")`; `"auto"` is the default.
- Internal `.mi_resolve_draws_auto()` in `R/multi_impute.R`. It applies the posterior route's input contract (`.multi_impute_posterior()` in `R/mi_posterior.R`) before any fitting: covariates, `species_col` or `multi_proportion_groups` given, or any trait type other than continuous after `preprocess_traits()`, give `"conformal"`; otherwise `"posterior"`. Missing-cell and tree checks stay with the chosen route so its own errors apply.
- Messages: the conformal fallback always prints its reason and that `with_imputations()` / `pool_mi()` refuse the draws; the posterior choice prints a line (when `verbose`) noting the sampler can take minutes. `result$draws_method` records the choice.
- Roxygen for `draws_method` documents `"auto"` and how to get the previous default.
- README, getting-started and gnn-architecture vignettes, NEWS.
- `tests/testthat/test-draws-auto.R`: default value; resolver to posterior (and silent when `verbose = FALSE`); fallback for a binary trait, covariates and `species_col` with the reason in the message; default `multi_impute()` on continuous data giving a `pigauto_posterior_mi` that pools through `with_imputations()` / `pool_mi()`; default on mixed data giving conformal draws that `with_imputations()` refuses.

## 3a. Decisions and Rejected Alternatives

- Shinichi: add `"auto"` as the default (proposed in the #199 audit).
- Mine: the fallback message is a `message()`, always shown, not a warning: the call does what was asked, and the message is about what the result can be used for. `multi_impute_trees()` is unchanged: it has no posterior route. No version bump here, so this PR does not depend on #200's 0.11.0.9001; the roxygen says "previously conformal" instead of naming a version. `"auto"` does not choose posterior for count or proportion traits, because the posterior route accepts only `continuous`.

## 4. Files Touched

Modified: `NEWS.md`, `R/multi_impute.R`, `README.md`, `man/multi_impute.Rd`, `vignettes/getting-started.Rmd`, `vignettes/gnn-architecture.Rmd`.
Created: `tests/testthat/test-draws-auto.R`, `docs/dev-log/after-task/2026-10-02-draws-auto.md`.

## 5. Checks Run

- `test-draws-auto.R` with `NOT_CRAN=true`: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 16 ]`.
- Before the full run, every test call to `multi_impute()` without `draws_method` was listed: the three real calls use factor traits or `species_col`, so they resolve to conformal as before (they now also print the fallback message).
- Full suite and `rcmdcheck --as-cran`: see the addendum.
- `slop_check.py` on all added prose lines: 0 findings.

## 6. Tests of the Tests

- The default-value test fails on origin/main (`"conformal"` is first there).
- The fallback tests match the reason text (`diet: binary`, `covariates`, `species_col`), so a resolver that falls back silently, or for the wrong reason, fails.
- The end-to-end test asserts the result class and that pooling succeeds, which the conformal route cannot do.

## 7a. Issue Ledger

- Runtime: the default call on continuous data now runs the MCMC sampler (about ten minutes per fit at 300 species and four traits on one core, branch pre-run timing) instead of a conformal pass. Stated in NEWS and the message; not otherwise mitigated.
- Count and proportion traits are not eligible for `"auto"` → posterior even though they are continuous-like; that follows the posterior route's current contract.

## 8. Consistency Audit

- `with_imputations()` and `pool_mi()` refusal messages already point to `draws_method = "posterior"`; still accurate.
- `multi_impute_trees()` keeps `"mc_dropout"`; its own roxygen is unchanged.
- PR #199's MI article and README section describe `"conformal"` as the current default; whichever of #199 and this PR merges second needs that wording updated (noted in the PR).

## 9. What Did Not Go Smoothly

Nothing material. The resolver calls `preprocess_traits()` a second time (wrapped in `suppressMessages()` so its type-detection messages are not printed twice); the cost is small next to either draws route.

## 10. Known Residuals

- Interaction with #200: under `gnn = FALSE` (the #200 default) the conformal fallback runs without the GNN; nothing in this PR depends on that.
- The posterior route's known limits on an unmerged branch (`arc/rubin-freq-bace`: under-coverage at n = 1000, lambda = 1; sampler failures at n = 100) now apply to more users by default.

## 11. Team Learning

- Make an `"auto"` option reuse the target route's own input contract, and leave the checks that need fitting to that route, so the resolver cannot drift into a second, different definition of eligibility.

## 12. Cross-Product Coverage

Covers: `multi_impute()` default resolution for continuous, binary, covariate and multi-observation inputs; the end-to-end pooled route.
Does NOT cover: `multi_impute_trees()`; count, proportion, ordinal or categorical traits through the posterior route (not implemented); `multi_proportion_groups` beyond the resolver branch; runtime at n > 1000.

## Section 5 addendum: final check lines

- Full suite, `TESTTHAT_MAX_FAILS=Inf NOT_CRAN=true devtools::test()`: FAIL 0, PASS 3104, SKIP 8 (all environmental).
- `rcmdcheck::rcmdcheck(args = c("--as-cran", "--no-manual"))`: 0 errors, 0 warnings, 1 note (development version number; pre-existing GitHub links in `common-pitfalls` returned HTTP 503 during the check).
