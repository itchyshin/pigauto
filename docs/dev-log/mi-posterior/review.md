# S5 review: posterior multiple imputation (2026-09-24)

## How the review ran

The review was fresh and adversarial, run as a workflow.
- Four lenses read the code at commit 5b0d6f9, which has the same R code as the campaign commit 69670d4:
  1. sampler mathematics;
  2. diagnostics and outputs;
  3. API, provenance and claims;
  4. test adequacy.
- Each finding went to two independent skeptics:
  - one traced it through the code;
  - one checked it against the spec and its practical impact.

  Each skeptic also said whether the finding changes the numbers of the running campaign.
- A fixer resolved the findings, a fresh checker re-tested them, and one repair round followed.

Raw verdicts are in the workflow journal
`subagents/workflows/wf_7ae0bdc8-8f1/journal.jsonl` (session directory). This file is the durable record.

## Result

- **Findings:** 28 in total.
  - 17 confirmed by both skeptics;
  - 8 upheld by one;
  - 3 refuted.
- **No error in the sampler mathematics.** The math lens re-derived every full conditional, the parameter
  expansion, the collapsed Metropolis moves (targets, Jacobians, prior terms) and the adaptation schedule.
  It raised no finding against them.
- **No finding changes the campaign numbers.** Every skeptic answered "affects campaign: false".
- **The algorithm is unchanged by the fixes.** Proof:
  - 33 of 36 functions in `R/mi_posterior.R` parse identically to 69670d4;
  - on a fixture, `.mip_fit()` output is bitwise identical to 69670d4;
  - the new O(n) tip depths are bit-identical to the old ones on 300 trees.

## Findings and dispositions

| Id | Verdict | Severity | Finding | Resolution (commit) |
|---|---|---|---|---|
| api#0 / math#0 | confirmed / plausible | required | `param_uncertainty = "none"` draws carried the proper provenance marker and could be pooled | New marker `pigauto_posterior_plugin_diagnostic`, refused by `with_imputations()` and `pool_mi()`; print says unsupported (dcc1174, tests a7a24f4) |
| api#1 | confirmed | required | NEWS/roxygen gave MC-dropout the conformal bias numbers | Stated per method: conformal -0.20 to -0.46 with 0 to 17% coverage; MC-dropout -0.03 to -0.38 (dcc1174) |
| api#2 | confirmed | required | Congeniality scope too broad | Narrowed to analyses linear in the imputed traits on the imputation scale; nonlinear terms, interactions, raw-scale analysis of logged traits and external covariates named as not covered (dcc1174; design.md 2.5) |
| diag#0 / tests#1 | confirmed | required | K = 1 crashed after the MCMC (diag of a scalar) | Diagonal indexing with cbind; K = 1 test (dcc1174, a7a24f4) |
| tests#0 | confirmed | required | No test could detect a wrong posterior | Kernel-agreement test: Metropolis-only vs Gibbs-only kernels within 3.5 combined MCSE. Catches 8 of 8 sampler mutants (a7a24f4) |
| math#1 | confirmed | minor | `build_henderson_S_inv()` formed a dense n x n vcv to read tip depths | Optional `tip_depths` argument; default path unchanged; posterior path passes O(n) depths (dcc1174) |
| math#3 / api#6 | confirmed | minor | Priors described as "MCMCglmm defaults" | Priors written out explicitly (dcc1174) |
| math#2 | plausible | minor | `keep_draws` described as a minimum | Now described as a target; behaviour unchanged (dcc1174) |
| api#3 / diag#1 | confirmed | minor | `species_col` error claimed duplicates when there were none | Message split on actual duplicates (dcc1174) |
| diag#2 / tests#8 | confirmed | minor | Fully observed input crashed after the MCMC | Clear error before any MCMC, counting only rows of `traits` (dcc1174, 20213e2) |
| api#7 | confirmed | minor | Error messages pointed users of the new path the wrong way | Messages point to `with_imputations()` / `draws_method = "posterior"` (dcc1174) |
| diag#4 / tests#7 | confirmed | minor | Folded R-hat untested | Scale-only test that only the folded half detects (a7a24f4) |
| tests#2 | confirmed | minor | Item-6 row-slicing test could not fail | Replaced by a one-sweep probe test that catches row-slicing mutants (a7a24f4) |
| tests#5 | confirmed | minor | No test that intervals use all kept draws or that the m datasets span chains | `draw_index` recorded; test checks both (a7a24f4) |
| tests#6 | plausible | minor | Decode check on the log scale only | Check on the original scale (a7a24f4) |
| diag#5 | plausible | minor | `n_chains` docs | Clarified (dcc1174) |
| diag#3 | plausible | minor | Same as api#0 | Resolved with api#0 |
| tests#3 | plausible | minor | Item-2 test partly tautological | Left: not a one-line fix; the kernel test catches the bare-Q mutant it targeted |
| api#4 | plausible | minor | ESS reported as Inf rather than NA when all ESS are NA | Left: cosmetic |
| api#5, api#8, tests#4 | refuted | minor | (refuted by both skeptics) | None |

One gap was found by the checker and fixed in the repair round. The "no missing cells" guard counted
padded tree tips that are not in `traits` (20213e2).

## Verification after the fixes (checker plus fixer)

- `test-mi-posterior.R`: 14 tests, 70 expectations, 0 failures.
- `test-multi-impute-posterior.R`: 15 tests, 130 expectations, 0 failures.
- Related suites (multi-impute, with-imputations, pool-mi, mi-provenance, mi-pool, henderson): 616
  expectations, 0 failures, 1 skip (pre-existing by design).
- `gate_exactness.R`: `EXACTNESS_OK`.

## Verdict

PROCEED, pending the final gate re-verification (G1, G2, G5a, G5c on the final commit) and the
campaign gates G6 to G8.
