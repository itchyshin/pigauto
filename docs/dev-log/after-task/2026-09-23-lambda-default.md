# After-task: Pagel's lambda in the joint baseline, estimated by default (brain D-278)

Lane feat/joint-lambda-default, worktree `pigauto-lambda-default`, draft PR #187 to main (not merged).
Platform Claude Code. Session ran 2026-09-22 16:15 to 2026-09-23 about 15:30 MDT, with one overnight pause
at the rate limit.

## 1. Goal

Make pigauto's in-house joint baseline estimate Pagel's lambda (a residual variance beside the
phylogenetic one) and flip `lambda_mode`'s default to `"estimate"`, covariates included. Evidence: the 18
core simulation cells against the committed #184 yardstick; draft PR separate from #184.

## 2. Implemented

- `fit_mvn_bm_inhouse()` and `fit_joint_solver()` take `lambda` ("fixed_1", "estimate", a numeric scalar
  or a per-column vector) and `lambda_cols`. Each continuous-family column gets its own lambda by profile
  REML (reusing `build_pagel_nll_cache()`); prediction uses Hadfield and Nakagawa's sparse precision on the
  Pagel-transformed tree, with the GLS mean subtracted first so the predictor matches the estimator. A
  block lambda feeds the Sigma M-step and the opt-in exact and refine paths. One eigendecomposition per
  column.
- `fit_baseline()`: default `"estimate"`; the joint and threshold-joint paths keep continuous columns
  under "estimate" (only "cv" and "bayes" still force the per-column path); new `lambda_fixed` argument
  for rebuilds; returns `lambda_per_trait`, `lambda_block`, `lambda_mode`. Discrete columns (binary,
  categorical, zi gate, ordinal) stay at lambda = 1 on every path.
- `bm_impute_col_with_cov()` estimates lambda with GLS coefficients inside the profile likelihood.
- `impute()`, `fit_pigauto()`, `multi_impute()` (new argument) and `multi_impute_trees()` default to
  "estimate"; `model_config` stores the three lambda fields; `multi_impute_trees()` records lambda per tree.
- `joint_solver = "rphylopars"` now asks phylopars for `model = "lambda"` and falls back to the in-house
  solver when phylopars returns tip predictions more than 10 times beyond the observed range.
- NEWS 0.11.0.9000 entry, roxygen, DESCRIPTION dev version.
- Benchmarks: `script/bench_lambda_datasets.R` (real data) and a Totoro launcher kept outside the repo.

## 3a. Decisions and Rejected Alternatives

- Per-trait lambda on continuous columns, block lambda only where a common R is required. Rejected: one
  shared lambda (a weak-signal trait would shrink a strong one) and per-trait lambda inside the Kronecker
  solve (no Kronecker form; needs a two-matrix REML, a spec non-goal).
- Discrete traits do not inherit lambda (Rose's plan review): the inherit would have broken the committed
  per-type guards and repeated a 19-point Trophic.Level loss recorded in NEWS.
- Recovery tests gate bias, not mean absolute error: at n = 300 the absolute error is about 0.13 by
  construction.
- Gate G12 failed at lambda 0.3 (17 to 28% of the gap closed, against 50% planned). Shinichi decided to
  keep "estimate" as the default (2026-09-23), since it was never worse than lambda = 1 in any simulation
  cell and never worse by more than 0.9% in 13 real-data cases.
- GNN arm capped at 100 seeds, and at 50 seeds for n = 1000 after another Totoro user took about 159 cores
  (Shinichi's call). The no-GNN arms, which carry the gate, kept 200 seeds.
- Kohaku (GPU) was offered mid-run; R was still being installed there, so the work stayed on Totoro.

## 4. Files Touched

R: `R/joint_mvn_solver.R`, `R/fit_baseline.R`, `R/bm_internal.R`, `R/pagel_lambda.R`,
`R/joint_mvn_baseline.R`, `R/joint_threshold_baseline.R`, `R/ovr_categorical.R` (comment only),
`R/fit_pigauto.R`, `R/impute.R`, `R/multi_impute.R`, `R/multi_impute_trees.R`.
Tests (new): `test-joint-lambda.R`, `test-lambda-covariates.R`, `test-lambda-dispatch.R`,
`test-lambda-default.R`, `fixtures/lambda_fixed1_reference_ab02e31.rds`.
Tests (edited): `test-lambda-per-type.R` (test 4 only), `test-pagel-lambda.R`, `test-joint-baseline.R`,
`test-joint-threshold-baseline.R` (pinned to fixed_1 where they compare against a lambda = 1 kernel).
Docs: `man/fit_baseline.Rd`, `man/fit_pigauto.Rd`, `man/impute.Rd`, `man/multi_impute.Rd`, `NEWS.md`,
`DESCRIPTION`. Scripts: `script/bench_lambda_datasets.R`, `dev/lambda_avonet_smoke.R`.
Evidence (force-added past the `docs/` ignore rule): `docs/dev-log/2026-09-22-joint-lambda-alignment.md`,
`docs/dev-log/lambda-default/*` (slice reports, benchmark, real data, reviews, logs, comparison scripts),
this file. Nothing under `script/campaign_*`, `BACE/` or `docs/dev-log/arc/` changed.

## 5. Checks Run

| check | result |
|---|---|
| full test suite, final build | FAIL 0, PASS 2596 |
| `rcmdcheck --as-cran`, final build | 0 errors, 0 warnings, 1 note (dev version) |
| AVONET300 smoke | +0.18% mean z-RMSE, lambda 0.99 to 0.995 |
| acceptance ledger, code leaves S2 to S6 | 14 of 14 gates met (run with a Python stand-in; see section 9) |
| pre-run cell, paired on 20 seeds | continuous z-RMSE -0.015 to -0.08 in every arm; discrete unchanged |
| wave 1, no-GNN arms, 18 cells x 200 seeds | in-house arm better in all 9 cells; details in `benchmark.md` |
| wave 2, GNN arm, 18 cells x 100 seeds (50 at n = 1000) | better in all 9 cells; lambda 0.7 closes 61 to 71% of the gap, lambda 0.3 13 to 29% |
| real data, 13 cases x 5 seeds x 4 modes, final build | 12 of 13 better or tied, worst +0.9% |

## 6. Tests of the Tests

- The fixed_1 reference was generated from a `git archive` of origin/main before any edit, so the 1e-12
  comparison cannot be self-consistent; Rose confirmed it matches on 8 setups.
- The recovery file carries a negative control: the estimator must separate a lambda = 1 DGP from 0.3.
- Rose's final review deliberately disabled the application of the estimated lambda and found only 2
  tests failed; the new `[lambda-default] estimated lambda is applied` test and the real `predict()`
  rebuild test were added so the default path itself catches that.
- The Rphylopars guard test stubs phylopars with explosive but finite output and asserts the fallback.
- The G5 oracle checks that the three discrete guards in `test-lambda-per-type.R` are byte-identical to
  origin/main, not just that the file passes.

## 7a. Issue Ledger

Fixed in this lane: ordinal columns picked up the lambda setting (Rose final review); a fully observed
continuous column reported lambda = 1 because the per-column kernel returns before estimating (found while
auditing section 12; predictions were unaffected); partial
`lambda_fixed` errored on joint paths; two eigendecompositions per column; phylopars lambda blow-ups;
stale roxygen and NEWS gaps; overclaims in the reports (discrete "identical", coverage range). On CI (2026-09-23 and 24):
the PMM acceptance test asserted the pooled median of 10 draws was a donor value, which fails when the
draws spread under the estimated lambda; it now asserts every draw is a donor value. The clamp no-op test
compared two separately trained fits bit for bit and failed on the macOS runner, where seeded training was
not repeatable; it now compares two predictions from one fit. Pinning torch to one thread was tried first
and dropped: torch ignores the setting once parallel work has started in the process.
Deferred: `suggest_next_observation()` still scores under lambda = 1 (NEWS known limitation); downward
bias of the lambda estimate at weak signal.

## 8. Consistency Audit

Swept every `lambda_mode` entry point (`fit_baseline`, `fit_pigauto`, `impute`, `multi_impute`,
`multi_impute_trees` fallback) and every joint delegate (joint MVN, threshold-joint, OVR, covariate
branch). `cross_validate()`, `compare_methods()` and `simulate_benchmark()` inherit the default through
`impute()` (named in NEWS). The simulation's arm library re-reads `model_config$lambda_mode`, which is now
stored. The one surface not moved is `suggest_next_observation()`.

## 9. What Did Not Go Smoothly

- The session ended mid-batch overnight (rate limit); two builders were resumed from saved transcripts.
- Totoro was shared from mid-morning (another user near 159 cores, load up to 665), which cut the GNN
  wave's share by half and forced the seed cap.
- My own `pkill -f` patterns twice matched the SSH session that ran them and killed it.
- The aggregation script takes one results folder; a two-folder call silently re-aggregated the wrong wave
  under a misleading name (renamed, nothing lost).
- Node.js disappeared from the Mac between days, so `gate-check.mjs` could not run; the ledger was
  re-verified with a Python stand-in applying the same rule (exit 0 and the expected text in output).
- Two of my own gate oracles were malformed (an EXPECT spanning two lines; an R escape error).

## 10. Known Residuals

- At true lambda 0.3 the default closes only 17 to 28% of the gap to the frequentist lambda arm; the
  per-column lambda estimate is biased low by about 0.05 at lambda 0.3 and 0.7 (n = 300).
- LepTraits flight duration at 2,000 species is worse under the estimate than at lambda = 1.
- The GNN arm has 100 seeds (50 at n = 1000), not 200.
- Rphylopars under lambda is about 26 times slower than under BM.
- Fit time at large n grows with the dense eigendecomposition (6.0 s at n = 3000 for three traits).
- The acceptance ledger was not re-run by the official Node checker.

## 11. Team Learning

- A gate that pins old behaviour (the per-type test 4) must be amended in the open, with the diff scoped
  and an oracle proving the rest of the file is untouched.
- On a shared Totoro, check other users' load before projecting wall time, not only free memory.
- Never `pkill -f` a pattern that also appears in the command running it.

## 12. Cross-Product Coverage

The default flip is cross-cutting: it changes every fit that does not pass `lambda_mode`.

Covers: continuous, count, proportion and zi-magnitude columns on the joint MVN, threshold-joint and
per-column paths; `multi_proportion` compositions (per-component lambda on the per-column path, measured
at the lower bound on an iid test composition); covariate fits; single-obs and multi-obs data; `impute`, `fit_pigauto`, `multi_impute`,
`multi_impute_trees` (lambda per tree); `predict()` (carried baseline); conformal intervals (recomputed
from the new baseline, coverage 0.957 to 0.969 at n >= 300); gnn = TRUE and gnn = FALSE.

Does NOT cover: binary, categorical, zi gate and ordinal traits (deliberately at lambda = 1);
`suggest_next_observation()`; the "cv" and "bayes" modes on the joint path (per-column only, and not with
covariates); separate phylogenetic and residual covariance matrices across traits; OU, EB, kappa and delta
models; GPU runs.
