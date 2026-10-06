# Lane plan: estimate Pagel's lambda in pigauto's joint baseline, and make it the default

Ordered by Shinichi 2026-09-22 (brain D-278): *"open a lane for spec decision 3 with lambda estimated as
the new default, covariates included; use the core slice as the benchmark; separate PR after #184."*
Parent spec: `specs/2026-05-18-pagel-lambda-baseline-design.md` (APPROVED 2026-05-18), decisions 1 and 3.
This file is the plan stub; the lane opens after PR #184 merges, in its own worktree and branch.

## Why

In the four-arm simulation every pigauto arm assumed lambda = 1 in the path mixed-type data takes:
`R/joint_mvn_solver.R` has no residual (non-phylogenetic) component, and `R/fit_baseline.R:219-224`
forces the per-column path whenever `lambda_mode` is not `"fixed_1"`, dropping the joint baseline, the
threshold-joint discrete path and covariates. The two arms that estimate the phylo / non-phylo split
(BACE; Rphylopars at `model = "lambda"`) were the only ones below the mean floor at lambda 0.3.

## Scope (files)

- `R/joint_mvn_solver.R`: add a residual variance beside the phylogenetic one in `fit_mvn_bm_inhouse()`;
  profile or ML over one lambda per trait block via the tree transform `R(lambda) = lambda R + (1 - lambda) I`
  (equivalently scaling internal branch lengths), reusing the existing sparse Henderson precision (O(n)).
  Store `lambda_per_trait` on the fit (spec 4.5).
- `R/fit_baseline.R`: `lambda_mode = "estimate"` no longer forces per-column; threshold-joint and OVR
  inherit (spec 8.4). Default flips to `"estimate"`; `"fixed_1"` stays available and bit-identical.
- `R/bm_internal.R` covariate-aware BM path: accept `lambda` (today it has none and silently fits at 1).
- `R/fit_pigauto.R`, `R/impute.R`, `R/multi_impute.R`: default `lambda_mode = "estimate"`; roxygen; NEWS.
- Out of scope: OU/EB/kappa/delta; multivariate lambda per off-diagonal; any change to the GNN.

## Tests (ship with the implementation)

- Recovery: simulated BM data with lambda in {0.3, 0.7, 1}, n = 300, 20 seeds: estimated lambda within
  0.05 of truth for each trait block.
- Back-compat: `lambda_mode = "fixed_1"` reproduces v0.11 outputs to 1e-8 on `tests/testthat` fixtures.
- Dispatcher: with lambda estimated, the joint path is taken (not per-column) when >= 2 BM-eligible columns.
- Covariates: lambda estimated with covariates present; no silent fallback to 1 (regression test for the
  2026-08-08 "Pagel lambda dropped with covariates" finding).
- Discrete: threshold-joint binary/ordinal and OVR categorical run under the estimated lambda.

## Benchmark (the acceptance gate)

Re-run the 18 core cells of the simulation for the pigauto arms only (`script/campaign_sim_totoro.sh core`
with `ARMS=gnn_off,gnn_off_rphylopars,gnn_on`, seeded resume-skip against a fresh OUT_DIR), aggregate
with `script/campaign_gnn_off_aggregate.R`, and compare to the committed `script/campaign_sim_results/summary.csv`.

- Target: at lambda 0.3, n = 100, pigauto GNN off z-RMSE from 1.031 toward `freq_lambda`'s 0.913 (close
  most of the gap); at lambda 0.7 from 0.916 toward 0.801.
- Guard: lambda = 1 cells move by less than 0.01 (0.487 / 0.426 / 0.379 at n = 100 / 300 / 1000).
- Guard: conformal coverage at n >= 300 stays at 0.95 to 0.96; the spec's AVONET check (< 5% RMSE change).

## Estimate

1 to 2 days of package work (D-139: state before running) plus about 4 h on Totoro for the benchmark.
Pre-run: one core cell (lambda 0.3, n = 300, 20 seeds) before the 18-cell run.
