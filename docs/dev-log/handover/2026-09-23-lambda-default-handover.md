# Handover: lambda-default lane (PR #187)

State on 2026-09-23 afternoon MDT. Written for whoever reviews or continues PR #187.

## Where things stand

- Branch `feat/joint-lambda-default`, worktree `/Users/z3437171/Dropbox/Github Local/pigauto-lambda-default`,
  pushed. Draft PR https://github.com/itchyshin/pigauto/pull/187 targets main. Not merged; merging is
  Shinichi's call.
- The implementation is complete. `lambda_mode = "estimate"` is the default, and the joint, covariate and
  per-column paths all estimate Pagel's lambda for continuous-family traits. Discrete traits stay at 1.
- Evidence is complete. Full suite and R CMD check were clean on the final build. The 18-cell simulation
  covers both waves, and the real-data check covers 13 cases. Two fresh reviews by Rose are on file (plan
  review and final review), and every finding from them is fixed or recorded.
- Shinichi made two decisions. He kept "estimate" as the default although gate G12 failed at lambda 0.3.
  He capped the GNN arm at 50 seeds for n = 1000.

## Read first

1. `docs/dev-log/after-task/2026-09-23-lambda-default.md`: what changed, checks, residuals, negative space.
2. `docs/dev-log/lambda-default/benchmark.md` and `real-data.md`: the numbers.
3. `docs/dev-log/2026-09-22-joint-lambda-alignment.md`: design, symbolic alignment, both reviews.

## Open items (none block review)

- The lambda estimate is biased low at weak signal, about -0.05 at true lambda 0.3 and 0.7 with n = 300.
  This is why lambda 0.3 closes only 13 to 29% of the gap. A follow-up lane could try a REML correction or
  a small-n CV choice.
- `suggest_next_observation()` still scores candidates at lambda = 1.
- The GNN arm has 100 seeds (50 at n = 1000). If the numbers are needed at full strength, rerun with
  `WAVE=gnn2 ARMS=gnn_on` from seed 51 once Totoro is quiet.
- Node.js is missing on the Mac, so the unlazy ledger was re-verified with a Python stand-in rather than
  `gate-check.mjs`.

## Compute left behind (all on Totoro, nothing running)

- `~/R/lib-lambda` holds the build used for the simulation waves. `~/R/lib-lambda2` holds the final build,
  used for the real-data rerun.
- `~/pigauto_lambda_src`, `~/pigauto_lambda_src2`: source copies, launchers, logs.
- `~/pigauto_sim/results/core_lambda_*` holds the simulation results. `core_lambda_core_gnn_summary.csv`
  was renamed `MISNAMED_fast_wave_summary.csv`, because it is a re-aggregation of the no-GNN wave under
  the wrong name.
- These can be deleted once the PR is merged. Nothing here is needed to reproduce the PR's tables.

## Resume

```sh
cd "/Users/z3437171/Dropbox/Github Local/pigauto-lambda-default"
git pull
NOT_CRAN=true Rscript -e 'devtools::load_all(); testthat::test_dir("tests/testthat")'
Rscript docs/dev-log/lambda-default/compare_benchmark.R docs/dev-log/lambda-default/core_lambda_fast_agg_summary.csv
```
