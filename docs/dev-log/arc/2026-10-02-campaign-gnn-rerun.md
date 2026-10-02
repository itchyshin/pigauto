# Rerun at current main: pigauto with and without the GNN, and the in-house vs Rphylopars joint solver

2026-10-02, Totoro, pigauto origin/main b565cad (after PR #187 joint lambda default, #191 REML lambda, #192 `predict_method = "auto"`, #195, #189), installed into a private library (`~/defaults-audit/Rlib`, torch 0.17.0 with its own libtorch) so the shared library used by another lane was not touched. Runners unchanged from 2026-09-19 except one line: the derived "GNN on, full baseline" arm now replays the held-out fit's per-trait `predict_method` route, as production `gnn = FALSE` does (`script/campaign_gnn_off_cell.R`). Purpose: re-measure, under the defaults now shipped, the two findings the defaults audit (`docs/dev-log/defaults-audit/2026-10-02-defaults-audit.md`, Decisions 2 and 3) leans on.

## Design

- Same DGPs, masks and metrics as `docs/dev-log/arc/2026-09-19-campaign-gnn-off-results.md`: `bm_mixed`, `ou_mixed`, `bace_dgp` at n in {100, 300, 1000}, AVONET300 at n = 300; 20 seeds per cell (200 cells); 30% MCAR on the user data, the same mask for every arm within a cell; z-RMSE averaged over continuous traits, accuracy over discrete traits, coverage of the production conformal interval (nominal 0.95). Monte Carlo SE over seeds in brackets.
- Default `lambda_mode = "estimate"`, `predict_method = "auto"`, `joint_solver = "inhouse"` unless an arm says otherwise. BACE was not rerun (its 2026-09-19 numbers stand; it was the costliest arm and is not a defaults question).
- Pass 1 (no torch): GNN off (default safety machinery), GNN off pure (`safety_floor = FALSE, phylo_signal_gate = FALSE`), raw `Rphylopars::phylopars(model = "BM")` on continuous columns, mean/mode floor; plus the solver runner (`script/campaign_solver_cell.R`): in-house default, `joint_solver = "rphylopars"` default and pure, raw Rphylopars. 48 concurrent cells x 2 threads; 400 cells, 0 errors, 12 min wall.
- Pass 2: GNN on as shipped (2000 epochs) and the derived GNN on predicting from the full baseline. 36 concurrent cells x 4 threads (144 cores). Merged per cell with `script/campaign_gnn_rerun_merge.R`, aggregated with `script/campaign_gnn_off_aggregate.R`.
- Aggregates: `script/campaign_gnn_rerun_results/` (`agg_gnn.rds`, `agg_solver.rds`, summary CSVs). Raw cell files on Totoro: `~/defaults-audit/results_{off,on,merged,solver}/`.

## Joint solver (pass 1)

Mean z-RMSE over continuous traits, GNN off, default safety machinery, 20 seeds.

| DGP | n | in-house (default) | `joint_solver = "rphylopars"` | paired difference (MCSE) | relative | raw Rphylopars |
|---|---:|---:|---:|---:|---:|---:|
| AVONET300 | 300 | 0.539 (0.056) | 0.440 (0.040) | -0.098 (0.025) | -14% | 0.547 (0.044) |
| bm_mixed | 100 | 0.274 | 0.272 | -0.002 (0.010) | +0.7% | 0.256 |
| bm_mixed | 300 | 0.156 | 0.153 | -0.002 (0.001) | -1.6% | 0.152 |
| bm_mixed | 1000 | 0.097 | 0.090 | -0.007 (0.001) | -7.7% | 0.090 |
| ou_mixed | 100 | 0.424 | 0.452 | +0.028 (0.010) | +6.5% | 0.417 |
| ou_mixed | 300 | 0.225 | 0.220 | -0.005 (0.002) | -2.3% | 0.217 |
| ou_mixed | 1000 | 0.143 | 0.134 | -0.010 (0.001) | -7.2% | 0.132 |
| bace_dgp | 100 | 0.985 | 0.991 | +0.007 (0.006) | +0.7% | 1.212 |
| bace_dgp | 300 | 0.927 | 0.923 | -0.004 (0.004) | -0.5% | 1.176 |
| bace_dgp | 1000 | 0.898 | 0.885 | -0.012 (0.002) | -1.4% | 1.144 |

Mean wall time per fit (s), in-house against Rphylopars solver: 4.8 against 104 on AVONET300; 1.4 against 17, 3.6 against 51, 24 against 174 on `bm_mixed` at n = 100, 300, 1000. Discrete accuracy and interval coverage do not differ between the solvers beyond one MCSE (accuracy 0.770 against 0.772 on AVONET300; coverage 0.97 to 0.98 for both).

What changed since 2026-09-19: the in-house default on AVONET300 moved from 0.790 to 0.539, now level with raw Rphylopars (0.547), so most of the real-data gap reported then has closed under the new lambda and prediction defaults. The Rphylopars solver inside pigauto is still 14% lower on AVONET300 and 7 to 8% lower on the BM and OU simulations at n = 1000, at 7 to 22 times the fit time, and is worse on `ou_mixed` at n = 100.

## With and without the GNN (pass 2)

200 cells, 0 errors, 23 min wall. Paired differences in mean continuous z-RMSE against the default GNN-off arm (negative = better), MCSE over 20 seeds in brackets.

| DGP | n | GNN off (default) | GNN on, as shipped | GNN on, full baseline | GNN off, pure |
|---|---:|---:|---:|---:|---:|
| AVONET300 | 300 | 0.539 | +0.101 (0.011), +23% | -0.018 (0.011), -1.4% | -0.002 (0.001) |
| bm_mixed | 100 | 0.274 | +0.040 (0.008), +15% | +0.003 (0.003) | -0.016 (0.010) |
| bm_mixed | 300 | 0.156 | +0.024 (0.002), +16% | +0.001 (0.001) | -0.000 (0.000) |
| bm_mixed | 1000 | 0.097 | +0.013 (0.002), +13% | +0.000 (0.000) | 0 |
| ou_mixed | 100 | 0.424 | +0.062 (0.014), +17% | +0.011 (0.008) | -0.014 (0.008) |
| ou_mixed | 300 | 0.225 | +0.027 (0.004), +13% | +0.001 (0.000) | -0.002 (0.001) |
| ou_mixed | 1000 | 0.143 | +0.017 (0.002), +12% | +0.000 (0.000) | 0 |
| bace_dgp | 100 | 0.985 | +0.008 (0.004), +0.9% | +0.001 (0.002) | -0.027 (0.006) |
| bace_dgp | 300 | 0.927 | +0.022 (0.007), +2.4% | +0.008 (0.007) | -0.023 (0.006) |
| bace_dgp | 1000 | 0.898 | +0.013 (0.006), +1.5% | -0.001 (0.006) | -0.008 (0.003) |

Discrete accuracy, GNN on as shipped against GNN off: 0.746 against 0.770 on AVONET300, 0.792 against 0.820 on `bm_mixed` at n = 100, equal within MCSE elsewhere. Production-interval coverage is 0.96 to 0.98 for GNN off at n >= 300 and 0.87 to 0.92 at n = 100; GNN on is 0.01 to 0.03 lower throughout. Mean wall time per fit on 4 threads: GNN off 1.4, 3.5 and 24 s at n = 100, 300 and 1000 on `bm_mixed` (4.7 s on AVONET300); GNN off pure 0.4, 1.0 and 4.1 s (1.9 s); GNN on 122, 170 and 496 s (181 s).

1. **As shipped, GNN on is worse than GNN off in every cell**, by 12 to 23% on the BM, OU and AVONET300 data and 1 to 2% on the low-signal DGP. The 2026-09-19 finding holds at current main, and on AVONET300 it is now larger: the GNN-off baseline improved (0.790 to 0.539) and GNN-on did not keep pace (0.775 to 0.640).
2. **Predicting from the full baseline removes the loss and nothing more.** GNN on from the full baseline is within two MCSEs of GNN off in all ten cells (AVONET300 -0.018, MCSE 0.011). At the current defaults the network adds no measurable accuracy on these data, at 20 to 90 times the fit time.
3. **New: the safety machinery now costs a little.** The pure baseline (`safety_floor = FALSE, phylo_signal_gate = FALSE`) is lower than the default GNN-off arm on `bace_dgp` at every n (-0.8 to -2.7%) and at n = 100 on the BM and OU data (-2.7 to -4.5%, about 1.7 MCSEs), and identical at n = 1000. On 2026-09-19, under `lambda_mode = "fixed_1"`, the safety machinery was what rescued the low-signal DGP; estimating lambda now does that shrinkage itself. The pure arm is also 3 to 6 times faster. This is a candidate question for the `safety_floor` and `phylo_signal_gate` defaults, measured on simulations only (AVONET300: -0.002).

## What this does not cover

- Multi-observation data, `multi_proportion`, `zi_count`, covariates and MNAR missingness.
- GNN-on arms with `joint_solver = "rphylopars"`.
- Data sizes above n = 1000; the in-house wall time at n = 1000 is dominated by gate calibration.
- 20 seeds give MCSEs of about 0.001 to 0.01 on the simulations and 0.03 to 0.06 on AVONET300; differences under about two MCSEs are not claims.
