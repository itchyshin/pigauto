# Arc C results: the AVONET300 continuous-trait gap is the in-house joint solver, and the Rphylopars solver inside pigauto beats raw Rphylopars

2026-09-19, Totoro, pigauto main ebbf63e. Runner `script/campaign_solver_cell.R`, aggregate `script/campaign_solver_aggregate.R`, rds `script/campaign_solver_results/solver_agg.rds`. 200 cells (AVONET300 x 20 seeds; bm_mixed, ou_mixed, bace_dgp x n {100, 300, 1000} x 20 seeds), the same masks as the with/without-GNN campaign, 11 GNN-off arms, 0 errors, 20 minutes wall at 100 concurrent single-thread cells.

## Result

| AVONET300, n = 300, 20 seeds | mean z-RMSE over 4 continuous traits (MCSE) | paired difference vs raw Rphylopars (MCSE) | wall s |
|---|---:|---:|---:|
| pigauto GNN off, `joint_solver = "rphylopars"`, pure | **0.433 (0.040)** | **-0.11 (0.02)** | 71 |
| pigauto GNN off, `joint_solver = "rphylopars"`, default safety machinery | 0.438 (0.040) | -0.11 (0.02) | 75 |
| raw Rphylopars, continuous columns only | 0.547 (0.044) | 0 | 0.7 |
| pigauto GNN off, `lambda_mode = "bayes"`, in-house | 0.781 (0.063) | +0.23 (0.03) | 2.6 |
| pigauto GNN off, in-house (current default), pure or default | 0.790 (0.064) | +0.24 (0.03) | 1.5 / 7.5 |
| pigauto GNN off, in-house, continuous columns only | 0.790 (0.064) | +0.24 (0.03) | 0.3 |
| pigauto GNN off, in-house `sigma_method = "fisher_ml"` | 0.790 (0.064) | +0.24 (0.03) | 176 |

Discrete accuracy on AVONET300 is 0.77 to 0.78 for every pigauto arm.

1. **The mixed-type path is not the cause.** Continuous-only pigauto equals the 7-trait pigauto to four decimals under the in-house solver: the threshold-joint liability step leaves the continuous columns untouched. The campaign note's sentence "all four traits go through the threshold-joint fit, which is why" was an over-strong inference and is withdrawn; the 2026-08-16 decomposition (on the old handover branch) had already exonerated this path.
2. **The in-house single-pass Sigma is the cause**, and `sigma_method = "fisher_ml"` does not repair it: its `optim()` fell back to single-pass in 200 of 200 cells (bit-identical results at 100x the cost).
3. **`joint_solver = "rphylopars"` inside pigauto closes the gap and goes past raw Rphylopars by 0.11 z-RMSE (about 20%)**: with a converged REML Sigma, the discrete traits add cross-trait information to the continuous ones, which is what the unified mixed-type design is for. This is on the same masks (paired), so it is not seed noise.
4. `lambda_mode = "bayes"` on the in-house solver moves AVONET by 1% only, but it is the best arm on the low-signal `bace_dgp` (0.91 vs 1.03 at n = 1000): the lambda shrinkage does what the safety floor does, a little better.

## Where the Rphylopars solver is not better

- BM and OU simulations: within one MCSE of the in-house solver at every n (paired differences of 0.001 to 0.03; at n = 100 the in-house pure arm is 0.01 to 0.02 better).
- Low-signal `bace_dgp`: the pure Rphylopars-solver arm sits at the raw-Rphylopars level (above the mean floor), because a BM fit at lambda = 1 with no shrinkage extrapolates signal that is not there; with the default safety machinery (`rphylo_def`) it equals the in-house default (0.96 vs 0.96 at n = 1000). The safety floor, not the solver, decides this regime.
- Cost: 71 to 75 s per fit at n = 300 on AVONET (K phylopars fits through the OVR categorical path, plus the threshold-joint fit) versus 1.5 s in-house; phylopars also emits `solve(): system is singular` warnings on some masks (it still returned finite output on all 20 seeds here).

## What this means for the default (Shinichi's decision, not taken here)

- Making `joint_solver = "rphylopars"` the default, with the safety machinery kept, gives about 30% better continuous accuracy on real data and no loss on the simulations, at 50x the fit time and with a `Suggests` dependency (a fallback to in-house already exists in `fit_joint_solver()`).
- The alternative that keeps pigauto self-contained is to make the in-house solver converge: a proper REML or ML Sigma with the cross-trait EM refinement that was disabled on 2026-05-17 after divergence. That is a methods slice with its own recovery study.
- Either way, the with/without-GNN campaign's absolute AVONET numbers should be re-stated with whichever solver becomes the default; the GNN-on vs GNN-off conclusion (the held-out-cell cost, not the network) is solver-independent.

## Not covered

Multi-observation, multi-proportion and zi_count traits; MNAR; the GNN-on arms with the Rphylopars solver (the 08-16 note measured those on 5 masks: the GNN then added net value on Mass); `pheno_correlated = FALSE` made no difference to raw Rphylopars on any dataset.
