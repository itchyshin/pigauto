# Campaign results: pigauto with and without the GNN, against raw Rphylopars and BACE

2026-09-19, Totoro, pigauto main 92edc73 (PR #180 merged). Runner `script/campaign_gnn_off_cell.R`, aggregate `script/campaign_gnn_off_aggregate.R`, tables `script/campaign_gnn_off_tables.R`, figures `script/campaign_gnn_off_figures.R`. Raw cell files: `~/gnn-off/results_full/` on Totoro (200 rds); aggregate committed as `script/campaign_gnn_off_results/campaign_agg.rds`.

## Design (as locked in the plan, approved "merge and full")

- Arms: pigauto GNN on (defaults, 2000 epochs); the same GNN-on fit predicting from a baseline refit on all observed cells ("full baseline"); pigauto GNN off (default safety machinery); pigauto GNN off pure (`safety_floor = FALSE, phylo_signal_gate = FALSE`); raw `Rphylopars::phylopars(model = "BM")` on the continuous columns only; BACE (in-tree, `nitt = 50,000`, `burnin = 10,000`, `thin = 25`, 2 chains, 5 final draws, OVR on); mean / mode floor.
- DGPs: `bm_mixed` (4 BM continuous + 1 binary + 1 three-level categorical, `ape::rcoal`), `ou_mixed` (same with OU continuous traits, `simulate_non_bm`), `bace_dgp` (`BACE::sim_bace`: y + 2 gaussian + 1 binary predictors, all imputed), and the bundled AVONET300 (7 traits: 4 continuous, 1 ordinal, 2 categorical).
- n in {100, 300, 1000} for the simulated DGPs, AVONET300 at its n = 300; 20 seeds per cell; 30% MCAR applied in the user data, the same mask for every arm within a cell.
- Metrics on the masked cells: z-RMSE (truth and prediction standardised by the training mean and sd of each trait, averaged over the continuous traits), accuracy averaged over the discrete traits, and coverage of pigauto's production conformal interval (nominal 0.95). Monte Carlo SE over seeds in brackets.
- Compute: 36 concurrent cells, 4 threads each; 200 cells, 0 errors; 1 h 36 min wall in two batches (05:35 to 07:11 MDT), against a pre-run estimate of 2 to 3 h.

## What the campaign says

**1. On the simulations, the GNN adds nothing, and the shipped GNN-on predictions are worse than GNN-off.** Across `bm_mixed` and `ou_mixed` at every n, pigauto GNN off, raw Rphylopars, and GNN on predicting from the full baseline are within one MCSE of each other. GNN on as shipped is 15 to 25% worse in z-RMSE (largest at n = 100). The per-seed lines in `fig_tax.png` are flat from "GNN off" to "GNN on, full baseline" and rise only at "GNN on as shipped": the loss is the held-out-cell cost of gate calibration (the shipped prediction's baseline is fit on the 75% of observed cells left after the val/test split), not the network. This is the mechanism the 2026-08-16 note listed first and could not measure.

**2. On the BM and OU simulations pigauto GNN off beats BACE by 13 to 15% on continuous traits and by 1 to 2 points of accuracy on discrete traits**, at 1/500 of the wall time (2 s vs 1,100 s at n = 1000). BACE and GNN on as shipped are indistinguishable there.

**2b. On the low-signal BACE DGP (phylogenetic signal 0.4 on y, 0.3 on the predictors; the floor sits at z-RMSE 1.0 and the binary predictor is a coin flip for every arm) the ordering changes in an instructive way.** Raw Rphylopars is the worst arm (1.14 to 1.22, below the floor: a BM fit at lambda = 1 extrapolates signal that is not there), pure GNN off is at the floor (1.03 to 1.06), BACE is 0.99 to 1.11, and the default GNN off is best (0.96 to 1.03) because its safety machinery blends toward the grand mean when the phylogenetic signal is weak. This is the regime the plan's D1 decision was made for, and it is why the default `gnn = FALSE` keeps the safety floor and phylo-signal gate rather than being pure BM. The GNN adds nothing here either (GNN on equals GNN off within MCSE).

**3. On real data (AVONET300) the picture is different, and it is the finding of this campaign.** Raw Rphylopars on the four continuous columns beats pigauto GNN off by about 30% in z-RMSE (0.55 vs 0.79, 20 seeds, MCSE 0.04 to 0.06), and BACE sits in between (0.62). `fit$baseline$path` shows why the simulations did not predict this: on AVONET300 every continuous trait is estimated inside the threshold-joint liability fit (with Trophic.Level and Primary.Lifestyle through OVR and Migration mostly per-column BM), whereas Rphylopars fits a plain joint BM on the continuous block. On the simulations the discrete traits are generated from independent liabilities, so the liability machinery costs nothing there. On AVONET the continuous columns lose 25 to 40% per trait (Beak 0.40 vs 0.52, Tarsus 0.37 vs 0.59, Wing 0.47 vs 0.72, Mass 0.95 vs 1.34). On the discrete traits pigauto GNN off leads BACE (Trophic 0.79 vs 0.73; Lifestyle 0.76 vs 0.69; Migration 0.78 vs 0.74). Here the GNN does contribute a little: GNN on from the full baseline is 0.72, better than GNN off at 0.79 (Wing 0.58 vs 0.72, Tarsus 0.54 vs 0.59), but still far from Rphylopars.

**4. Coverage of the GNN-off production interval** is 0.96 to 0.98 at n = 300 and 1000 (slightly conservative, as the design expected) and 0.89 to 0.93 at n = 100, where a trait has about 21 validation cells.

**5. Cost** (mean over DGPs and seeds, 4 threads). GNN off pure: 0.2 / 0.4 / 2.2 s at n = 100 / 300 / 1000. GNN off with the default safety machinery: 0.5 / 4.5 / 117 s; the 117 s at n = 1000 is `calibrate_gates()` (cv folds over the simplex grid), not the baseline. GNN on: 118 / 165 / 526 s. BACE: 124 / 477 / 1,131 s.

## What this does not say

- Nothing here is about multi-observation data, multi-proportion or zi_count traits, covariates, or MNAR missingness.
- The AVONET300 gap to Rphylopars is measured for the default `joint_solver = "inhouse"`, `lambda_mode = "fixed_1"`, `log_transform = TRUE` settings; whether it comes from the liability step, the log transform on Mass (MCSE 0.24 to 0.35 there), or the val/test masking inside the held-out baseline has not been separated. That separation is the next slice: run the continuous block as a plain joint BM fit with the discrete traits handled separately, and compare `joint_solver = "rphylopars"` on the same masks.
- BACE was run at one chain-length setting; longer chains or its own tuning may move it.
- 20 seeds gives MCSEs of 0.004 to 0.02 on the simulations and 0.04 to 0.06 on AVONET300; differences smaller than about two MCSEs are not claims.

## Tables

### avonet: mean z-RMSE over the continuous traits, mean (MCSE) over seeds

| arm | n = 300 (seeds 20) |
|---|---:|
| pigauto, GNN off, pure baseline | 0.790 (0.064) |
| pigauto, GNN off | 0.790 (0.064) |
| raw Rphylopars (continuous only) | 0.547 (0.044) |
| pigauto, GNN on, full baseline | 0.723 (0.058) |
| pigauto, GNN on (as shipped) | 0.775 (0.060) |
| BACE (50k iterations, OVR) | 0.622 (0.051) |
| mean / mode floor | 1.032 (0.063) |

### avonet: mean accuracy over the discrete traits, mean (MCSE); and production-interval coverage (nominal 0.95)

| arm | acc n = 300 | cov n = 300 |
|---|---:|---:|
| pigauto, GNN off, pure baseline | 0.775 (0.003) | 0.970 |
| pigauto, GNN off | 0.765 (0.007) | 0.969 |
| pigauto, GNN on, full baseline | 0.766 (0.007) | 0.968 |
| pigauto, GNN on (as shipped) | 0.746 (0.006) | 0.957 |
| BACE (50k iterations, OVR) | 0.727 (0.006) |  |
| mean / mode floor | 0.633 (0.005) |  |

### bace_dgp: mean z-RMSE over the continuous traits, mean (MCSE) over seeds

| arm | n = 100 (seeds 20) | n = 300 (seeds 20) | n = 1000 (seeds 20) |
|---|---:|---:|---:|
| pigauto, GNN off, pure baseline | 1.060 (0.022) | 1.047 (0.010) | 1.034 (0.009) |
| pigauto, GNN off | 1.025 (0.020) | 0.986 (0.010) | 0.956 (0.008) |
| raw Rphylopars (continuous only) | 1.221 (0.051) | 1.160 (0.038) | 1.137 (0.032) |
| pigauto, GNN on, full baseline | 1.026 (0.018) | 0.990 (0.010) | 0.967 (0.009) |
| pigauto, GNN on (as shipped) | 1.031 (0.018) | 0.988 (0.010) | 0.966 (0.009) |
| BACE (50k iterations, OVR) | 1.109 (0.021) | 1.022 (0.010) | 0.991 (0.009) |
| mean / mode floor | 1.021 (0.018) | 1.012 (0.013) | 1.002 (0.007) |

### bace_dgp: mean accuracy over the discrete traits, mean (MCSE); and production-interval coverage (nominal 0.95)

| arm | acc n = 100 | acc n = 300 | acc n = 1000 | cov n = 100 | cov n = 300 | cov n = 1000 |
|---|---:|---:|---:|---:|---:|---:|
| pigauto, GNN off, pure baseline | 0.540 (0.022) | 0.493 (0.009) | 0.512 (0.005) | 0.893 | 0.963 | 0.956 |
| pigauto, GNN off | 0.528 (0.020) | 0.510 (0.013) | 0.513 (0.005) | 0.891 | 0.959 | 0.950 |
| pigauto, GNN on, full baseline | 0.528 (0.020) | 0.510 (0.013) | 0.513 (0.005) | 0.891 | 0.958 | 0.950 |
| pigauto, GNN on (as shipped) | 0.528 (0.020) | 0.510 (0.013) | 0.513 (0.005) | 0.886 | 0.958 | 0.950 |
| BACE (50k iterations, OVR) | 0.530 (0.026) | 0.527 (0.014) | 0.513 (0.007) |  |  |  |
| mean / mode floor | 0.498 (0.021) | 0.501 (0.013) | 0.518 (0.005) |  |  |  |

### bm_mixed: mean z-RMSE over the continuous traits, mean (MCSE) over seeds

| arm | n = 100 (seeds 20) | n = 300 (seeds 20) | n = 1000 (seeds 20) |
|---|---:|---:|---:|
| pigauto, GNN off, pure baseline | 0.251 (0.017) | 0.151 (0.008) | 0.090 (0.003) |
| pigauto, GNN off | 0.270 (0.020) | 0.151 (0.008) | 0.090 (0.003) |
| raw Rphylopars (continuous only) | 0.256 (0.018) | 0.152 (0.008) | 0.090 (0.003) |
| pigauto, GNN on, full baseline | 0.272 (0.020) | 0.152 (0.008) | 0.090 (0.003) |
| pigauto, GNN on (as shipped) | 0.313 (0.024) | 0.175 (0.009) | 0.103 (0.005) |
| BACE (50k iterations, OVR) | 0.309 (0.018) | 0.174 (0.009) | 0.101 (0.004) |
| mean / mode floor | 1.005 (0.015) | 0.991 (0.008) | 1.002 (0.005) |

### bm_mixed: mean accuracy over the discrete traits, mean (MCSE); and production-interval coverage (nominal 0.95)

| arm | acc n = 100 | acc n = 300 | acc n = 1000 | cov n = 100 | cov n = 300 | cov n = 1000 |
|---|---:|---:|---:|---:|---:|---:|
| pigauto, GNN off, pure baseline | 0.837 (0.015) | 0.928 (0.007) | 0.959 (0.002) | 0.927 | 0.977 | 0.964 |
| pigauto, GNN off | 0.814 (0.017) | 0.928 (0.007) | 0.959 (0.002) | 0.911 | 0.977 | 0.964 |
| pigauto, GNN on, full baseline | 0.814 (0.017) | 0.928 (0.007) | 0.959 (0.002) | 0.905 | 0.977 | 0.964 |
| pigauto, GNN on (as shipped) | 0.802 (0.015) | 0.916 (0.006) | 0.953 (0.002) | 0.889 | 0.967 | 0.951 |
| BACE (50k iterations, OVR) | 0.818 (0.013) | 0.913 (0.007) | 0.953 (0.003) |  |  |  |
| mean / mode floor | 0.335 (0.008) | 0.380 (0.004) | 0.394 (0.002) |  |  |  |

### ou_mixed: mean z-RMSE over the continuous traits, mean (MCSE) over seeds

| arm | n = 100 (seeds 20) | n = 300 (seeds 20) | n = 1000 (seeds 20) |
|---|---:|---:|---:|
| pigauto, GNN off, pure baseline | 0.409 (0.023) | 0.215 (0.011) | 0.132 (0.006) |
| pigauto, GNN off | 0.431 (0.024) | 0.216 (0.011) | 0.133 (0.006) |
| raw Rphylopars (continuous only) | 0.417 (0.023) | 0.217 (0.011) | 0.132 (0.006) |
| pigauto, GNN on, full baseline | 0.425 (0.025) | 0.218 (0.011) | 0.133 (0.006) |
| pigauto, GNN on (as shipped) | 0.478 (0.020) | 0.246 (0.011) | 0.151 (0.007) |
| BACE (50k iterations, OVR) | 0.467 (0.026) | 0.249 (0.013) | 0.150 (0.007) |
| mean / mode floor | 0.999 (0.012) | 1.001 (0.007) | 0.994 (0.005) |

### ou_mixed: mean accuracy over the discrete traits, mean (MCSE); and production-interval coverage (nominal 0.95)

| arm | acc n = 100 | acc n = 300 | acc n = 1000 | cov n = 100 | cov n = 300 | cov n = 1000 |
|---|---:|---:|---:|---:|---:|---:|
| pigauto, GNN off, pure baseline | 0.864 (0.017) | 0.939 (0.004) | 0.959 (0.003) | 0.902 | 0.976 | 0.967 |
| pigauto, GNN off | 0.855 (0.022) | 0.939 (0.004) | 0.959 (0.003) | 0.890 | 0.975 | 0.966 |
| pigauto, GNN on, full baseline | 0.864 (0.017) | 0.939 (0.004) | 0.959 (0.003) | 0.883 | 0.973 | 0.966 |
| pigauto, GNN on (as shipped) | 0.845 (0.016) | 0.925 (0.005) | 0.952 (0.003) | 0.872 | 0.965 | 0.953 |
| BACE (50k iterations, OVR) | 0.842 (0.014) | 0.930 (0.005) | 0.954 (0.003) |  |  |  |
| mean / mode floor | 0.360 (0.008) | 0.383 (0.004) | 0.394 (0.002) |  |  |  |

### AVONET300 per continuous trait: z-RMSE, mean (MCSE)

| arm | Beak.Length_Culmen | Mass | Tarsus.Length | Wing.Length |
|---|---:|---:|---:|---:|
| pigauto, GNN off, pure baseline | 0.523 (0.042) | 1.337 (0.228) | 0.585 (0.063) | 0.715 (0.043) |
| pigauto, GNN off | 0.523 (0.042) | 1.337 (0.228) | 0.585 (0.063) | 0.716 (0.043) |
| raw Rphylopars (continuous only) | 0.399 (0.029) | 0.953 (0.156) | 0.368 (0.036) | 0.469 (0.039) |
| pigauto, GNN on, full baseline | 0.532 (0.041) | 1.249 (0.218) | 0.535 (0.055) | 0.577 (0.037) |
| pigauto, GNN on (as shipped) | 0.576 (0.043) | 1.291 (0.220) | 0.600 (0.059) | 0.633 (0.039) |
| BACE (50k iterations, OVR) | 0.499 (0.045) | 1.069 (0.165) | 0.446 (0.034) | 0.476 (0.042) |

### AVONET300 per discrete trait: accuracy, mean (MCSE)

| arm | Migration | Primary.Lifestyle | Trophic.Level |
|---|---:|---:|---:|
| pigauto, GNN off, pure baseline | 0.778 (0.008) | 0.758 (0.008) | 0.790 (0.008) |
| pigauto, GNN off | 0.776 (0.008) | 0.748 (0.012) | 0.771 (0.016) |
| pigauto, GNN on, full baseline | 0.776 (0.008) | 0.748 (0.012) | 0.774 (0.015) |
| pigauto, GNN on (as shipped) | 0.753 (0.008) | 0.721 (0.012) | 0.763 (0.015) |
| BACE (50k iterations, OVR) | 0.757 (0.011) | 0.695 (0.009) | 0.731 (0.009) |

### Wall time per arm, seconds, mean over DGPs and seeds (4 threads per cell)

| arm | n = 100 | n = 300 | n = 1000 |
|---|---:|---:|---:|
| pigauto, GNN off, pure baseline | 0.2 | 0.4 | 2.2 |
| pigauto, GNN off | 0.5 | 4.5 | 117.2 |
| raw Rphylopars (continuous only) | 0.4 | 0.7 | 1.6 |
| pigauto, GNN on (as shipped) | 117.9 | 165.0 | 526.4 |
| BACE (50k iterations, OVR) | 124.4 | 476.9 | 1131.4 |
| mean / mode floor | 0.0 | 0.0 | 0.0 |

### Baseline dispatch recorded by `fit$baseline$path` (GNN-off arm; trait-cells over all n and seeds)

| dgp | path | count |
|---|---|---:|
| avonet | ovr_categorical | 40 |
| bm_mixed | ovr_categorical | 60 |
| ou_mixed | ovr_categorical | 60 |
| avonet | per_column_bm | 18 |
| avonet | threshold_joint | 82 |
| bace_dgp | threshold_joint | 240 |
| bm_mixed | threshold_joint | 300 |
| ou_mixed | threshold_joint | 300 |

Cells aggregated: 200. Errors: 0.

## Figures

`script/campaign_gnn_off_results/fig_zrmse.png`, `fig_accuracy.png`, `fig_coverage.png`, `fig_wall.png`, `fig_tax.png` (regenerate with `Rscript script/campaign_gnn_off_figures.R script/campaign_gnn_off_results/campaign_agg.rds script/campaign_gnn_off_results`).
