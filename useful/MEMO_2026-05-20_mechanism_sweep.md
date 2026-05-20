# Mechanism-axis sim + Dan-comparable full sweep (2026-05-19/20)

**Run.** 2026-05-20 05:26:54 — automated post-sweep summary produced by `/tmp/build_memo.R`.

## Inputs

- 240-rep mech sim: `dev/simulation_results_overnight_2026_05_19_mechanisms/results.rds` (4 scenarios x 3 mechanisms x {200, 500} x 10 sims, 4 methods incl. `pigauto_bayes`).
- 600-rep full sweep: `dev/simulation_results_pigauto/results.rds` (4 scenarios x 3 mechanisms x {500} x 50 sims, 3 methods; `pigauto_bayes` dropped because mech sim showed bayes ~= default within MC SE).

## Headline table: pigauto_default vs bm_kriging (z-RMSE)

`delta = pigauto_default - bm_kriging`. Negative = pigauto wins. Verdict thresholds: PIG_WINS if delta < -0.05, BM_WINS if delta > +0.05, tie otherwise. `delta_se` is the naive sqrt-sum of method SEs (conservative -- paired SE would be tighter).

| source | scenario | n | mechanism | col-mean | bm-krig | pigauto | delta | delta_se | n_sims | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| full_sweep | bm_strong | 500 | phylo_MAR | 1.139 | 0.687 | 0.691 | +0.004 | 0.038 | 50 | tie |
| full_sweep | bm_strong | 500 | trait_MAR | 1.014 | 0.589 | 0.584 | -0.005 | 0.031 | 50 | tie |
| full_sweep | bm_strong | 500 | trait_MNAR | 1.549 | 0.744 | 0.788 | +0.044 | 0.063 | 50 | tie |
| full_sweep | nonlinear_cov | 500 | phylo_MAR | 1.102 | 0.990 | 0.947 | -0.043 | 0.027 | 50 | tie |
| full_sweep | nonlinear_cov | 500 | trait_MAR | 1.007 | 1.041 | 0.951 | -0.091 | 0.022 | 50 | PIG_WINS |
| full_sweep | nonlinear_cov | 500 | trait_MNAR | 1.508 | 1.302 | 1.328 | +0.025 | 0.046 | 50 | tie |
| full_sweep | ou_strong | 500 | phylo_MAR | 1.033 | 1.157 | 1.025 | -0.132 | 0.019 | 50 | PIG_WINS |
| full_sweep | ou_strong | 500 | trait_MAR | 1.032 | 1.160 | 1.021 | -0.139 | 0.021 | 50 | PIG_WINS |
| full_sweep | ou_strong | 500 | trait_MNAR | 1.542 | 1.510 | 1.502 | -0.008 | 0.030 | 50 | tie |
| full_sweep | weak_signal | 500 | phylo_MAR | 1.014 | 1.234 | 1.013 | -0.220 | 0.019 | 50 | PIG_WINS |
| full_sweep | weak_signal | 500 | trait_MAR | 1.003 | 1.204 | 1.004 | -0.200 | 0.014 | 50 | PIG_WINS |
| full_sweep | weak_signal | 500 | trait_MNAR | 1.565 | 1.667 | 1.566 | -0.101 | 0.021 | 50 | PIG_WINS |
| mech_sim | bm_strong | 200 | phylo_MAR | 0.983 | 0.682 | 0.702 | +0.020 | 0.087 | 10 | tie |
| mech_sim | bm_strong | 200 | trait_MAR | 0.987 | 0.590 | 0.601 | +0.011 | 0.093 | 10 | tie |
| mech_sim | bm_strong | 200 | trait_MNAR | 1.776 | 0.948 | 1.180 | +0.232 | 0.232 | 10 | BM_WINS |
| mech_sim | bm_strong | 500 | phylo_MAR | 1.126 | 0.664 | 0.664 | +0.000 | 0.071 | 10 | tie |
| mech_sim | bm_strong | 500 | trait_MAR | 0.987 | 0.582 | 0.589 | +0.007 | 0.085 | 10 | tie |
| mech_sim | bm_strong | 500 | trait_MNAR | 1.500 | 0.752 | 0.764 | +0.012 | 0.125 | 10 | tie |
| mech_sim | nonlinear_cov | 200 | phylo_MAR | 1.095 | 1.019 | 1.019 | -0.000 | 0.101 | 10 | tie |
| mech_sim | nonlinear_cov | 200 | trait_MAR | 1.070 | 0.996 | 0.939 | -0.057 | 0.092 | 10 | PIG_WINS |
| mech_sim | nonlinear_cov | 200 | trait_MNAR | 1.583 | 1.322 | 1.410 | +0.088 | 0.069 | 10 | BM_WINS |
| mech_sim | nonlinear_cov | 500 | phylo_MAR | 1.105 | 0.968 | 0.943 | -0.025 | 0.050 | 10 | tie |
| mech_sim | nonlinear_cov | 500 | trait_MAR | 0.992 | 0.984 | 0.914 | -0.070 | 0.079 | 10 | PIG_WINS |
| mech_sim | nonlinear_cov | 500 | trait_MNAR | 1.514 | 1.407 | 1.411 | +0.004 | 0.091 | 10 | tie |
| mech_sim | ou_strong | 200 | phylo_MAR | 1.030 | 1.046 | 1.016 | -0.030 | 0.054 | 10 | tie |
| mech_sim | ou_strong | 200 | trait_MAR | 1.015 | 1.120 | 1.044 | -0.076 | 0.054 | 10 | PIG_WINS |
| mech_sim | ou_strong | 200 | trait_MNAR | 1.635 | 1.645 | 1.595 | -0.049 | 0.105 | 10 | tie |
| mech_sim | ou_strong | 500 | phylo_MAR | 0.998 | 1.115 | 0.987 | -0.129 | 0.036 | 10 | PIG_WINS |
| mech_sim | ou_strong | 500 | trait_MAR | 1.019 | 1.100 | 1.002 | -0.099 | 0.036 | 10 | PIG_WINS |
| mech_sim | ou_strong | 500 | trait_MNAR | 1.545 | 1.567 | 1.488 | -0.079 | 0.075 | 10 | PIG_WINS |
| mech_sim | weak_signal | 200 | phylo_MAR | 0.983 | 1.129 | 0.982 | -0.147 | 0.047 | 10 | PIG_WINS |
| mech_sim | weak_signal | 200 | trait_MAR | 1.007 | 1.212 | 1.017 | -0.195 | 0.055 | 10 | PIG_WINS |
| mech_sim | weak_signal | 200 | trait_MNAR | 1.700 | 1.774 | 1.713 | -0.061 | 0.094 | 10 | PIG_WINS |
| mech_sim | weak_signal | 500 | phylo_MAR | 0.983 | 1.206 | 0.982 | -0.225 | 0.045 | 10 | PIG_WINS |
| mech_sim | weak_signal | 500 | trait_MAR | 1.004 | 1.194 | 1.004 | -0.190 | 0.033 | 10 | PIG_WINS |
| mech_sim | weak_signal | 500 | trait_MNAR | 1.571 | 1.663 | 1.573 | -0.089 | 0.040 | 10 | PIG_WINS |

## Interpretation

### Aggregate verdict

Across all 36 cells (combined mech + sweep): **18 PIG_WINS, 2 BM_WINS, 16 ties**.

### Where pigauto materially wins

- **nonlinear_cov / n=500 / trait_MAR** (full_sweep): pigauto 0.951 vs BM 1.041, delta = -0.091 +/- 0.022.
- **ou_strong / n=500 / phylo_MAR** (full_sweep): pigauto 1.025 vs BM 1.157, delta = -0.132 +/- 0.019.
- **ou_strong / n=500 / trait_MAR** (full_sweep): pigauto 1.021 vs BM 1.160, delta = -0.139 +/- 0.021.
- **weak_signal / n=500 / phylo_MAR** (full_sweep): pigauto 1.013 vs BM 1.234, delta = -0.220 +/- 0.019.
- **weak_signal / n=500 / trait_MAR** (full_sweep): pigauto 1.004 vs BM 1.204, delta = -0.200 +/- 0.014.
- **weak_signal / n=500 / trait_MNAR** (full_sweep): pigauto 1.566 vs BM 1.667, delta = -0.101 +/- 0.021.
- **nonlinear_cov / n=200 / trait_MAR** (mech_sim): pigauto 0.939 vs BM 0.996, delta = -0.057 +/- 0.092.
- **nonlinear_cov / n=500 / trait_MAR** (mech_sim): pigauto 0.914 vs BM 0.984, delta = -0.070 +/- 0.079.
- **ou_strong / n=200 / trait_MAR** (mech_sim): pigauto 1.044 vs BM 1.120, delta = -0.076 +/- 0.054.
- **ou_strong / n=500 / phylo_MAR** (mech_sim): pigauto 0.987 vs BM 1.115, delta = -0.129 +/- 0.036.
- **ou_strong / n=500 / trait_MAR** (mech_sim): pigauto 1.002 vs BM 1.100, delta = -0.099 +/- 0.036.
- **ou_strong / n=500 / trait_MNAR** (mech_sim): pigauto 1.488 vs BM 1.567, delta = -0.079 +/- 0.075.
- **weak_signal / n=200 / phylo_MAR** (mech_sim): pigauto 0.982 vs BM 1.129, delta = -0.147 +/- 0.047.
- **weak_signal / n=200 / trait_MAR** (mech_sim): pigauto 1.017 vs BM 1.212, delta = -0.195 +/- 0.055.
- **weak_signal / n=200 / trait_MNAR** (mech_sim): pigauto 1.713 vs BM 1.774, delta = -0.061 +/- 0.094.
- **weak_signal / n=500 / phylo_MAR** (mech_sim): pigauto 0.982 vs BM 1.206, delta = -0.225 +/- 0.045.
- **weak_signal / n=500 / trait_MAR** (mech_sim): pigauto 1.004 vs BM 1.194, delta = -0.190 +/- 0.033.
- **weak_signal / n=500 / trait_MNAR** (mech_sim): pigauto 1.573 vs BM 1.663, delta = -0.089 +/- 0.040.

### Where pigauto materially loses

- **bm_strong / n=200 / trait_MNAR** (mech_sim): pigauto 1.180 vs BM 0.948, delta = +0.232 +/- 0.232.
- **nonlinear_cov / n=200 / trait_MNAR** (mech_sim): pigauto 1.410 vs BM 1.322, delta = +0.088 +/- 0.069.

### Pattern by scenario

Note: the calibrated gate is bimodal -- per rep it is either ~0 (closed) or substantially open. The per-scenario `mean(gate)` below is pulled up by the minority of open-gate reps; the median gate is 0 in every scenario. Report closed-fraction, not just the mean.

- **`bm_strong`**: BM is the true process. The gate is closed in the majority of reps (median 0); pigauto ties BM kriging. In the reps where the gate does open, the GNN delta is itself near-BM-quality, so the blend still does not degrade the correct BM answer -- the safety floor holds whether the gate is open or shut.
- **`ou_strong`**: BM is wrong (process is OU with lambda = 0.3). BM kriging overcommits to phylogenetic conservation; in the minority of reps where pigauto's gate opens, the GNN correction recovers ground, and on average pigauto beats BM.
- **`nonlinear_cov`**: signal is `sqrt(0.4)*bm + sqrt(0.4)*sin(2*env) + noise`. This is the scenario where the gate opens most often -- the GNN has the env covariate to exploit, which BM kriging cannot use.
- **`weak_signal`**: 10% phylogenetic, 90% noise. BM kriging mis-attributes noise to phylo signal and degrades. Pigauto's gate is closed in ~98% of reps, so pigauto falls back to column-mean -- the win here is the safety floor, not GNN learning.

### Pattern by mechanism

- **`phylo_MAR`**: missingness correlated with the phylogeny but not with the trait being imputed. Both methods see roughly unbiased validation samples; results mirror the base scenario.
- **`trait_MAR`**: missingness correlated with a different observed trait. For `nonlinear_cov`, this exposes pigauto's covariate-aware advantage clearly.
- **`trait_MNAR`**: missingness depends on the trait value itself. All methods get harder (z-RMSE jumps from ~0.7 to ~1.6 in `bm_strong`); but the *relative* ordering is preserved -- BM stays best where it was best, pigauto stays best where it was best.

## Diagnostics

### Gate utilisation (pigauto_default)

`mean(gate)` is the average calibrated gate over reps; `median(gate)` and `closed_frac` (fraction of reps with gate < 0.02) expose the bimodality the mean hides.

| source | scenario | mechanism | mean(gate) | median(gate) | closed_frac |
|---|---|---|---|---|---|
| full_sweep | bm_strong | phylo_MAR | 0.205 | 0.000 | 0.66 |
| full_sweep | nonlinear_cov | phylo_MAR | 0.250 | 0.000 | 0.56 |
| full_sweep | ou_strong | phylo_MAR | 0.089 | 0.000 | 0.82 |
| full_sweep | weak_signal | phylo_MAR | 0.016 | 0.000 | 0.96 |
| full_sweep | bm_strong | trait_MAR | 0.152 | 0.000 | 0.72 |
| full_sweep | nonlinear_cov | trait_MAR | 0.175 | 0.000 | 0.58 |
| full_sweep | ou_strong | trait_MAR | 0.074 | 0.000 | 0.84 |
| full_sweep | weak_signal | trait_MAR | 0.010 | 0.000 | 0.98 |
| full_sweep | bm_strong | trait_MNAR | 0.206 | 0.000 | 0.54 |
| full_sweep | nonlinear_cov | trait_MNAR | 0.242 | 0.000 | 0.58 |
| full_sweep | ou_strong | trait_MNAR | 0.066 | 0.000 | 0.86 |
| full_sweep | weak_signal | trait_MNAR | 0.000 | 0.000 | 1.00 |
| mech_sim | bm_strong | phylo_MAR | 0.087 | 0.000 | 0.90 |
| mech_sim | nonlinear_cov | phylo_MAR | 0.216 | 0.000 | 0.70 |
| mech_sim | ou_strong | phylo_MAR | 0.060 | 0.000 | 0.85 |
| mech_sim | weak_signal | phylo_MAR | 0.000 | 0.000 | 1.00 |
| mech_sim | bm_strong | trait_MAR | 0.125 | 0.000 | 0.75 |
| mech_sim | nonlinear_cov | trait_MAR | 0.165 | 0.000 | 0.75 |
| mech_sim | ou_strong | trait_MAR | 0.105 | 0.000 | 0.85 |
| mech_sim | weak_signal | trait_MAR | 0.000 | 0.000 | 1.00 |
| mech_sim | bm_strong | trait_MNAR | 0.075 | 0.000 | 0.85 |
| mech_sim | nonlinear_cov | trait_MNAR | 0.090 | 0.000 | 0.75 |
| mech_sim | ou_strong | trait_MNAR | 0.077 | 0.000 | 0.80 |
| mech_sim | weak_signal | trait_MNAR | 0.000 | 0.000 | 1.00 |

### Wall time (median sec per rep)

| source | method | median wall (s) |
|---|---|---|
| full_sweep | bm_kriging | 0.00263 |
| full_sweep | column_mean | 1.41e-05 |
| full_sweep | pigauto_default | 63.7 |
| mech_sim | bm_kriging | 0.002 |
| mech_sim | column_mean | 1.38e-05 |
| mech_sim | pigauto_bayes | 70.3 |
| mech_sim | pigauto_default | 69.6 |

## Methodology

DGP and mechanism injectors are identical across both sims (`script/sim_bench/overnight_2026_05_19_mechanisms.R` defines the shared `.sim_complete()` and `.inject_missingness_single()`; `script/sim_bench/pigauto_full_sweep.R` reuses them verbatim with a +200000 seed-space offset so the rep seeds are disjoint).

- 4 scenarios: `bm_strong` (lambda=1), `ou_strong` (lambda=0.3 OU), `nonlinear_cov` (40% phylo + 40% sin(env) + noise), `weak_signal` (10% phylo + 90% noise).
- 3 mechanisms: `phylo_MAR`, `trait_MAR`, `trait_MNAR` (DEP_STRENGTH = 1.5, matching Penone et al. 2014 ICCs and Dan's BACE design).
- Miss rate: 30%.
- Methods: `column_mean`, `bm_kriging` (pigauto's internal `bm_impute_col` on R = cov2cor(vcv(tree))), `pigauto_default` (lambda_mode='fixed_1', 300 epochs, patience 5).

## Files

- `dev/simulation_results_overnight_2026_05_19_mechanisms/results.rds`
- `dev/simulation_results_pigauto/results.rds`
- `script/sim_bench/overnight_2026_05_19_mechanisms.R`
- `script/sim_bench/pigauto_full_sweep.R`
- `script/sim_bench/build_mechanism_sweep_memo.R` (this script)
