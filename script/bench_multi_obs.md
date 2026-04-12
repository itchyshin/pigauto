# Multi-observation-per-species benchmark

Run on: 2026-04-11 10:31:40
Machine: Darwin 25.4.0 (arm64), R 4.5.2
Species: 200, reps: 2, epochs: 200
Total wall: 37.1 min

## Data-generating process

```
CTmax_ij = mu + phylo_i + beta * acclim_temp_ij + epsilon_ij
mu = 30.0, sigma_residual = 0.50
acclim_temp ~ N(20.0, 5.0)
obs_per_species ~ Pois(5), clamped [1, 20]
```

## Methods

- **species_mean**: mean of observed CTmax per species; grand mean for unobserved species.
- **pigauto_no_cov**: `impute(traits, tree, species_col = "species")`
- **pigauto_cov**: `impute(traits, tree, species_col = "species", covariates = acclim_temp)`

## Observation-level RMSE (lower is better)

| method | lambda | beta | sp_miss | obs_RMSE |
|--------|--------|------|---------|----------|
| species_mean | 0.5 | 0.0 | 0.5 | 1.1232 |
| pigauto_no_cov | 0.5 | 0.0 | 0.5 | 0.8914 |
| pigauto_cov | 0.5 | 0.0 | 0.5 | 0.8894 |
| species_mean | 0.5 | 0.0 | 0.8 | 1.4126 |
| pigauto_no_cov | 0.5 | 0.0 | 0.8 | 1.0683 |
| pigauto_cov | 0.5 | 0.0 | 0.8 | 1.0713 |
| species_mean | 0.5 | 0.5 | 0.5 | 2.9272 |
| pigauto_no_cov | 0.5 | 0.5 | 0.5 | 2.8668 |
| pigauto_cov | 0.5 | 0.5 | 0.5 | 2.4848 |
| species_mean | 0.5 | 0.5 | 0.8 | 2.9247 |
| pigauto_no_cov | 0.5 | 0.5 | 0.8 | 2.9400 |
| pigauto_cov | 0.5 | 0.5 | 0.8 | 2.6059 |
| species_mean | 0.5 | 1.0 | 0.5 | 5.1071 |
| pigauto_no_cov | 0.5 | 1.0 | 0.5 | 5.4290 |
| pigauto_cov | 0.5 | 1.0 | 0.5 | 4.6016 |
| species_mean | 0.5 | 1.0 | 0.8 | 5.2278 |
| pigauto_no_cov | 0.5 | 1.0 | 0.8 | 5.2843 |
| pigauto_cov | 0.5 | 1.0 | 0.8 | 4.4811 |
| species_mean | 0.9 | 0.0 | 0.5 | 1.3816 |
| pigauto_no_cov | 0.9 | 0.0 | 0.5 | 1.0907 |
| pigauto_cov | 0.9 | 0.0 | 0.5 | 1.0878 |
| species_mean | 0.9 | 0.0 | 0.8 | 1.9755 |
| pigauto_no_cov | 0.9 | 0.0 | 0.8 | 1.5500 |
| pigauto_cov | 0.9 | 0.0 | 0.8 | 1.5517 |
| species_mean | 0.9 | 0.5 | 0.5 | 3.0603 |
| pigauto_no_cov | 0.9 | 0.5 | 0.5 | 2.9719 |
| pigauto_cov | 0.9 | 0.5 | 0.5 | 2.6218 |
| species_mean | 0.9 | 0.5 | 0.8 | 3.0626 |
| pigauto_no_cov | 0.9 | 0.5 | 0.8 | 3.1912 |
| pigauto_cov | 0.9 | 0.5 | 0.8 | 2.8671 |
| species_mean | 0.9 | 1.0 | 0.5 | 5.5368 |
| pigauto_no_cov | 0.9 | 1.0 | 0.5 | 5.7449 |
| pigauto_cov | 0.9 | 1.0 | 0.5 | 4.9300 |
| species_mean | 0.9 | 1.0 | 0.8 | 5.4523 |
| pigauto_no_cov | 0.9 | 1.0 | 0.8 | 5.4962 |
| pigauto_cov | 0.9 | 1.0 | 0.8 | 4.7032 |

## Species-level RMSE (unobserved species, lower is better)

| method | lambda | beta | sp_miss | sp_RMSE |
|--------|--------|------|---------|---------|
| species_mean | 0.5 | 0.0 | 0.5 | 1.1091 |
| pigauto_no_cov | 0.5 | 0.0 | 0.5 | 0.7804 |
| pigauto_cov | 0.5 | 0.0 | 0.5 | 0.7775 |
| species_mean | 0.5 | 0.0 | 0.8 | 1.3434 |
| pigauto_no_cov | 0.5 | 0.0 | 0.8 | 0.9338 |
| pigauto_cov | 0.5 | 0.0 | 0.8 | 0.9374 |
| species_mean | 0.5 | 0.5 | 0.5 | 1.2797 |
| pigauto_no_cov | 0.5 | 0.5 | 0.5 | 1.1324 |
| pigauto_cov | 0.5 | 0.5 | 0.5 | 1.0915 |
| species_mean | 0.5 | 0.5 | 0.8 | 1.4960 |
| pigauto_no_cov | 0.5 | 0.5 | 0.8 | 1.4235 |
| pigauto_cov | 0.5 | 0.5 | 0.8 | 1.3971 |
| species_mean | 0.5 | 1.0 | 0.5 | 1.3076 |
| pigauto_no_cov | 0.5 | 1.0 | 0.5 | 2.1506 |
| pigauto_cov | 0.5 | 1.0 | 0.5 | 1.9491 |
| species_mean | 0.5 | 1.0 | 0.8 | 1.2900 |
| pigauto_no_cov | 0.5 | 1.0 | 0.8 | 1.9825 |
| pigauto_cov | 0.5 | 1.0 | 0.8 | 1.9005 |
| species_mean | 0.9 | 0.0 | 0.5 | 1.4088 |
| pigauto_no_cov | 0.9 | 0.0 | 0.5 | 1.0291 |
| pigauto_cov | 0.9 | 0.0 | 0.5 | 1.0265 |
| species_mean | 0.9 | 0.0 | 0.8 | 1.9242 |
| pigauto_no_cov | 0.9 | 0.0 | 0.8 | 1.3973 |
| pigauto_cov | 0.9 | 0.0 | 0.8 | 1.3987 |
| species_mean | 0.9 | 0.5 | 0.5 | 1.6372 |
| pigauto_no_cov | 0.9 | 0.5 | 0.5 | 1.4668 |
| pigauto_cov | 0.9 | 0.5 | 0.5 | 1.4294 |
| species_mean | 0.9 | 0.5 | 0.8 | 1.7644 |
| pigauto_no_cov | 0.9 | 0.5 | 0.8 | 1.8639 |
| pigauto_cov | 0.9 | 0.5 | 0.8 | 1.8123 |
| species_mean | 0.9 | 1.0 | 0.5 | 1.8208 |
| pigauto_no_cov | 0.9 | 1.0 | 0.5 | 2.2964 |
| pigauto_cov | 0.9 | 1.0 | 0.5 | 2.1426 |
| species_mean | 0.9 | 1.0 | 0.8 | 1.8283 |
| pigauto_no_cov | 0.9 | 1.0 | 0.8 | 1.9941 |
| pigauto_cov | 0.9 | 1.0 | 0.8 | 1.8396 |

## Observation-level Pearson r (higher is better)

| method | lambda | beta | sp_miss | pearson_r |
|--------|--------|------|---------|-----------|
| species_mean | 0.5 | 0.0 | 0.5 | 0.3781 |
| pigauto_no_cov | 0.5 | 0.0 | 0.5 | 0.6754 |
| pigauto_cov | 0.5 | 0.0 | 0.5 | 0.6753 |
| species_mean | 0.5 | 0.0 | 0.8 | 0.2043 |
| pigauto_no_cov | 0.5 | 0.0 | 0.8 | 0.7062 |
| pigauto_cov | 0.5 | 0.0 | 0.8 | 0.7069 |
| species_mean | 0.5 | 0.5 | 0.5 | 0.1267 |
| pigauto_no_cov | 0.5 | 0.5 | 0.5 | 0.2764 |
| pigauto_cov | 0.5 | 0.5 | 0.5 | 0.5356 |
| species_mean | 0.5 | 0.5 | 0.8 | 0.0599 |
| pigauto_no_cov | 0.5 | 0.5 | 0.8 | 0.2114 |
| pigauto_cov | 0.5 | 0.5 | 0.8 | 0.4482 |
| species_mean | 0.5 | 1.0 | 0.5 | 0.0729 |
| pigauto_no_cov | 0.5 | 1.0 | 0.5 | 0.0565 |
| pigauto_cov | 0.5 | 1.0 | 0.5 | 0.4225 |
| species_mean | 0.5 | 1.0 | 0.8 | -0.0051 |
| pigauto_no_cov | 0.5 | 1.0 | 0.8 | 0.1096 |
| pigauto_cov | 0.5 | 1.0 | 0.8 | 0.5574 |
| species_mean | 0.9 | 0.0 | 0.5 | 0.3523 |
| pigauto_no_cov | 0.9 | 0.0 | 0.5 | 0.6847 |
| pigauto_cov | 0.9 | 0.0 | 0.5 | 0.6858 |
| species_mean | 0.9 | 0.0 | 0.8 | 0.2048 |
| pigauto_no_cov | 0.9 | 0.0 | 0.8 | 0.6514 |
| pigauto_cov | 0.9 | 0.0 | 0.8 | 0.6503 |
| species_mean | 0.9 | 0.5 | 0.5 | 0.1570 |
| pigauto_no_cov | 0.9 | 0.5 | 0.5 | 0.3530 |
| pigauto_cov | 0.9 | 0.5 | 0.5 | 0.5282 |
| species_mean | 0.9 | 0.5 | 0.8 | 0.0676 |
| pigauto_no_cov | 0.9 | 0.5 | 0.8 | 0.2218 |
| pigauto_cov | 0.9 | 0.5 | 0.8 | 0.3976 |
| species_mean | 0.9 | 1.0 | 0.5 | 0.0644 |
| pigauto_no_cov | 0.9 | 1.0 | 0.5 | 0.1069 |
| pigauto_cov | 0.9 | 1.0 | 0.5 | 0.4326 |
| species_mean | 0.9 | 1.0 | 0.8 | 0.0448 |
| pigauto_no_cov | 0.9 | 1.0 | 0.8 | 0.1541 |
| pigauto_cov | 0.9 | 1.0 | 0.8 | 0.5500 |

## Covariate lift (pigauto_cov / pigauto_no_cov RMSE ratio)

Ratio < 1 means covariates help; ratio > 1 means they hurt.

| lambda | beta | sp_miss | RMSE_nocov | RMSE_cov | ratio |
|--------|------|---------|------------|----------|-------|
| 0.5 | 0.0 | 0.5 | 0.8914 | 0.8894 | 0.998 |
| 0.5 | 0.0 | 0.8 | 1.0683 | 1.0713 | 1.003 |
| 0.5 | 0.5 | 0.5 | 2.8668 | 2.4848 | 0.867 |
| 0.5 | 0.5 | 0.8 | 2.9400 | 2.6059 | 0.886 |
| 0.5 | 1.0 | 0.5 | 5.4290 | 4.6016 | 0.848 |
| 0.5 | 1.0 | 0.8 | 5.2843 | 4.4811 | 0.848 |
| 0.9 | 0.0 | 0.5 | 1.0907 | 1.0878 | 0.997 |
| 0.9 | 0.0 | 0.8 | 1.5500 | 1.5517 | 1.001 |
| 0.9 | 0.5 | 0.5 | 2.9719 | 2.6218 | 0.882 |
| 0.9 | 0.5 | 0.8 | 3.1912 | 2.8671 | 0.898 |
| 0.9 | 1.0 | 0.5 | 5.7449 | 4.9300 | 0.858 |
| 0.9 | 1.0 | 0.8 | 5.4962 | 4.7032 | 0.856 |

---
Generated: 2026-04-11 10:31
