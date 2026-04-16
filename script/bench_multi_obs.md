# Multi-observation-per-species benchmark

Run on: 2026-04-16 11:31:46
Machine: Darwin 25.4.0 (arm64), R 4.5.2
Species: 200, reps: 2, epochs: 200
Total wall: 54.0 min

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
| pigauto_no_cov | 0.5 | 0.0 | 0.5 | 0.8941 |
| pigauto_cov | 0.5 | 0.0 | 0.5 | 0.8876 |
| species_mean | 0.5 | 0.0 | 0.8 | 1.4126 |
| pigauto_no_cov | 0.5 | 0.0 | 0.8 | 1.0684 |
| pigauto_cov | 0.5 | 0.0 | 0.8 | 1.0663 |
| species_mean | 0.5 | 0.5 | 0.5 | 2.9272 |
| pigauto_no_cov | 0.5 | 0.5 | 0.5 | 2.9220 |
| pigauto_cov | 0.5 | 0.5 | 0.5 | 2.4830 |
| species_mean | 0.5 | 0.5 | 0.8 | 2.9247 |
| pigauto_no_cov | 0.5 | 0.5 | 0.8 | 2.9385 |
| pigauto_cov | 0.5 | 0.5 | 0.8 | 2.5868 |
| species_mean | 0.5 | 1.0 | 0.5 | 5.1071 |
| pigauto_no_cov | 0.5 | 1.0 | 0.5 | 5.4872 |
| pigauto_cov | 0.5 | 1.0 | 0.5 | 4.6799 |
| species_mean | 0.5 | 1.0 | 0.8 | 5.2278 |
| pigauto_no_cov | 0.5 | 1.0 | 0.8 | 5.3110 |
| pigauto_cov | 0.5 | 1.0 | 0.8 | 4.4848 |
| species_mean | 0.9 | 0.0 | 0.5 | 1.3816 |
| pigauto_no_cov | 0.9 | 0.0 | 0.5 | 1.0914 |
| pigauto_cov | 0.9 | 0.0 | 0.5 | 1.0886 |
| species_mean | 0.9 | 0.0 | 0.8 | 1.9755 |
| pigauto_no_cov | 0.9 | 0.0 | 0.8 | 1.5510 |
| pigauto_cov | 0.9 | 0.0 | 0.8 | 1.5505 |
| species_mean | 0.9 | 0.5 | 0.5 | 3.0603 |
| pigauto_no_cov | 0.9 | 0.5 | 0.5 | 3.0138 |
| pigauto_cov | 0.9 | 0.5 | 0.5 | 2.6206 |
| species_mean | 0.9 | 0.5 | 0.8 | 3.0626 |
| pigauto_no_cov | 0.9 | 0.5 | 0.8 | 3.1889 |
| pigauto_cov | 0.9 | 0.5 | 0.8 | 2.8739 |
| species_mean | 0.9 | 1.0 | 0.5 | 5.5368 |
| pigauto_no_cov | 0.9 | 1.0 | 0.5 | 6.0190 |
| pigauto_cov | 0.9 | 1.0 | 0.5 | 4.8947 |
| species_mean | 0.9 | 1.0 | 0.8 | 5.4523 |
| pigauto_no_cov | 0.9 | 1.0 | 0.8 | 5.5354 |
| pigauto_cov | 0.9 | 1.0 | 0.8 | 4.6893 |

## Species-level RMSE (unobserved species, lower is better)

| method | lambda | beta | sp_miss | sp_RMSE |
|--------|--------|------|---------|---------|
| species_mean | 0.5 | 0.0 | 0.5 | 1.1091 |
| pigauto_no_cov | 0.5 | 0.0 | 0.5 | 0.7838 |
| pigauto_cov | 0.5 | 0.0 | 0.5 | 0.7761 |
| species_mean | 0.5 | 0.0 | 0.8 | 1.3434 |
| pigauto_no_cov | 0.5 | 0.0 | 0.8 | 0.9338 |
| pigauto_cov | 0.5 | 0.0 | 0.8 | 0.9322 |
| species_mean | 0.5 | 0.5 | 0.5 | 1.2797 |
| pigauto_no_cov | 0.5 | 0.5 | 0.5 | 1.2326 |
| pigauto_cov | 0.5 | 0.5 | 0.5 | 1.0930 |
| species_mean | 0.5 | 0.5 | 0.8 | 1.4960 |
| pigauto_no_cov | 0.5 | 0.5 | 0.8 | 1.4148 |
| pigauto_cov | 0.5 | 0.5 | 0.8 | 1.3862 |
| species_mean | 0.5 | 1.0 | 0.5 | 1.3076 |
| pigauto_no_cov | 0.5 | 1.0 | 0.5 | 2.2932 |
| pigauto_cov | 0.5 | 1.0 | 0.5 | 1.9749 |
| species_mean | 0.5 | 1.0 | 0.8 | 1.2900 |
| pigauto_no_cov | 0.5 | 1.0 | 0.8 | 2.1172 |
| pigauto_cov | 0.5 | 1.0 | 0.8 | 1.8723 |
| species_mean | 0.9 | 0.0 | 0.5 | 1.4088 |
| pigauto_no_cov | 0.9 | 0.0 | 0.5 | 1.0289 |
| pigauto_cov | 0.9 | 0.0 | 0.5 | 1.0260 |
| species_mean | 0.9 | 0.0 | 0.8 | 1.9242 |
| pigauto_no_cov | 0.9 | 0.0 | 0.8 | 1.3988 |
| pigauto_cov | 0.9 | 0.0 | 0.8 | 1.3973 |
| species_mean | 0.9 | 0.5 | 0.5 | 1.6372 |
| pigauto_no_cov | 0.9 | 0.5 | 0.5 | 1.5484 |
| pigauto_cov | 0.9 | 0.5 | 0.5 | 1.4426 |
| species_mean | 0.9 | 0.5 | 0.8 | 1.7644 |
| pigauto_no_cov | 0.9 | 0.5 | 0.8 | 1.8577 |
| pigauto_cov | 0.9 | 0.5 | 0.8 | 1.8367 |
| species_mean | 0.9 | 1.0 | 0.5 | 1.8208 |
| pigauto_no_cov | 0.9 | 1.0 | 0.5 | 2.8135 |
| pigauto_cov | 0.9 | 1.0 | 0.5 | 2.1230 |
| species_mean | 0.9 | 1.0 | 0.8 | 1.8283 |
| pigauto_no_cov | 0.9 | 1.0 | 0.8 | 2.1001 |
| pigauto_cov | 0.9 | 1.0 | 0.8 | 1.8515 |

## Observation-level Pearson r (higher is better)

| method | lambda | beta | sp_miss | pearson_r |
|--------|--------|------|---------|-----------|
| species_mean | 0.5 | 0.0 | 0.5 | 0.3781 |
| pigauto_no_cov | 0.5 | 0.0 | 0.5 | 0.6758 |
| pigauto_cov | 0.5 | 0.0 | 0.5 | 0.6761 |
| species_mean | 0.5 | 0.0 | 0.8 | 0.2043 |
| pigauto_no_cov | 0.5 | 0.0 | 0.8 | 0.7050 |
| pigauto_cov | 0.5 | 0.0 | 0.8 | 0.7062 |
| species_mean | 0.5 | 0.5 | 0.5 | 0.1267 |
| pigauto_no_cov | 0.5 | 0.5 | 0.5 | 0.2319 |
| pigauto_cov | 0.5 | 0.5 | 0.5 | 0.5377 |
| species_mean | 0.5 | 0.5 | 0.8 | 0.0599 |
| pigauto_no_cov | 0.5 | 0.5 | 0.8 | 0.2118 |
| pigauto_cov | 0.5 | 0.5 | 0.8 | 0.4599 |
| species_mean | 0.5 | 1.0 | 0.5 | 0.0729 |
| pigauto_no_cov | 0.5 | 1.0 | 0.5 | 0.0476 |
| pigauto_cov | 0.5 | 1.0 | 0.5 | 0.3913 |
| species_mean | 0.5 | 1.0 | 0.8 | -0.0051 |
| pigauto_no_cov | 0.5 | 1.0 | 0.8 | 0.1115 |
| pigauto_cov | 0.5 | 1.0 | 0.8 | 0.5514 |
| species_mean | 0.9 | 0.0 | 0.5 | 0.3523 |
| pigauto_no_cov | 0.9 | 0.0 | 0.5 | 0.6839 |
| pigauto_cov | 0.9 | 0.0 | 0.5 | 0.6854 |
| species_mean | 0.9 | 0.0 | 0.8 | 0.2048 |
| pigauto_no_cov | 0.9 | 0.0 | 0.8 | 0.6501 |
| pigauto_cov | 0.9 | 0.0 | 0.8 | 0.6500 |
| species_mean | 0.9 | 0.5 | 0.5 | 0.1570 |
| pigauto_no_cov | 0.9 | 0.5 | 0.5 | 0.3184 |
| pigauto_cov | 0.9 | 0.5 | 0.5 | 0.5338 |
| species_mean | 0.9 | 0.5 | 0.8 | 0.0676 |
| pigauto_no_cov | 0.9 | 0.5 | 0.8 | 0.2216 |
| pigauto_cov | 0.9 | 0.5 | 0.8 | 0.3997 |
| species_mean | 0.9 | 1.0 | 0.5 | 0.0644 |
| pigauto_no_cov | 0.9 | 1.0 | 0.5 | 0.0710 |
| pigauto_cov | 0.9 | 1.0 | 0.5 | 0.4441 |
| species_mean | 0.9 | 1.0 | 0.8 | 0.0448 |
| pigauto_no_cov | 0.9 | 1.0 | 0.8 | 0.1556 |
| pigauto_cov | 0.9 | 1.0 | 0.8 | 0.5605 |

## Covariate lift (pigauto_cov / pigauto_no_cov RMSE ratio)

Ratio < 1 means covariates help; ratio > 1 means they hurt.

| lambda | beta | sp_miss | RMSE_nocov | RMSE_cov | ratio |
|--------|------|---------|------------|----------|-------|
| 0.5 | 0.0 | 0.5 | 0.8941 | 0.8876 | 0.993 |
| 0.5 | 0.0 | 0.8 | 1.0684 | 1.0663 | 0.998 |
| 0.5 | 0.5 | 0.5 | 2.9220 | 2.4830 | 0.850 |
| 0.5 | 0.5 | 0.8 | 2.9385 | 2.5868 | 0.880 |
| 0.5 | 1.0 | 0.5 | 5.4872 | 4.6799 | 0.853 |
| 0.5 | 1.0 | 0.8 | 5.3110 | 4.4848 | 0.844 |
| 0.9 | 0.0 | 0.5 | 1.0914 | 1.0886 | 0.997 |
| 0.9 | 0.0 | 0.8 | 1.5510 | 1.5505 | 1.000 |
| 0.9 | 0.5 | 0.5 | 3.0138 | 2.6206 | 0.870 |
| 0.9 | 0.5 | 0.8 | 3.1889 | 2.8739 | 0.901 |
| 0.9 | 1.0 | 0.5 | 6.0190 | 4.8947 | 0.813 |
| 0.9 | 1.0 | 0.8 | 5.5354 | 4.6893 | 0.847 |

---
Generated: 2026-04-16 11:31
