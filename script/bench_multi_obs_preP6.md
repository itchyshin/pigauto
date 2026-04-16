# Multi-observation-per-species benchmark

Run on: 2026-04-14 06:37:34
Machine: Darwin 25.4.0 (arm64), R 4.5.2
Species: 200, reps: 2, epochs: 200
Total wall: 43.8 min

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
| pigauto_no_cov | 0.5 | 0.0 | 0.5 | 0.8916 |
| pigauto_cov | 0.5 | 0.0 | 0.5 | 0.8899 |
| species_mean | 0.5 | 0.0 | 0.8 | 1.4126 |
| pigauto_no_cov | 0.5 | 0.0 | 0.8 | 1.0670 |
| pigauto_cov | 0.5 | 0.0 | 0.8 | 1.0743 |
| species_mean | 0.5 | 0.5 | 0.5 | 2.9272 |
| pigauto_no_cov | 0.5 | 0.5 | 0.5 | 2.8634 |
| pigauto_cov | 0.5 | 0.5 | 0.5 | 2.4807 |
| species_mean | 0.5 | 0.5 | 0.8 | 2.9247 |
| pigauto_no_cov | 0.5 | 0.5 | 0.8 | 2.9355 |
| pigauto_cov | 0.5 | 0.5 | 0.8 | 2.6046 |
| species_mean | 0.5 | 1.0 | 0.5 | 5.1071 |
| pigauto_no_cov | 0.5 | 1.0 | 0.5 | 5.4333 |
| pigauto_cov | 0.5 | 1.0 | 0.5 | 4.5837 |
| species_mean | 0.5 | 1.0 | 0.8 | 5.2278 |
| pigauto_no_cov | 0.5 | 1.0 | 0.8 | 5.2851 |
| pigauto_cov | 0.5 | 1.0 | 0.8 | 4.4910 |
| species_mean | 0.9 | 0.0 | 0.5 | 1.3816 |
| pigauto_no_cov | 0.9 | 0.0 | 0.5 | 1.0903 |
| pigauto_cov | 0.9 | 0.0 | 0.5 | 1.0887 |
| species_mean | 0.9 | 0.0 | 0.8 | 1.9755 |
| pigauto_no_cov | 0.9 | 0.0 | 0.8 | 1.5485 |
| pigauto_cov | 0.9 | 0.0 | 0.8 | 1.5498 |
| species_mean | 0.9 | 0.5 | 0.5 | 3.0603 |
| pigauto_no_cov | 0.9 | 0.5 | 0.5 | 2.9719 |
| pigauto_cov | 0.9 | 0.5 | 0.5 | 2.6190 |
| species_mean | 0.9 | 0.5 | 0.8 | 3.0626 |
| pigauto_no_cov | 0.9 | 0.5 | 0.8 | 3.1822 |
| pigauto_cov | 0.9 | 0.5 | 0.8 | 2.8723 |
| species_mean | 0.9 | 1.0 | 0.5 | 5.5368 |
| pigauto_no_cov | 0.9 | 1.0 | 0.5 | 5.7517 |
| pigauto_cov | 0.9 | 1.0 | 0.5 | 4.9428 |
| species_mean | 0.9 | 1.0 | 0.8 | 5.4523 |
| pigauto_no_cov | 0.9 | 1.0 | 0.8 | 5.4978 |
| pigauto_cov | 0.9 | 1.0 | 0.8 | 4.7126 |

## Species-level RMSE (unobserved species, lower is better)

| method | lambda | beta | sp_miss | sp_RMSE |
|--------|--------|------|---------|---------|
| species_mean | 0.5 | 0.0 | 0.5 | 1.1091 |
| pigauto_no_cov | 0.5 | 0.0 | 0.5 | 0.7806 |
| pigauto_cov | 0.5 | 0.0 | 0.5 | 0.7783 |
| species_mean | 0.5 | 0.0 | 0.8 | 1.3434 |
| pigauto_no_cov | 0.5 | 0.0 | 0.8 | 0.9321 |
| pigauto_cov | 0.5 | 0.0 | 0.8 | 0.9409 |
| species_mean | 0.5 | 0.5 | 0.5 | 1.2797 |
| pigauto_no_cov | 0.5 | 0.5 | 0.5 | 1.1267 |
| pigauto_cov | 0.5 | 0.5 | 0.5 | 1.0864 |
| species_mean | 0.5 | 0.5 | 0.8 | 1.4960 |
| pigauto_no_cov | 0.5 | 0.5 | 0.8 | 1.4140 |
| pigauto_cov | 0.5 | 0.5 | 0.8 | 1.3946 |
| species_mean | 0.5 | 1.0 | 0.5 | 1.3076 |
| pigauto_no_cov | 0.5 | 1.0 | 0.5 | 2.1586 |
| pigauto_cov | 0.5 | 1.0 | 0.5 | 1.9306 |
| species_mean | 0.5 | 1.0 | 0.8 | 1.2900 |
| pigauto_no_cov | 0.5 | 1.0 | 0.8 | 1.9882 |
| pigauto_cov | 0.5 | 1.0 | 0.8 | 1.9025 |
| species_mean | 0.9 | 0.0 | 0.5 | 1.4088 |
| pigauto_no_cov | 0.9 | 0.0 | 0.5 | 1.0287 |
| pigauto_cov | 0.9 | 0.0 | 0.5 | 1.0273 |
| species_mean | 0.9 | 0.0 | 0.8 | 1.9242 |
| pigauto_no_cov | 0.9 | 0.0 | 0.8 | 1.3960 |
| pigauto_cov | 0.9 | 0.0 | 0.8 | 1.3961 |
| species_mean | 0.9 | 0.5 | 0.5 | 1.6372 |
| pigauto_no_cov | 0.9 | 0.5 | 0.5 | 1.4662 |
| pigauto_cov | 0.9 | 0.5 | 0.5 | 1.4273 |
| species_mean | 0.9 | 0.5 | 0.8 | 1.7644 |
| pigauto_no_cov | 0.9 | 0.5 | 0.8 | 1.8506 |
| pigauto_cov | 0.9 | 0.5 | 0.8 | 1.8194 |
| species_mean | 0.9 | 1.0 | 0.5 | 1.8208 |
| pigauto_no_cov | 0.9 | 1.0 | 0.5 | 2.3095 |
| pigauto_cov | 0.9 | 1.0 | 0.5 | 2.1465 |
| species_mean | 0.9 | 1.0 | 0.8 | 1.8283 |
| pigauto_no_cov | 0.9 | 1.0 | 0.8 | 1.9981 |
| pigauto_cov | 0.9 | 1.0 | 0.8 | 1.8830 |

## Observation-level Pearson r (higher is better)

| method | lambda | beta | sp_miss | pearson_r |
|--------|--------|------|---------|-----------|
| species_mean | 0.5 | 0.0 | 0.5 | 0.3781 |
| pigauto_no_cov | 0.5 | 0.0 | 0.5 | 0.6753 |
| pigauto_cov | 0.5 | 0.0 | 0.5 | 0.6754 |
| species_mean | 0.5 | 0.0 | 0.8 | 0.2043 |
| pigauto_no_cov | 0.5 | 0.0 | 0.8 | 0.7062 |
| pigauto_cov | 0.5 | 0.0 | 0.8 | 0.7068 |
| species_mean | 0.5 | 0.5 | 0.5 | 0.1267 |
| pigauto_no_cov | 0.5 | 0.5 | 0.5 | 0.2765 |
| pigauto_cov | 0.5 | 0.5 | 0.5 | 0.5375 |
| species_mean | 0.5 | 0.5 | 0.8 | 0.0599 |
| pigauto_no_cov | 0.5 | 0.5 | 0.8 | 0.2113 |
| pigauto_cov | 0.5 | 0.5 | 0.8 | 0.4485 |
| species_mean | 0.5 | 1.0 | 0.5 | 0.0729 |
| pigauto_no_cov | 0.5 | 1.0 | 0.5 | 0.0565 |
| pigauto_cov | 0.5 | 1.0 | 0.5 | 0.4310 |
| species_mean | 0.5 | 1.0 | 0.8 | -0.0051 |
| pigauto_no_cov | 0.5 | 1.0 | 0.8 | 0.1097 |
| pigauto_cov | 0.5 | 1.0 | 0.8 | 0.5535 |
| species_mean | 0.9 | 0.0 | 0.5 | 0.3523 |
| pigauto_no_cov | 0.9 | 0.0 | 0.5 | 0.6846 |
| pigauto_cov | 0.9 | 0.0 | 0.5 | 0.6855 |
| species_mean | 0.9 | 0.0 | 0.8 | 0.2048 |
| pigauto_no_cov | 0.9 | 0.0 | 0.8 | 0.6514 |
| pigauto_cov | 0.9 | 0.0 | 0.8 | 0.6509 |
| species_mean | 0.9 | 0.5 | 0.5 | 0.1570 |
| pigauto_no_cov | 0.9 | 0.5 | 0.5 | 0.3531 |
| pigauto_cov | 0.9 | 0.5 | 0.5 | 0.5290 |
| species_mean | 0.9 | 0.5 | 0.8 | 0.0676 |
| pigauto_no_cov | 0.9 | 0.5 | 0.8 | 0.2221 |
| pigauto_cov | 0.9 | 0.5 | 0.8 | 0.3972 |
| species_mean | 0.9 | 1.0 | 0.5 | 0.0644 |
| pigauto_no_cov | 0.9 | 1.0 | 0.5 | 0.1069 |
| pigauto_cov | 0.9 | 1.0 | 0.5 | 0.4274 |
| species_mean | 0.9 | 1.0 | 0.8 | 0.0448 |
| pigauto_no_cov | 0.9 | 1.0 | 0.8 | 0.1544 |
| pigauto_cov | 0.9 | 1.0 | 0.8 | 0.5524 |

## Covariate lift (pigauto_cov / pigauto_no_cov RMSE ratio)

Ratio < 1 means covariates help; ratio > 1 means they hurt.

| lambda | beta | sp_miss | RMSE_nocov | RMSE_cov | ratio |
|--------|------|---------|------------|----------|-------|
| 0.5 | 0.0 | 0.5 | 0.8916 | 0.8899 | 0.998 |
| 0.5 | 0.0 | 0.8 | 1.0670 | 1.0743 | 1.007 |
| 0.5 | 0.5 | 0.5 | 2.8634 | 2.4807 | 0.866 |
| 0.5 | 0.5 | 0.8 | 2.9355 | 2.6046 | 0.887 |
| 0.5 | 1.0 | 0.5 | 5.4333 | 4.5837 | 0.844 |
| 0.5 | 1.0 | 0.8 | 5.2851 | 4.4910 | 0.850 |
| 0.9 | 0.0 | 0.5 | 1.0903 | 1.0887 | 0.998 |
| 0.9 | 0.0 | 0.8 | 1.5485 | 1.5498 | 1.001 |
| 0.9 | 0.5 | 0.5 | 2.9719 | 2.6190 | 0.881 |
| 0.9 | 0.5 | 0.8 | 3.1822 | 2.8723 | 0.903 |
| 0.9 | 1.0 | 0.5 | 5.7517 | 4.9428 | 0.859 |
| 0.9 | 1.0 | 0.8 | 5.4978 | 4.7126 | 0.857 |

---
Generated: 2026-04-14 06:37
