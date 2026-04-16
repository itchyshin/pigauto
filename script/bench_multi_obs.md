# Multi-observation-per-species benchmark

Run on: 2026-04-16 08:32:43
Machine: Darwin 25.4.0 (arm64), R 4.5.2
Species: 200, reps: 2, epochs: 200
Total wall: 51.1 min

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
| pigauto_no_cov | 0.5 | 0.0 | 0.5 | 0.8897 |
| pigauto_cov | 0.5 | 0.0 | 0.5 | 0.8930 |
| species_mean | 0.5 | 0.0 | 0.8 | 1.4126 |
| pigauto_no_cov | 0.5 | 0.0 | 0.8 | 1.0756 |
| pigauto_cov | 0.5 | 0.0 | 0.8 | 1.0663 |
| species_mean | 0.5 | 0.5 | 0.5 | 2.9272 |
| pigauto_no_cov | 0.5 | 0.5 | 0.5 | 2.9458 |
| pigauto_cov | 0.5 | 0.5 | 0.5 | 2.4665 |
| species_mean | 0.5 | 0.5 | 0.8 | 2.9247 |
| pigauto_no_cov | 0.5 | 0.5 | 0.8 | 2.9210 |
| pigauto_cov | 0.5 | 0.5 | 0.8 | 2.6101 |
| species_mean | 0.5 | 1.0 | 0.5 | 5.1071 |
| pigauto_no_cov | 0.5 | 1.0 | 0.5 | 5.4940 |
| pigauto_cov | 0.5 | 1.0 | 0.5 | 4.6165 |
| species_mean | 0.5 | 1.0 | 0.8 | 5.2278 |
| pigauto_no_cov | 0.5 | 1.0 | 0.8 | 5.3326 |
| pigauto_cov | 0.5 | 1.0 | 0.8 | 4.5332 |
| species_mean | 0.9 | 0.0 | 0.5 | 1.3816 |
| pigauto_no_cov | 0.9 | 0.0 | 0.5 | 1.0909 |
| pigauto_cov | 0.9 | 0.0 | 0.5 | 1.0851 |
| species_mean | 0.9 | 0.0 | 0.8 | 1.9755 |
| pigauto_no_cov | 0.9 | 0.0 | 0.8 | 1.5475 |
| pigauto_cov | 0.9 | 0.0 | 0.8 | 1.5569 |
| species_mean | 0.9 | 0.5 | 0.5 | 3.0603 |
| pigauto_no_cov | 0.9 | 0.5 | 0.5 | 2.9619 |
| pigauto_cov | 0.9 | 0.5 | 0.5 | 2.6099 |
| species_mean | 0.9 | 0.5 | 0.8 | 3.0626 |
| pigauto_no_cov | 0.9 | 0.5 | 0.8 | 3.1975 |
| pigauto_cov | 0.9 | 0.5 | 0.8 | 2.8575 |
| species_mean | 0.9 | 1.0 | 0.5 | 5.5368 |
| pigauto_no_cov | 0.9 | 1.0 | 0.5 | 5.8925 |
| pigauto_cov | 0.9 | 1.0 | 0.5 | 4.8986 |
| species_mean | 0.9 | 1.0 | 0.8 | 5.4523 |
| pigauto_no_cov | 0.9 | 1.0 | 0.8 | 5.5057 |
| pigauto_cov | 0.9 | 1.0 | 0.8 | 4.7083 |

## Species-level RMSE (unobserved species, lower is better)

| method | lambda | beta | sp_miss | sp_RMSE |
|--------|--------|------|---------|---------|
| species_mean | 0.5 | 0.0 | 0.5 | 1.1091 |
| pigauto_no_cov | 0.5 | 0.0 | 0.5 | 0.7790 |
| pigauto_cov | 0.5 | 0.0 | 0.5 | 0.7826 |
| species_mean | 0.5 | 0.0 | 0.8 | 1.3434 |
| pigauto_no_cov | 0.5 | 0.0 | 0.8 | 0.9418 |
| pigauto_cov | 0.5 | 0.0 | 0.8 | 0.9318 |
| species_mean | 0.5 | 0.5 | 0.5 | 1.2797 |
| pigauto_no_cov | 0.5 | 0.5 | 0.5 | 1.2879 |
| pigauto_cov | 0.5 | 0.5 | 0.5 | 1.0864 |
| species_mean | 0.5 | 0.5 | 0.8 | 1.4960 |
| pigauto_no_cov | 0.5 | 0.5 | 0.8 | 1.3945 |
| pigauto_cov | 0.5 | 0.5 | 0.8 | 1.4040 |
| species_mean | 0.5 | 1.0 | 0.5 | 1.3076 |
| pigauto_no_cov | 0.5 | 1.0 | 0.5 | 2.3073 |
| pigauto_cov | 0.5 | 1.0 | 0.5 | 1.9697 |
| species_mean | 0.5 | 1.0 | 0.8 | 1.2900 |
| pigauto_no_cov | 0.5 | 1.0 | 0.8 | 2.0638 |
| pigauto_cov | 0.5 | 1.0 | 0.8 | 1.9612 |
| species_mean | 0.9 | 0.0 | 0.5 | 1.4088 |
| pigauto_no_cov | 0.9 | 0.0 | 0.5 | 1.0280 |
| pigauto_cov | 0.9 | 0.0 | 0.5 | 1.0226 |
| species_mean | 0.9 | 0.0 | 0.8 | 1.9242 |
| pigauto_no_cov | 0.9 | 0.0 | 0.8 | 1.3945 |
| pigauto_cov | 0.9 | 0.0 | 0.8 | 1.4026 |
| species_mean | 0.9 | 0.5 | 0.5 | 1.6372 |
| pigauto_no_cov | 0.9 | 0.5 | 0.5 | 1.4687 |
| pigauto_cov | 0.9 | 0.5 | 0.5 | 1.4177 |
| species_mean | 0.9 | 0.5 | 0.8 | 1.7644 |
| pigauto_no_cov | 0.9 | 0.5 | 0.8 | 1.8658 |
| pigauto_cov | 0.9 | 0.5 | 0.8 | 1.8333 |
| species_mean | 0.9 | 1.0 | 0.5 | 1.8208 |
| pigauto_no_cov | 0.9 | 1.0 | 0.5 | 2.6505 |
| pigauto_cov | 0.9 | 1.0 | 0.5 | 2.1212 |
| species_mean | 0.9 | 1.0 | 0.8 | 1.8283 |
| pigauto_no_cov | 0.9 | 1.0 | 0.8 | 1.9914 |
| pigauto_cov | 0.9 | 1.0 | 0.8 | 1.8952 |

## Observation-level Pearson r (higher is better)

| method | lambda | beta | sp_miss | pearson_r |
|--------|--------|------|---------|-----------|
| species_mean | 0.5 | 0.0 | 0.5 | 0.3781 |
| pigauto_no_cov | 0.5 | 0.0 | 0.5 | 0.6749 |
| pigauto_cov | 0.5 | 0.0 | 0.5 | 0.6753 |
| species_mean | 0.5 | 0.0 | 0.8 | 0.2043 |
| pigauto_no_cov | 0.5 | 0.0 | 0.8 | 0.7049 |
| pigauto_cov | 0.5 | 0.0 | 0.8 | 0.7060 |
| species_mean | 0.5 | 0.5 | 0.5 | 0.1267 |
| pigauto_no_cov | 0.5 | 0.5 | 0.5 | 0.2194 |
| pigauto_cov | 0.5 | 0.5 | 0.5 | 0.5447 |
| species_mean | 0.5 | 0.5 | 0.8 | 0.0599 |
| pigauto_no_cov | 0.5 | 0.5 | 0.8 | 0.2214 |
| pigauto_cov | 0.5 | 0.5 | 0.8 | 0.4519 |
| species_mean | 0.5 | 1.0 | 0.5 | 0.0729 |
| pigauto_no_cov | 0.5 | 1.0 | 0.5 | 0.0641 |
| pigauto_cov | 0.5 | 1.0 | 0.5 | 0.4196 |
| species_mean | 0.5 | 1.0 | 0.8 | -0.0051 |
| pigauto_no_cov | 0.5 | 1.0 | 0.8 | 0.1004 |
| pigauto_cov | 0.5 | 1.0 | 0.8 | 0.5403 |
| species_mean | 0.9 | 0.0 | 0.5 | 0.3523 |
| pigauto_no_cov | 0.9 | 0.0 | 0.5 | 0.6831 |
| pigauto_cov | 0.9 | 0.0 | 0.5 | 0.6856 |
| species_mean | 0.9 | 0.0 | 0.8 | 0.2048 |
| pigauto_no_cov | 0.9 | 0.0 | 0.8 | 0.6508 |
| pigauto_cov | 0.9 | 0.0 | 0.8 | 0.6487 |
| species_mean | 0.9 | 0.5 | 0.5 | 0.1570 |
| pigauto_no_cov | 0.9 | 0.5 | 0.5 | 0.3555 |
| pigauto_cov | 0.9 | 0.5 | 0.5 | 0.5306 |
| species_mean | 0.9 | 0.5 | 0.8 | 0.0676 |
| pigauto_no_cov | 0.9 | 0.5 | 0.8 | 0.2226 |
| pigauto_cov | 0.9 | 0.5 | 0.8 | 0.4065 |
| species_mean | 0.9 | 1.0 | 0.5 | 0.0644 |
| pigauto_no_cov | 0.9 | 1.0 | 0.5 | 0.0875 |
| pigauto_cov | 0.9 | 1.0 | 0.5 | 0.4437 |
| species_mean | 0.9 | 1.0 | 0.8 | 0.0448 |
| pigauto_no_cov | 0.9 | 1.0 | 0.8 | 0.1399 |
| pigauto_cov | 0.9 | 1.0 | 0.8 | 0.5462 |

## Covariate lift (pigauto_cov / pigauto_no_cov RMSE ratio)

Ratio < 1 means covariates help; ratio > 1 means they hurt.

| lambda | beta | sp_miss | RMSE_nocov | RMSE_cov | ratio |
|--------|------|---------|------------|----------|-------|
| 0.5 | 0.0 | 0.5 | 0.8897 | 0.8930 | 1.004 |
| 0.5 | 0.0 | 0.8 | 1.0756 | 1.0663 | 0.991 |
| 0.5 | 0.5 | 0.5 | 2.9458 | 2.4665 | 0.837 |
| 0.5 | 0.5 | 0.8 | 2.9210 | 2.6101 | 0.894 |
| 0.5 | 1.0 | 0.5 | 5.4940 | 4.6165 | 0.840 |
| 0.5 | 1.0 | 0.8 | 5.3326 | 4.5332 | 0.850 |
| 0.9 | 0.0 | 0.5 | 1.0909 | 1.0851 | 0.995 |
| 0.9 | 0.0 | 0.8 | 1.5475 | 1.5569 | 1.006 |
| 0.9 | 0.5 | 0.5 | 2.9619 | 2.6099 | 0.881 |
| 0.9 | 0.5 | 0.8 | 3.1975 | 2.8575 | 0.894 |
| 0.9 | 1.0 | 0.5 | 5.8925 | 4.8986 | 0.831 |
| 0.9 | 1.0 | 0.8 | 5.5057 | 4.7083 | 0.855 |

---
Generated: 2026-04-16 08:32
