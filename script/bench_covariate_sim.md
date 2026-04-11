# Covariate effectiveness simulation benchmark

- Species: 300 (tree300)
- Traits: 3 continuous (simulated)
- Covariates: 3 (uncorrelated normal)
- Scenarios: 4 (lambda x beta combinations)
- Missingness: 0.2, 0.4, 0.6
- Replicates: 3
- Total wall time: 91.3 min

## Mean RMSE by scenario and method

| scenario | lambda | beta | miss_frac | method | mean_RMSE | mean_r |
|----------|--------|------|-----------|--------|-----------|--------|
| low_phylo_strong_env | 0.1 | 1.5 | 0.20 | pigauto | 3.1543 | 0.0414 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.20 | pigauto_covs | 2.9053 | 0.1940 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.40 | pigauto | 3.2042 | 0.0172 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.40 | pigauto_covs | 2.8244 | 0.2841 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.60 | pigauto | 3.1681 | 0.0174 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.60 | pigauto_covs | 2.7102 | 0.3352 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.20 | pigauto | 1.5618 | 0.0667 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.20 | pigauto_covs | 1.5080 | 0.1291 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.40 | pigauto | 1.5627 | 0.0813 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.40 | pigauto_covs | 1.5006 | 0.1628 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.60 | pigauto | 1.6277 | 0.0379 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.60 | pigauto_covs | 1.5551 | 0.1214 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.20 | pigauto | 1.4569 | 0.2193 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.20 | pigauto_covs | 1.4261 | 0.2605 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.40 | pigauto | 1.4036 | 0.2507 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.40 | pigauto_covs | 1.3626 | 0.2938 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.60 | pigauto | 1.4211 | 0.1460 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.60 | pigauto_covs | 1.3676 | 0.2203 |
| high_phylo_no_env | 0.9 | 0.0 | 0.20 | pigauto | 0.5658 | 0.5853 |
| high_phylo_no_env | 0.9 | 0.0 | 0.20 | pigauto_covs | 0.5650 | 0.5842 |
| high_phylo_no_env | 0.9 | 0.0 | 0.40 | pigauto | 0.6092 | 0.6276 |
| high_phylo_no_env | 0.9 | 0.0 | 0.40 | pigauto_covs | 0.6097 | 0.6280 |
| high_phylo_no_env | 0.9 | 0.0 | 0.60 | pigauto | 0.6400 | 0.5974 |
| high_phylo_no_env | 0.9 | 0.0 | 0.60 | pigauto_covs | 0.6392 | 0.5981 |

## Covariate lift (pigauto_covs / pigauto RMSE ratio)

| scenario | lambda | beta | miss_frac | RMSE_no_cov | RMSE_covs | ratio |
|----------|--------|------|-----------|-------------|-----------|-------|
| low_phylo_strong_env | 0.1 | 1.5 | 0.20 | 3.1543 | 2.9053 | 0.921 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.40 | 3.2042 | 2.8244 | 0.881 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.60 | 3.1681 | 2.7102 | 0.855 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.20 | 1.5618 | 1.5080 | 0.966 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.40 | 1.5627 | 1.5006 | 0.960 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.60 | 1.6277 | 1.5551 | 0.955 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.20 | 1.4569 | 1.4261 | 0.979 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.40 | 1.4036 | 1.3626 | 0.971 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.60 | 1.4211 | 1.3676 | 0.962 |
| high_phylo_no_env | 0.9 | 0.0 | 0.20 | 0.5658 | 0.5650 | 0.999 |
| high_phylo_no_env | 0.9 | 0.0 | 0.40 | 0.6092 | 0.6097 | 1.001 |
| high_phylo_no_env | 0.9 | 0.0 | 0.60 | 0.6400 | 0.6392 | 0.999 |

---
Generated:
2026-04-11 07:10
