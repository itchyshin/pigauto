# Covariate effectiveness simulation benchmark

- Species: 300 (tree300)
- Traits: 3 continuous (simulated)
- Covariates: 3 (uncorrelated normal)
- Scenarios: 4 (lambda x beta combinations)
- Missingness: 0.2, 0.4, 0.6
- Replicates: 3
- Total wall time: 86.8 min

## Mean RMSE by scenario and method

| scenario | lambda | beta | miss_frac | method | mean_RMSE | mean_r |
|----------|--------|------|-----------|--------|-----------|--------|
| low_phylo_strong_env | 0.1 | 1.5 | 0.20 | pigauto | 3.1334 | 0.0537 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.20 | pigauto_covs | 2.8698 | 0.2037 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.40 | pigauto | 3.2004 | 0.0197 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.40 | pigauto_covs | 2.8101 | 0.2948 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.60 | pigauto | 3.1707 | 0.0166 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.60 | pigauto_covs | 2.7015 | 0.3392 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.20 | pigauto | 1.5617 | 0.0685 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.20 | pigauto_covs | 1.5092 | 0.1285 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.40 | pigauto | 1.5606 | 0.0834 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.40 | pigauto_covs | 1.4984 | 0.1648 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.60 | pigauto | 1.6245 | 0.0397 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.60 | pigauto_covs | 1.5568 | 0.1201 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.20 | pigauto | 1.4573 | 0.2186 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.20 | pigauto_covs | 1.4213 | 0.2639 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.40 | pigauto | 1.4049 | 0.2500 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.40 | pigauto_covs | 1.3622 | 0.2940 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.60 | pigauto | 1.4203 | 0.1466 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.60 | pigauto_covs | 1.3680 | 0.2175 |
| high_phylo_no_env | 0.9 | 0.0 | 0.20 | pigauto | 0.5660 | 0.5854 |
| high_phylo_no_env | 0.9 | 0.0 | 0.20 | pigauto_covs | 0.5652 | 0.5842 |
| high_phylo_no_env | 0.9 | 0.0 | 0.40 | pigauto | 0.6095 | 0.6273 |
| high_phylo_no_env | 0.9 | 0.0 | 0.40 | pigauto_covs | 0.6092 | 0.6283 |
| high_phylo_no_env | 0.9 | 0.0 | 0.60 | pigauto | 0.6399 | 0.5974 |
| high_phylo_no_env | 0.9 | 0.0 | 0.60 | pigauto_covs | 0.6392 | 0.5981 |

## Covariate lift (pigauto_covs / pigauto RMSE ratio)

| scenario | lambda | beta | miss_frac | RMSE_no_cov | RMSE_covs | ratio |
|----------|--------|------|-----------|-------------|-----------|-------|
| low_phylo_strong_env | 0.1 | 1.5 | 0.20 | 3.1334 | 2.8698 | 0.916 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.40 | 3.2004 | 2.8101 | 0.878 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.60 | 3.1707 | 2.7015 | 0.852 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.20 | 1.5617 | 1.5092 | 0.966 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.40 | 1.5606 | 1.4984 | 0.960 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.60 | 1.6245 | 1.5568 | 0.958 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.20 | 1.4573 | 1.4213 | 0.975 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.40 | 1.4049 | 1.3622 | 0.970 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.60 | 1.4203 | 1.3680 | 0.963 |
| high_phylo_no_env | 0.9 | 0.0 | 0.20 | 0.5660 | 0.5652 | 0.999 |
| high_phylo_no_env | 0.9 | 0.0 | 0.40 | 0.6095 | 0.6092 | 1.000 |
| high_phylo_no_env | 0.9 | 0.0 | 0.60 | 0.6399 | 0.6392 | 0.999 |

---
Generated:
2026-04-13 16:34
