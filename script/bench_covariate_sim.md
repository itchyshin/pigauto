# Covariate effectiveness simulation benchmark

- Species: 300 (tree300)
- Traits: 3 continuous (simulated)
- Covariates: 3 (uncorrelated normal)
- Scenarios: 4 (lambda x beta combinations)
- Missingness: 0.2, 0.4, 0.6
- Replicates: 3
- Total wall time: 110.8 min

## Mean RMSE by scenario and method

| scenario | lambda | beta | miss_frac | method | mean_RMSE | mean_r |
|----------|--------|------|-----------|--------|-----------|--------|
| low_phylo_strong_env | 0.1 | 1.5 | 0.20 | pigauto | 2.7143 | 0.4209 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.20 | pigauto_covs | 2.4950 | 0.4897 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.40 | pigauto | 2.8274 | 0.3367 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.40 | pigauto_covs | 2.4769 | 0.4927 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.60 | pigauto | 2.8496 | 0.3287 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.60 | pigauto_covs | 2.6084 | 0.4551 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.20 | pigauto | 1.4934 | 0.2491 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.20 | pigauto_covs | 1.4614 | 0.2895 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.40 | pigauto | 1.5037 | 0.2245 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.40 | pigauto_covs | 1.4675 | 0.2699 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.60 | pigauto | 1.6106 | 0.1441 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.60 | pigauto_covs | 1.5678 | 0.2170 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.20 | pigauto | 1.4701 | 0.2079 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.20 | pigauto_covs | 1.4522 | 0.2531 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.40 | pigauto | 1.4458 | 0.2191 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.40 | pigauto_covs | 1.3817 | 0.2740 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.60 | pigauto | 1.4700 | 0.1340 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.60 | pigauto_covs | 1.3996 | 0.2238 |
| high_phylo_no_env | 0.9 | 0.0 | 0.20 | pigauto | 0.5678 | 0.5804 |
| high_phylo_no_env | 0.9 | 0.0 | 0.20 | pigauto_covs | 0.5660 | 0.5815 |
| high_phylo_no_env | 0.9 | 0.0 | 0.40 | pigauto | 0.6105 | 0.6272 |
| high_phylo_no_env | 0.9 | 0.0 | 0.40 | pigauto_covs | 0.6090 | 0.6274 |
| high_phylo_no_env | 0.9 | 0.0 | 0.60 | pigauto | 0.6460 | 0.5884 |
| high_phylo_no_env | 0.9 | 0.0 | 0.60 | pigauto_covs | 0.6461 | 0.5900 |

## Covariate lift (pigauto_covs / pigauto RMSE ratio)

| scenario | lambda | beta | miss_frac | RMSE_no_cov | RMSE_covs | ratio |
|----------|--------|------|-----------|-------------|-----------|-------|
| low_phylo_strong_env | 0.1 | 1.5 | 0.20 | 2.7143 | 2.4950 | 0.919 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.40 | 2.8274 | 2.4769 | 0.876 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.60 | 2.8496 | 2.6084 | 0.915 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.20 | 1.4934 | 1.4614 | 0.979 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.40 | 1.5037 | 1.4675 | 0.976 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.60 | 1.6106 | 1.5678 | 0.973 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.20 | 1.4701 | 1.4522 | 0.988 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.40 | 1.4458 | 1.3817 | 0.956 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.60 | 1.4700 | 1.3996 | 0.952 |
| high_phylo_no_env | 0.9 | 0.0 | 0.20 | 0.5678 | 0.5660 | 0.997 |
| high_phylo_no_env | 0.9 | 0.0 | 0.40 | 0.6105 | 0.6090 | 0.998 |
| high_phylo_no_env | 0.9 | 0.0 | 0.60 | 0.6460 | 0.6461 | 1.000 |

---
Generated:
2026-04-16 09:32
