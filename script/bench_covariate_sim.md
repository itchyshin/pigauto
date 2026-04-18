# Covariate effectiveness simulation benchmark

- Species: 300 (tree300)
- Traits: 3 continuous (simulated)
- Covariates: 3 (uncorrelated normal)
- Scenarios: 4 (lambda x beta combinations)
- Missingness: 0.2, 0.4, 0.6
- Replicates: 3
- Total wall time: 119.6 min

## Mean RMSE by scenario and method

| scenario | lambda | beta | miss_frac | method | mean_RMSE | mean_r |
|----------|--------|------|-----------|--------|-----------|--------|
| low_phylo_strong_env | 0.1 | 1.5 | 0.20 | pigauto | 2.7136 | 0.4208 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.20 | pigauto_covs | 2.5046 | 0.4964 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.40 | pigauto | 2.8370 | 0.3365 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.40 | pigauto_covs | 2.5128 | 0.4887 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.60 | pigauto | 2.8860 | 0.3274 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.60 | pigauto_covs | 2.5846 | 0.4782 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.20 | pigauto | 1.5049 | 0.2497 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.20 | pigauto_covs | 1.4635 | 0.3010 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.40 | pigauto | 1.5298 | 0.2268 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.40 | pigauto_covs | 1.4739 | 0.2782 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.60 | pigauto | 1.6484 | 0.1464 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.60 | pigauto_covs | 1.6004 | 0.2008 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.20 | pigauto | 1.4887 | 0.2104 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.20 | pigauto_covs | 1.4369 | 0.2773 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.40 | pigauto | 1.4657 | 0.2189 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.40 | pigauto_covs | 1.3757 | 0.2858 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.60 | pigauto | 1.4837 | 0.1329 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.60 | pigauto_covs | 1.4317 | 0.2012 |
| high_phylo_no_env | 0.9 | 0.0 | 0.20 | pigauto | 0.5716 | 0.5804 |
| high_phylo_no_env | 0.9 | 0.0 | 0.20 | pigauto_covs | 0.5695 | 0.5809 |
| high_phylo_no_env | 0.9 | 0.0 | 0.40 | pigauto | 0.6105 | 0.6280 |
| high_phylo_no_env | 0.9 | 0.0 | 0.40 | pigauto_covs | 0.6085 | 0.6284 |
| high_phylo_no_env | 0.9 | 0.0 | 0.60 | pigauto | 0.6476 | 0.5886 |
| high_phylo_no_env | 0.9 | 0.0 | 0.60 | pigauto_covs | 0.6459 | 0.5888 |

## Covariate lift (pigauto_covs / pigauto RMSE ratio)

| scenario | lambda | beta | miss_frac | RMSE_no_cov | RMSE_covs | ratio |
|----------|--------|------|-----------|-------------|-----------|-------|
| low_phylo_strong_env | 0.1 | 1.5 | 0.20 | 2.7136 | 2.5046 | 0.923 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.40 | 2.8370 | 2.5128 | 0.886 |
| low_phylo_strong_env | 0.1 | 1.5 | 0.60 | 2.8860 | 2.5846 | 0.896 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.20 | 1.5049 | 1.4635 | 0.972 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.40 | 1.5298 | 1.4739 | 0.963 |
| mod_phylo_strong_env | 0.3 | 1.0 | 0.60 | 1.6484 | 1.6004 | 0.971 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.20 | 1.4887 | 1.4369 | 0.965 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.40 | 1.4657 | 1.3757 | 0.939 |
| high_phylo_mod_env | 0.7 | 0.5 | 0.60 | 1.4837 | 1.4317 | 0.965 |
| high_phylo_no_env | 0.9 | 0.0 | 0.20 | 0.5716 | 0.5695 | 0.996 |
| high_phylo_no_env | 0.9 | 0.0 | 0.40 | 0.6105 | 0.6085 | 0.997 |
| high_phylo_no_env | 0.9 | 0.0 | 0.60 | 0.6476 | 0.6459 | 0.997 |

---
Generated:
2026-04-18 08:01
