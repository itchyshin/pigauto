# Delhey plumage lightness benchmark (covariates)

- Species: 5809 (Delhey et al. 2019)
- Traits: lightness_male, lightness_female
- Covariates: annual_mean_temperature, annual_precipitation, percent_tree_cover, mean_temperature_of_warmest_quarter, precipitation_of_warmest_quarter, midLatitude
- Methods: mean, BM_only, pigauto, pigauto_covs
- Missingness: 0.2, 0.4, 0.6
- Replicates: 3
- Total wall time: 107.7 min

## Test-set RMSE

| method | miss_frac | trait | mean_RMSE |
|--------|-----------|-------|-----------|
| BM_only | 0.20 | lightness_female | 0.6259 |
| mean | 0.20 | lightness_female | 1.0139 |
| pigauto | 0.20 | lightness_female | 0.6250 |
| pigauto_covs | 0.20 | lightness_female | 0.6252 |
| BM_only | 0.20 | lightness_male | 0.5961 |
| mean | 0.20 | lightness_male | 1.0062 |
| pigauto | 0.20 | lightness_male | 0.5955 |
| pigauto_covs | 0.20 | lightness_male | 0.5955 |
| BM_only | 0.40 | lightness_female | 0.7063 |
| mean | 0.40 | lightness_female | 0.9911 |
| pigauto | 0.40 | lightness_female | 0.7056 |
| pigauto_covs | 0.40 | lightness_female | 0.7063 |
| BM_only | 0.40 | lightness_male | 0.6473 |
| mean | 0.40 | lightness_male | 1.0239 |
| pigauto | 0.40 | lightness_male | 0.6498 |
| pigauto_covs | 0.40 | lightness_male | 0.6473 |
| BM_only | 0.60 | lightness_female | 0.7273 |
| mean | 0.60 | lightness_female | 0.9757 |
| pigauto | 0.60 | lightness_female | 0.7263 |
| pigauto_covs | 0.60 | lightness_female | 0.7264 |
| BM_only | 0.60 | lightness_male | 0.6943 |
| mean | 0.60 | lightness_male | 0.9807 |
| pigauto | 0.60 | lightness_male | 0.6943 |
| pigauto_covs | 0.60 | lightness_male | 0.6943 |

## Test-set Pearson r

| method | miss_frac | trait | mean_r |
|--------|-----------|-------|--------|
| BM_only | 0.20 | lightness_female | 0.7890 |
| pigauto | 0.20 | lightness_female | 0.7905 |
| pigauto_covs | 0.20 | lightness_female | 0.7884 |
| BM_only | 0.20 | lightness_male | 0.8067 |
| pigauto | 0.20 | lightness_male | 0.8069 |
| pigauto_covs | 0.20 | lightness_male | 0.8064 |
| BM_only | 0.40 | lightness_female | 0.7120 |
| pigauto | 0.40 | lightness_female | 0.7113 |
| pigauto_covs | 0.40 | lightness_female | 0.7120 |
| BM_only | 0.40 | lightness_male | 0.7759 |
| pigauto | 0.40 | lightness_male | 0.7752 |
| pigauto_covs | 0.40 | lightness_male | 0.7759 |
| BM_only | 0.60 | lightness_female | 0.6842 |
| pigauto | 0.60 | lightness_female | 0.6840 |
| pigauto_covs | 0.60 | lightness_female | 0.6836 |
| BM_only | 0.60 | lightness_male | 0.7151 |
| pigauto | 0.60 | lightness_male | 0.7151 |
| pigauto_covs | 0.60 | lightness_male | 0.7151 |

## Covariate lift (pigauto_covs / pigauto RMSE ratio)

| miss_frac | trait | RMSE_pigauto | RMSE_covs | ratio |
|-----------|-------|--------------|-----------|-------|
| 0.20 | lightness_female | 0.6250 | 0.6252 | 1.000 |
| 0.20 | lightness_male | 0.5955 | 0.5955 | 1.000 |
| 0.40 | lightness_female | 0.7056 | 0.7063 | 1.001 |
| 0.40 | lightness_male | 0.6498 | 0.6473 | 0.996 |
| 0.60 | lightness_female | 0.7263 | 0.7264 | 1.000 |
| 0.60 | lightness_male | 0.6943 | 0.6943 | 1.000 |

---
Generated:
2026-04-10 12:15
