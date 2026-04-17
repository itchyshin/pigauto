# Delhey plumage lightness benchmark (covariates)

- Species: 5809 (Delhey et al. 2019)
- Traits: lightness_male, lightness_female
- Covariates: annual_mean_temperature, annual_precipitation, percent_tree_cover, mean_temperature_of_warmest_quarter, precipitation_of_warmest_quarter, midLatitude
- Methods: mean, BM_only, pigauto, pigauto_covs
- Missingness: 0.2, 0.4, 0.6
- Replicates: 3
- Total wall time: 131.4 min

## Test-set RMSE

| method | miss_frac | trait | mean_RMSE |
|--------|-----------|-------|-----------|
| BM_only | 0.20 | lightness_female | 0.7654 |
| mean | 0.20 | lightness_female | 1.0139 |
| pigauto | 0.20 | lightness_female | 0.7534 |
| pigauto_covs | 0.20 | lightness_female | 0.7514 |
| BM_only | 0.20 | lightness_male | 0.7557 |
| mean | 0.20 | lightness_male | 1.0062 |
| pigauto | 0.20 | lightness_male | 0.7473 |
| pigauto_covs | 0.20 | lightness_male | 0.7505 |
| BM_only | 0.40 | lightness_female | 0.7964 |
| mean | 0.40 | lightness_female | 0.9911 |
| pigauto | 0.40 | lightness_female | 0.7893 |
| pigauto_covs | 0.40 | lightness_female | 0.7964 |
| BM_only | 0.40 | lightness_male | 0.7747 |
| mean | 0.40 | lightness_male | 1.0239 |
| pigauto | 0.40 | lightness_male | 0.7692 |
| pigauto_covs | 0.40 | lightness_male | 0.7707 |
| BM_only | 0.60 | lightness_female | 0.7924 |
| mean | 0.60 | lightness_female | 0.9757 |
| pigauto | 0.60 | lightness_female | 0.7924 |
| pigauto_covs | 0.60 | lightness_female | 0.7828 |
| BM_only | 0.60 | lightness_male | 0.7878 |
| mean | 0.60 | lightness_male | 0.9807 |
| pigauto | 0.60 | lightness_male | 0.7807 |
| pigauto_covs | 0.60 | lightness_male | 0.7831 |

## Test-set Pearson r

| method | miss_frac | trait | mean_r |
|--------|-----------|-------|--------|
| BM_only | 0.20 | lightness_female | 0.6646 |
| pigauto | 0.20 | lightness_female | 0.6744 |
| pigauto_covs | 0.20 | lightness_female | 0.6726 |
| BM_only | 0.20 | lightness_male | 0.6692 |
| pigauto | 0.20 | lightness_male | 0.6760 |
| pigauto_covs | 0.20 | lightness_male | 0.6749 |
| BM_only | 0.40 | lightness_female | 0.6130 |
| pigauto | 0.40 | lightness_female | 0.6163 |
| pigauto_covs | 0.40 | lightness_female | 0.6130 |
| BM_only | 0.40 | lightness_male | 0.6594 |
| pigauto | 0.40 | lightness_male | 0.6625 |
| pigauto_covs | 0.40 | lightness_male | 0.6623 |
| BM_only | 0.60 | lightness_female | 0.6040 |
| pigauto | 0.60 | lightness_female | 0.6040 |
| pigauto_covs | 0.60 | lightness_female | 0.6048 |
| BM_only | 0.60 | lightness_male | 0.6154 |
| pigauto | 0.60 | lightness_male | 0.6209 |
| pigauto_covs | 0.60 | lightness_male | 0.6177 |

## Covariate lift (pigauto_covs / pigauto RMSE ratio)

| miss_frac | trait | RMSE_pigauto | RMSE_covs | ratio |
|-----------|-------|--------------|-----------|-------|
| 0.20 | lightness_female | 0.7534 | 0.7514 | 0.997 |
| 0.20 | lightness_male | 0.7473 | 0.7505 | 1.004 |
| 0.40 | lightness_female | 0.7893 | 0.7964 | 1.009 |
| 0.40 | lightness_male | 0.7692 | 0.7707 | 1.002 |
| 0.60 | lightness_female | 0.7924 | 0.7828 | 0.988 |
| 0.60 | lightness_male | 0.7807 | 0.7831 | 1.003 |

---
Generated:
2026-04-14 08:05
