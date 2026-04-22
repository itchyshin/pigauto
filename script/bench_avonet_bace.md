# AVONET x pigauto + BACE head-to-head

Mode: bundled AVONET 300. n = 300 species x 7 traits.
Seed = 2026, miss_frac = 0.30, n_imputations = 20.
**BACE skipped** (not installed or failed).

## Per-trait metrics

```
          method              trait    metric        value n_cells  wall_s
   mean_baseline               Mass      rmse 2521.4896370      90   0.006
   mean_baseline               Mass pearson_r           NA      90   0.006
   mean_baseline Beak.Length_Culmen      rmse   38.0406157      90   0.006
   mean_baseline Beak.Length_Culmen pearson_r           NA      90   0.006
   mean_baseline      Tarsus.Length      rmse   29.8727964      90   0.006
   mean_baseline      Tarsus.Length pearson_r           NA      90   0.006
   mean_baseline        Wing.Length      rmse  120.2149687      90   0.006
   mean_baseline        Wing.Length pearson_r           NA      90   0.006
   mean_baseline      Trophic.Level  accuracy    0.4666667      90   0.006
   mean_baseline  Primary.Lifestyle  accuracy    0.5777778      90   0.006
   mean_baseline          Migration  accuracy    0.7666667      90   0.006
 pigauto_default               Mass      rmse 1716.2190157      90 220.420
 pigauto_default               Mass pearson_r    0.9766346      90 220.420
 pigauto_default Beak.Length_Culmen      rmse   18.1003632      90 220.420
 pigauto_default Beak.Length_Culmen pearson_r    0.8812052      90 220.420
 pigauto_default      Tarsus.Length      rmse   14.8444306      90 220.420
 pigauto_default      Tarsus.Length pearson_r    0.8838990      90 220.420
 pigauto_default        Wing.Length      rmse   37.5069502      90 220.420
 pigauto_default        Wing.Length pearson_r    0.9501292      90 220.420
 pigauto_default      Trophic.Level  accuracy    0.7222222      90 220.420
 pigauto_default  Primary.Lifestyle  accuracy    0.7555556      90 220.420
 pigauto_default          Migration  accuracy    0.7666667      90 220.420
```
