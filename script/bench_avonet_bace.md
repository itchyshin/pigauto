# AVONET x pigauto + BACE head-to-head

Mode: bundled AVONET 300. n = 300 species x 7 traits.
Seed = 2026, miss_frac = 0.30, n_imputations = 20.


## Per-trait metrics

```
          method              trait    metric        value n_cells  wall_s
   mean_baseline               Mass      rmse 2621.8216613      90   0.014
   mean_baseline               Mass pearson_r           NA      90   0.014
   mean_baseline Beak.Length_Culmen      rmse   34.5080932      90   0.014
   mean_baseline Beak.Length_Culmen pearson_r           NA      90   0.014
   mean_baseline      Tarsus.Length      rmse   41.2449180      90   0.014
   mean_baseline      Tarsus.Length pearson_r           NA      90   0.014
   mean_baseline        Wing.Length      rmse  122.5179080      90   0.014
   mean_baseline        Wing.Length pearson_r           NA      90   0.014
   mean_baseline      Trophic.Level  accuracy    0.5666667      90   0.014
   mean_baseline  Primary.Lifestyle  accuracy    0.5777778      90   0.014
   mean_baseline          Migration  accuracy    0.7555556      90   0.014
 pigauto_default               Mass      rmse 2434.9191088      90 362.003
 pigauto_default               Mass pearson_r    0.6917224      90 362.003
 pigauto_default Beak.Length_Culmen      rmse   22.2669802      90 362.003
 pigauto_default Beak.Length_Culmen pearson_r    0.8226840      90 362.003
 pigauto_default      Tarsus.Length      rmse   23.2998135      90 362.003
 pigauto_default      Tarsus.Length pearson_r    0.9072575      90 362.003
 pigauto_default        Wing.Length      rmse   75.5859578      90 362.003
 pigauto_default        Wing.Length pearson_r    0.8581673      90 362.003
 pigauto_default      Trophic.Level  accuracy    0.7888889      90 362.003
 pigauto_default  Primary.Lifestyle  accuracy    0.6888889      90 362.003
 pigauto_default          Migration  accuracy    0.5777778      90 362.003
    bace_default               Mass      rmse 1921.1302664      90  59.360
    bace_default               Mass pearson_r    0.7827951      90  59.360
    bace_default Beak.Length_Culmen      rmse   14.9053081      90  59.360
    bace_default Beak.Length_Culmen pearson_r    0.9035781      90  59.360
    bace_default      Tarsus.Length      rmse   11.8420530      90  59.360
    bace_default      Tarsus.Length pearson_r    0.9578286      90  59.360
    bace_default        Wing.Length      rmse   48.7803437      90  59.360
    bace_default        Wing.Length pearson_r    0.9189314      90  59.360
    bace_default      Trophic.Level  accuracy    0.5000000      90  59.360
    bace_default  Primary.Lifestyle  accuracy    0.3888889      90  59.360
    bace_default          Migration  accuracy    0.5777778      90  59.360
```
