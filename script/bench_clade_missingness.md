# Phase 8.3: clade-correlated missingness (realistic MAR)

n=300, mixed types, n_reps=3, focal trait=c1
Total wall: 16.7 min

## Focal-trait metrics per (method, target_frac)

```
          method target_frac trait    metric     value
 pigauto_default        0.10    c1 pearson_r 0.9727866
     pigauto_em5        0.10    c1 pearson_r 0.9727866
 pigauto_default        0.25    c1 pearson_r 0.9739005
     pigauto_em5        0.25    c1 pearson_r 0.9739005
 pigauto_default        0.40    c1 pearson_r 0.9648520
     pigauto_em5        0.40    c1 pearson_r 0.9648520
   mean_baseline        0.10    c1      rmse 0.9948995
 pigauto_default        0.10    c1      rmse 0.1799307
     pigauto_em5        0.10    c1      rmse 0.1799307
   mean_baseline        0.25    c1      rmse 1.1111784
 pigauto_default        0.25    c1      rmse 0.1428676
     pigauto_em5        0.25    c1      rmse 0.1428676
   mean_baseline        0.40    c1      rmse 1.1578697
 pigauto_default        0.40    c1      rmse 0.1596219
     pigauto_em5        0.40    c1      rmse 0.1596219
```

## All traits summary

```
          method target_frac trait    metric      value
   mean_baseline        0.10    b1  accuracy 0.78888889
 pigauto_default        0.10    b1  accuracy 1.00000000
     pigauto_em5        0.10    b1  accuracy 1.00000000
   mean_baseline        0.25    b1  accuracy 0.77777778
 pigauto_default        0.25    b1  accuracy 1.00000000
     pigauto_em5        0.25    b1  accuracy 1.00000000
   mean_baseline        0.40    b1  accuracy 0.76666667
 pigauto_default        0.40    b1  accuracy 0.90000000
     pigauto_em5        0.40    b1  accuracy 0.90000000
   mean_baseline        0.10  cat1  accuracy 0.60000000
 pigauto_default        0.10  cat1  accuracy 0.93333333
     pigauto_em5        0.10  cat1  accuracy 0.93333333
   mean_baseline        0.25  cat1  accuracy 0.60000000
 pigauto_default        0.25  cat1  accuracy 0.86666667
     pigauto_em5        0.25  cat1  accuracy 0.86666667
   mean_baseline        0.40  cat1  accuracy 0.53333333
 pigauto_default        0.40  cat1  accuracy 0.95000000
     pigauto_em5        0.40  cat1  accuracy 0.95000000
 pigauto_default        0.10    c1 pearson_r 0.97278657
     pigauto_em5        0.10    c1 pearson_r 0.97278657
 pigauto_default        0.25    c1 pearson_r 0.97390047
     pigauto_em5        0.25    c1 pearson_r 0.97390047
 pigauto_default        0.40    c1 pearson_r 0.96485203
     pigauto_em5        0.40    c1 pearson_r 0.96485203
 pigauto_default        0.10    c2 pearson_r 0.94795179
     pigauto_em5        0.10    c2 pearson_r 0.94795179
 pigauto_default        0.25    c2 pearson_r 0.96859484
     pigauto_em5        0.25    c2 pearson_r 0.97159844
 pigauto_default        0.40    c2 pearson_r 0.96403003
     pigauto_em5        0.40    c2 pearson_r 0.96370931
 pigauto_default        0.10    c3 pearson_r 0.92847806
     pigauto_em5        0.10    c3 pearson_r 0.92847806
 pigauto_default        0.25    c3 pearson_r 0.93335990
     pigauto_em5        0.25    c3 pearson_r 0.93335990
 pigauto_default        0.40    c3 pearson_r 0.98212068
     pigauto_em5        0.40    c3 pearson_r 0.98212068
   mean_baseline        0.10    c1      rmse 0.99489953
 pigauto_default        0.10    c1      rmse 0.17993069
     pigauto_em5        0.10    c1      rmse 0.17993069
   mean_baseline        0.25    c1      rmse 1.11117838
 pigauto_default        0.25    c1      rmse 0.14286756
     pigauto_em5        0.25    c1      rmse 0.14286756
   mean_baseline        0.40    c1      rmse 1.15786969
 pigauto_default        0.40    c1      rmse 0.15962195
     pigauto_em5        0.40    c1      rmse 0.15962195
   mean_baseline        0.10    c2      rmse 0.59537920
 pigauto_default        0.10    c2      rmse 0.16296841
     pigauto_em5        0.10    c2      rmse 0.16296841
   mean_baseline        0.25    c2      rmse 0.61078611
 pigauto_default        0.25    c2      rmse 0.16740555
     pigauto_em5        0.25    c2      rmse 0.15247326
   mean_baseline        0.40    c2      rmse 0.58748484
 pigauto_default        0.40    c2      rmse 0.14163093
     pigauto_em5        0.40    c2      rmse 0.14043682
   mean_baseline        0.10    c3      rmse 0.54474468
 pigauto_default        0.10    c3      rmse 0.16328236
     pigauto_em5        0.10    c3      rmse 0.16328236
   mean_baseline        0.25    c3      rmse 0.52737924
 pigauto_default        0.25    c3      rmse 0.12823525
     pigauto_em5        0.25    c3      rmse 0.12823525
   mean_baseline        0.40    c3      rmse 0.57153835
 pigauto_default        0.40    c3      rmse 0.08214369
     pigauto_em5        0.40    c3      rmse 0.08214369
```
