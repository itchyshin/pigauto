# Phase 8.3: clade-correlated missingness (realistic MAR)

n=300, mixed types, n_reps=3, focal trait=c1
Total wall: 17.9 min

## Focal-trait metrics per (method, target_frac)

```
          method target_frac trait               metric     value
 pigauto_default        0.10    c1 coverage95_conformal 0.9666667
     pigauto_em5        0.10    c1 coverage95_conformal 0.9666667
 pigauto_default        0.25    c1 coverage95_conformal 0.9066667
     pigauto_em5        0.25    c1 coverage95_conformal 0.9000000
 pigauto_default        0.40    c1 coverage95_conformal 0.9166667
     pigauto_em5        0.40    c1 coverage95_conformal 0.9166667
 pigauto_default        0.10    c1 coverage95_mcdropout 0.8000000
     pigauto_em5        0.10    c1 coverage95_mcdropout 0.7833333
 pigauto_default        0.25    c1 coverage95_mcdropout 0.9266667
     pigauto_em5        0.25    c1 coverage95_mcdropout 0.9000000
 pigauto_default        0.40    c1 coverage95_mcdropout 0.9000000
     pigauto_em5        0.40    c1 coverage95_mcdropout 0.8791667
 pigauto_default        0.10    c1            pearson_r 0.9692203
     pigauto_em5        0.10    c1            pearson_r 0.9747650
 pigauto_default        0.25    c1            pearson_r 0.9721579
     pigauto_em5        0.25    c1            pearson_r 0.9711420
 pigauto_default        0.40    c1            pearson_r 0.9626111
     pigauto_em5        0.40    c1            pearson_r 0.9600209
   mean_baseline        0.10    c1                 rmse 0.9948995
 pigauto_default        0.10    c1                 rmse 0.1810246
     pigauto_em5        0.10    c1                 rmse 0.1778391
   mean_baseline        0.25    c1                 rmse 1.1111784
 pigauto_default        0.25    c1                 rmse 0.1471566
     pigauto_em5        0.25    c1                 rmse 0.1503378
   mean_baseline        0.40    c1                 rmse 1.1578697
 pigauto_default        0.40    c1                 rmse 0.1660802
     pigauto_em5        0.40    c1                 rmse 0.1696286
```

## All traits summary

```
          method target_frac trait               metric      value
   mean_baseline        0.10    b1             accuracy 0.78888889
 pigauto_default        0.10    b1             accuracy 1.00000000
     pigauto_em5        0.10    b1             accuracy 1.00000000
   mean_baseline        0.25    b1             accuracy 0.77777778
 pigauto_default        0.25    b1             accuracy 1.00000000
     pigauto_em5        0.25    b1             accuracy 1.00000000
   mean_baseline        0.40    b1             accuracy 0.76666667
 pigauto_default        0.40    b1             accuracy 0.90000000
     pigauto_em5        0.40    b1             accuracy 0.90000000
   mean_baseline        0.10  cat1             accuracy 0.60000000
 pigauto_default        0.10  cat1             accuracy 0.93333333
     pigauto_em5        0.10  cat1             accuracy 0.93333333
   mean_baseline        0.25  cat1             accuracy 0.60000000
 pigauto_default        0.25  cat1             accuracy 0.86666667
     pigauto_em5        0.25  cat1             accuracy 0.86666667
   mean_baseline        0.40  cat1             accuracy 0.53333333
 pigauto_default        0.40  cat1             accuracy 0.95000000
     pigauto_em5        0.40  cat1             accuracy 0.95000000
 pigauto_default        0.10    c1 coverage95_conformal 0.96666667
     pigauto_em5        0.10    c1 coverage95_conformal 0.96666667
 pigauto_default        0.25    c1 coverage95_conformal 0.90666667
     pigauto_em5        0.25    c1 coverage95_conformal 0.90000000
 pigauto_default        0.40    c1 coverage95_conformal 0.91666667
     pigauto_em5        0.40    c1 coverage95_conformal 0.91666667
 pigauto_default        0.10    c2 coverage95_conformal 0.96666667
     pigauto_em5        0.10    c2 coverage95_conformal 0.96666667
 pigauto_default        0.25    c2 coverage95_conformal 0.96666667
     pigauto_em5        0.25    c2 coverage95_conformal 1.00000000
 pigauto_default        0.40    c2 coverage95_conformal 0.90000000
     pigauto_em5        0.40    c2 coverage95_conformal 0.86666667
 pigauto_default        0.10    c3 coverage95_conformal 0.98333333
     pigauto_em5        0.10    c3 coverage95_conformal 0.98333333
 pigauto_default        0.25    c3 coverage95_conformal 0.96666667
     pigauto_em5        0.25    c3 coverage95_conformal 0.98333333
 pigauto_default        0.40    c3 coverage95_conformal 1.00000000
     pigauto_em5        0.40    c3 coverage95_conformal 1.00000000
 pigauto_default        0.10    c1 coverage95_mcdropout 0.80000000
     pigauto_em5        0.10    c1 coverage95_mcdropout 0.78333333
 pigauto_default        0.25    c1 coverage95_mcdropout 0.92666667
     pigauto_em5        0.25    c1 coverage95_mcdropout 0.90000000
 pigauto_default        0.40    c1 coverage95_mcdropout 0.90000000
     pigauto_em5        0.40    c1 coverage95_mcdropout 0.87916667
 pigauto_default        0.10    c2 coverage95_mcdropout 0.81666667
     pigauto_em5        0.10    c2 coverage95_mcdropout 0.83333333
 pigauto_default        0.25    c2 coverage95_mcdropout 0.90000000
     pigauto_em5        0.25    c2 coverage95_mcdropout 0.86666667
 pigauto_default        0.40    c2 coverage95_mcdropout 0.80000000
     pigauto_em5        0.40    c2 coverage95_mcdropout 0.86666667
 pigauto_default        0.10    c3 coverage95_mcdropout 0.75000000
     pigauto_em5        0.10    c3 coverage95_mcdropout 0.81666667
 pigauto_default        0.25    c3 coverage95_mcdropout 0.88333333
     pigauto_em5        0.25    c3 coverage95_mcdropout 0.90000000
 pigauto_default        0.40    c3 coverage95_mcdropout 0.90000000
     pigauto_em5        0.40    c3 coverage95_mcdropout 0.90000000
 pigauto_default        0.10    c1            pearson_r 0.96922027
     pigauto_em5        0.10    c1            pearson_r 0.97476504
 pigauto_default        0.25    c1            pearson_r 0.97215793
     pigauto_em5        0.25    c1            pearson_r 0.97114196
 pigauto_default        0.40    c1            pearson_r 0.96261109
     pigauto_em5        0.40    c1            pearson_r 0.96002087
 pigauto_default        0.10    c2            pearson_r 0.93671435
     pigauto_em5        0.10    c2            pearson_r 0.93103007
 pigauto_default        0.25    c2            pearson_r 0.96799825
     pigauto_em5        0.25    c2            pearson_r 0.97724327
 pigauto_default        0.40    c2            pearson_r 0.96027986
     pigauto_em5        0.40    c2            pearson_r 0.96202012
 pigauto_default        0.10    c3            pearson_r 0.94146470
     pigauto_em5        0.10    c3            pearson_r 0.92321592
 pigauto_default        0.25    c3            pearson_r 0.91857139
     pigauto_em5        0.25    c3            pearson_r 0.94119207
 pigauto_default        0.40    c3            pearson_r 0.98186256
     pigauto_em5        0.40    c3            pearson_r 0.98115640
   mean_baseline        0.10    c1                 rmse 0.99489953
 pigauto_default        0.10    c1                 rmse 0.18102464
     pigauto_em5        0.10    c1                 rmse 0.17783907
   mean_baseline        0.25    c1                 rmse 1.11117838
 pigauto_default        0.25    c1                 rmse 0.14715661
     pigauto_em5        0.25    c1                 rmse 0.15033776
   mean_baseline        0.40    c1                 rmse 1.15786969
 pigauto_default        0.40    c1                 rmse 0.16608021
     pigauto_em5        0.40    c1                 rmse 0.16962863
   mean_baseline        0.10    c2                 rmse 0.59537920
 pigauto_default        0.10    c2                 rmse 0.17212871
     pigauto_em5        0.10    c2                 rmse 0.18144597
   mean_baseline        0.25    c2                 rmse 0.61078611
 pigauto_default        0.25    c2                 rmse 0.16055338
     pigauto_em5        0.25    c2                 rmse 0.13975372
   mean_baseline        0.40    c2                 rmse 0.58748484
 pigauto_default        0.40    c2                 rmse 0.14984367
     pigauto_em5        0.40    c2                 rmse 0.14725711
   mean_baseline        0.10    c3                 rmse 0.54474468
 pigauto_default        0.10    c3                 rmse 0.14840690
     pigauto_em5        0.10    c3                 rmse 0.16480326
   mean_baseline        0.25    c3                 rmse 0.52737924
 pigauto_default        0.25    c3                 rmse 0.14058209
     pigauto_em5        0.25    c3                 rmse 0.12001811
   mean_baseline        0.40    c3                 rmse 0.57153835
 pigauto_default        0.40    c3                 rmse 0.08018079
     pigauto_em5        0.40    c3                 rmse 0.08508044
```
