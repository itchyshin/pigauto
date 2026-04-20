# Phase 8.2: evolutionary-model sweep

n=300, K=4 continuous traits, miss_frac=0.30, n_reps=3
Total wall: 23.7 min

## Means per (method, model, metric), averaged over traits and reps

```
          method        model    metric      value
 pigauto_default           BM pearson_r 0.97446045
     pigauto_em5           BM pearson_r 0.97460255
 pigauto_default    nonlinear pearson_r 0.69394955
     pigauto_em5    nonlinear pearson_r 0.69283117
 pigauto_default           OU pearson_r 0.94362183
     pigauto_em5           OU pearson_r 0.94364050
 pigauto_default regime_shift pearson_r 0.99732355
     pigauto_em5 regime_shift pearson_r 0.99731414
   mean_baseline           BM      rmse 0.76377155
 pigauto_default           BM      rmse 0.12657521
     pigauto_em5           BM      rmse 0.12646105
   mean_baseline    nonlinear      rmse 0.99720920
 pigauto_default    nonlinear      rmse 0.62560534
     pigauto_em5    nonlinear      rmse 0.62710943
   mean_baseline           OU      rmse 0.98611772
 pigauto_default           OU      rmse 0.31788618
     pigauto_em5           OU      rmse 0.31814343
   mean_baseline regime_shift      rmse 1.01404776
 pigauto_default regime_shift      rmse 0.06957536
     pigauto_em5 regime_shift      rmse 0.06936114
```
