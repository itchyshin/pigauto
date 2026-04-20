# Phase 8.2: evolutionary-model sweep

n=300, K=4 continuous traits, miss_frac=0.30, n_reps=3
Total wall: 22.8 min

## Means per (method, model, metric), averaged over traits and reps

```
          method        model               metric      value
 pigauto_default           BM coverage95_conformal 0.96944444
     pigauto_em5           BM coverage95_conformal 0.97129630
 pigauto_default    nonlinear coverage95_conformal 0.96481481
     pigauto_em5    nonlinear coverage95_conformal 0.96018519
 pigauto_default           OU coverage95_conformal 0.96851852
     pigauto_em5           OU coverage95_conformal 0.96759259
 pigauto_default regime_shift coverage95_conformal 0.95833333
     pigauto_em5 regime_shift coverage95_conformal 0.95092593
 pigauto_default           BM coverage95_mcdropout 0.84166667
     pigauto_em5           BM coverage95_mcdropout 0.86203704
 pigauto_default    nonlinear coverage95_mcdropout 0.57777778
     pigauto_em5    nonlinear coverage95_mcdropout 0.70925926
 pigauto_default           OU coverage95_mcdropout 0.83240741
     pigauto_em5           OU coverage95_mcdropout 0.83796296
 pigauto_default regime_shift coverage95_mcdropout 0.86944444
     pigauto_em5 regime_shift coverage95_mcdropout 0.85740741
 pigauto_default           BM            pearson_r 0.97325786
     pigauto_em5           BM            pearson_r 0.97371030
 pigauto_default    nonlinear            pearson_r 0.68921020
     pigauto_em5    nonlinear            pearson_r 0.67097040
 pigauto_default           OU            pearson_r 0.94155078
     pigauto_em5           OU            pearson_r 0.94210398
 pigauto_default regime_shift            pearson_r 0.99721615
     pigauto_em5 regime_shift            pearson_r 0.99726158
   mean_baseline           BM                 rmse 0.76377155
 pigauto_default           BM                 rmse 0.12944808
     pigauto_em5           BM                 rmse 0.12991972
   mean_baseline    nonlinear                 rmse 0.99720920
 pigauto_default    nonlinear                 rmse 0.61805490
     pigauto_em5    nonlinear                 rmse 0.66321809
   mean_baseline           OU                 rmse 0.98611772
 pigauto_default           OU                 rmse 0.32393724
     pigauto_em5           OU                 rmse 0.32269340
   mean_baseline regime_shift                 rmse 1.01404776
 pigauto_default regime_shift                 rmse 0.07130286
     pigauto_em5 regime_shift                 rmse 0.07044602
```
