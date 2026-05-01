# Multi-obs covariate-lift on real bird phylogeny (tree300)

Tree: 300 tips (AVONET 300 real bird phylogeny)
Reps: 2, epochs: 200

Story: bench_multi_obs.R uses ape::rtree() Yule sim tree at n=200.
This script is the same DGP on the REAL bird phylogeny tree300.
Confirms the architectural lift survives realistic phylogenetic structure.

## Per-cell results (mean across reps, observation-level RMSE)

```
 lambda beta sp_missing_frac         method  obs_rmse   sp_rmse obs_pearson_r
    0.5  0.0             0.5    pigauto_cov  3.455902  3.740638     0.7742663
    0.5  0.0             0.5 pigauto_no_cov  3.459234  3.744138     0.7736494
    0.5  0.0             0.5   species_mean  5.039046  5.562925     0.3819110
    0.5  0.0             0.8    pigauto_cov  4.142678  4.202468     0.5118587
    0.5  0.0             0.8 pigauto_no_cov  4.140201  4.199232     0.5118387
    0.5  0.0             0.8   species_mean  4.713208  4.829396     0.2090858
    0.5  0.5             0.5    pigauto_cov  4.554960  4.082586     0.7125781
    0.5  0.5             0.5 pigauto_no_cov  4.741945  4.085739     0.6848099
    0.5  0.5             0.5   species_mean  6.150708  5.973083     0.3561700
    0.5  0.5             0.8    pigauto_cov  5.471255  5.092852     0.5990901
    0.5  0.5             0.8 pigauto_no_cov  5.628496  5.118410     0.5764017
    0.5  0.5             0.8   species_mean  6.805463  6.464736     0.2385914
    0.5  1.0             0.5    pigauto_cov  6.037733  4.507409     0.6950427
    0.5  1.0             0.5 pigauto_no_cov  6.647976  4.664206     0.6259007
    0.5  1.0             0.5   species_mean  8.016506  6.353169     0.3120333
    0.5  1.0             0.8    pigauto_cov  6.411539  4.809998     0.5805100
    0.5  1.0             0.8 pigauto_no_cov  6.923501  4.803569     0.4720641
    0.5  1.0             0.8   species_mean  7.718648  5.901182     0.1543793
    0.9  0.0             0.5    pigauto_cov  4.331124  4.869675     0.8248721
    0.9  0.0             0.5 pigauto_no_cov  4.333593  4.869082     0.8246016
    0.9  0.0             0.5   species_mean  7.006796  7.849298     0.4060096
    0.9  0.0             0.8    pigauto_cov  6.324988  6.407826     0.5805051
    0.9  0.0             0.8 pigauto_no_cov  6.325640  6.405235     0.5793276
    0.9  0.0             0.8   species_mean  7.675158  7.933005     0.2035968
    0.9  0.5             0.5    pigauto_cov  5.273962  5.148424     0.7760590
    0.9  0.5             0.5 pigauto_no_cov  5.405248  5.137673     0.7619549
    0.9  0.5             0.5   species_mean  7.833238  8.020699     0.3737272
    0.9  0.5             0.8    pigauto_cov  6.499223  6.428835     0.8054702
    0.9  0.5             0.8 pigauto_no_cov  6.600528  6.450133     0.7985474
    0.9  0.5             0.8   species_mean 10.730721 10.848151     0.2173854
    0.9  1.0             0.5    pigauto_cov  6.687078  5.599496     0.6225137
    0.9  1.0             0.5 pigauto_no_cov  7.221730  5.634646     0.5477142
    0.9  1.0             0.5   species_mean  8.157640  6.903883     0.2839769
    0.9  1.0             0.8    pigauto_cov  7.397728  5.889718     0.6446607
    0.9  1.0             0.8 pigauto_no_cov  7.832687  5.859125     0.5847585
    0.9  1.0             0.8   species_mean  9.506671  7.818611     0.1692544
```

## RMSE-lift summary (pigauto_cov / pigauto_no_cov)

```
 lambda beta sp_missing_frac obs_rmse.pigauto_no_cov obs_rmse.pigauto_cov
    0.5  0.0             0.5                3.459234             3.455902
    0.9  0.0             0.5                4.333593             4.331124
    0.5  0.5             0.5                4.741945             4.554960
    0.9  0.5             0.5                5.405248             5.273962
    0.5  1.0             0.5                6.647976             6.037733
    0.9  1.0             0.5                7.221730             6.687078
    0.5  0.0             0.8                4.140201             4.142678
    0.9  0.0             0.8                6.325640             6.324988
    0.5  0.5             0.8                5.628496             5.471255
    0.9  0.5             0.8                6.600528             6.499223
    0.5  1.0             0.8                6.923501             6.411539
    0.9  1.0             0.8                7.832687             7.397728
  cov_lift_obs
  9.458007e-04
  6.295189e-04
  3.943323e-02
  2.434307e-02
  9.167402e-02
  7.472788e-02
 -5.304107e-04
  4.430618e-05
  2.796055e-02
  1.539609e-02
  7.349296e-02
  5.576150e-02
```
