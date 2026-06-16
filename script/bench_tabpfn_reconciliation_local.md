# TabPFN Reconciliation Benchmark

- Generated: 2026-06-10 22:38:19 MDT
- Git commit: `6cec7178407200a08199d1efeeaa894af22010f1`
- Scales: `50, 75, 300`
- Regimes: `same_row, phylo_only, same_row_lappe, shuffled_same_row_lappe, same_row_lappe_nfa`
- Replicates: `3`
- Dry run: `FALSE`

## Aim

This benchmark tests whether Russell-style TabPFN gains are mostly explained by same-row cross-trait features, phylogenetic features, or their combination. The shuffled regime keeps the same feature distribution but breaks row-level cross-trait alignment.

## Status Counts

                         method status Freq
                       baseline     ok   36
                        pigauto     ok   36
              tabpfn_phylo_only     ok   36
                tabpfn_same_row     ok   36
          tabpfn_same_row_lappe     ok   36
      tabpfn_same_row_lappe_nfa     ok   36
 tabpfn_shuffled_same_row_lappe     ok   36

## Mean RMSE By Method

                         method scale_n   rmse pearson_r coverage_95 wall_sec
          tabpfn_same_row_lappe      50 0.3920    0.8252      0.9298   7.1436
                tabpfn_same_row      50 0.4007    0.8162      0.9087   7.0026
      tabpfn_same_row_lappe_nfa      50 0.4975    0.6910      0.8390   7.1498
                        pigauto      50 0.5364    0.6601      0.9917 147.1530
                       baseline      50 0.6165    0.6212          NA   0.3597
              tabpfn_phylo_only      50 0.7275    0.4754      0.8928   7.3033
 tabpfn_shuffled_same_row_lappe      50 0.7653    0.3221      0.8854   7.1212
          tabpfn_same_row_lappe      75 0.4347    0.8820      0.8776   7.7444
                tabpfn_same_row      75 0.4401    0.8856      0.8521   7.1591
      tabpfn_same_row_lappe_nfa      75 0.5101    0.8290      0.8717   7.1378
                        pigauto      75 0.5819    0.8053      0.9917 153.0263
                       baseline      75 0.7981    0.6524          NA   0.0670
              tabpfn_phylo_only      75 1.0122    0.2534      0.8362   7.8743
 tabpfn_shuffled_same_row_lappe      75 1.0507    0.1934      0.8258   7.3662
          tabpfn_same_row_lappe     300 0.3345    0.9416      0.9742   9.3677
                tabpfn_same_row     300 0.3611    0.9312      0.9760   8.9992
      tabpfn_same_row_lappe_nfa     300 0.4065    0.9011      0.9541   9.3391
                        pigauto     300 0.4713    0.8906      0.9919 157.8227
                       baseline     300 0.5452    0.8581          NA   0.1630
 tabpfn_shuffled_same_row_lappe     300 0.8052    0.6475      0.9302   9.3160
              tabpfn_phylo_only     300 0.8082    0.6430      0.9240   8.7858

## Mean RMSE By Trait

                         method scale_n              trait   rmse pearson_r
                tabpfn_same_row      50 Beak.Length_Culmen 0.5361  0.835884
          tabpfn_same_row_lappe      50 Beak.Length_Culmen 0.5464  0.852152
                       baseline      50 Beak.Length_Culmen 0.6775  0.657620
                        pigauto      50 Beak.Length_Culmen 0.7305  0.629571
      tabpfn_same_row_lappe_nfa      50 Beak.Length_Culmen 0.7776  0.581646
 tabpfn_shuffled_same_row_lappe      50 Beak.Length_Culmen 0.9400  0.034914
              tabpfn_phylo_only      50 Beak.Length_Culmen 0.9654  0.162091
      tabpfn_same_row_lappe_nfa      50               Mass 0.2366  0.975535
          tabpfn_same_row_lappe      50               Mass 0.2711  0.949791
                tabpfn_same_row      50               Mass 0.3222  0.938982
                        pigauto      50               Mass 0.4949  0.875377
                       baseline      50               Mass 0.6166  0.843199
 tabpfn_shuffled_same_row_lappe      50               Mass 0.6207  0.495290
              tabpfn_phylo_only      50               Mass 0.6464  0.401712
                tabpfn_same_row      50      Tarsus.Length 0.4829  0.601294
          tabpfn_same_row_lappe      50      Tarsus.Length 0.4936  0.583490
                        pigauto      50      Tarsus.Length 0.6083  0.583445
                       baseline      50      Tarsus.Length 0.6322  0.566598
              tabpfn_phylo_only      50      Tarsus.Length 0.6374  0.797860
      tabpfn_same_row_lappe_nfa      50      Tarsus.Length 0.6376  0.295515
 tabpfn_shuffled_same_row_lappe      50      Tarsus.Length 0.7178  0.644383
          tabpfn_same_row_lappe      50        Wing.Length 0.2570  0.915440
                tabpfn_same_row      50        Wing.Length 0.2616  0.888490
                        pigauto      50        Wing.Length 0.3118  0.551881
      tabpfn_same_row_lappe_nfa      50        Wing.Length 0.3381  0.911379
                       baseline      50        Wing.Length 0.5395  0.417429
              tabpfn_phylo_only      50        Wing.Length 0.6608  0.540131
 tabpfn_shuffled_same_row_lappe      50        Wing.Length 0.7828  0.113634
          tabpfn_same_row_lappe      75 Beak.Length_Culmen 0.6461  0.808033
                        pigauto      75 Beak.Length_Culmen 0.6642  0.867086
                tabpfn_same_row      75 Beak.Length_Culmen 0.6731  0.811031
      tabpfn_same_row_lappe_nfa      75 Beak.Length_Culmen 0.7420  0.758588
                       baseline      75 Beak.Length_Culmen 0.8291  0.820829
              tabpfn_phylo_only      75 Beak.Length_Culmen 1.0781  0.355413
 tabpfn_shuffled_same_row_lappe      75 Beak.Length_Culmen 1.1588  0.232544
                tabpfn_same_row      75               Mass 0.2727  0.980373
          tabpfn_same_row_lappe      75               Mass 0.2966  0.978946
      tabpfn_same_row_lappe_nfa      75               Mass 0.2986  0.971387
                        pigauto      75               Mass 0.3950  0.971059
                       baseline      75               Mass 0.9590  0.773127
              tabpfn_phylo_only      75               Mass 1.2007  0.427559
 tabpfn_shuffled_same_row_lappe      75               Mass 1.2186  0.374700
                tabpfn_same_row      75      Tarsus.Length 0.4142  0.843828
          tabpfn_same_row_lappe      75      Tarsus.Length 0.4252  0.823869
      tabpfn_same_row_lappe_nfa      75      Tarsus.Length 0.5183  0.718045
                        pigauto      75      Tarsus.Length 0.5734  0.623161
                       baseline      75      Tarsus.Length 0.6831  0.420166
 tabpfn_shuffled_same_row_lappe      75      Tarsus.Length 0.8535  0.007664
              tabpfn_phylo_only      75      Tarsus.Length 0.8593 -0.006252
          tabpfn_same_row_lappe      75        Wing.Length 0.3710  0.917009
                tabpfn_same_row      75        Wing.Length 0.4001  0.907266
      tabpfn_same_row_lappe_nfa      75        Wing.Length 0.4814  0.868075
                        pigauto      75        Wing.Length 0.6951  0.737078
                       baseline      75        Wing.Length 0.7211  0.595521
              tabpfn_phylo_only      75        Wing.Length 0.9105  0.236793
 tabpfn_shuffled_same_row_lappe      75        Wing.Length 0.9719  0.158612
          tabpfn_same_row_lappe     300 Beak.Length_Culmen 0.4772  0.883930
                tabpfn_same_row     300 Beak.Length_Culmen 0.5529  0.854334
                        pigauto     300 Beak.Length_Culmen 0.6583  0.788484
                       baseline     300 Beak.Length_Culmen 0.6711  0.780288
      tabpfn_same_row_lappe_nfa     300 Beak.Length_Culmen 0.7715  0.719726
 tabpfn_shuffled_same_row_lappe     300 Beak.Length_Culmen 0.8927  0.550860
              tabpfn_phylo_only     300 Beak.Length_Culmen 0.8997  0.539248
                tabpfn_same_row     300               Mass 0.2086  0.975232
      tabpfn_same_row_lappe_nfa     300               Mass 0.2127  0.977891
          tabpfn_same_row_lappe     300               Mass 0.2203  0.974073
                        pigauto     300               Mass 0.3046  0.953973
                       baseline     300               Mass 0.4227  0.918350
              tabpfn_phylo_only     300               Mass 0.7539  0.670499
 tabpfn_shuffled_same_row_lappe     300               Mass 0.7585  0.663949
          tabpfn_same_row_lappe     300      Tarsus.Length 0.3424  0.945532
                tabpfn_same_row     300      Tarsus.Length 0.3936  0.930738
      tabpfn_same_row_lappe_nfa     300      Tarsus.Length 0.4030  0.930725
                        pigauto     300      Tarsus.Length 0.4899  0.893151
                       baseline     300      Tarsus.Length 0.5119  0.877737
 tabpfn_shuffled_same_row_lappe     300      Tarsus.Length 0.7882  0.657550
              tabpfn_phylo_only     300      Tarsus.Length 0.7959  0.648295
      tabpfn_same_row_lappe_nfa     300        Wing.Length 0.2388  0.976095
                tabpfn_same_row     300        Wing.Length 0.2892  0.964685
          tabpfn_same_row_lappe     300        Wing.Length 0.2980  0.962916
                        pigauto     300        Wing.Length 0.4324  0.926946
                       baseline     300        Wing.Length 0.5751  0.856061
 tabpfn_shuffled_same_row_lappe     300        Wing.Length 0.7813  0.717561
              tabpfn_phylo_only     300        Wing.Length 0.7833  0.714156
 coverage_95
      0.8741
      1.0000
          NA
      0.9667
      0.6963
      0.8630
      0.8926
      0.8857
      0.8857
      0.8857
      1.0000
          NA
      0.8571
      0.8095
      0.9167
      0.9167
      1.0000
          NA
      0.8690
      0.8571
      0.8214
      0.9167
      0.9583
      1.0000
      0.9167
          NA
      1.0000
      1.0000
      0.8561
      1.0000
      0.7727
      0.8084
          NA
      0.7911
      0.7495
      0.9444
      0.8889
      0.9444
      1.0000
          NA
      0.8889
      0.8333
      0.9394
      0.9394
      0.8783
      1.0000
          NA
      0.9722
      0.9167
      0.8259
      0.7519
      0.8556
      0.9667
          NA
      0.7481
      0.7481
      0.9829
      0.9567
      0.9756
          NA
      0.8800
      0.8620
      0.8457
      0.9737
      0.9561
      0.9649
      1.0000
          NA
      0.9561
      0.9474
      0.9810
      0.9810
      0.9881
      0.9921
          NA
      0.9587
      0.9492
      0.9922
      0.9926
      0.9679
      1.0000
          NA
      0.9528
      0.9450

## Pigauto Gate Audit

 scale_n              trait   rmse r_cal_bm r_cal_gnn r_cal_mean
      50 Beak.Length_Culmen 0.7305   0.5833    0.2103    0.20635
      50               Mass 0.4949   0.4500    0.5500    0.00000
      50      Tarsus.Length 0.6083   0.4833    0.5167    0.00000
      50        Wing.Length 0.3118   0.2833    0.7167    0.00000
      75 Beak.Length_Culmen 0.6642   0.3500    0.5759    0.07407
      75               Mass 0.3950   0.0000    0.9833    0.01667
      75      Tarsus.Length 0.5734   0.4670    0.5330    0.00000
      75        Wing.Length 0.6951   0.3833    0.2833    0.33333
     300 Beak.Length_Culmen 0.6583   0.7500    0.2500    0.00000
     300               Mass 0.3046   0.4000    0.6000    0.00000
     300      Tarsus.Length 0.4899   0.8140    0.1860    0.00000
     300        Wing.Length 0.4324   0.2500    0.6807    0.06930

