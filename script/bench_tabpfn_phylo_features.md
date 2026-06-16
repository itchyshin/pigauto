# TabPFN Phylo-Feature Benchmark

- Generated: 2026-06-10 10:53:36 MDT
- Git commit: `a4c09b4ed2317cfa8d68d9a60035c405a13d176c`
- Scales: `50, 75, 300`
- Variants: `plain, lappe, lappe_nfa, knn`
- Replicates: `3`
- Dry run: `FALSE`

## Design

Continuous AVONET targets are held out with pigauto's `make_missing_splits()` and scored on the same test cells used for the BM baseline and pigauto. TabPFN feature sets are built from same-row traits, optional Laplacian phylogenetic eigenvectors, and optional nearest-training-target aggregate features from cophenetic distances. Split conformal intervals use validation residuals from the same split.

## Status Counts

           method status Freq
         baseline     ok   36
          pigauto     ok   36
       tabpfn_knn     ok   36
     tabpfn_lappe     ok   36
 tabpfn_lappe_nfa     ok   36
     tabpfn_plain     ok   36

## Test Summary

           method scale_n              trait   rmse pearson_r coverage_95
         baseline      50 Beak.Length_Culmen 0.8032    0.6299          NA
          pigauto      50 Beak.Length_Culmen 0.6153    0.7193      0.9667
       tabpfn_knn      50 Beak.Length_Culmen 0.8531    0.7969      0.6741
     tabpfn_lappe      50 Beak.Length_Culmen 0.5541    0.9324      0.9444
 tabpfn_lappe_nfa      50 Beak.Length_Culmen 0.7162    0.8649      0.6407
     tabpfn_plain      50 Beak.Length_Culmen 0.4830    0.9372      0.8741
         baseline      50               Mass 0.8020    0.6771          NA
          pigauto      50               Mass 0.3937    0.9577      1.0000
       tabpfn_knn      50               Mass 0.3249    0.9679      0.8381
     tabpfn_lappe      50               Mass 0.3345    0.9479      0.7905
 tabpfn_lappe_nfa      50               Mass 0.2962    0.9682      0.7905
     tabpfn_plain      50               Mass 0.3516    0.9434      0.7905
         baseline      50      Tarsus.Length 0.4495    0.7973          NA
          pigauto      50      Tarsus.Length 0.4429    0.8627      1.0000
       tabpfn_knn      50      Tarsus.Length 0.4671    0.4158      0.9167
     tabpfn_lappe      50      Tarsus.Length 0.4035    0.6742      0.9167
 tabpfn_lappe_nfa      50      Tarsus.Length 0.5024    0.4182      0.9167
     tabpfn_plain      50      Tarsus.Length 0.3754    0.7150      0.9167
         baseline      50        Wing.Length 0.7404    0.8227          NA
          pigauto      50        Wing.Length 0.4602    0.9460      1.0000
       tabpfn_knn      50        Wing.Length 0.3185    0.9558      0.8611
     tabpfn_lappe      50        Wing.Length 0.3694    0.9443      0.8611
 tabpfn_lappe_nfa      50        Wing.Length 0.3467    0.9460      0.8611
     tabpfn_plain      50        Wing.Length 0.3587    0.9474      0.8472
         baseline      75 Beak.Length_Culmen 0.7903    0.7733          NA
          pigauto      75 Beak.Length_Culmen 0.6441    0.8308      0.9722
       tabpfn_knn      75 Beak.Length_Culmen 0.6881    0.7741      0.7866
     tabpfn_lappe      75 Beak.Length_Culmen 0.6148    0.8168      0.7449
 tabpfn_lappe_nfa      75 Beak.Length_Culmen 0.6180    0.7893      0.8283
     tabpfn_plain      75 Beak.Length_Culmen 0.6633    0.8097      0.6616
         baseline      75               Mass 0.8085    0.7585          NA
          pigauto      75               Mass 0.3185    0.9790      1.0000
       tabpfn_knn      75               Mass 0.2394    0.9820      1.0000
     tabpfn_lappe      75               Mass 0.2555    0.9830      0.8889
 tabpfn_lappe_nfa      75               Mass 0.2543    0.9779      0.9444
     tabpfn_plain      75               Mass 0.2451    0.9827      0.9444
         baseline      75      Tarsus.Length 0.6800    0.6071          NA
          pigauto      75      Tarsus.Length 0.5266    0.7503      1.0000
       tabpfn_knn      75      Tarsus.Length 0.5260    0.7982      0.8753
     tabpfn_lappe      75      Tarsus.Length 0.4026    0.8397      0.9394
 tabpfn_lappe_nfa      75      Tarsus.Length 0.5244    0.8142      0.8449
     tabpfn_plain      75      Tarsus.Length 0.4164    0.8494      0.9394
         baseline      75        Wing.Length 0.5087    0.7462          NA
          pigauto      75        Wing.Length 0.3962    0.8108      1.0000
       tabpfn_knn      75        Wing.Length 0.3653    0.8996      0.9259
     tabpfn_lappe      75        Wing.Length 0.3221    0.9241      0.9259
 tabpfn_lappe_nfa      75        Wing.Length 0.3873    0.8765      0.8889
     tabpfn_plain      75        Wing.Length 0.3063    0.9258      0.8519
         baseline     300 Beak.Length_Culmen 0.5732    0.8129          NA
          pigauto     300 Beak.Length_Culmen 0.5725    0.8151      1.0000
       tabpfn_knn     300 Beak.Length_Culmen 0.6836    0.7286      0.9218
     tabpfn_lappe     300 Beak.Length_Culmen 0.4777    0.8734      0.9834
 tabpfn_lappe_nfa     300 Beak.Length_Culmen 0.6556    0.7570      0.9294
     tabpfn_plain     300 Beak.Length_Culmen 0.5224    0.8496      0.9730
         baseline     300               Mass 0.4763    0.8560          NA
          pigauto     300               Mass 0.3719    0.9226      1.0000
       tabpfn_knn     300               Mass 0.2004    0.9758      0.9649
     tabpfn_lappe     300               Mass 0.2279    0.9702      0.9649
 tabpfn_lappe_nfa     300               Mass 0.2057    0.9762      0.9649
     tabpfn_plain     300               Mass 0.2190    0.9698      0.9737
         baseline     300      Tarsus.Length 0.5458    0.8629          NA
          pigauto     300      Tarsus.Length 0.5218    0.8782      0.9921
       tabpfn_knn     300      Tarsus.Length 0.4411    0.9198      0.9734
     tabpfn_lappe     300      Tarsus.Length 0.3882    0.9313      0.9639
 tabpfn_lappe_nfa     300      Tarsus.Length 0.4172    0.9253      0.9829
     tabpfn_plain     300      Tarsus.Length 0.4338    0.9178      0.9553
         baseline     300        Wing.Length 0.5687    0.8507          NA
          pigauto     300        Wing.Length 0.4908    0.8903      0.9926
       tabpfn_knn     300        Wing.Length 0.2467    0.9720      0.9640
     tabpfn_lappe     300        Wing.Length 0.2990    0.9598      0.9396
 tabpfn_lappe_nfa     300        Wing.Length 0.2311    0.9756      0.9730
     tabpfn_plain     300        Wing.Length 0.3018    0.9576      0.9746
