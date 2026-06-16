# TabPFN AVONET Discrete Benchmark

- Generated: 2026-06-10 18:17:37 MDT
- Git commit: `4414ac16297fe8687c4cfacce15027c462afcaca`
- Scales: `50, 75, 300`
- Variants: `plain, lappe, lappe_nfa, knn`
- Replicates: `3`
- Dry run: `FALSE`

## Design

AVONET categorical/ordinal targets are held out with pigauto's `make_missing_splits()` and scored as classification tasks on the same test rows used for baseline and pigauto. TabPFN feature sets use same-row latent traits, optional Laplacian phylogenetic eigenvectors, and optional nearest-training-target class-proportion features. This benchmark reports accuracy and balanced accuracy only; it does not attempt classification prediction sets.

## Status Counts

           method status Freq
         baseline     ok   27
         majority     ok   27
          pigauto     ok   27
       tabpfn_knn     ok   27
     tabpfn_lappe     ok   27
 tabpfn_lappe_nfa     ok   27
     tabpfn_plain     ok   27

## Test Summary

           method scale_n             trait        type accuracy
         majority      50         Migration     ordinal   0.9333
     tabpfn_lappe      50         Migration     ordinal   0.9333
     tabpfn_plain      50         Migration     ordinal   0.9333
         baseline      50         Migration     ordinal   0.8381
          pigauto      50         Migration     ordinal   0.8381
       tabpfn_knn      50         Migration     ordinal   0.7111
 tabpfn_lappe_nfa      50         Migration     ordinal   0.7111
         baseline      50 Primary.Lifestyle categorical   0.6595
          pigauto      50 Primary.Lifestyle categorical   0.6595
         majority      50 Primary.Lifestyle categorical   0.6119
     tabpfn_lappe      50 Primary.Lifestyle categorical   0.6119
       tabpfn_knn      50 Primary.Lifestyle categorical   0.5738
     tabpfn_plain      50 Primary.Lifestyle categorical   0.5643
 tabpfn_lappe_nfa      50 Primary.Lifestyle categorical   0.5262
         baseline      50     Trophic.Level categorical   0.7027
     tabpfn_lappe      50     Trophic.Level categorical   0.5866
     tabpfn_plain      50     Trophic.Level categorical   0.5866
          pigauto      50     Trophic.Level categorical   0.5519
       tabpfn_knn      50     Trophic.Level categorical   0.4993
         majority      50     Trophic.Level categorical   0.4913
 tabpfn_lappe_nfa      50     Trophic.Level categorical   0.4343
         baseline      75         Migration     ordinal   0.8364
         majority      75         Migration     ordinal   0.7758
          pigauto      75         Migration     ordinal   0.7758
     tabpfn_lappe      75         Migration     ordinal   0.7758
     tabpfn_plain      75         Migration     ordinal   0.7758
       tabpfn_knn      75         Migration     ordinal   0.7152
 tabpfn_lappe_nfa      75         Migration     ordinal   0.7152
         baseline      75 Primary.Lifestyle categorical   0.6667
     tabpfn_lappe      75 Primary.Lifestyle categorical   0.5889
     tabpfn_plain      75 Primary.Lifestyle categorical   0.5667
          pigauto      75 Primary.Lifestyle categorical   0.5111
       tabpfn_knn      75 Primary.Lifestyle categorical   0.4667
 tabpfn_lappe_nfa      75 Primary.Lifestyle categorical   0.4667
         majority      75 Primary.Lifestyle categorical   0.4444
         baseline      75     Trophic.Level categorical   0.6699
          pigauto      75     Trophic.Level categorical   0.6093
     tabpfn_plain      75     Trophic.Level categorical   0.5945
         majority      75     Trophic.Level categorical   0.5339
     tabpfn_lappe      75     Trophic.Level categorical   0.5339
 tabpfn_lappe_nfa      75     Trophic.Level categorical   0.5314
       tabpfn_knn      75     Trophic.Level categorical   0.5115
         majority     300         Migration     ordinal   0.7246
     tabpfn_lappe     300         Migration     ordinal   0.7246
         baseline     300         Migration     ordinal   0.7217
          pigauto     300         Migration     ordinal   0.7217
     tabpfn_plain     300         Migration     ordinal   0.7148
       tabpfn_knn     300         Migration     ordinal   0.6520
 tabpfn_lappe_nfa     300         Migration     ordinal   0.6123
          pigauto     300 Primary.Lifestyle categorical   0.7653
         baseline     300 Primary.Lifestyle categorical   0.7478
     tabpfn_lappe     300 Primary.Lifestyle categorical   0.7457
     tabpfn_plain     300 Primary.Lifestyle categorical   0.7238
         majority     300 Primary.Lifestyle categorical   0.6024
       tabpfn_knn     300 Primary.Lifestyle categorical   0.5953
 tabpfn_lappe_nfa     300 Primary.Lifestyle categorical   0.5679
         baseline     300     Trophic.Level categorical   0.8007
          pigauto     300     Trophic.Level categorical   0.8007
     tabpfn_lappe     300     Trophic.Level categorical   0.6481
     tabpfn_plain     300     Trophic.Level categorical   0.5504
         majority     300     Trophic.Level categorical   0.4409
 tabpfn_lappe_nfa     300     Trophic.Level categorical   0.4226
       tabpfn_knn     300     Trophic.Level categorical   0.4057
 balanced_accuracy
            0.8333
            0.8333
            0.8333
            0.7381
            0.7381
            0.6111
            0.6111
            0.5417
            0.5417
            0.3889
            0.3889
            0.4861
            0.3333
            0.4583
            0.5926
            0.4907
            0.4583
            0.4630
            0.3750
            0.3889
            0.3148
            0.6111
            0.4444
            0.4444
            0.4444
            0.4444
            0.4583
            0.4583
            0.5833
            0.3750
            0.4083
            0.3750
            0.2333
            0.2333
            0.2500
            0.5602
            0.4861
            0.4630
            0.3889
            0.3889
            0.4236
            0.4375
            0.3333
            0.3449
            0.4764
            0.4764
            0.3291
            0.3467
            0.3178
            0.5526
            0.5530
            0.4918
            0.4827
            0.2167
            0.2955
            0.2926
            0.6677
            0.6677
            0.5061
            0.3936
            0.3056
            0.3369
            0.2962
