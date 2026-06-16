# TabPFN Phylo-Feature Benchmark

- Generated: 2026-06-10 11:14:43 MDT
- Git commit: `caf2388f057c383c480bc616da166539420bbc6e`
- Scales: `2000`
- Variants: `plain, lappe, lappe_nfa`
- Replicates: `1`
- Dry run: `FALSE`

## Design

Continuous AVONET targets are held out with pigauto's `make_missing_splits()` and scored on the same test cells used for the BM baseline and pigauto. TabPFN feature sets are built from same-row traits, optional Laplacian phylogenetic eigenvectors, and optional nearest-training-target aggregate features from cophenetic distances. Split conformal intervals use validation residuals from the same split.

## Status Counts

           method status Freq
         baseline     ok    4
          pigauto     ok    4
     tabpfn_lappe     ok    4
 tabpfn_lappe_nfa     ok    4
     tabpfn_plain     ok    4

## Test Summary

           method scale_n              trait   rmse pearson_r coverage_95
         baseline    2000 Beak.Length_Culmen 0.4048    0.9144          NA
          pigauto    2000 Beak.Length_Culmen 0.4048    0.9144      1.0000
     tabpfn_lappe    2000 Beak.Length_Culmen 0.2844    0.9584      0.9502
 tabpfn_lappe_nfa    2000 Beak.Length_Culmen 0.3340    0.9424      0.9544
     tabpfn_plain    2000 Beak.Length_Culmen 0.5027    0.8686      0.9461
         baseline    2000               Mass 0.3357    0.9424          NA
          pigauto    2000               Mass 0.3428    0.9407      1.0000
     tabpfn_lappe    2000               Mass 0.2289    0.9735      0.9567
 tabpfn_lappe_nfa    2000               Mass 0.2291    0.9734      0.9437
     tabpfn_plain    2000               Mass 0.2797    0.9605      0.9437
         baseline    2000      Tarsus.Length 0.3267    0.9377          NA
          pigauto    2000      Tarsus.Length 0.3267    0.9377      0.9962
     tabpfn_lappe    2000      Tarsus.Length 0.2896    0.9510      0.9624
 tabpfn_lappe_nfa    2000      Tarsus.Length 0.2872    0.9519      0.9586
     tabpfn_plain    2000      Tarsus.Length 0.4236    0.8922      0.9511
         baseline    2000        Wing.Length 0.8750    0.7381          NA
          pigauto    2000        Wing.Length 0.8750    0.7381      0.9961
     tabpfn_lappe    2000        Wing.Length 0.8558    0.7520      0.9225
 tabpfn_lappe_nfa    2000        Wing.Length 0.8417    0.7605      0.9380
     tabpfn_plain    2000        Wing.Length 0.8809    0.7335      0.9612
