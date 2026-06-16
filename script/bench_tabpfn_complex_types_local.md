# TabPFN Complex-Type Benchmark

- Generated: 2026-06-10 20:21:30 MDT
- Git commit: `b6dd2d752960db2abd63e3e9e0aad6978c73d034`
- Scales: `150, 300`
- Scenarios: `zi_moderate, zi_sparse, zi_many_zeros, multi_moderate, multi_high_phylo, multi_K8`
- Variants: `lappe, lappe_nfa`
- Replicates: `1`

## Design

ZI-count uses one TabPFN classifier for the non-zero gate and one TabPFN regressor for the conditional magnitude. Multi-proportion uses independent TabPFN regressions on the z-scored CLR latent columns and is scored by pigauto's existing Aitchison, CLR RMSE, simplex MAE, and dominant-component accuracy metrics.

## Status Counts

           method status Freq
         baseline     ok   18
             mean     ok   18
          pigauto     ok   18
     tabpfn_lappe     ok   18
 tabpfn_lappe_nfa     ok   18

## Test Summary

           method         scenario scale_n trait             type   rmse    mae
         baseline multi_high_phylo     150  comp multi_proportion     NA     NA
             mean multi_high_phylo     150  comp multi_proportion     NA     NA
          pigauto multi_high_phylo     150  comp multi_proportion     NA     NA
     tabpfn_lappe multi_high_phylo     150  comp multi_proportion     NA     NA
 tabpfn_lappe_nfa multi_high_phylo     150  comp multi_proportion     NA     NA
         baseline multi_high_phylo     300  comp multi_proportion     NA     NA
             mean multi_high_phylo     300  comp multi_proportion     NA     NA
          pigauto multi_high_phylo     300  comp multi_proportion     NA     NA
     tabpfn_lappe multi_high_phylo     300  comp multi_proportion     NA     NA
 tabpfn_lappe_nfa multi_high_phylo     300  comp multi_proportion     NA     NA
         baseline         multi_K8     150  comp multi_proportion     NA     NA
             mean         multi_K8     150  comp multi_proportion     NA     NA
          pigauto         multi_K8     150  comp multi_proportion     NA     NA
     tabpfn_lappe         multi_K8     150  comp multi_proportion     NA     NA
 tabpfn_lappe_nfa         multi_K8     150  comp multi_proportion     NA     NA
         baseline         multi_K8     300  comp multi_proportion     NA     NA
             mean         multi_K8     300  comp multi_proportion     NA     NA
          pigauto         multi_K8     300  comp multi_proportion     NA     NA
     tabpfn_lappe         multi_K8     300  comp multi_proportion     NA     NA
 tabpfn_lappe_nfa         multi_K8     300  comp multi_proportion     NA     NA
         baseline   multi_moderate     150  comp multi_proportion     NA     NA
             mean   multi_moderate     150  comp multi_proportion     NA     NA
          pigauto   multi_moderate     150  comp multi_proportion     NA     NA
     tabpfn_lappe   multi_moderate     150  comp multi_proportion     NA     NA
 tabpfn_lappe_nfa   multi_moderate     150  comp multi_proportion     NA     NA
         baseline   multi_moderate     300  comp multi_proportion     NA     NA
             mean   multi_moderate     300  comp multi_proportion     NA     NA
          pigauto   multi_moderate     300  comp multi_proportion     NA     NA
     tabpfn_lappe   multi_moderate     300  comp multi_proportion     NA     NA
 tabpfn_lappe_nfa   multi_moderate     300  comp multi_proportion     NA     NA
         baseline    zi_many_zeros     150   zi1         zi_count 17.675 17.331
             mean    zi_many_zeros     150   zi1         zi_count 11.971 10.293
          pigauto    zi_many_zeros     150   zi1         zi_count 17.675 17.331
     tabpfn_lappe    zi_many_zeros     150   zi1         zi_count 18.212 16.595
 tabpfn_lappe_nfa    zi_many_zeros     150   zi1         zi_count 22.136 19.051
         baseline    zi_many_zeros     150   zi2         zi_count 10.818  8.576
             mean    zi_many_zeros     150   zi2         zi_count 19.587 13.932
          pigauto    zi_many_zeros     150   zi2         zi_count 10.818  8.576
     tabpfn_lappe    zi_many_zeros     150   zi2         zi_count 12.614 10.343
 tabpfn_lappe_nfa    zi_many_zeros     150   zi2         zi_count 12.477  9.450
         baseline    zi_many_zeros     300   zi1         zi_count  9.741  8.173
             mean    zi_many_zeros     300   zi1         zi_count 13.695 10.458
          pigauto    zi_many_zeros     300   zi1         zi_count  9.741  8.173
     tabpfn_lappe    zi_many_zeros     300   zi1         zi_count  9.478  8.349
 tabpfn_lappe_nfa    zi_many_zeros     300   zi1         zi_count 19.067 15.804
         baseline    zi_many_zeros     300   zi2         zi_count 23.434 18.276
             mean    zi_many_zeros     300   zi2         zi_count 20.573 16.581
          pigauto    zi_many_zeros     300   zi2         zi_count 27.877 22.420
     tabpfn_lappe    zi_many_zeros     300   zi2         zi_count 24.891 18.670
 tabpfn_lappe_nfa    zi_many_zeros     300   zi2         zi_count 33.727 28.474
         baseline      zi_moderate     150   zi1         zi_count 14.065  9.745
             mean      zi_moderate     150   zi1         zi_count 17.641 12.812
          pigauto      zi_moderate     150   zi1         zi_count 14.065  9.745
     tabpfn_lappe      zi_moderate     150   zi1         zi_count 15.271 10.859
 tabpfn_lappe_nfa      zi_moderate     150   zi1         zi_count 14.353 10.418
         baseline      zi_moderate     150   zi2         zi_count  9.012  7.582
             mean      zi_moderate     150   zi2         zi_count  8.823  7.560
          pigauto      zi_moderate     150   zi2         zi_count  9.433  8.173
     tabpfn_lappe      zi_moderate     150   zi2         zi_count  9.765  8.165
 tabpfn_lappe_nfa      zi_moderate     150   zi2         zi_count 12.100 11.204
         baseline      zi_moderate     300   zi1         zi_count  8.996  5.462
             mean      zi_moderate     300   zi1         zi_count  8.928  5.392
          pigauto      zi_moderate     300   zi1         zi_count  9.153  5.516
     tabpfn_lappe      zi_moderate     300   zi1         zi_count  9.529  6.082
 tabpfn_lappe_nfa      zi_moderate     300   zi1         zi_count 10.425  6.535
         baseline      zi_moderate     300   zi2         zi_count 22.895 14.050
             mean      zi_moderate     300   zi2         zi_count 24.525 14.644
          pigauto      zi_moderate     300   zi2         zi_count 27.259 18.536
     tabpfn_lappe      zi_moderate     300   zi2         zi_count 26.784 18.168
 tabpfn_lappe_nfa      zi_moderate     300   zi2         zi_count 30.063 23.694
         baseline        zi_sparse     150   zi1         zi_count  3.519  2.617
             mean        zi_sparse     150   zi1         zi_count  4.921  3.800
          pigauto        zi_sparse     150   zi1         zi_count  7.236  4.471
     tabpfn_lappe        zi_sparse     150   zi1         zi_count  4.951  3.199
 tabpfn_lappe_nfa        zi_sparse     150   zi1         zi_count  5.597  3.958
         baseline        zi_sparse     150   zi2         zi_count 14.954 10.493
             mean        zi_sparse     150   zi2         zi_count  7.922  5.333
          pigauto        zi_sparse     150   zi2         zi_count 14.954 10.493
     tabpfn_lappe        zi_sparse     150   zi2         zi_count  7.877  5.083
 tabpfn_lappe_nfa        zi_sparse     150   zi2         zi_count  7.270  4.940
         baseline        zi_sparse     300   zi1         zi_count  9.002  5.434
             mean        zi_sparse     300   zi1         zi_count 10.514  6.852
          pigauto        zi_sparse     300   zi1         zi_count  9.002  5.434
     tabpfn_lappe        zi_sparse     300   zi1         zi_count 10.522  6.629
 tabpfn_lappe_nfa        zi_sparse     300   zi1         zi_count 10.923  7.082
         baseline        zi_sparse     300   zi2         zi_count 19.564 14.165
             mean        zi_sparse     300   zi2         zi_count 21.759 15.981
          pigauto        zi_sparse     300   zi2         zi_count 19.564 14.165
     tabpfn_lappe        zi_sparse     300   zi2         zi_count 19.231 14.187
 tabpfn_lappe_nfa        zi_sparse     300   zi2         zi_count 28.782 23.276
 zero_accuracy   brier aitchison rmse_clr simplex_mae accuracy
            NA      NA     2.455   0.7092     0.11439   0.5789
            NA      NA     3.773   1.0987     0.17132   0.3684
            NA      NA     2.455   0.7092     0.11439   0.5789
            NA      NA     2.857   0.8469     0.14500   0.3158
            NA      NA     3.983   1.1062     0.17276   0.5263
            NA      NA     2.432   0.5907     0.11200   0.6316
            NA      NA     4.091   0.9956     0.18935   0.2895
            NA      NA     2.432   0.5907     0.11200   0.6316
            NA      NA     2.666   0.6572     0.12400   0.6579
            NA      NA     3.325   0.8096     0.14915   0.5000
            NA      NA     5.657   0.8371     0.09739   0.6316
            NA      NA     7.087   1.0737     0.17095   0.2105
            NA      NA     5.657   0.8371     0.09739   0.6316
            NA      NA     5.962   0.8954     0.13674   0.4211
            NA      NA     7.573   1.1423     0.16421   0.2105
            NA      NA     7.089   0.9162     0.13467   0.3684
            NA      NA     7.919   0.9909     0.15389   0.2105
            NA      NA     7.089   0.9162     0.13467   0.3684
            NA      NA     6.995   0.8878     0.13181   0.3947
            NA      NA     8.094   1.0456     0.16390   0.2368
            NA      NA     4.359   0.9214     0.18314   0.5263
            NA      NA     4.147   0.9310     0.18343   0.5263
            NA      NA     4.359   0.9214     0.18314   0.5263
            NA      NA     4.093   0.9114     0.18551   0.5789
            NA      NA     4.989   1.1077     0.20653   0.4211
            NA      NA     4.470   0.8773     0.17025   0.5000
            NA      NA     5.137   0.9868     0.20756   0.4211
            NA      NA     4.470   0.8773     0.17025   0.5000
            NA      NA     4.388   0.8672     0.18175   0.4211
            NA      NA     5.781   1.1163     0.23894   0.3421
        0.8000 0.17496        NA       NA          NA       NA
        1.0000 0.20935        NA       NA          NA       NA
        0.8000 0.17496        NA       NA          NA       NA
        0.2000 0.52213        NA       NA          NA       NA
        0.4000 0.59594        NA       NA          NA       NA
        1.0000 0.09752        NA       NA          NA       NA
        1.0000 0.20515        NA       NA          NA       NA
        1.0000 0.09752        NA       NA          NA       NA
        0.6667 0.17089        NA       NA          NA       NA
        0.6667 0.24575        NA       NA          NA       NA
        0.6667 0.19247        NA       NA          NA       NA
        1.0000 0.20869        NA       NA          NA       NA
        0.6667 0.19247        NA       NA          NA       NA
        0.4444 0.33194        NA       NA          NA       NA
        0.1111 0.73938        NA       NA          NA       NA
        0.5000 0.27764        NA       NA          NA       NA
        1.0000 0.20564        NA       NA          NA       NA
        0.3750 0.46452        NA       NA          NA       NA
        0.6250 0.34988        NA       NA          NA       NA
        0.5000 0.50463        NA       NA          NA       NA
        1.0000 0.10877        NA       NA          NA       NA
        1.0000 0.14667        NA       NA          NA       NA
        1.0000 0.10877        NA       NA          NA       NA
        0.9000 0.17155        NA       NA          NA       NA
        0.8000 0.19972        NA       NA          NA       NA
        0.5556 0.24214        NA       NA          NA       NA
        1.0000 0.16222        NA       NA          NA       NA
        0.5556 0.27434        NA       NA          NA       NA
        0.5556 0.31606        NA       NA          NA       NA
        0.2222 0.68053        NA       NA          NA       NA
        0.4762 0.25488        NA       NA          NA       NA
        1.0000 0.16090        NA       NA          NA       NA
        0.4762 0.27079        NA       NA          NA       NA
        0.4762 0.34687        NA       NA          NA       NA
        0.5714 0.34474        NA       NA          NA       NA
        0.7143 0.22551        NA       NA          NA       NA
        1.0000 0.16192        NA       NA          NA       NA
        0.5000 0.31509        NA       NA          NA       NA
        0.2857 0.41312        NA       NA          NA       NA
        0.2857 0.59998        NA       NA          NA       NA
        0.8000 0.15879        NA       NA          NA       NA
        1.0000 0.16162        NA       NA          NA       NA
        0.8000 0.15426        NA       NA          NA       NA
        0.7000 0.30380        NA       NA          NA       NA
        0.6000 0.39324        NA       NA          NA       NA
        0.5000 0.30176        NA       NA          NA       NA
        1.0000 0.17734        NA       NA          NA       NA
        0.5000 0.30176        NA       NA          NA       NA
        0.3333 0.43768        NA       NA          NA       NA
        0.3333 0.56538        NA       NA          NA       NA
        0.8421 0.14826        NA       NA          NA       NA
        1.0000 0.15070        NA       NA          NA       NA
        0.8421 0.14826        NA       NA          NA       NA
        0.6842 0.15172        NA       NA          NA       NA
        0.5263 0.27141        NA       NA          NA       NA
        0.8824 0.12202        NA       NA          NA       NA
        1.0000 0.16837        NA       NA          NA       NA
        0.8824 0.12202        NA       NA          NA       NA
        0.8824 0.09175        NA       NA          NA       NA
        0.8235 0.13468        NA       NA          NA       NA
