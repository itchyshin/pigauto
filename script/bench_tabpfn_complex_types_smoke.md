# TabPFN Complex-Type Benchmark

- Generated: 2026-06-10 19:52:02 MDT
- Git commit: `b6dd2d752960db2abd63e3e9e0aad6978c73d034`
- Scales: `100`
- Scenarios: `zi_moderate, multi_moderate`
- Variants: `lappe`
- Replicates: `1`

## Design

ZI-count uses one TabPFN classifier for the non-zero gate and one TabPFN regressor for the conditional magnitude. Multi-proportion uses independent TabPFN regressions on the z-scored CLR latent columns and is scored by pigauto's existing Aitchison, CLR RMSE, simplex MAE, and dominant-component accuracy metrics.

## Status Counts

       method status Freq
     baseline     ok    3
         mean     ok    3
 tabpfn_lappe     ok    3

## Test Summary

       method       scenario scale_n trait             type  rmse   mae
     baseline multi_moderate     100  comp multi_proportion    NA    NA
         mean multi_moderate     100  comp multi_proportion    NA    NA
 tabpfn_lappe multi_moderate     100  comp multi_proportion    NA    NA
     baseline    zi_moderate     100   zi1         zi_count 23.75 19.25
         mean    zi_moderate     100   zi1         zi_count 42.05 27.59
 tabpfn_lappe    zi_moderate     100   zi1         zi_count 34.36 24.67
     baseline    zi_moderate     100   zi2         zi_count 15.33 14.73
         mean    zi_moderate     100   zi2         zi_count 20.41 14.65
 tabpfn_lappe    zi_moderate     100   zi2         zi_count 14.61 12.96
 zero_accuracy   brier aitchison rmse_clr simplex_mae accuracy
            NA      NA     3.930   0.8109      0.1520   0.6154
            NA      NA     4.955   1.0086      0.2134   0.3077
            NA      NA     3.896   0.8213      0.1289   0.5385
         0.375 0.28424        NA       NA          NA       NA
         1.000 0.17192        NA       NA          NA       NA
         0.375 0.31942        NA       NA          NA       NA
         1.000 0.08673        NA       NA          NA       NA
         1.000 0.15471        NA       NA          NA       NA
         1.000 0.01244        NA       NA          NA       NA
