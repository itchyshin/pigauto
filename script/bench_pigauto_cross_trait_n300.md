# Pigauto Cross-Trait Diagnostic

- Generated: 2026-06-11 05:29:04 MDT
- Git commit: `77c04a760bf678d2b4d90cfd119189a4888cf0ab`
- Scales: `300`
- Replicates: `3`
- Pigauto configs: `default, trait_attention, relaxed_gate`

## Aim

This benchmark tests whether pigauto-side cross-trait settings can close the continuous AVONET same-row TabPFN gap, and whether a simple linear same-row model already captures much of the signal.

## Status Counts

                  method status Freq
                baseline     ok   12
         linear_same_row     ok   12
   linear_same_row_lappe     ok   12
         pigauto_default     ok   12
    pigauto_relaxed_gate     ok   12
 pigauto_trait_attention     ok   12

## Mean RMSE By Method

                  method scale_n      rmse pearson_r coverage_95  wall_sec
    pigauto_relaxed_gate     300 4.682e-01    0.8914      0.9919 1.948e+02
 pigauto_trait_attention     300 4.932e-01    0.8765      0.9921 2.365e+02
         pigauto_default     300 5.023e-01    0.8774      0.9901 2.178e+02
         linear_same_row     300 5.054e-01    0.8650          NA 1.583e-03
                baseline     300 5.452e-01    0.8581          NA 5.230e-01
   linear_same_row_lappe     300 3.888e+13    0.6920          NA 3.000e-03

## Mean RMSE By Trait

                  method scale_n              trait      rmse pearson_r
    pigauto_relaxed_gate     300 Beak.Length_Culmen 6.629e-01    0.7843
         pigauto_default     300 Beak.Length_Culmen 6.689e-01    0.7801
                baseline     300 Beak.Length_Culmen 6.711e-01    0.7803
         linear_same_row     300 Beak.Length_Culmen 6.754e-01    0.7654
 pigauto_trait_attention     300 Beak.Length_Culmen 6.949e-01    0.7592
   linear_same_row_lappe     300 Beak.Length_Culmen 7.077e-01    0.7299
         pigauto_default     300               Mass 3.201e-01    0.9525
    pigauto_relaxed_gate     300               Mass 3.226e-01    0.9523
 pigauto_trait_attention     300               Mass 3.365e-01    0.9527
         linear_same_row     300               Mass 4.081e-01    0.9027
                baseline     300               Mass 4.227e-01    0.9183
   linear_same_row_lappe     300               Mass 1.417e+14    0.6756
    pigauto_relaxed_gate     300      Tarsus.Length 4.982e-01    0.8877
         pigauto_default     300      Tarsus.Length 5.065e-01    0.8828
         linear_same_row     300      Tarsus.Length 5.115e-01    0.8764
 pigauto_trait_attention     300      Tarsus.Length 5.117e-01    0.8781
                baseline     300      Tarsus.Length 5.119e-01    0.8777
   linear_same_row_lappe     300      Tarsus.Length 1.379e+13    0.4496
    pigauto_relaxed_gate     300        Wing.Length 3.893e-01    0.9412
         linear_same_row     300        Wing.Length 4.264e-01    0.9155
 pigauto_trait_attention     300        Wing.Length 4.296e-01    0.9160
   linear_same_row_lappe     300        Wing.Length 4.381e-01    0.9128
         pigauto_default     300        Wing.Length 5.138e-01    0.8942
                baseline     300        Wing.Length 5.751e-01    0.8561
 coverage_95
      0.9756
      0.9756
          NA
          NA
      0.9837
          NA
      1.0000
      1.0000
      1.0000
          NA
          NA
          NA
      0.9921
      0.9921
          NA
      0.9921
          NA
          NA
      1.0000
          NA
      0.9926
          NA
      0.9926
          NA

## Gate Audit

                  method scale_n              trait   rmse r_cal_bm r_cal_gnn
         pigauto_default     300 Beak.Length_Culmen 0.6689   0.6852   0.31481
    pigauto_relaxed_gate     300 Beak.Length_Culmen 0.6629   0.6667   0.33333
 pigauto_trait_attention     300 Beak.Length_Culmen 0.6949   0.7963   0.16667
         pigauto_default     300               Mass 0.3201   0.4167   0.56667
    pigauto_relaxed_gate     300               Mass 0.3226   0.4833   0.51667
 pigauto_trait_attention     300               Mass 0.3365   0.5000   0.45000
         pigauto_default     300      Tarsus.Length 0.5065   0.6500   0.35000
    pigauto_relaxed_gate     300      Tarsus.Length 0.4982   0.8500   0.15000
 pigauto_trait_attention     300      Tarsus.Length 0.5117   0.9333   0.06667
         pigauto_default     300        Wing.Length 0.5138   0.4833   0.43333
    pigauto_relaxed_gate     300        Wing.Length 0.3893   0.1000   0.86491
 pigauto_trait_attention     300        Wing.Length 0.4296   0.2833   0.69912
 r_cal_mean
    0.00000
    0.00000
    0.03704
    0.01667
    0.00000
    0.05000
    0.00000
    0.00000
    0.00000
    0.08333
    0.03509
    0.01754
