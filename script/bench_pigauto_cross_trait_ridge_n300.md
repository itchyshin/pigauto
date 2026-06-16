# Pigauto Cross-Trait Diagnostic

- Generated: 2026-06-11 05:33:51 MDT
- Git commit: `77c04a760bf678d2b4d90cfd119189a4888cf0ab`
- Scales: `300`
- Replicates: `3`
- Run pigauto configs: `FALSE`
- Ridge lambda: `0.001`
- Pigauto configs: `default`

## Aim

This benchmark tests whether pigauto-side cross-trait settings can close the continuous AVONET same-row TabPFN gap, and whether a ridge-stabilized same-row model already captures much of the signal.

## Status Counts

               method status Freq
             baseline     ok   12
       ridge_same_row     ok   12
 ridge_same_row_lappe     ok   12

## Mean RMSE By Method

               method scale_n      rmse pearson_r coverage_95  wall_sec
       ridge_same_row     300 5.054e-01    0.8650          NA 0.0006667
             baseline     300 5.452e-01    0.8581          NA 0.5313333
 ridge_same_row_lappe     300 3.543e+09    0.7082          NA 0.0004167

## Mean RMSE By Trait

               method scale_n              trait      rmse pearson_r
             baseline     300 Beak.Length_Culmen 6.711e-01    0.7803
       ridge_same_row     300 Beak.Length_Culmen 6.754e-01    0.7654
 ridge_same_row_lappe     300 Beak.Length_Culmen 7.851e-01    0.7139
       ridge_same_row     300               Mass 4.081e-01    0.9027
             baseline     300               Mass 4.227e-01    0.9183
 ridge_same_row_lappe     300               Mass 3.097e+07    0.6714
       ridge_same_row     300      Tarsus.Length 5.115e-01    0.8764
             baseline     300      Tarsus.Length 5.119e-01    0.8777
 ridge_same_row_lappe     300      Tarsus.Length 1.414e+10    0.5326
       ridge_same_row     300        Wing.Length 4.264e-01    0.9155
 ridge_same_row_lappe     300        Wing.Length 4.316e-01    0.9149
             baseline     300        Wing.Length 5.751e-01    0.8561
 coverage_95
          NA
          NA
          NA
          NA
          NA
          NA
          NA
          NA
          NA
          NA
          NA
          NA

