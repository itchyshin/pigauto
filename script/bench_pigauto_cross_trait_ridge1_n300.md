# Pigauto Cross-Trait Diagnostic

- Generated: 2026-06-11 05:34:46 MDT
- Git commit: `77c04a760bf678d2b4d90cfd119189a4888cf0ab`
- Scales: `300`
- Replicates: `3`
- Run pigauto configs: `FALSE`
- Ridge lambda: `1`
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
       ridge_same_row     300 5.052e-01    0.8651          NA 0.0009167
             baseline     300 5.452e-01    0.8581          NA 0.5280000
 ridge_same_row_lappe     300 3.164e+07    0.7285          NA 0.0008333

## Mean RMSE By Trait

               method scale_n              trait      rmse pearson_r
             baseline     300 Beak.Length_Culmen 6.711e-01    0.7803
       ridge_same_row     300 Beak.Length_Culmen 6.749e-01    0.7658
 ridge_same_row_lappe     300 Beak.Length_Culmen 7.829e-01    0.7147
       ridge_same_row     300               Mass 4.078e-01    0.9029
             baseline     300               Mass 4.227e-01    0.9183
 ridge_same_row_lappe     300               Mass 2.789e+07    0.6716
       ridge_same_row     300      Tarsus.Length 5.118e-01    0.8765
             baseline     300      Tarsus.Length 5.119e-01    0.8777
 ridge_same_row_lappe     300      Tarsus.Length 9.868e+07    0.6129
       ridge_same_row     300        Wing.Length 4.264e-01    0.9155
 ridge_same_row_lappe     300        Wing.Length 4.311e-01    0.9149
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

