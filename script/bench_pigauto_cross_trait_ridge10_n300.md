# Pigauto Cross-Trait Diagnostic

- Generated: 2026-06-11 05:34:48 MDT
- Git commit: `77c04a760bf678d2b4d90cfd119189a4888cf0ab`
- Scales: `300`
- Replicates: `3`
- Run pigauto configs: `FALSE`
- Ridge lambda: `10`
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
       ridge_same_row     300 5.048e-01    0.8660          NA 0.0008333
             baseline     300 5.452e-01    0.8581          NA 0.5376667
 ridge_same_row_lappe     300 3.795e+07    0.7303          NA 0.0005000

## Mean RMSE By Trait

               method scale_n              trait      rmse pearson_r
             baseline     300 Beak.Length_Culmen 6.711e-01    0.7803
       ridge_same_row     300 Beak.Length_Culmen 6.712e-01    0.7682
 ridge_same_row_lappe     300 Beak.Length_Culmen 7.660e-01    0.7206
       ridge_same_row     300               Mass 4.061e-01    0.9041
             baseline     300               Mass 4.227e-01    0.9183
 ridge_same_row_lappe     300               Mass 3.150e+07    0.6726
             baseline     300      Tarsus.Length 5.119e-01    0.8777
       ridge_same_row     300      Tarsus.Length 5.149e-01    0.8765
 ridge_same_row_lappe     300      Tarsus.Length 1.203e+08    0.6133
       ridge_same_row     300        Wing.Length 4.269e-01    0.9151
 ridge_same_row_lappe     300        Wing.Length 4.292e-01    0.9146
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

