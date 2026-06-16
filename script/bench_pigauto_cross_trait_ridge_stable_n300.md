# Pigauto Cross-Trait Diagnostic

- Generated: 2026-06-11 05:38:13 MDT
- Git commit: `77c04a760bf678d2b4d90cfd119189a4888cf0ab`
- Scales: `300`
- Replicates: `3`
- Run pigauto configs: `FALSE`
- Ridge lambda: `1`
- Ridge min SD: `1e-06`
- Pigauto configs: `default`

## Aim

This benchmark tests whether pigauto-side cross-trait settings can close the continuous AVONET same-row TabPFN gap, and whether a ridge-stabilized same-row model already captures much of the signal.

## Status Counts

               method status Freq
             baseline     ok   12
       ridge_same_row     ok   12
 ridge_same_row_lappe     ok   12

## Mean RMSE By Method

               method scale_n   rmse pearson_r coverage_95  wall_sec
       ridge_same_row     300 0.5052    0.8651          NA 0.0010000
 ridge_same_row_lappe     300 0.5378    0.8508          NA 0.0005833
             baseline     300 0.5452    0.8581          NA 0.5206667

## Mean RMSE By Trait

               method scale_n              trait   rmse pearson_r coverage_95
             baseline     300 Beak.Length_Culmen 0.6711    0.7803          NA
       ridge_same_row     300 Beak.Length_Culmen 0.6749    0.7658          NA
 ridge_same_row_lappe     300 Beak.Length_Culmen 0.7829    0.7147          NA
       ridge_same_row     300               Mass 0.4078    0.9029          NA
 ridge_same_row_lappe     300               Mass 0.4162    0.9031          NA
             baseline     300               Mass 0.4227    0.9183          NA
       ridge_same_row     300      Tarsus.Length 0.5118    0.8765          NA
             baseline     300      Tarsus.Length 0.5119    0.8777          NA
 ridge_same_row_lappe     300      Tarsus.Length 0.5210    0.8705          NA
       ridge_same_row     300        Wing.Length 0.4264    0.9155          NA
 ridge_same_row_lappe     300        Wing.Length 0.4311    0.9149          NA
             baseline     300        Wing.Length 0.5751    0.8561          NA

