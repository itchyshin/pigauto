# Pigauto Cross-Trait Diagnostic

- Generated: 2026-06-11 04:55:46 MDT
- Git commit: `77c04a760bf678d2b4d90cfd119189a4888cf0ab`
- Scales: `50`
- Replicates: `1`
- Pigauto configs: `default, trait_attention`

## Aim

This benchmark tests whether pigauto-side cross-trait settings can close the continuous AVONET same-row TabPFN gap, and whether a simple linear same-row model already captures much of the signal.

## Status Counts

                  method status Freq
                baseline     ok    1
         linear_same_row     ok    1
   linear_same_row_lappe     ok    1
         pigauto_default     ok    1
 pigauto_trait_attention     ok    1

## Mean RMSE By Method

                  method scale_n   rmse pearson_r coverage_95 wall_sec
         linear_same_row      50 0.5930    0.9955          NA    0.003
   linear_same_row_lappe      50 0.6898    0.9948          NA    0.002
                baseline      50 0.8898    0.9092          NA    0.729
         pigauto_default      50 1.1820        NA         0.8   42.900
 pigauto_trait_attention      50 1.1820        NA         0.8   47.845

## Mean RMSE By Trait

                  method scale_n trait   rmse pearson_r coverage_95
         linear_same_row      50  Mass 0.5930    0.9955          NA
   linear_same_row_lappe      50  Mass 0.6898    0.9948          NA
                baseline      50  Mass 0.8898    0.9092          NA
         pigauto_default      50  Mass 1.1820        NA         0.8
 pigauto_trait_attention      50  Mass 1.1820        NA         0.8

## Gate Audit

                  method scale_n trait  rmse r_cal_bm r_cal_gnn r_cal_mean
         pigauto_default      50  Mass 1.182        0         0          1
 pigauto_trait_attention      50  Mass 1.182        0         0          1
