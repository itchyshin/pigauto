# TabPFN Simulated Mixed Scalar-Type Benchmark

- Generated: 2026-06-10 18:34:28 MDT
- Git commit: `f6c939033e5c97e11131a25671f49c4378d2f4e9`
- Scales: `50`
- Scenarios: `mixed_moderate`
- Variants: `lappe`
- Replicates: `1`
- Dry run: `FALSE`

## Design

Each cell simulates one mixed scalar dataset with two continuous traits and one trait each for count, proportion, binary, ordinal, and categorical data. Missing cells are held out with pigauto's `make_missing_splits()`. TabPFN uses regression for continuous/count/ proportion latent targets and classification for binary/ordinal/ categorical targets. This branch-local check does not cover `zi_count`, `multi_proportion`, classification prediction sets, multi-tree Rubin pooling, multi-observation covariates, or active-imputation guidance.

## Status Counts

       method status Freq
     baseline     ok    7
    mean_mode     ok    7
 tabpfn_lappe     ok    7

## Test Summary

       method       scenario scale_n         trait        type           task
     baseline mixed_moderate      50      binary_1      binary classification
    mean_mode mixed_moderate      50      binary_1      binary classification
 tabpfn_lappe mixed_moderate      50      binary_1      binary classification
     baseline mixed_moderate      50 categorical_1 categorical classification
    mean_mode mixed_moderate      50 categorical_1 categorical classification
 tabpfn_lappe mixed_moderate      50 categorical_1 categorical classification
     baseline mixed_moderate      50  continuous_1  continuous     regression
    mean_mode mixed_moderate      50  continuous_1  continuous     regression
 tabpfn_lappe mixed_moderate      50  continuous_1  continuous     regression
     baseline mixed_moderate      50  continuous_2  continuous     regression
    mean_mode mixed_moderate      50  continuous_2  continuous     regression
 tabpfn_lappe mixed_moderate      50  continuous_2  continuous     regression
     baseline mixed_moderate      50       count_1       count     regression
    mean_mode mixed_moderate      50       count_1       count     regression
 tabpfn_lappe mixed_moderate      50       count_1       count     regression
     baseline mixed_moderate      50     ordinal_1     ordinal classification
    mean_mode mixed_moderate      50     ordinal_1     ordinal classification
 tabpfn_lappe mixed_moderate      50     ordinal_1     ordinal classification
     baseline mixed_moderate      50  proportion_1  proportion     regression
    mean_mode mixed_moderate      50  proportion_1  proportion     regression
 tabpfn_lappe mixed_moderate      50  proportion_1  proportion     regression
   rmse pearson_r accuracy balanced_accuracy
     NA        NA   0.5000           0.50000
     NA        NA   1.0000           1.00000
     NA        NA   0.0000           0.00000
     NA        NA   0.6667           0.55556
     NA        NA   0.3333           0.33333
     NA        NA   0.8333           0.66667
 0.7189   0.81810       NA                NA
 1.2912        NA       NA                NA
 1.2254   0.50412       NA                NA
 0.5476   0.77646       NA                NA
 0.9061        NA       NA                NA
 0.6327   0.70186       NA                NA
 0.9734   0.30459       NA                NA
 1.3979        NA       NA                NA
 1.0172   0.46985       NA                NA
     NA        NA   0.0000           0.00000
     NA        NA   0.0000           0.00000
     NA        NA   0.1667           0.08333
 1.0845  -0.05877       NA                NA
 1.0477        NA       NA                NA
 1.1214  -0.83247       NA                NA
