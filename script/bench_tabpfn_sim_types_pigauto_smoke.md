# TabPFN Simulated Mixed Scalar-Type Benchmark

- Generated: 2026-06-10 18:37:58 MDT
- Git commit: `f6c939033e5c97e11131a25671f49c4378d2f4e9`
- Scales: `75`
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
      pigauto     ok    7
 tabpfn_lappe     ok    7

## Test Summary

       method       scenario scale_n         trait        type           task
     baseline mixed_moderate      75      binary_1      binary classification
    mean_mode mixed_moderate      75      binary_1      binary classification
      pigauto mixed_moderate      75      binary_1      binary classification
 tabpfn_lappe mixed_moderate      75      binary_1      binary classification
     baseline mixed_moderate      75 categorical_1 categorical classification
    mean_mode mixed_moderate      75 categorical_1 categorical classification
      pigauto mixed_moderate      75 categorical_1 categorical classification
 tabpfn_lappe mixed_moderate      75 categorical_1 categorical classification
     baseline mixed_moderate      75  continuous_1  continuous     regression
    mean_mode mixed_moderate      75  continuous_1  continuous     regression
      pigauto mixed_moderate      75  continuous_1  continuous     regression
 tabpfn_lappe mixed_moderate      75  continuous_1  continuous     regression
     baseline mixed_moderate      75  continuous_2  continuous     regression
    mean_mode mixed_moderate      75  continuous_2  continuous     regression
      pigauto mixed_moderate      75  continuous_2  continuous     regression
 tabpfn_lappe mixed_moderate      75  continuous_2  continuous     regression
     baseline mixed_moderate      75       count_1       count     regression
    mean_mode mixed_moderate      75       count_1       count     regression
      pigauto mixed_moderate      75       count_1       count     regression
 tabpfn_lappe mixed_moderate      75       count_1       count     regression
     baseline mixed_moderate      75     ordinal_1     ordinal classification
    mean_mode mixed_moderate      75     ordinal_1     ordinal classification
      pigauto mixed_moderate      75     ordinal_1     ordinal classification
 tabpfn_lappe mixed_moderate      75     ordinal_1     ordinal classification
     baseline mixed_moderate      75  proportion_1  proportion     regression
    mean_mode mixed_moderate      75  proportion_1  proportion     regression
      pigauto mixed_moderate      75  proportion_1  proportion     regression
 tabpfn_lappe mixed_moderate      75  proportion_1  proportion     regression
   rmse pearson_r accuracy balanced_accuracy
     NA        NA   1.0000            1.0000
     NA        NA   0.6429            0.5000
     NA        NA   0.9286            0.9444
     NA        NA   0.6429            0.6333
     NA        NA   0.5000            0.5000
     NA        NA   0.5000            0.5000
     NA        NA   0.5000            0.5000
     NA        NA   0.7500            0.7500
 0.6383    0.8586       NA                NA
 1.1883        NA       NA                NA
 0.6383    0.8586       NA                NA
 0.8787    0.7549       NA                NA
 0.5599    0.7745       NA                NA
 1.0456        NA       NA                NA
 0.5599    0.7745       NA                NA
 0.6724    0.7432       NA                NA
 0.5844    0.8200       NA                NA
 1.1986        NA       NA                NA
 0.5844    0.8200       NA                NA
 0.7116    0.8704       NA                NA
     NA        NA   0.3846            0.4000
     NA        NA   0.1538            0.2000
     NA        NA   0.3846            0.4000
     NA        NA   0.2308            0.2333
 0.7826    0.7618       NA                NA
 0.8754        NA       NA                NA
 0.7826    0.7618       NA                NA
 0.6733    0.6406       NA                NA
