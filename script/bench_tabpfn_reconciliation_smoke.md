# TabPFN Reconciliation Benchmark

- Generated: 2026-06-10 21:51:06 MDT
- Git commit: `6cec7178407200a08199d1efeeaa894af22010f1`
- Scales: `50`
- Regimes: `same_row, phylo_only, shuffled_same_row_lappe`
- Replicates: `1`
- Dry run: `FALSE`

## Aim

This benchmark tests whether Russell-style TabPFN gains are mostly explained by same-row cross-trait features, phylogenetic features, or their combination. The shuffled regime keeps the same feature distribution but breaks row-level cross-trait alignment.

## Status Counts

                         method status Freq
              tabpfn_phylo_only     ok    1
                tabpfn_same_row     ok    1
 tabpfn_shuffled_same_row_lappe     ok    1

## Mean RMSE By Method

                         method scale_n   rmse pearson_r coverage_95 wall_sec
                tabpfn_same_row      50 0.3429    0.9811         0.8    9.216
              tabpfn_phylo_only      50 0.7045    0.9901         1.0    8.126
 tabpfn_shuffled_same_row_lappe      50 0.7301    0.9924         1.0    7.493

## Mean RMSE By Trait

                         method scale_n trait   rmse pearson_r coverage_95
                tabpfn_same_row      50  Mass 0.3429    0.9811         0.8
              tabpfn_phylo_only      50  Mass 0.7045    0.9901         1.0
 tabpfn_shuffled_same_row_lappe      50  Mass 0.7301    0.9924         1.0

