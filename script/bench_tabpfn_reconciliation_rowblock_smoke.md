# TabPFN Reconciliation Benchmark

- Generated: 2026-06-11 06:25:44 MDT
- Git commit: `e0329d8503f2cbc8a79997a8182d8914f6db8dba`
- Scales: `50`
- Regimes: `same_row, same_row_lappe, phylo_only`
- Replicates: `1`
- Split mode: `row_block_all`
- Dry run: `TRUE`
- Pigauto configs: `default`

## Aim

This benchmark uses row-block holdouts so validation/test species rows have their selected trait cells hidden together. In `row_block_all`, every latent column for held-out rows is masked, so same-row TabPFN features cannot borrow observed cells from the same species at prediction time.

## Status Counts

                method  status Freq
     tabpfn_phylo_only dry_run    1
       tabpfn_same_row dry_run    1
 tabpfn_same_row_lappe dry_run    1

## Test Summary

No scored rows were produced. Inspect `Status Counts`.
