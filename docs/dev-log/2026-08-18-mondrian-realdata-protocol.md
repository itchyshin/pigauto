# Stage C: real-data Mondrian confirmation protocol

Independently of Stage B, FishBase and PanTHERIA will mask only cells that were
originally observed.  For each trait, compare `conformal_method = "split"` and
`"mondrian"`, retaining masked-cell count, calibration/validation count,
activation or fallback, coverage, interval width, MCSE, runtime, and errors.

The only promotable result is empirical calibration improvement with transparent
fallback and a stated width cost.  This study cannot establish an unconditional
95% guarantee or justify changing the default.

## Execution contract

`script/mondrian_confirmation/00_prepare_realdata_input.R` prepares the
canonical PanTHERIA or FishBase input from the existing source caches. It fails
closed when source data or a tree are unavailable, labels cannot align exactly,
or a retained trait has fewer than 20 observed cells; it never substitutes a
synthetic data set or fallback tree. The runner
`script/mondrian_confirmation/01_run_masked_confirmation.R` accepts that
prepared `list(data, tree, dataset)` RDS input, masks only originally observed cells,
and retains a `mask_receipt.rds`.  It runs `conformal_method = "split"` and
`"mondrian"` with the identical seed and mask, retaining method-specific RDS
receipts containing trait-wise masked-cell counts, interval counts, coverage,
width, elapsed time, errors, and Mondrian activation/fallback metadata.
`02_summarise_masked_confirmation.R` fails closed if either method errored and
otherwise writes the paired reporting table.  Prepare FishBase and PanTHERIA
inputs using their existing canonical loaders; do not use naturally missing
cells as truth.  Execution remains separately compute-gated.

## Structural smoke

The runner was exercised locally on a 30-tip, one-trait AVONET subset with a
shared six-cell observed-only mask and five epochs.  It wrote the mask receipt,
both method receipts, and the combined table successfully.  The expected
three-cell validation warnings make this a harness check only: it is neither
a FishBase/PanTHERIA result nor evidence about interval coverage.
