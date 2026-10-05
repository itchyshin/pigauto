# PanTHERIA: `log_transform = FALSE` sensitivity (2026-10-05)

Requested by Shinichi ("run the log_transform = FALSE sensitivity on PanTHERIA"). Follow-up named in
`../results.md` (real-data section): PanTHERIA's columns are already log values, and pigauto's default
`log_transform = TRUE` logs every all-positive continuous column again, so the imputation model was linear on a
log(log) scale while the pre-registered analysis (`script/mi_realdata/pairs.R`) is on the log scale.

## Run

- Code: this branch at dbe3304 (main 45011fc plus the opt-in `MI_REALDATA_LOG_TRANSFORM` in
  `script/mi_realdata/01_run.R`; unset keeps pigauto's default and the receipt records the value).
- The 6 PanTHERIA cells of G8 (4,027 species, 4 traits; 3 MCAR and 3 clade-structured masks), the same masks as
  before, pigauto defaults (residual prior "sep"), `MI_REALDATA_LOG_TRANSFORM=FALSE`. Launcher `launch_logtf.sh`.
- Totoro, 6 cores, 05:26 to 06:18 (52 min; estimate 1.5 h). 6 of 6 ok, all converged (max R-hat 1.012, min ESS 409).
- Smoke first: one cell with a 40-sweep sampler confirmed the switch reaches `multi_impute()` and the receipt.

Comparison: the same 6 cells under the default (`log_transform = TRUE`), same code otherwise (main 45011fc, "sep"),
run 2026-10-04 (`script/rubin_study/data/realdata_sep/` on `arc/rubin-freq-bace`). Only `log_transform` differs.
`compare_logtf.R` reproduces `compare_logtf.txt` from the two output folders.

## Result: the double log caused the slope shifts

Pooled slope minus the complete-row reference, in reference SEs, 3 masks per row:

| pair | mask | log_transform = TRUE | FALSE |
|---|---|---|---|
| body mass ~ head-body length | MCAR | +2.6 to +3.9 | -0.9 to +0.6 |
| body mass ~ head-body length | structured | -3.2 to -0.4 | -1.1 to +0.9 |
| gestation ~ body mass | MCAR | -2.1 to -0.9 | +0.2 to +1.0 |
| gestation ~ body mass | structured | 0.0 to +1.1 | +1.0 to +1.8 |
| longevity ~ body mass | MCAR | -3.1 to -0.7 | -0.1 to +0.7 |
| longevity ~ body mass | structured | -1.1 to -0.1 | -0.7 to +0.8 |

- Mean absolute shift 0.76 reference SEs (TRUE: 1.63); largest 1.75 (TRUE: 3.92); 14 of 18 pair-cells within 5% of
  the reference (TRUE: 11).
- The sign pattern by mask type for body mass ~ head-body length is gone. Its pooled SE is also smaller (0.025 to
  0.030 against 0.034 to 0.043).
- Remaining: gestation ~ body mass sits above the reference in all 6 cells (+0.2 to +1.8 SEs), within 2 SEs.

Per-trait model coverage (mean of 6 cells): body mass 0.963 (TRUE 0.916), head-body length 0.946 (0.941), gestation
0.934 (0.933), longevity 0.937 (0.938); overall 0.945 (0.932) against split conformal 0.948. Model coverage is below
split conformal in 11 of 24 trait-cells (TRUE: 18).

## What this means for users

`log_transform = TRUE` is pigauto's default and logs every all-positive continuous column. Data that are already on
a log scale are then imputed on log(log), which here shifted downstream slopes by up to 4 SEs in a direction that
depended on the missingness pattern. Users with already-logged columns should pass `log_transform = FALSE`. Whether
pigauto should detect or warn about this is a separate package decision; not done here.

## Does NOT cover

AVONET and FishBase (raw measurements; the default is the right scale there); a coverage study (one mask set per
dataset); the G8 gate (6 cells only, and G8 is defined at pigauto's defaults).
