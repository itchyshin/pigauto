# Real-data preview (9 of 10 cells; G8 pending the FishBase cell)

All nine finished cells ran at commit 69670d4 on Totoro (`receipts.csv`, column `code_sha`). The
FishBase cell (`fishbase-structured-m20260818`) had not written a final receipt when these tables were
made, so it is absent from every table. G8 is not evaluated here.

| File | What it holds |
|---|---|
| `coverage_table.csv` | Per-trait coverage of the masked cells: model-based, split conformal and Mondrian conformal (36 trait-cells: PanTHERIA 6 masks x 4 traits, AVONET 3 masks x 4 traits). |
| `slope_table.csv` | Per cell and analysis pair: complete-row reference PGLS slope and SE, MI-pooled slope and SE, `diff_ref_se` = (MI - reference) / reference SE, `rel_diff`. |
| `pair_summary.csv` | The same per analysis pair over cells (reported 5% criterion, not gated). |
| `convergence_table.csv` | Max split R-hat and min bulk ESS per cell. |
| `summary.md` | The printed summary of all four tables. |
| `receipts.csv` | One row per planned cell from its receipt (`mi_posterior.rds`): status, species, kept traits, code SHA. |

Provenance: `script/mi_realdata/02_summarise.R` (at the commit that adds this folder) run on the ten
receipt folders reproduces the first five files byte for byte (checked 2026-09-24). `receipts.csv` was
extracted from the same receipts with `readRDS()`. The receipts themselves are not committed.
