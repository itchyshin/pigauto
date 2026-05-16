# BACE results snapshot

Pinned BACE benchmark outputs used by pigauto's
`.github/workflows/bace-headtohead.yml` CI workflow as the
comparison baseline.

## Provenance

- Source: `daniel1noble/BACE` Actions run **25329857467**
  (manual `workflow_dispatch` on 2026-05-13).
- BACE version: see Dan's commit at the time of the run.
- Datasets: AVONET, PanTHERIA, AmphiBIO, BIEN, GlobTherm, LepTraits.
- Config: BACE default (matches BACE's `dev/0[03-7]_benchmark_*.R`).
- Snapshot date: 2026-05-16.

## Schema

Each `<dataset>.rds` is a tibble with columns:

| col | type | notes |
|---|---|---|
| `dataset` | chr | one of {avonet, pantheria, amphibio, bien, globtherm, leptraits} |
| `trait` | chr | trait column name in BACE's bench output |
| `type` | chr | continuous / count / binary / ordinal / categorical / proportion / zi_count |
| `method` | chr | always `"BACE_snapshot"` |
| `imputation_idx` | int | 1..N_IMP (matches BACE's pool of imputations) |
| `rmse`, `mae`, `pearson_r` | dbl | continuous-family metrics, NA otherwise |
| `accuracy`, `brier` | dbl | discrete-family metrics, NA otherwise |
| `time_sec` | dbl | per-dataset wall time (replicated across rows) |

## Re-snapshotting

Manual, out-of-band. See `script/gha/snapshot_bace.R` and
`script/gha/README.md` for the procedure.
