# #189 acceptance sweep under residual_prior = "sep" (PR #204)

2026-10-03. Approved by Shinichi ("go validation"). Totoro, 40 cores, 05:21 to 18:38.

- Code: `script/mi_gls/` at ebcd466 (this branch), pigauto cd3a4cb, env `MI_POST_RESIDUAL_PRIOR=sep`; each cell
  records `residual_prior` in `sampler_control`.
- Grid: the same 40 regimes x 200 reps (8,000 cells), seeds and masks as the "iw" sweep of record
  (`../sim_summary.csv`, `../cell_coverage.csv`). `CAMPAIGN_DONE ran=8000 ok=8000 failed=0`.
- Raw files on Totoro: `~/hsq_work/mi_post/out_sep_ebcd4666b0fd4f46bc2429b858d540c7f1351392/`.

Files here:
- `sim_summary.csv`, `sim_summary.md`, `cell_coverage.csv`: `03_summarise_v2.R` output (`summarise.log`:
  `SETTINGS 0 of 8000`, one code sha).
- `g6.log`: `04_acceptance.R` on `sim_summary.csv`, ends `SIM_ACCEPT_PASS`.
- `g7.log`: `05_cell_coverage.R` on `cell_coverage.csv`, ends `CELL_COVERAGE_PASS`.
- `compare_iw_sep.R` / `.txt`: posterior_full rows, "iw" against "sep". Run from the repo root with a directory
  holding `sim_summary_iw.csv` and `cell_coverage_iw.csv` (copies of `../sim_summary.csv`, `../cell_coverage.csv`) and
  this folder's two CSVs.

## Gate outcomes, "iw" (of record) against "sep"

| | iw | sep |
|---|---|---|
| G6 | pass | pass |
| pooled relative SE ratio, 24 gated phylolm rows, band [0.95, 1.10] | 1.059 | 1.077 |
| gated rows outside [0.90, 1.15] (reported, not gated) | 35, 36, 38 | 21, 35, 38 (1.19, 1.18, 1.17) |
| gated mean abs paired bias / max | 0.0052 / 0.0142 | 0.0051 / 0.0161 |
| gated coverage range (complete minus MI, worst) | 0.845 to 0.975 (0.025) | 0.865 to 0.970 (0.030) |
| converged fits, gated + stress | 15,998 of 16,000 | 16,000 of 16,000 |
| G7 | pass | pass |
| gated MCAR per-cell coverage range, band [0.92, 0.98] | 0.933 to 0.955 | 0.945 to 0.956 |
| stress rule violations (G6, not gated) | 2 | 2 (bias, regimes 1 and 9 phylolm) |

The small negative slope bias in the 16 lambda = 1 twin rows (the caveat in `../results.md`) remains under "sep": -0.002
to -0.016 (iw -0.004 to -0.014). It shrinks at n = 1000 (by about 30%) and grows slightly at n = 300 (by about 15%).
This is a different and much smaller effect than the Rubin study's lambda = 1, rho = 0.5 bias (-0.048 at n = 1000 under
"iw"), which "sep" removes (see `docs/dev-log/arc/2026-10-02-rubin-pigauto-campaign.md` on `arc/rubin-freq-bace`).
