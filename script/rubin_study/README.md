# Rubin simulation study: one folder for the publication

Everything needed to report the freq-vs-BACE multiple-imputation simulation: the write-up and every result table.

- `study.qmd`: design (ADEMP), computing (clusters, nodes, wall times, versions), results, failures, limitations,
  and the planned discrete extension. Render with `quarto render script/rubin_study/study.qmd`; `study.html` is
  the rendered copy (self-contained).
- `data/`: per-fit tables (`fits.csv`, `fit_cells.csv.gz`, `fit_estimands.csv.gz`, `bace_diagnostics.csv`,
  `failure_records.csv`) and the aggregated tables in `data/agg/`.
- `data/discrete/`: the discrete-trait results go here, from the same datasets (same cells and seeds) as the
  continuous results, so the two parts pair dataset by dataset. `study.qmd` picks them up when present.
- `build_data.R`: rebuilds `data/` from the result pool (`~/pigauto_rubin_pool`) after `script/rubin_pool.sh`
  and `script/rubin_campaign_aggregate.R`.
- `script/rubin_discrete_aggregate.R` writes the tables in `data/discrete/`; `build_data.R [POOL] [DISC_POOL]` adds their compute record.

The raw per-fit `rds` files (556 MB) are not in git. They stay in `~/pigauto_rubin_pool` on the Mac and in the
results folders on nibi, fir, rorqual and Totoro. Never delete them.
