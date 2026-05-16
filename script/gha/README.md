# pigauto CI bench harness (`script/gha/`)

Scripts driving the `.github/workflows/bace-headtohead.yml` workflow.
Never edit these in isolation — they form a fixed set keyed by the
6-dataset matrix.

## Files

| file | role |
|---|---|
| `_ci_config.R` | shared `PIGAUTO_CI_CONFIG` + helpers (`.normalize_eval`, `.eval_per_imputation`, `.run_bench`, `.build_tax_tree`, `.download_to_cache`, `.write_skip_marker`) |
| `run_bench_<dataset>.R` × 6 | matrix-job entry points (avonet / pantheria / amphibio / bien / globtherm / leptraits) |
| `make_headtohead_report.R` | aggregate job's report generator |
| `stage_ci_run.R` | stages per-run summary into `useful/ci_runs/<date>_<run-id>/` |
| `snapshot_bace.R` | one-off Mac script for re-snapshotting BACE |
| `results/` | per-dataset output dir (gitignored, kept by `.gitkeep`) |

## Config

The locked-in CI config (in `_ci_config.R`):

```r
PIGAUTO_CI_CONFIG <- list(
  subset_n          = 2000L,     # cap; smaller datasets use all available
  n_imputations     = 10L,
  missing_frac      = 0.30,
  seed              = 2026L,
  pool_method       = "median",
  clamp_outliers    = FALSE,
  phylo_signal_gate = TRUE,
  mc_cores          = 4L
)
```

This differs from each laptop bench's historical defaults; the CI
config is the single source of truth when interpreting CI artifacts.

## How to run

Manual only:

1. Go to <https://github.com/itchyshin/pigauto/actions>.
2. Pick **"pigauto vs BACE cross-dataset benchmark"**.
3. **Run workflow** → choose `main` (or the feature branch).
4. Wait ~3-6 hours for the 6-job matrix to finish.
5. Review the draft PR titled **"CI run <run-id> — pigauto vs BACE"**.

## How to re-snapshot BACE

When Dan reruns BACE's `benchmark.yml`:

```bash
gh run download <bace-run-id> -R daniel1noble/BACE -D /tmp/bace_art
Rscript script/gha/snapshot_bace.R /tmp/bace_art
# Update useful/bace_results_snapshot/README.md provenance
git add useful/bace_results_snapshot/
git commit -m "snapshot: BACE run <bace-run-id>"
```

`snapshot_bace.R` reads BACE's `summary_metrics.csv` + `run_info.csv`
from each `bench-<dataset>/run_*/<date>/` subtree, drops the
`mean_baseline` / `column_mean` rows BACE includes for reference,
projects `(correlation, mae_fit, accuracy, brier)` to pigauto's
canonical schema, and writes 6 `.rds` files under
`useful/bace_results_snapshot/`.

## Known partial-overlap notes

The CI wrappers target a subset of traits per dataset; BACE's
snapshot covers a different (sometimes broader) subset. The
aggregator merges on `(dataset, trait, type)` and shows NA where one
side is missing. Aligning trait names is a follow-up polish task.

The BIEN wrapper expects `BIEN` + `V.PhyloMaker2` packages installed;
on first runs without those it writes a skip-marker and the
aggregator handles the missing dataset gracefully.

## Spec

`specs/2026-05-16-bace-headtohead-ci-design.md` (committed
`6311a6e`).
