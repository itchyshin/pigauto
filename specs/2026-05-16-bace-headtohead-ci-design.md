# BACE head-to-head CI workflow — design spec

**Status:** approved (brainstorming 2026-05-16).
**Branch:** `feature/bace-headtohead-ci` (to be created).
**Plan:** `plans/2026-05-16-bace-headtohead-ci.md` (next).
**Mirrors:** `daniel1noble/BACE` `.github/workflows/benchmark.yml`
(Actions run [25329857467](https://github.com/daniel1noble/BACE/actions/runs/25329857467)).

---

## 1. Context

Dan ran BACE through a 6-dataset cross-package benchmark on GitHub
Actions ([Actions run
25329857467](https://github.com/daniel1noble/BACE/actions/runs/25329857467)),
using AVONET / PanTHERIA / AmphiBIO / BIEN / GlobTherm / LepTraits. The
workflow file is `daniel1noble/BACE` →
`.github/workflows/benchmark.yml` — a 6-job matrix on `ubuntu-latest`
with a `workflow_dispatch` trigger, `timeout-minutes: 345`, a heartbeat
shell trick to keep the runner watchdog awake during long parallel-R
phases, per-dataset artifacts, and an `aggregate` report job that posts
to the workflow summary.

Pigauto already has equivalent per-dataset bench scripts —
`script/bench_amphibio.R`, `bench_pantheria_full.R`, `bench_bien.R`,
`bench_globtherm_covariates.R`, `bench_leptraits_covariates.R`,
`bench_avonet_full_local.R` — but they have only ever been run on the
author's Mac. There has been no equivalent CI run for pigauto, and no
mechanical apples-to-apples head-to-head against Dan's BACE numbers in
this repo. The existing one-off head-to-head scripts
(`bench_bace_avonet_head_to_head.R`,
`bench_pantheria_bace_head_to_head.R`) cover only 2 of the 6 datasets
and re-fit BACE locally inside pigauto, which is slow and not what we
want at the CI bar.

This spec defines a sibling workflow that runs pigauto's 6 dataset
benches on GitHub Actions with the same matrix-fan-out + aggregate
pattern Dan used, snapshots BACE's published per-dataset summary into
pigauto's repo for the comparison, and produces a head-to-head table
per CI run. Manual trigger only (no auto-runs on push). Hardware
budget: 4 CPU / 16 GB / 14 GB SSD per public-repo runner; jobs cap at
5h45m wall.

## 2. Goals

1. Reproduce a pigauto vs BACE head-to-head — same 6 datasets, same
   30% MCAR held-out, same seed (2026) — without touching either
   author's laptop. CI-as-reproduction.
2. Make every pigauto bench in this set use the same documented
   default config (subset_n=2000, N_IMP=10, missing_frac=0.30,
   seed=2026, pool_method="median", clamp_outliers=FALSE,
   phylo_signal_gate=TRUE) so the head-to-head is interpretable.
3. Persist per-run artifacts and a short summary into the repo at
   `useful/ci_runs/<date>_<run-id>/` so the audit trail survives the
   90-day GHA retention window.
4. Snapshot BACE's published per-dataset results into
   `useful/bace_results_snapshot/` so the head-to-head report is
   self-contained — no live cross-repo fetch at CI time, no
   re-running BACE from inside pigauto's workflow.

## 3. Non-goals

This spec does **not** cover:

- Auto-runs on push or pull-request. Manual `workflow_dispatch` only.
  The bench is too long (~hours per dataset) and the math doesn't move
  on every commit.
- Cross-OS testing. Linux-only on standard hosted runners. macOS /
  Windows runners cost 10× and 2× more and add nothing for a
  numerical benchmark.
- Re-running BACE inside this workflow. BACE results come from a
  pinned snapshot, not a live fit. (Re-snapshotting is a manual,
  separate step — see §4.6.)
- Changing pigauto's bench config defaults or its existing
  per-dataset bench scripts beyond minimal CLI plumbing (see §4.3).
- Replacing the two existing local head-to-head scripts
  (`bench_bace_avonet_head_to_head.R`,
  `bench_pantheria_bace_head_to_head.R`). Those continue to exist for
  laptop runs that DO re-fit BACE.
- Modifying `R/`, `DESCRIPTION`, `NAMESPACE`, `tests/`, or `NEWS.md`.
  The CI workflow is a benchmarking harness, not a package change.
- Pushing CI results back to `main` from inside the workflow. The
  per-run summary commits go on a branch / PR, never directly on
  `main`.

## 4. Architecture

### 4.1 Workflow shape

Mirrors Dan's `.github/workflows/benchmark.yml` line-for-line on
structure, with pigauto-specific script paths and a few small
additions for the snapshot-based head-to-head.

```yaml
name: pigauto vs BACE cross-dataset benchmark
on:
  workflow_dispatch:

jobs:
  benchmark:
    name: ${{ matrix.dataset }}
    runs-on: ubuntu-latest
    timeout-minutes: 345    # 5h 45m, ~15 min margin below GHA 6h ceiling
    strategy:
      fail-fast: false
      matrix:
        include:
          - dataset: avonet
            script: script/gha/run_bench_avonet.R
          - dataset: pantheria
            script: script/gha/run_bench_pantheria.R
          - dataset: amphibio
            script: script/gha/run_bench_amphibio.R
          - dataset: bien
            script: script/gha/run_bench_bien.R
          - dataset: globtherm
            script: script/gha/run_bench_globtherm.R
          - dataset: leptraits
            script: script/gha/run_bench_leptraits.R

    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}
      R_KEEP_PKG_SOURCE: yes

    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: 'release'
          use-public-rspm: true
      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          extra-packages: any::devtools
          needs: ''
      - name: Run ${{ matrix.dataset }} benchmark
        run: |
          (while sleep 60; do echo "[heartbeat $(date -u +%H:%M:%S)] still running ${{ matrix.dataset }}"; done) &
          HEARTBEAT_PID=$!
          trap 'kill $HEARTBEAT_PID 2>/dev/null || true' EXIT
          Rscript ${{ matrix.script }}
      - name: Upload per-dataset artifact
        if: always()
        uses: actions/upload-artifact@v4
        with:
          name: bench-${{ matrix.dataset }}
          path: script/gha/results/${{ matrix.dataset }}/
          if-no-files-found: warn

  aggregate:
    name: Aggregate + head-to-head report
    needs: benchmark
    if: always()
    runs-on: ubuntu-latest
    timeout-minutes: 30
    permissions:
      contents: write
      pull-requests: write
    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: 'release'
          use-public-rspm: true
      - name: Install report-time R packages
        run: |
          install.packages(c("knitr", "rmarkdown"),
            repos = "https://packagemanager.posit.co/cran/__linux__/jammy/latest")
        shell: Rscript {0}
      - uses: r-lib/actions/setup-pandoc@v2
      - name: Download all per-dataset artifacts
        uses: actions/download-artifact@v4
        with:
          path: script/gha/results/_artifacts
          pattern: bench-*
      - name: Generate cross-dataset head-to-head report
        run: Rscript script/gha/make_headtohead_report.R
      - name: Publish report inline to workflow summary
        run: |
          REPORT=$(find script/gha/results/cross_dataset \
                        -name "report.md" -type f | sort -r | head -1)
          if [ -n "$REPORT" ] && [ -f "$REPORT" ]; then
            cat "$REPORT" >> "$GITHUB_STEP_SUMMARY"
          else
            echo "## Report not generated" >> "$GITHUB_STEP_SUMMARY"
            echo "Check the matrix job logs for failures." >> "$GITHUB_STEP_SUMMARY"
          fi
      - name: Upload cross-dataset summary
        if: always()
        uses: actions/upload-artifact@v4
        with:
          name: cross-dataset-summary
          path: script/gha/results/cross_dataset/
          if-no-files-found: warn
      - name: Stage per-run summary into useful/ci_runs/
        if: always()
        run: Rscript script/gha/stage_ci_run.R "${{ github.run_id }}"
      - name: Open draft PR with per-run summary
        if: always()
        uses: peter-evans/create-pull-request@v6
        with:
          branch: ci-run/${{ github.run_id }}
          commit-message: "ci: archive CI run ${{ github.run_id }}"
          title: "CI run ${{ github.run_id }} (${{ github.event.head_commit.message != '' && github.sha || 'manual' }})"
          body: |
            Auto-staged per-run summary from pigauto-vs-BACE head-to-head CI.
            See `useful/ci_runs/` for the new directory.
            Workflow run: ${{ github.server_url }}/${{ github.repository }}/actions/runs/${{ github.run_id }}
          draft: true
          add-paths: |
            useful/ci_runs/**
```

### 4.2 Why the matrix at all

A single sequential job would exceed the 6h ceiling on the largest
datasets (PanTHERIA + LepTraits each push 30-60 min wall on the
author's Mac; the public-repo Ubuntu runner is closer in spec but
slower than a 16-core Mac because `mc.cores = 4` not 16). Six parallel
matrix entries fit comfortably under both the per-job 6h ceiling and
the Free-plan 20-concurrent-job limit. `fail-fast: false` so one
dataset's failure doesn't kill the others.

### 4.3 Pigauto's documented-default config

Each `script/gha/run_bench_<dataset>.R` is a thin wrapper that loads
`devtools::load_all()` and calls into the existing
`script/bench_<dataset>.R` (or its `_covariates` / `_full` variant)
with one shared config object:

```r
PIGAUTO_CI_CONFIG <- list(
  subset_n         = 2000L,   # cap, not target: datasets smaller than 2000 use all available species
  n_imputations    = 10L,
  missing_frac     = 0.30,
  seed             = 2026L,
  pool_method      = "median",
  clamp_outliers   = FALSE,
  phylo_signal_gate = TRUE,
  mc_cores         = 4L
)
```

`subset_n` is a cap: a dataset with fewer than `subset_n` species
runs on all of them (AVONET300 = 300 species → uses all 300; BIEN
~10k species → subset to 2000).

This is the locked-in CI config — different from each laptop bench's
historical defaults (some used `n_imputations = 20`, some used
`mc.cores = 16`). The CI config is documented in
`script/gha/README.md` and is the source of truth when interpreting
CI artifacts. The existing `script/bench_*.R` scripts are NOT edited;
the GHA wrappers re-implement the bench loop on top of pigauto's
public API (`impute()`, `multi_impute()`, `evaluate_imputation()`)
using the shared CI config.

Per-dataset output: `script/gha/results/<dataset>/results.rds` +
`results.md` + `timings.json`. The `.rds` carries the full per-trait
per-method per-imputation tidy data.frame; the `.md` is a short
human-readable summary; `timings.json` records wall time per stage
(data load / fit / predict / evaluate) for cost auditing.

### 4.4 BACE-snapshot-based comparison

BACE's results are NOT recomputed inside this workflow. They come
from a curated snapshot at:

```
useful/bace_results_snapshot/
├── README.md           # how / when / what version snapshotted
├── avonet.rds          # tidy per-trait per-imputation results
├── pantheria.rds
├── amphibio.rds
├── bien.rds
├── globtherm.rds
└── leptraits.rds
```

Each `.rds` is a `tibble` with columns
`dataset, trait, type, method, imputation_idx, rmse, mae, pearson_r,
accuracy, brier, time_sec` (NA where not applicable to the trait
type), normalised to pigauto's `evaluate_imputation()` output
schema. For the snapshot files, `method = "BACE_snapshot"`; the
fresh pigauto outputs use `method = "pigauto_ci"`. The aggregate
job's
`make_headtohead_report.R` reads these alongside pigauto's six
fresh `.rds` outputs and produces:

- A summary table per dataset: pigauto vs BACE RMSE (continuous /
  count / ordinal / proportion / zi_count magnitude) and accuracy /
  Brier (binary / categorical), median + IQR across imputations and
  traits.
- A wins/ties/losses count per dataset.
- A cross-dataset bar chart (`png`) per metric.
- A single combined `report.md`.

### 4.5 Per-run summary persistence

After the aggregate job finishes, a final lightweight step commits
the per-run summary onto a new branch:

```
useful/ci_runs/2026-05-16_<run-id>/
├── report.md
├── pigauto_per_dataset/
│   ├── avonet.md
│   ├── pantheria.md
│   ├── amphibio.md
│   ├── bien.md
│   ├── globtherm.md
│   └── leptraits.md
└── timings.json
```

The full per-imputation `.rds` payloads stay in GHA artifacts (90-day
retention is fine for those). What lands in the repo is small:
`.md` summaries + `timings.json` only.

The commit is made on a `ci-run/<date>-<run-id>` branch, not on
`main`. The workflow opens a draft PR titled "CI run YYYY-MM-DD
<run-id>" for the maintainer to review and decide whether to merge,
relabel, or discard. This keeps `main` clean and gives the
maintainer veto over what becomes part of the audit trail.

### 4.6 Re-snapshotting BACE (manual, out-of-band)

When BACE has a new release or Dan reruns the benchmark, the
maintainer manually:

1. Downloads Dan's latest `bench-<dataset>` artifacts from BACE's
   GHA run via `gh run download <run-id> -R daniel1noble/BACE`.
2. Runs `script/gha/snapshot_bace.R <path-to-downloaded-artifacts>`
   locally on the Mac (this script is part of the deliverable — see
   §5). The script reads BACE's per-dataset `.rds` payloads (which
   carry BACE's native long-format results) and projects them onto
   pigauto's evaluate_imputation output schema — adding `method =
   "BACE_snapshot"`, dropping BACE-specific columns, NA-filling
   metrics that BACE doesn't compute. The projection mapping is
   defined in `script/gha/snapshot_bace.R` and is the single point
   where BACE's column names get translated.
3. Updates `useful/bace_results_snapshot/README.md` with the new
   provenance (BACE version + BACE run-id + snapshot date).
4. Commits to `main` via a normal PR.

This is NOT a CI step. Re-snapshotting is rare and deliberate.

## 5. Files to create

| Path | Role | LoC ballpark |
|---|---|---|
| `.github/workflows/bace-headtohead.yml` | The workflow itself; mirrors Dan's structure | ~110 |
| `script/gha/README.md` | CI config, file map, how to re-run, how to re-snapshot BACE | ~80 |
| `script/gha/_ci_config.R` | The shared `PIGAUTO_CI_CONFIG` list | ~25 |
| `script/gha/run_bench_avonet.R` | Wrapper around AVONET workflow | ~60 |
| `script/gha/run_bench_pantheria.R` | Wrapper around PanTHERIA workflow | ~60 |
| `script/gha/run_bench_amphibio.R` | Wrapper around AmphiBIO workflow | ~60 |
| `script/gha/run_bench_bien.R` | Wrapper around BIEN workflow | ~60 |
| `script/gha/run_bench_globtherm.R` | Wrapper around GlobTherm workflow | ~60 |
| `script/gha/run_bench_leptraits.R` | Wrapper around LepTraits workflow | ~60 |
| `script/gha/make_headtohead_report.R` | Aggregate job's report generator | ~200 |
| `script/gha/stage_ci_run.R` | Stage per-run summary into `useful/ci_runs/<date>_<run-id>/` for PR | ~60 |
| `script/gha/snapshot_bace.R` | One-off Mac script for re-snapshotting | ~80 |
| `script/gha/results/.gitkeep` | Empty marker (rest of dir is gitignored) | 0 |
| `useful/bace_results_snapshot/README.md` | Provenance for the snapshot | ~40 |
| `useful/bace_results_snapshot/<dataset>.rds` × 6 | One per dataset (binary, ignored by Rbuild) | n/a |
| `useful/ci_runs/.gitkeep` | Empty marker | 0 |

### 5.1 `.gitignore` additions

```
script/gha/results/*
!script/gha/results/.gitkeep
useful/ci_runs/*/_artifacts/
```

`useful/bace_results_snapshot/*.rds` are NOT gitignored — they're
the pinned snapshot, committed to the repo deliberately. Size budget:
~50 KB each.

### 5.2 `.Rbuildignore` additions

```
^script/gha$
^useful/bace_results_snapshot$
^useful/ci_runs$
^\\.github/workflows/bace-headtohead\\.yml$
```

Nothing in this spec ships in the installed package.

## 6. Verification

### 6.1 Local pre-flight (on the Mac, before pushing)

```r
# 1. Sanity-run the smallest wrapper end-to-end locally
Rscript script/gha/run_bench_globtherm.R     # ~10 min on Mac

# 2. Aggregate locally with one dataset's worth of fresh output +
#    six BACE snapshots — confirms make_headtohead_report.R runs.
Rscript script/gha/make_headtohead_report.R

# 3. Confirm package still loads and tests pass.
devtools::load_all()
devtools::test()
devtools::check()
```

### 6.2 GHA acceptance

After the first workflow_dispatch run completes:

- 6 `bench-<dataset>` artifacts present, each containing
  `results.rds`, `results.md`, `timings.json`.
- 1 `cross-dataset-summary` artifact containing `report.md` +
  per-trait-type bar charts.
- Workflow summary tab in the GHA UI shows the full `report.md`
  inline.
- A draft PR exists on a `ci-run/<date>-<run-id>` branch carrying
  the `useful/ci_runs/<date>_<run-id>/` summary.
- All 6 jobs are green OR the failure is a known data-availability
  issue (BIEN download flake, etc.) and the per-job log makes that
  clear.

### 6.3 Documented head-to-head result format

`report.md` MUST include for each dataset:

- N species used (after subsetting to 2000 if applicable)
- Trait types covered
- Per-trait pigauto vs BACE median RMSE / accuracy / Brier with
  IQR across `n_imputations` imputations
- Wins / ties / losses count
- Wall time per method

Format must be machine-readable (markdown tables with consistent
column headers) so future runs can be diffed.

## 7. Files NOT changed

Out of scope for this spec; will fail PR review if touched:

- `DESCRIPTION`, `NAMESPACE`, `R/*`, `tests/*`, `NEWS.md`
- `script/bench_*.R` (existing per-dataset Mac drivers stay as-is)
- `script/bench_*_bace_head_to_head.R` (existing local
  head-to-heads continue to exist for laptop runs)
- `_pkgdown.yml` (CI output is repo audit trail, not a website
  navbar entry — at least not in this spec)
- `BACE/` (the vendored sister package — never touched by pigauto
  work)

---

## Future-work pointers (not in this spec)

- Adding a 7th matrix entry once `phyloTraitData` ships and the
  Delhey-5809 bench joins the cross-dataset set.
- Replacing `useful/bace_results_snapshot/` with a live BACE-as-a-
  submodule pin if Dan ever stabilises a `bace` package on CRAN /
  GitHub-installable.
- Promoting the workflow_dispatch trigger to a monthly
  `schedule:` cron once the matrix has been validated for 3+
  successive manual runs without regression.
- Wiring the per-run summary into the pkgdown site under a new
  `Methodology > Cross-package CI` section, replacing the current
  static `bench_bace_avonet_head_to_head.html` entry. Out of scope
  here because pkgdown rebuild + navbar surgery is its own piece
  of work.
