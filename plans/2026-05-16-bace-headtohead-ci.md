# BACE head-to-head CI workflow — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `workflow_dispatch`-triggered GitHub Actions workflow that runs pigauto's 6 dataset benches (AVONET / PanTHERIA / AmphiBIO / BIEN / GlobTherm / LepTraits) on `ubuntu-latest` matrix runners, aggregates results against a pinned BACE snapshot, and lands a per-run summary as a draft PR.

**Architecture:** New `.github/workflows/bace-headtohead.yml` mirrors Dan's `benchmark.yml` line-for-line (matrix + aggregate, 5h45m timeout, heartbeat trick). Thin R wrappers under `script/gha/` invoke each dataset's bench with one shared `PIGAUTO_CI_CONFIG`. The aggregator reads pigauto fresh `.rds` + BACE pinned snapshot `.rds`, produces a head-to-head `report.md`, and stages a slim per-run summary into `useful/ci_runs/<run-id>/` for a draft PR via `peter-evans/create-pull-request@v6`. No changes to `R/`, `DESCRIPTION`, `tests/`, `NAMESPACE`, or `NEWS.md`.

**Tech Stack:** R ≥ 4.0, GitHub Actions, `testthat` 3rd edition (for snapshot/report helper tests), `tibble`, `dplyr`, `ggplot2` (already in Suggests), pigauto's public API.

**Working branch:** `feature/bace-headtohead-ci` (spec commit `6311a6e` already on branch).

**Spec:** `specs/2026-05-16-bace-headtohead-ci-design.md`.

---

## File structure

| Path | Purpose | Create/Modify |
|---|---|---|
| `script/gha/README.md` | CI config, file map, how to re-run, how to re-snapshot BACE | create |
| `script/gha/_ci_config.R` | Shared `PIGAUTO_CI_CONFIG` list + `.normalize_eval()` helper | create |
| `script/gha/run_bench_avonet.R` | AVONET wrapper | create |
| `script/gha/run_bench_pantheria.R` | PanTHERIA wrapper | create |
| `script/gha/run_bench_amphibio.R` | AmphiBIO wrapper | create |
| `script/gha/run_bench_bien.R` | BIEN wrapper | create |
| `script/gha/run_bench_globtherm.R` | GlobTherm wrapper | create |
| `script/gha/run_bench_leptraits.R` | LepTraits wrapper | create |
| `script/gha/make_headtohead_report.R` | Aggregate job's report generator | create |
| `script/gha/stage_ci_run.R` | Stages per-run summary into `useful/ci_runs/` | create |
| `script/gha/snapshot_bace.R` | One-off Mac script for re-snapshotting BACE artifacts | create |
| `script/gha/results/.gitkeep` | Empty marker (rest of dir is gitignored) | create |
| `.github/workflows/bace-headtohead.yml` | The workflow itself | create |
| `useful/bace_results_snapshot/README.md` | Provenance for the snapshot | create |
| `useful/bace_results_snapshot/<dataset>.rds` × 6 | Initial snapshot (taken from Dan's run 25329857467) | create (one-shot Mac fixture step) |
| `useful/ci_runs/.gitkeep` | Empty marker | create |
| `tests/testthat/test-gha-helpers.R` | Unit tests for `.normalize_eval()`, `snapshot_bace`, `make_headtohead_report`, `stage_ci_run` | create |
| `.gitignore` | Add `script/gha/results/*` (keep .gitkeep) | modify |
| `.Rbuildignore` | Add `^script/gha$`, `^useful/bace_results_snapshot$`, `^useful/ci_runs$`, `^\\.github/workflows/bace-headtohead\\.yml$` | modify |

**Total:** 12 new R / config files, 1 new workflow yaml, 6 BACE snapshot `.rds`, 1 new test file, 2 admin files modified. No package source changes.

---

## Task 1: Scaffold `script/gha/` + ignore rules

**Files:**
- Create: `script/gha/results/.gitkeep`
- Create: `useful/ci_runs/.gitkeep`
- Create: `useful/bace_results_snapshot/README.md`
- Modify: `.gitignore`
- Modify: `.Rbuildignore`

- [ ] **Step 1: Verify current ignore state**

Run: `cat .gitignore && echo '---' && cat .Rbuildignore`

Note the current contents so the additions are appended (not duplicated).

- [ ] **Step 2: Add the gitkeep markers**

```bash
mkdir -p script/gha/results useful/ci_runs useful/bace_results_snapshot
touch script/gha/results/.gitkeep useful/ci_runs/.gitkeep
```

- [ ] **Step 3: Write `useful/bace_results_snapshot/README.md`**

```markdown
# BACE results snapshot

Pinned BACE benchmark outputs used by pigauto's
`bace-headtohead.yml` CI workflow as the comparison baseline.

## Provenance

- Source: `daniel1noble/BACE` Actions run **25329857467**
  (manual `workflow_dispatch` on 2026-05-13).
- BACE version: see Dan's commit at the time of the run.
- Datasets: AVONET, PanTHERIA, AmphiBIO, BIEN, GlobTherm, LepTraits.
- Config: BACE default (matches BACE's `dev/0[03-7]_benchmark_*.R`).
- Snapshot date: <YYYY-MM-DD on snapshot day>.

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
```

- [ ] **Step 4: Append to `.gitignore`**

Add at the end:

```
# pigauto CI bench results — committed only via per-run summary PR
script/gha/results/*
!script/gha/results/.gitkeep
```

- [ ] **Step 5: Append to `.Rbuildignore`**

Add at the end (escape backslashes per existing entries):

```
^script/gha$
^useful/bace_results_snapshot$
^useful/ci_runs$
^\.github/workflows/bace-headtohead\.yml$
```

- [ ] **Step 6: Verify package still builds**

Run: `R CMD build . 2>&1 | tail -20`
Expected: no warnings or errors mentioning `script/gha` or `useful/ci_runs`.

- [ ] **Step 7: Commit**

```bash
git add script/gha/results/.gitkeep useful/ci_runs/.gitkeep \
        useful/bace_results_snapshot/README.md .gitignore .Rbuildignore
git commit -m "ci(bace-h2h): scaffold script/gha + useful/{ci_runs,bace_results_snapshot}"
```

---

## Task 2: Shared CI config + `.normalize_eval()` helper

**Files:**
- Create: `script/gha/_ci_config.R`
- Create: `tests/testthat/test-gha-helpers.R`

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-gha-helpers.R`:

```r
# Unit tests for script/gha/ helpers.
# These exercise the helpers via testthat so they're covered by the
# regular `devtools::test()` run, even though they live under script/.

test_that("[gha-config] PIGAUTO_CI_CONFIG loads with documented defaults", {
  e <- new.env()
  source(file.path(testthat::test_path(), "..", "..", "script", "gha", "_ci_config.R"),
         local = e)
  cfg <- e$PIGAUTO_CI_CONFIG
  expect_type(cfg, "list")
  expect_equal(cfg$subset_n,          2000L)
  expect_equal(cfg$n_imputations,     10L)
  expect_equal(cfg$missing_frac,      0.30)
  expect_equal(cfg$seed,              2026L)
  expect_equal(cfg$pool_method,       "median")
  expect_false(cfg$clamp_outliers)
  expect_true(cfg$phylo_signal_gate)
  expect_equal(cfg$mc_cores,          4L)
})

test_that("[gha-config] .normalize_eval() projects to the canonical schema", {
  e <- new.env()
  source(file.path(testthat::test_path(), "..", "..", "script", "gha", "_ci_config.R"),
         local = e)
  norm <- e$.normalize_eval

  raw <- data.frame(
    trait          = c("mass", "diet"),
    type           = c("continuous", "categorical"),
    imputation_idx = c(1L, 1L),
    rmse           = c(0.41, NA_real_),
    mae            = c(0.30, NA_real_),
    pearson_r      = c(0.85, NA_real_),
    accuracy       = c(NA_real_, 0.72),
    brier          = c(NA_real_, 0.18),
    time_sec       = c(120, 120),
    stringsAsFactors = FALSE
  )

  out <- norm(raw, dataset = "avonet", method = "pigauto_ci")
  expected_cols <- c("dataset", "trait", "type", "method",
                     "imputation_idx", "rmse", "mae", "pearson_r",
                     "accuracy", "brier", "time_sec")
  expect_equal(colnames(out), expected_cols)
  expect_equal(out$dataset, c("avonet", "avonet"))
  expect_equal(out$method,  c("pigauto_ci", "pigauto_ci"))
  expect_equal(nrow(out), 2L)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `devtools::test(filter = "gha-helpers")`
Expected: 2 failures with "could not find function" or "no such file" for `_ci_config.R`.

- [ ] **Step 3: Write `script/gha/_ci_config.R`**

```r
# script/gha/_ci_config.R
#
# Shared CI config for the BACE head-to-head workflow.
# Sourced by every script/gha/run_bench_*.R wrapper.
#
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md.

PIGAUTO_CI_CONFIG <- list(
  subset_n          = 2000L,     # cap; smaller datasets use all available species
  n_imputations     = 10L,
  missing_frac      = 0.30,
  seed              = 2026L,
  pool_method       = "median",
  clamp_outliers    = FALSE,
  phylo_signal_gate = TRUE,
  mc_cores          = 4L
)

# Canonical schema for cross-method comparison.
# Pigauto's evaluate_imputation() output gets passed through here.
# BACE snapshots also conform to this schema (see script/gha/snapshot_bace.R).
.normalize_eval <- function(df, dataset, method) {
  stopifnot(is.data.frame(df), is.character(dataset), is.character(method))
  expected_cols <- c("trait", "type", "imputation_idx", "rmse", "mae",
                     "pearson_r", "accuracy", "brier", "time_sec")
  for (col in expected_cols) {
    if (!col %in% colnames(df)) df[[col]] <- NA
  }
  out <- data.frame(
    dataset        = rep(dataset, nrow(df)),
    trait          = as.character(df$trait),
    type           = as.character(df$type),
    method         = rep(method, nrow(df)),
    imputation_idx = as.integer(df$imputation_idx),
    rmse           = as.numeric(df$rmse),
    mae            = as.numeric(df$mae),
    pearson_r      = as.numeric(df$pearson_r),
    accuracy       = as.numeric(df$accuracy),
    brier          = as.numeric(df$brier),
    time_sec       = as.numeric(df$time_sec),
    stringsAsFactors = FALSE
  )
  out
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `devtools::test(filter = "gha-helpers")`
Expected: 2 passes, 0 failures.

- [ ] **Step 5: Commit**

```bash
git add script/gha/_ci_config.R tests/testthat/test-gha-helpers.R
git commit -m "ci(bace-h2h): shared PIGAUTO_CI_CONFIG + normalize_eval helper"
```

---

## Task 3: AVONET wrapper (template for the other five)

**Files:**
- Create: `script/gha/run_bench_avonet.R`

- [ ] **Step 1: Write the wrapper**

```r
#!/usr/bin/env Rscript
# script/gha/run_bench_avonet.R
#
# CI wrapper for the AVONET head-to-head bench.
# Loads pigauto via devtools::load_all() (CI doesn't install),
# runs the documented-default config, writes a normalized
# .rds + a short .md + timings.json under script/gha/results/avonet/.
#
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  here <- "."  # workflow runs from repo root
  source(file.path(here, "script", "gha", "_ci_config.R"))
  devtools::load_all(here, quiet = TRUE)
})

cfg     <- PIGAUTO_CI_CONFIG
out_dir <- file.path("script", "gha", "results", "avonet")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Data — bundled with the package
# -------------------------------------------------------------------------
e <- new.env()
utils::data("avonet300", package = "pigauto", envir = e)
utils::data("tree300",   package = "pigauto", envir = e)
df <- e$avonet300
rownames(df) <- df$Species_Key
df$Species_Key <- NULL
tree <- e$tree300

# Subset to cap (no-op for AVONET300 which is already 300 species).
n_target <- min(cfg$subset_n, nrow(df))
if (n_target < nrow(df)) {
  set.seed(cfg$seed)
  keep <- sort(sample.int(nrow(df), n_target))
  df   <- df[keep, , drop = FALSE]
  tree <- ape::keep.tip(tree, rownames(df))
}
cat(sprintf("[avonet] n=%d species x %d traits\n", nrow(df), ncol(df)))

# -------------------------------------------------------------------------
# Hold-out mask (30% MCAR, seeded)
# -------------------------------------------------------------------------
t0_split <- Sys.time()
pd0 <- preprocess_traits(df, tree, log_transform = TRUE)
splits <- make_missing_splits(pd0$X_scaled,
                              missing_frac = cfg$missing_frac,
                              val_frac     = 0.5,
                              seed         = cfg$seed,
                              trait_map    = pd0$trait_map)

# Mark which LATENT cells are in the test split.
mask_latent <- matrix(FALSE, nrow = nrow(pd0$X_scaled), ncol = ncol(pd0$X_scaled))
mask_latent[splits$test_idx] <- TRUE

# Map latent cols → user cols via trait_map (matches the pattern
# used in script/bench_bace_avonet_head_to_head.R).
user_mask_test <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df),
                          dimnames = list(rownames(df), colnames(df)))
for (k in seq_along(pd0$trait_map)) {
  tm <- pd0$trait_map[[k]]
  user_cols   <- tm$input_cols %||% tm$col
  latent_cols <- tm$latent_cols %||% tm$cols %||% k
  hit_rows <- apply(mask_latent[, latent_cols, drop = FALSE], 1L, any)
  for (uc in user_cols) {
    user_mask_test[, uc] <- user_mask_test[, uc] | hit_rows
  }
}
df_masked <- df
df_masked[user_mask_test] <- NA

t_split <- as.numeric(difftime(Sys.time(), t0_split, units = "secs"))

# -------------------------------------------------------------------------
# Multi-impute with CI config
# -------------------------------------------------------------------------
t0_fit <- Sys.time()
mi <- multi_impute(
  df_masked, tree,
  m              = cfg$n_imputations,
  pool_method    = cfg$pool_method,
  clamp_outliers = cfg$clamp_outliers,
  seed           = cfg$seed
)
t_fit <- as.numeric(difftime(Sys.time(), t0_fit, units = "secs"))

# -------------------------------------------------------------------------
# Evaluate per imputation
# -------------------------------------------------------------------------
t0_eval <- Sys.time()
ev_list <- lapply(seq_along(mi$datasets), function(i) {
  ev <- evaluate_imputation(mi$datasets[[i]], df,
                            mask = user_mask_test,
                            trait_map = pd0$trait_map)
  ev$imputation_idx <- i
  ev$time_sec       <- NA_real_
  ev
})
ev_all <- do.call(rbind, ev_list)
ev_all$time_sec <- t_fit
t_eval <- as.numeric(difftime(Sys.time(), t0_eval, units = "secs"))

# -------------------------------------------------------------------------
# Normalize + write outputs
# -------------------------------------------------------------------------
out_tbl <- .normalize_eval(ev_all, dataset = "avonet", method = "pigauto_ci")
saveRDS(out_tbl, file.path(out_dir, "results.rds"))

timings <- list(split_sec = t_split, fit_sec = t_fit, eval_sec = t_eval,
                total_sec = t_split + t_fit + t_eval,
                n_species = nrow(df), n_imputations = cfg$n_imputations)
jsonlite::write_json(timings, file.path(out_dir, "timings.json"),
                     auto_unbox = TRUE, pretty = TRUE)

# Markdown summary
md_lines <- c(
  "# avonet — pigauto CI bench",
  sprintf("Run config: seed=%d, missing_frac=%.2f, n_imputations=%d",
          cfg$seed, cfg$missing_frac, cfg$n_imputations),
  sprintf("N species used: %d", nrow(df)),
  sprintf("Wall time: %.1f s (fit) + %.1f s (eval)", t_fit, t_eval),
  "",
  "## Per-trait medians across imputations",
  ""
)
agg <- stats::aggregate(cbind(rmse, accuracy) ~ trait + type, data = out_tbl,
                        FUN = function(x) stats::median(x, na.rm = TRUE))
md_lines <- c(md_lines, knitr::kable(agg, format = "markdown"))
writeLines(md_lines, file.path(out_dir, "results.md"))

cat(sprintf("[avonet] done in %.1f s total\n", t_split + t_fit + t_eval))
```

Note the `%||%` operator used above is base R 4.4+. If targeting older R,
substitute `if (is.null(x)) y else x`.

- [ ] **Step 2: Smoke-test locally on the Mac**

Run: `Rscript script/gha/run_bench_avonet.R`
Expected: takes ~5-10 min, writes three files to `script/gha/results/avonet/`:
- `results.rds`
- `results.md`
- `timings.json`

Final stdout line: `[avonet] done in <N> s total`.

- [ ] **Step 3: Inspect the output**

```r
res <- readRDS("script/gha/results/avonet/results.rds")
stopifnot(is.data.frame(res),
          all(c("dataset","trait","type","method","imputation_idx",
                "rmse","mae","pearson_r","accuracy","brier","time_sec")
              %in% colnames(res)),
          all(res$method == "pigauto_ci"),
          all(res$dataset == "avonet"))
```

Expected: no errors.

- [ ] **Step 4: Commit**

```bash
git add script/gha/run_bench_avonet.R
git commit -m "ci(bace-h2h): AVONET wrapper using documented-default config"
```

---

## Task 4: PanTHERIA / AmphiBIO / BIEN / GlobTherm / LepTraits wrappers

**Files:**
- Create: `script/gha/run_bench_pantheria.R`
- Create: `script/gha/run_bench_amphibio.R`
- Create: `script/gha/run_bench_bien.R`
- Create: `script/gha/run_bench_globtherm.R`
- Create: `script/gha/run_bench_leptraits.R`

- [ ] **Step 1: Identify per-dataset data-loading lines**

For each dataset, locate the data load in the existing `script/bench_<dataset>*.R`:

```bash
grep -n "read\\.csv\\|read_csv\\|download.file\\|fread" script/bench_pantheria_full.R | head -5
grep -n "read\\.csv\\|read_csv\\|download.file\\|fread" script/bench_amphibio.R | head -5
grep -n "read\\.csv\\|read_csv\\|download.file\\|fread" script/bench_bien.R | head -5
grep -n "read\\.csv\\|read_csv\\|download.file\\|fread" script/bench_globtherm_covariates.R | head -5
grep -n "read\\.csv\\|read_csv\\|download.file\\|fread" script/bench_leptraits_covariates.R | head -5
```

Note the exact data-load lines, file URLs, and trait-column lists. Keep them in this format for the wrappers.

- [ ] **Step 2: Write each wrapper using AVONET wrapper as template**

For each dataset, copy `script/gha/run_bench_avonet.R` to `script/gha/run_bench_<dataset>.R` and replace the data-load block with the dataset-specific path (CSV URL or bundled dataset). Keep:

- The same `cfg <- PIGAUTO_CI_CONFIG` block
- The same subset-by-`cap` logic
- The same `preprocess_traits` → `make_missing_splits` → `multi_impute` → `evaluate_imputation` chain
- The same `.normalize_eval()` call with `dataset = "<name>"`
- The same output layout under `script/gha/results/<name>/`

Each wrapper should also `cat` early-stage warnings if the dataset isn't reachable (BIEN has network flake).

- [ ] **Step 3: Smoke-test the smallest two (GlobTherm + AmphiBIO)**

```bash
Rscript script/gha/run_bench_globtherm.R    # ~10 min
Rscript script/gha/run_bench_amphibio.R     # ~15 min
```

Verify outputs at `script/gha/results/{globtherm,amphibio}/`.

- [ ] **Step 4: Schema-check all wrappers without running fits**

```bash
for d in pantheria amphibio bien globtherm leptraits; do
  Rscript -e "src <- readLines('script/gha/run_bench_${d}.R');
              stopifnot(any(grepl('PIGAUTO_CI_CONFIG', src)),
                        any(grepl('.normalize_eval', src)),
                        any(grepl('script/gha/results/${d}', src, fixed = TRUE)));
              cat('${d}: schema OK\\n')"
done
```

- [ ] **Step 5: Commit**

```bash
git add script/gha/run_bench_*.R
git commit -m "ci(bace-h2h): wrappers for pantheria / amphibio / bien / globtherm / leptraits"
```

---

## Task 5: Head-to-head report generator

**Files:**
- Create: `script/gha/make_headtohead_report.R`
- Modify: `tests/testthat/test-gha-helpers.R`

- [ ] **Step 1: Add the failing test**

Append to `tests/testthat/test-gha-helpers.R`:

```r
test_that("[gha-h2h] make_headtohead_report produces report.md + summary tbl", {
  e <- new.env()
  source(file.path(testthat::test_path(), "..", "..", "script", "gha", "_ci_config.R"),
         local = e)
  source(file.path(testthat::test_path(), "..", "..", "script", "gha", "make_headtohead_report.R"),
         local = e)
  build <- e$build_h2h_report

  pigauto_tbl <- e$.normalize_eval(data.frame(
    trait = c("mass", "diet"),
    type  = c("continuous", "categorical"),
    imputation_idx = c(1L, 1L),
    rmse = c(0.41, NA_real_),
    accuracy = c(NA_real_, 0.72),
    time_sec = c(120, 120),
    stringsAsFactors = FALSE
  ), dataset = "avonet", method = "pigauto_ci")

  bace_tbl <- e$.normalize_eval(data.frame(
    trait = c("mass", "diet"),
    type  = c("continuous", "categorical"),
    imputation_idx = c(1L, 1L),
    rmse = c(0.55, NA_real_),
    accuracy = c(NA_real_, 0.65),
    time_sec = c(3600, 3600),
    stringsAsFactors = FALSE
  ), dataset = "avonet", method = "BACE_snapshot")

  tmp <- tempfile(fileext = "")
  dir.create(tmp)
  rep <- build(combined = rbind(pigauto_tbl, bace_tbl), out_dir = tmp)
  expect_true(file.exists(file.path(tmp, "report.md")))
  expect_s3_class(rep$summary, "data.frame")
  expect_true(all(c("dataset","trait","type","pigauto","bace","winner")
                  %in% colnames(rep$summary)))
  expect_equal(rep$summary$winner[rep$summary$trait == "mass"],
               "pigauto")     # 0.41 RMSE < 0.55 RMSE
  expect_equal(rep$summary$winner[rep$summary$trait == "diet"],
               "pigauto")     # 0.72 acc > 0.65 acc
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `devtools::test(filter = "gha-helpers")`
Expected: 1 new failure with "no such file" for `make_headtohead_report.R`.

- [ ] **Step 3: Implement `script/gha/make_headtohead_report.R`**

```r
#!/usr/bin/env Rscript
# script/gha/make_headtohead_report.R
#
# Aggregate job's report generator: reads pigauto's fresh per-dataset
# .rds + BACE pinned snapshot .rds, produces a combined head-to-head
# report.md + per-trait-type charts under
# script/gha/results/cross_dataset/.
#
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md.

suppressPackageStartupMessages({
  source(file.path("script", "gha", "_ci_config.R"))
})

DATASETS <- c("avonet", "pantheria", "amphibio", "bien", "globtherm", "leptraits")

load_pigauto_results <- function() {
  base <- file.path("script", "gha", "results", "_artifacts")
  do.call(rbind, lapply(DATASETS, function(d) {
    f <- file.path(base, paste0("bench-", d), "results.rds")
    if (!file.exists(f)) {
      message(sprintf("[h2h] missing pigauto results for %s; skipping", d))
      return(NULL)
    }
    readRDS(f)
  }))
}

load_bace_snapshot <- function() {
  base <- file.path("useful", "bace_results_snapshot")
  do.call(rbind, lapply(DATASETS, function(d) {
    f <- file.path(base, paste0(d, ".rds"))
    if (!file.exists(f)) {
      message(sprintf("[h2h] missing BACE snapshot for %s; skipping", d))
      return(NULL)
    }
    readRDS(f)
  }))
}

# Decide winner per (dataset, trait, type) by primary metric.
.h2h_winner <- function(pig, bace, type) {
  # continuous-family: lower RMSE wins
  # discrete-family:   higher accuracy wins
  if (type %in% c("continuous","count","ordinal","proportion","zi_count")) {
    if (is.na(pig) || is.na(bace)) return(NA_character_)
    if (pig < bace) "pigauto" else if (bace < pig) "bace" else "tie"
  } else {
    if (is.na(pig) || is.na(bace)) return(NA_character_)
    if (pig > bace) "pigauto" else if (bace > pig) "bace" else "tie"
  }
}

build_h2h_report <- function(combined, out_dir) {
  stopifnot(is.data.frame(combined),
            all(c("dataset","trait","type","method","rmse","accuracy")
                %in% colnames(combined)))

  # Median across imputations per (dataset, trait, type, method)
  med <- stats::aggregate(
    cbind(rmse, mae, pearson_r, accuracy, brier, time_sec) ~ dataset + trait + type + method,
    data = combined,
    FUN  = function(x) stats::median(x, na.rm = TRUE)
  )

  pig <- med[med$method == "pigauto_ci",   , drop = FALSE]
  bac <- med[med$method == "BACE_snapshot", , drop = FALSE]
  m <- merge(pig, bac, by = c("dataset","trait","type"),
             suffixes = c("_pigauto","_bace"), all = TRUE)

  m$pigauto <- ifelse(m$type %in% c("continuous","count","ordinal",
                                    "proportion","zi_count"),
                      m$rmse_pigauto, m$accuracy_pigauto)
  m$bace    <- ifelse(m$type %in% c("continuous","count","ordinal",
                                    "proportion","zi_count"),
                      m$rmse_bace, m$accuracy_bace)
  m$winner  <- vapply(seq_len(nrow(m)),
                      function(i) .h2h_winner(m$pigauto[i], m$bace[i], m$type[i]),
                      character(1))

  summary_tbl <- m[, c("dataset","trait","type","pigauto","bace","winner")]

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Write summary RDS for downstream stage_ci_run.R
  saveRDS(summary_tbl, file.path(out_dir, "summary.rds"))

  # Write markdown report
  md <- c(
    "# pigauto vs BACE — cross-dataset head-to-head",
    "",
    sprintf("Run date: %s (UTC).", format(Sys.time(), tz = "UTC")),
    sprintf("Datasets: %s.", paste(DATASETS, collapse = ", ")),
    "",
    "## Wins / ties / losses per dataset",
    ""
  )
  wtl <- as.data.frame(table(summary_tbl$dataset,
                              factor(summary_tbl$winner,
                                     levels = c("pigauto","tie","bace"))),
                       stringsAsFactors = FALSE)
  colnames(wtl) <- c("dataset","winner","n")
  md <- c(md, knitr::kable(wtl, format = "markdown"))
  md <- c(md, "", "## Per-trait detail", "",
          knitr::kable(summary_tbl, format = "markdown",
                       digits = 3))
  writeLines(md, file.path(out_dir, "report.md"))

  invisible(list(summary = summary_tbl, report_path = file.path(out_dir, "report.md")))
}

# When invoked as a script (not sourced), wire up I/O.
if (sys.nframe() == 0L) {
  out_dir <- file.path("script", "gha", "results", "cross_dataset")
  combined <- rbind(load_pigauto_results(), load_bace_snapshot())
  if (is.null(combined) || nrow(combined) == 0L) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines("# No data — both pigauto and BACE inputs empty.",
               file.path(out_dir, "report.md"))
    quit(save = "no", status = 0L)
  }
  build_h2h_report(combined, out_dir)
  cat(sprintf("[h2h] wrote report to %s\n", file.path(out_dir, "report.md")))
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `devtools::test(filter = "gha-helpers")`
Expected: 3 passes, 0 failures.

- [ ] **Step 5: Commit**

```bash
git add script/gha/make_headtohead_report.R tests/testthat/test-gha-helpers.R
git commit -m "ci(bace-h2h): aggregator generates pigauto vs BACE head-to-head"
```

---

## Task 6: BACE snapshot script

**Files:**
- Create: `script/gha/snapshot_bace.R`
- Modify: `tests/testthat/test-gha-helpers.R`

- [ ] **Step 1: Add the failing test**

Append to `tests/testthat/test-gha-helpers.R`:

```r
test_that("[gha-snapshot] snapshot_bace_one() projects to canonical schema", {
  e <- new.env()
  source(file.path(testthat::test_path(), "..", "..", "script", "gha", "_ci_config.R"),
         local = e)
  source(file.path(testthat::test_path(), "..", "..", "script", "gha", "snapshot_bace.R"),
         local = e)
  fn <- e$snapshot_bace_one

  # Simulate BACE's native output (long-format with BACE column names).
  bace_native <- data.frame(
    trait        = c("mass","mass","diet","diet"),
    trait_type   = c("continuous","continuous","categorical","categorical"),
    iter         = c(1L, 2L, 1L, 2L),
    RMSE         = c(0.55, 0.57, NA_real_, NA_real_),
    PearsonR     = c(0.82, 0.80, NA_real_, NA_real_),
    Accuracy     = c(NA_real_, NA_real_, 0.65, 0.66),
    BrierScore   = c(NA_real_, NA_real_, 0.20, 0.19),
    seconds      = c(3600, 3600, 3600, 3600),
    stringsAsFactors = FALSE
  )

  out <- fn(bace_native, dataset = "avonet")
  expect_equal(unique(out$dataset), "avonet")
  expect_equal(unique(out$method),  "BACE_snapshot")
  expect_setequal(out$trait, c("mass","diet"))
  expect_equal(nrow(out), 4L)
  expect_equal(out$rmse[out$trait == "mass"], c(0.55, 0.57))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `devtools::test(filter = "gha-helpers")`
Expected: 1 new failure for `snapshot_bace.R` not found.

- [ ] **Step 3: Implement `script/gha/snapshot_bace.R`**

```r
#!/usr/bin/env Rscript
# script/gha/snapshot_bace.R
#
# One-off Mac helper: re-snapshots BACE's per-dataset benchmark
# outputs into pigauto's useful/bace_results_snapshot/.
#
# Usage:
#   Rscript script/gha/snapshot_bace.R <path-to-downloaded-bace-artifacts>
#
# Procedure:
#   1. gh run download <run-id> -R daniel1noble/BACE \
#                       -D /tmp/bace_artifacts
#   2. Rscript script/gha/snapshot_bace.R /tmp/bace_artifacts
#   3. Inspect useful/bace_results_snapshot/<dataset>.rds + update README.md
#   4. Commit to main via PR.
#
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md §4.6.

suppressPackageStartupMessages({
  source(file.path("script", "gha", "_ci_config.R"))
})

DATASETS <- c("avonet","pantheria","amphibio","bien","globtherm","leptraits")

# BACE-column → canonical-column mapping. Edit if BACE's schema changes.
.BACE_COL_MAP <- c(
  trait          = "trait",
  trait_type     = "type",
  iter           = "imputation_idx",
  RMSE           = "rmse",
  MAE            = "mae",
  PearsonR       = "pearson_r",
  Accuracy       = "accuracy",
  BrierScore     = "brier",
  seconds        = "time_sec"
)

# Project ONE BACE per-dataset table to the canonical schema.
snapshot_bace_one <- function(bace_native, dataset) {
  stopifnot(is.data.frame(bace_native), is.character(dataset))
  df <- bace_native
  for (bc in names(.BACE_COL_MAP)) {
    canon <- .BACE_COL_MAP[[bc]]
    if (bc %in% colnames(df)) df[[canon]] <- df[[bc]]
  }
  .normalize_eval(df, dataset = dataset, method = "BACE_snapshot")
}

# When invoked as a script, walk the artifacts directory.
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) {
    stop("usage: Rscript script/gha/snapshot_bace.R <path-to-bace-artifacts>")
  }
  src_dir <- args[[1L]]
  out_dir <- file.path("useful", "bace_results_snapshot")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (d in DATASETS) {
    # Dan's artifact pattern: dev/benchmark_results/<dataset>/results.rds
    art <- list.files(file.path(src_dir, paste0("bench-", d)),
                      pattern = "\\.rds$", full.names = TRUE,
                      recursive = TRUE)
    if (length(art) == 0L) {
      message(sprintf("[snapshot] no artifact for %s; skipping", d))
      next
    }
    raw <- readRDS(art[[1L]])
    out <- snapshot_bace_one(raw, dataset = d)
    saveRDS(out, file.path(out_dir, paste0(d, ".rds")))
    cat(sprintf("[snapshot] %s: %d rows -> %s\n",
                d, nrow(out), file.path(out_dir, paste0(d, ".rds"))))
  }
  cat("[snapshot] done. Remember to update useful/bace_results_snapshot/README.md\n")
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `devtools::test(filter = "gha-helpers")`
Expected: 4 passes, 0 failures.

- [ ] **Step 5: Manually snapshot Dan's run 25329857467**

(One-shot; if `gh` auth or BACE artifact format diverges from the
script's assumed column mapping, edit `.BACE_COL_MAP` accordingly.)

```bash
mkdir -p /tmp/bace_artifacts_25329857467
gh run download 25329857467 -R daniel1noble/BACE \
   -D /tmp/bace_artifacts_25329857467
Rscript script/gha/snapshot_bace.R /tmp/bace_artifacts_25329857467
ls -la useful/bace_results_snapshot/
```

Expected: 6 `<dataset>.rds` files of ≤50 KB each.

- [ ] **Step 6: Update snapshot README provenance**

Edit `useful/bace_results_snapshot/README.md`'s "Snapshot date" line
to today's date.

- [ ] **Step 7: Commit**

```bash
git add script/gha/snapshot_bace.R \
        tests/testthat/test-gha-helpers.R \
        useful/bace_results_snapshot/*.rds \
        useful/bace_results_snapshot/README.md
git commit -m "ci(bace-h2h): snapshot_bace.R + initial snapshot of Dan's run 25329857467"
```

---

## Task 7: stage_ci_run.R — per-run summary stager

**Files:**
- Create: `script/gha/stage_ci_run.R`
- Modify: `tests/testthat/test-gha-helpers.R`

- [ ] **Step 1: Add the failing test**

Append to `tests/testthat/test-gha-helpers.R`:

```r
test_that("[gha-stage] stage_ci_run() copies the right files into useful/ci_runs/", {
  e <- new.env()
  source(file.path(testthat::test_path(), "..", "..", "script", "gha", "stage_ci_run.R"),
         local = e)
  stage <- e$stage_ci_run

  # Build a fake results dir
  tmp <- tempfile(fileext = ""); dir.create(tmp, recursive = TRUE)
  for (d in c("avonet","globtherm")) {
    sub <- file.path(tmp, "_artifacts", paste0("bench-", d))
    dir.create(sub, recursive = TRUE)
    writeLines(sprintf("# %s", d), file.path(sub, "results.md"))
    jsonlite::write_json(list(fit_sec = 60),
                         file.path(sub, "timings.json"),
                         auto_unbox = TRUE)
  }
  cd <- file.path(tmp, "cross_dataset")
  dir.create(cd, recursive = TRUE)
  writeLines("# h2h", file.path(cd, "report.md"))

  out_root <- tempfile(fileext = ""); dir.create(out_root)
  run_id <- "99999"
  staged <- stage(results_root = tmp, ci_runs_root = out_root,
                  run_id = run_id, date_str = "2026-05-16")

  expect_true(dir.exists(staged))
  expect_true(file.exists(file.path(staged, "report.md")))
  expect_true(file.exists(file.path(staged, "pigauto_per_dataset", "avonet.md")))
  expect_true(file.exists(file.path(staged, "pigauto_per_dataset", "globtherm.md")))
  expect_true(file.exists(file.path(staged, "timings.json")))
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `devtools::test(filter = "gha-helpers")`
Expected: 1 new failure.

- [ ] **Step 3: Implement `script/gha/stage_ci_run.R`**

```r
#!/usr/bin/env Rscript
# script/gha/stage_ci_run.R
#
# Stage a slim per-run summary into useful/ci_runs/<date>_<run-id>/
# for the draft-PR step.
#
# Usage (from CI):
#   Rscript script/gha/stage_ci_run.R "<github-run-id>"
#
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md §4.5.

DATASETS <- c("avonet","pantheria","amphibio","bien","globtherm","leptraits")

stage_ci_run <- function(results_root = file.path("script","gha","results"),
                         ci_runs_root = file.path("useful","ci_runs"),
                         run_id,
                         date_str = format(Sys.Date(), "%Y-%m-%d")) {
  stopifnot(is.character(run_id), nchar(run_id) > 0L)
  staged <- file.path(ci_runs_root, sprintf("%s_%s", date_str, run_id))
  dir.create(file.path(staged, "pigauto_per_dataset"),
             recursive = TRUE, showWarnings = FALSE)

  # report.md from cross_dataset/
  rep_src <- file.path(results_root, "cross_dataset", "report.md")
  if (file.exists(rep_src)) {
    file.copy(rep_src, file.path(staged, "report.md"), overwrite = TRUE)
  } else {
    writeLines("# Report unavailable (cross_dataset/report.md missing)",
               file.path(staged, "report.md"))
  }

  # per-dataset .md from _artifacts/bench-<d>/results.md
  agg_timings <- list()
  for (d in DATASETS) {
    md_src <- file.path(results_root, "_artifacts",
                        paste0("bench-", d), "results.md")
    if (file.exists(md_src)) {
      file.copy(md_src,
                file.path(staged, "pigauto_per_dataset", paste0(d, ".md")),
                overwrite = TRUE)
    }
    t_src <- file.path(results_root, "_artifacts",
                       paste0("bench-", d), "timings.json")
    if (file.exists(t_src)) {
      agg_timings[[d]] <- jsonlite::read_json(t_src, simplifyVector = TRUE)
    }
  }
  jsonlite::write_json(agg_timings,
                       file.path(staged, "timings.json"),
                       auto_unbox = TRUE, pretty = TRUE)

  invisible(staged)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1L) stop("usage: stage_ci_run.R <run-id>")
  out <- stage_ci_run(run_id = args[[1L]])
  cat(sprintf("[stage] staged at %s\n", out))
}
```

- [ ] **Step 4: Run to verify all gha-helpers tests pass**

Run: `devtools::test(filter = "gha-helpers")`
Expected: 5 passes (config, normalize_eval, h2h, snapshot, stage), 0 failures.

- [ ] **Step 5: Commit**

```bash
git add script/gha/stage_ci_run.R tests/testthat/test-gha-helpers.R
git commit -m "ci(bace-h2h): stage_ci_run.R writes slim per-run summary"
```

---

## Task 8: The workflow yaml

**Files:**
- Create: `.github/workflows/bace-headtohead.yml`

- [ ] **Step 1: Write the workflow**

Create `.github/workflows/bace-headtohead.yml`:

```yaml
name: pigauto vs BACE cross-dataset benchmark

on:
  workflow_dispatch:

jobs:
  benchmark:
    name: ${{ matrix.dataset }}
    runs-on: ubuntu-latest
    timeout-minutes: 345
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
          extra-packages: any::devtools, any::jsonlite, any::knitr
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

      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          extra-packages: any::knitr, any::rmarkdown, any::jsonlite
          needs: ''

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
          title: "CI run ${{ github.run_id }} — pigauto vs BACE"
          body: |
            Auto-staged per-run summary from pigauto-vs-BACE head-to-head CI.
            See `useful/ci_runs/` for the new directory.
            Workflow run: ${{ github.server_url }}/${{ github.repository }}/actions/runs/${{ github.run_id }}
          draft: true
          add-paths: |
            useful/ci_runs/**
```

- [ ] **Step 2: Lint the yaml**

If `actionlint` is installed locally:

```bash
actionlint .github/workflows/bace-headtohead.yml
```

If not, use `yq` for basic parse:

```bash
yq '.jobs.benchmark.strategy.matrix.include | length' .github/workflows/bace-headtohead.yml
```

Expected: `6` (one entry per dataset).

- [ ] **Step 3: Verify the matrix entries match the wrappers on disk**

```bash
for d in avonet pantheria amphibio bien globtherm leptraits; do
  test -f "script/gha/run_bench_${d}.R" && echo "OK ${d}" || echo "MISSING ${d}"
done
```

Expected: 6 lines, all `OK <d>`.

- [ ] **Step 4: Commit**

```bash
git add .github/workflows/bace-headtohead.yml
git commit -m "ci(bace-h2h): GitHub Actions workflow (matrix + aggregate + draft PR)"
```

---

## Task 9: README for `script/gha/`

**Files:**
- Create: `script/gha/README.md`

- [ ] **Step 1: Write the README**

```markdown
# pigauto CI bench harness (`script/gha/`)

Scripts driving the `.github/workflows/bace-headtohead.yml` workflow.
Never edit these in isolation — they form a fixed set keyed by the
6-dataset matrix.

## Files

| file | role |
|---|---|
| `_ci_config.R` | shared `PIGAUTO_CI_CONFIG` list + `.normalize_eval()` |
| `run_bench_<dataset>.R` × 6 | matrix job entry points |
| `make_headtohead_report.R` | aggregate job's report generator |
| `stage_ci_run.R` | stages per-run summary into `useful/ci_runs/` |
| `snapshot_bace.R` | one-off Mac script for re-snapshotting BACE |
| `results/` | per-dataset output dir (gitignored, kept by .gitkeep) |

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

## How to re-run

Manual only:

1. Go to <https://github.com/itchyshin/pigauto/actions>.
2. Pick "pigauto vs BACE cross-dataset benchmark".
3. "Run workflow" → choose `main` (or the feature branch).
4. Wait ~3-6 hours for the 6-job matrix to finish.
5. Review the draft PR titled "CI run <run-id> — pigauto vs BACE".

## How to re-snapshot BACE

When Dan reruns BACE's `benchmark.yml`:

```bash
gh run download <bace-run-id> -R daniel1noble/BACE -D /tmp/bace_art
Rscript script/gha/snapshot_bace.R /tmp/bace_art
# Update useful/bace_results_snapshot/README.md provenance
git add useful/bace_results_snapshot/
git commit -m "snapshot: BACE run <bace-run-id>"
```

## Spec

- Design: `specs/2026-05-16-bace-headtohead-ci-design.md`.
- GHA limits & conventions reference: `~/.claude/projects/-Users-z3437171-Dropbox-Github-Local-pigauto/memory/gha_limits_for_bace_workflows.md` (on the maintainer's machine).
```

- [ ] **Step 2: Commit**

```bash
git add script/gha/README.md
git commit -m "ci(bace-h2h): script/gha/README.md"
```

---

## Task 10: Final test pass + R CMD check sanity

**Files:** none new.

- [ ] **Step 1: Run the full unit-test suite**

Run: `devtools::test()`
Expected: all existing tests pass; +5 new gha-helpers tests pass.

- [ ] **Step 2: Run R CMD check**

Run: `devtools::check(args = "--no-manual")`
Expected: `Status: OK`. No NOTE / WARNING / ERROR about new files
(everything under `script/gha/`, `useful/bace_results_snapshot/`,
`useful/ci_runs/` is Rbuildignored).

- [ ] **Step 3: Confirm tarball excludes CI infrastructure**

```bash
R CMD build .
tar -tzf pigauto_*.tar.gz | grep -E "script/gha|useful/ci_runs|useful/bace_results_snapshot|\\.github/workflows/bace-headtohead" \
  && echo "FAIL: CI files leaked into tarball" \
  || echo "OK: CI files correctly Rbuildignored"
rm pigauto_*.tar.gz
```

Expected: `OK: CI files correctly Rbuildignored`.

- [ ] **Step 4: Final commit (if any tweaks)**

If steps 1-3 surface anything that needs fixing, fix and:

```bash
git add -p
git commit -m "ci(bace-h2h): fixes from final R CMD check pass"
```

If nothing needs fixing, skip this step.

---

## Task 11: Push branch + open PR

**Files:** none new.

- [ ] **Step 1: Push the feature branch**

```bash
git push -u origin feature/bace-headtohead-ci
```

- [ ] **Step 2: Open the PR**

```bash
gh pr create --title "ci: pigauto vs BACE cross-dataset benchmark workflow" --body "$(cat <<'EOF'
## Summary
- New `.github/workflows/bace-headtohead.yml` mirrors `daniel1noble/BACE`'s
  `benchmark.yml`: 6-dataset matrix on `ubuntu-latest`, manual
  `workflow_dispatch` trigger, 5h45m timeout, heartbeat trick, per-dataset
  artifacts, aggregate job that produces a head-to-head report.
- Documented-default pigauto CI config in `script/gha/_ci_config.R`.
- Pinned BACE snapshot under `useful/bace_results_snapshot/` (sourced
  from Dan's Actions run 25329857467).
- Per-run summaries land as draft PRs under
  `useful/ci_runs/<date>_<run-id>/`.
- Zero changes to `R/`, `DESCRIPTION`, `tests/` outside the new
  `test-gha-helpers.R` file, `NAMESPACE`, or `NEWS.md`.

## Spec
`specs/2026-05-16-bace-headtohead-ci-design.md` (committed `6311a6e`).

## Test plan
- [x] `devtools::test()` — all green incl. new `test-gha-helpers.R`
- [x] `devtools::check(args = "--no-manual")` — Status: OK
- [x] `R CMD build` tarball does NOT include CI infra
- [x] Locally smoke-tested `run_bench_globtherm.R` end-to-end
- [ ] (post-merge) trigger workflow_dispatch once and inspect the
  resulting draft PR

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

- [ ] **Step 3: Verify CI checks on the PR**

```bash
gh pr checks
```

Expected: existing `R-CMD-check` workflow passes (the new
`bace-headtohead.yml` does NOT auto-run because it's
`workflow_dispatch`-only).

---

## Verification (end-state assertions)

After all tasks complete:

1. `devtools::test()` — full suite passes; `test-gha-helpers.R` shows
   5 tests / ≥5 expectations.
2. `devtools::check()` — `Status: OK`.
3. `R CMD build` tarball excludes `script/gha/`, `useful/ci_runs/`,
   `useful/bace_results_snapshot/`, and `.github/workflows/bace-headtohead.yml`.
4. `.github/workflows/bace-headtohead.yml` parses (yq / actionlint).
5. `useful/bace_results_snapshot/` contains 6 `.rds` files + README.
6. `script/gha/` contains: 1 config, 6 wrappers, 1 aggregator, 1
   stager, 1 snapshot tool, 1 README, 1 `.gitkeep` in `results/`.
7. One manual `workflow_dispatch` produces:
   - 6 `bench-<dataset>` artifacts
   - 1 `cross-dataset-summary` artifact
   - Workflow summary tab with the full report.md
   - A draft PR titled "CI run <run-id> — pigauto vs BACE" carrying
     a new `useful/ci_runs/<date>_<run-id>/` directory.

## Out of scope (deferred)

- Adding the workflow to a regular `schedule:` cron (gated on 3+
  clean manual runs first)
- Wiring the per-run report into the pkgdown site
- Adding `phyloTraitData` Delhey-5809 as a 7th matrix entry
- macOS / Windows runners
- Replacing snapshot-mode BACE with a live BACE fit inside CI
