#!/usr/bin/env Rscript
# script/gha/snapshot_bace.R
#
# One-off Mac helper: re-snapshots BACE's per-dataset benchmark
# outputs into pigauto's useful/bace_results_snapshot/.
#
# Usage:
#   gh run download <bace-run-id> -R daniel1noble/BACE -D /tmp/bace_art
#   Rscript script/gha/snapshot_bace.R /tmp/bace_art
#
# BACE's artifact format (as of run 25329857467, 2026-05-04):
#   <art_root>/bench-<dataset>/run_*_<date>/summary_metrics.csv
#   <art_root>/bench-<dataset>/run_*_<date>/run_info.csv
#
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md §4.6.

DATASETS <- c("avonet","pantheria","amphibio","bien","globtherm","leptraits")

# Project ONE BACE summary table to the canonical schema.
# `bace_native` is the contents of summary_metrics.csv (wide format,
# one row per (method, trait)). `run_info` is run_info.csv (typically
# 1 row per dataset with runtime_min).
snapshot_bace_one <- function(bace_native, run_info, dataset) {
  stopifnot(is.data.frame(bace_native), is.character(dataset))
  # Drop non-bace methods (mean_baseline / column_mean from BACE's own bench)
  bace_only <- bace_native[tolower(bace_native$method) == "bace", , drop = FALSE]
  if (nrow(bace_only) == 0L) {
    return(.normalize_eval(
      data.frame(trait = character(0), type = character(0),
                 imputation_idx = integer(0), rmse = numeric(0),
                 mae = numeric(0), pearson_r = numeric(0),
                 accuracy = numeric(0), brier = numeric(0),
                 time_sec = numeric(0)),
      dataset = dataset, method = "BACE_snapshot"
    ))
  }
  rt_sec <- if (!is.null(run_info) && "runtime_min" %in% colnames(run_info)) {
    as.numeric(run_info$runtime_min[1]) * 60
  } else NA_real_

  out <- data.frame(
    trait          = as.character(bace_only$trait),
    type           = as.character(bace_only$type),
    imputation_idx = 1L,
    rmse           = suppressWarnings(as.numeric(bace_only$rmse)),
    mae            = suppressWarnings(as.numeric(bace_only$mae_fit)),
    pearson_r      = suppressWarnings(as.numeric(bace_only$correlation)),
    accuracy       = suppressWarnings(as.numeric(bace_only$accuracy)),
    brier          = suppressWarnings(as.numeric(bace_only$brier)),
    time_sec       = rt_sec,
    stringsAsFactors = FALSE
  )
  .normalize_eval(out, dataset = dataset, method = "BACE_snapshot")
}

# Find <art_root>/bench-<dataset>/run_*_<date>/summary_metrics.csv.
.find_bace_csv <- function(art_root, dataset, name) {
  cand <- list.files(file.path(art_root, paste0("bench-", dataset)),
                      pattern = paste0("^", name, "$"),
                      recursive = TRUE, full.names = TRUE)
  if (length(cand) == 0L) NA_character_ else cand[[1L]]
}

# When invoked as a script, walk the artifacts directory.
if (sys.nframe() == 0L) {
  source(file.path("script", "gha", "_ci_config.R"))
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) {
    stop("usage: Rscript script/gha/snapshot_bace.R <path-to-bace-artifacts>")
  }
  src_dir <- args[[1L]]
  out_dir <- file.path("useful", "bace_results_snapshot")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (d in DATASETS) {
    sm_path <- .find_bace_csv(src_dir, d, "summary_metrics.csv")
    ri_path <- .find_bace_csv(src_dir, d, "run_info.csv")
    if (is.na(sm_path)) {
      message(sprintf("[snapshot] no summary_metrics.csv for %s; skipping", d))
      next
    }
    sm <- utils::read.csv(sm_path, stringsAsFactors = FALSE)
    ri <- if (!is.na(ri_path)) utils::read.csv(ri_path, stringsAsFactors = FALSE) else NULL
    out <- snapshot_bace_one(sm, run_info = ri, dataset = d)
    saveRDS(out, file.path(out_dir, paste0(d, ".rds")))
    cat(sprintf("[snapshot] %s: %d rows -> %s\n",
                d, nrow(out), file.path(out_dir, paste0(d, ".rds"))))
  }
  cat("[snapshot] done. Remember to update useful/bace_results_snapshot/README.md\n")
}
