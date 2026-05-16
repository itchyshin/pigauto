#!/usr/bin/env Rscript
# script/gha/run_bench_pantheria.R
#
# CI wrapper for the PanTHERIA mammals bench (BACE-compat).
# Loads pantheria_traits + pantheria_tree from BACE's snapshot.
#
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  here <- "."
  source(file.path(here, "script", "gha", "_ci_config.R"))
  source(file.path(here, "script", "gha", "_bace_compat.R"))
  devtools::load_all(here, quiet = TRUE)
})

DATASET <- "pantheria"
OUT_DIR <- file.path("script", "gha", "results", DATASET)

tryCatch({
  load(file.path("useful", "bace_data_snapshot", "data", "pantheria_traits.rda"))
  load(file.path("useful", "bace_data_snapshot", "data", "pantheria_tree.rda"))
  .run_bench_bace_compat(traits_df = pantheria_traits,
                          tree      = pantheria_tree,
                          dataset   = DATASET,
                          out_dir   = OUT_DIR,
                          cfg       = PIGAUTO_CI_CONFIG)
}, error = function(e) {
  .write_skip_marker(DATASET, OUT_DIR, conditionMessage(e))
})
