#!/usr/bin/env Rscript
# script/gha/run_bench_avonet.R
#
# CI wrapper for the AVONET head-to-head bench (BACE-compat).
# Loads avonet_traits + avonet_tree from useful/bace_data_snapshot/data/
# — the same .rda files BACE used in run 25329857467. Then applies
# BACE's exact subset+mask procedure (seed=2026, n=2000, 30% MCAR).
#
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  here <- "."
  source(file.path(here, "script", "gha", "_ci_config.R"))
  source(file.path(here, "script", "gha", "_bace_compat.R"))
  devtools::load_all(here, quiet = TRUE)
})

DATASET <- "avonet"
OUT_DIR <- file.path("script", "gha", "results", DATASET)

tryCatch({
  load(file.path("useful", "bace_data_snapshot", "data", "avonet_traits.rda"))
  load(file.path("useful", "bace_data_snapshot", "data", "avonet_tree.rda"))
  .run_bench_bace_compat(traits_df = avonet_traits,
                          tree      = avonet_tree,
                          dataset   = DATASET,
                          out_dir   = OUT_DIR,
                          cfg       = PIGAUTO_CI_CONFIG)
}, error = function(e) {
  .write_skip_marker(DATASET, OUT_DIR, conditionMessage(e))
})
