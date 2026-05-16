#!/usr/bin/env Rscript
# script/gha/run_bench_bien.R
#
# CI wrapper for the BIEN plants bench (BACE-compat).
# Uses BACE's pre-built bien_traits + bien_tree, so no BIEN or
# V.PhyloMaker2 install is needed in CI.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  here <- "."
  source(file.path(here, "script", "gha", "_ci_config.R"))
  source(file.path(here, "script", "gha", "_bace_compat.R"))
  devtools::load_all(here, quiet = TRUE)
})

DATASET <- "bien"
OUT_DIR <- file.path("script", "gha", "results", DATASET)

tryCatch({
  load(file.path("useful", "bace_data_snapshot", "data", "bien_traits.rda"))
  load(file.path("useful", "bace_data_snapshot", "data", "bien_tree.rda"))
  .run_bench_bace_compat(traits_df = bien_traits,
                          tree      = bien_tree,
                          dataset   = DATASET,
                          out_dir   = OUT_DIR,
                          cfg       = PIGAUTO_CI_CONFIG)
}, error = function(e) {
  .write_skip_marker(DATASET, OUT_DIR, conditionMessage(e))
})
