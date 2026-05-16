#!/usr/bin/env Rscript
# script/gha/run_bench_globtherm.R
#
# CI wrapper for the GlobTherm thermal-tolerance bench (BACE-compat).

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  here <- "."
  source(file.path(here, "script", "gha", "_ci_config.R"))
  source(file.path(here, "script", "gha", "_bace_compat.R"))
  devtools::load_all(here, quiet = TRUE)
})

DATASET <- "globtherm"
OUT_DIR <- file.path("script", "gha", "results", DATASET)

tryCatch({
  load(file.path("useful", "bace_data_snapshot", "data", "globtherm_traits.rda"))
  load(file.path("useful", "bace_data_snapshot", "data", "globtherm_tree.rda"))
  .run_bench_bace_compat(traits_df = globtherm_traits,
                          tree      = globtherm_tree,
                          dataset   = DATASET,
                          out_dir   = OUT_DIR,
                          cfg       = PIGAUTO_CI_CONFIG)
}, error = function(e) {
  .write_skip_marker(DATASET, OUT_DIR, conditionMessage(e))
})
