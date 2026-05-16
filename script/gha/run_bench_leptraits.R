#!/usr/bin/env Rscript
# script/gha/run_bench_leptraits.R
#
# CI wrapper for the LepTraits butterflies bench (BACE-compat).

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  here <- "."
  source(file.path(here, "script", "gha", "_ci_config.R"))
  source(file.path(here, "script", "gha", "_bace_compat.R"))
  devtools::load_all(here, quiet = TRUE)
})

DATASET <- "leptraits"
OUT_DIR <- file.path("script", "gha", "results", DATASET)

tryCatch({
  load(file.path("useful", "bace_data_snapshot", "data", "leptraits_traits.rda"))
  load(file.path("useful", "bace_data_snapshot", "data", "leptraits_tree.rda"))
  .run_bench_bace_compat(traits_df = leptraits_traits,
                          tree      = leptraits_tree,
                          dataset   = DATASET,
                          out_dir   = OUT_DIR,
                          cfg       = PIGAUTO_CI_CONFIG)
}, error = function(e) {
  .write_skip_marker(DATASET, OUT_DIR, conditionMessage(e))
})
