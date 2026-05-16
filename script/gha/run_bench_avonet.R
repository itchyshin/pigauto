#!/usr/bin/env Rscript
# script/gha/run_bench_avonet.R
#
# CI wrapper for the AVONET head-to-head bench.
# Bundled data, so this should never need network access.
#
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  here <- "."
  source(file.path(here, "script", "gha", "_ci_config.R"))
  devtools::load_all(here, quiet = TRUE)
})

DATASET <- "avonet"
OUT_DIR <- file.path("script", "gha", "results", DATASET)

tryCatch({
  e <- new.env()
  utils::data("avonet300", package = "pigauto", envir = e)
  utils::data("tree300",   package = "pigauto", envir = e)
  df <- e$avonet300
  rownames(df) <- df$Species_Key
  df$Species_Key <- NULL
  tree <- e$tree300
  .run_bench(df, tree, dataset = DATASET, out_dir = OUT_DIR,
             cfg = PIGAUTO_CI_CONFIG, log_transform = TRUE)
}, error = function(e) {
  .write_skip_marker(DATASET, OUT_DIR, conditionMessage(e))
})
