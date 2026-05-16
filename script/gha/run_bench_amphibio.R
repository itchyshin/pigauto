#!/usr/bin/env Rscript
# script/gha/run_bench_amphibio.R
#
# CI wrapper for the AmphiBIO bench.
# Downloads AmphiBIO_v1.zip (1.4 MB, figshare 4644424) then builds a
# taxonomic tree from Order/Family/Genus/Species (Grafen brlens).
#
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  here <- "."
  source(file.path(here, "script", "gha", "_ci_config.R"))
  devtools::load_all(here, quiet = TRUE)
})

DATASET <- "amphibio"
OUT_DIR <- file.path("script", "gha", "results", DATASET)

tryCatch({
  zip_path <- .download_to_cache(
    "https://ndownloader.figshare.com/files/8828578",
    "AmphiBIO_v1.zip"
  )
  csv_name <- "AmphiBIO_v1.csv"
  csv_dir  <- .cache_path("amphibio")
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  csv_path <- file.path(csv_dir, csv_name)
  if (!file.exists(csv_path)) {
    utils::unzip(zip_path, files = csv_name, exdir = csv_dir)
  }
  amphi <- utils::read.csv(csv_path, stringsAsFactors = FALSE)

  tax_cols  <- c("Order", "Family", "Genus", "Species")
  pick      <- c("Body_size_mm", "Body_mass_g")
  pick      <- intersect(pick, names(amphi))
  stopifnot(all(tax_cols %in% names(amphi)), length(pick) > 0L)

  df <- amphi[, c(tax_cols, pick), drop = FALSE]
  df <- df[nzchar(df$Species) & !is.na(df$Species), , drop = FALSE]
  df <- df[!duplicated(df$Species), , drop = FALSE]
  for (v in pick) df[[v]] <- suppressWarnings(as.numeric(df[[v]]))
  # Drop rows with all picked-trait NA
  df <- df[rowSums(!is.na(df[, pick, drop = FALSE])) > 0L, , drop = FALSE]

  rownames(df) <- df$Species
  tree <- .build_tax_tree(df[, tax_cols, drop = FALSE],
                           seed = PIGAUTO_CI_CONFIG$seed)
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(df)))
  df   <- df[tree$tip.label, , drop = FALSE]
  # Hand .run_bench only the trait columns; taxonomy columns are GNN-irrelevant.
  df_traits <- df[, pick, drop = FALSE]

  .run_bench(df_traits, tree, dataset = DATASET, out_dir = OUT_DIR,
             cfg = PIGAUTO_CI_CONFIG, log_transform = TRUE)
}, error = function(e) {
  .write_skip_marker(DATASET, OUT_DIR, conditionMessage(e))
})
