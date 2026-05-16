#!/usr/bin/env Rscript
# script/gha/run_bench_globtherm.R
#
# CI wrapper for the GlobTherm thermal-tolerance bench.
# Bennett et al. 2018, Sci. Data 5: 180022. Figshare DOI 10.6084/m9.figshare.5749274.
# Tree built taxonomically from Class/Order/Family/Genus/Species.
#
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  here <- "."
  source(file.path(here, "script", "gha", "_ci_config.R"))
  devtools::load_all(here, quiet = TRUE)
})

DATASET <- "globtherm"
OUT_DIR <- file.path("script", "gha", "results", DATASET)

tryCatch({
  # Figshare CSV — Bennett 2018 GlobTherm
  csv_path <- .download_to_cache(
    "https://ndownloader.figshare.com/files/9885149",
    "globTherm.csv"
  )
  gt <- utils::read.csv(csv_path, stringsAsFactors = FALSE)

  ECTO_CLASSES <- c("Lepidosauria", "Insecta", "Amphibia", "Actinopteri",
                     "Arachnida", "Malacostraca", "Bivalvia", "Gastropoda")
  if ("Class" %in% colnames(gt)) {
    gt <- gt[gt$Class %in% ECTO_CLASSES, , drop = FALSE]
  }
  trait_cols <- intersect(c("Tmax", "tmin"), colnames(gt))
  stopifnot(length(trait_cols) > 0L)
  ok <- rowSums(!is.na(gt[, trait_cols, drop = FALSE])) > 0L &
        !is.na(gt$Class) & !is.na(gt$Order) & !is.na(gt$Family) &
        nzchar(gt$Class) & nzchar(gt$Order) & nzchar(gt$Family) &
        !is.na(gt$Genus) & nzchar(gt$Genus)
  gt <- gt[ok, , drop = FALSE]
  sp_col <- if ("scientificNameStd" %in% colnames(gt)) "scientificNameStd" else "Species"
  gt$Species_Key <- gsub(" ", "_", gt[[sp_col]])
  gt <- gt[!duplicated(gt$Species_Key) & nzchar(gt$Species_Key), , drop = FALSE]
  rownames(gt) <- gt$Species_Key

  df_traits <- gt[, trait_cols, drop = FALSE]
  for (v in trait_cols) df_traits[[v]] <- suppressWarnings(as.numeric(df_traits[[v]]))

  tax <- data.frame(
    Order = gt$Order, Family = gt$Family, Genus = gt$Genus,
    Species = gt$Species_Key, row.names = gt$Species_Key,
    stringsAsFactors = FALSE
  )
  tree <- .build_tax_tree(tax, seed = PIGAUTO_CI_CONFIG$seed)
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(df_traits)))
  df_traits <- df_traits[tree$tip.label, , drop = FALSE]

  .run_bench(df_traits, tree, dataset = DATASET, out_dir = OUT_DIR,
             cfg = PIGAUTO_CI_CONFIG, log_transform = FALSE)
}, error = function(e) {
  .write_skip_marker(DATASET, OUT_DIR, conditionMessage(e))
})
