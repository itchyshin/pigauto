#!/usr/bin/env Rscript
# script/gha/run_bench_pantheria.R
#
# CI wrapper for the PanTHERIA mammals bench.
# Downloads PanTHERIA_1-0 (~1 MB tab-separated) from esapubs.
# Tree is built taxonomically from MSW93 Order/Family/Genus/Species.
#
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  here <- "."
  source(file.path(here, "script", "gha", "_ci_config.R"))
  devtools::load_all(here, quiet = TRUE)
})

DATASET <- "pantheria"
OUT_DIR <- file.path("script", "gha", "results", DATASET)

tryCatch({
  txt_path <- .download_to_cache(
    "https://esapubs.org/archive/ecol/E090/184/PanTHERIA_1-0_WR93_Aug2008.txt",
    "pantheria.txt"
  )
  pan <- utils::read.table(txt_path, header = TRUE, sep = "\t",
                            na.strings = "-999", stringsAsFactors = FALSE,
                            quote = "", comment.char = "")
  pan$Species_Key <- gsub("[^A-Za-z0-9_]", "",
                          paste(pan$MSW93_Genus, pan$MSW93_Species, sep = "_"))
  ok <- complete.cases(pan[, c("MSW93_Order","MSW93_Family","MSW93_Genus","MSW93_Species")]) &
        nzchar(pan$MSW93_Order) & nzchar(pan$MSW93_Family) &
        nzchar(pan$MSW93_Genus) & nzchar(pan$MSW93_Species)
  pan <- pan[ok, , drop = FALSE]
  pan <- pan[!duplicated(pan$Species_Key), , drop = FALSE]
  rownames(pan) <- pan$Species_Key

  # 8 canonical traits (mix of continuous + ordinal); pigauto auto-detects types.
  pick <- c("X5.1_AdultBodyMass_g", "X9.1_GestationLen_d",
             "X15.1_LitterSize", "X16.1_LittersPerYear",
             "X17.1_MaxLongevity_m", "X22.1_HomeRange_km2",
             "X23.1_SexualMaturityAge_d", "X26.4_GR_MRLat_dd")
  pick <- intersect(pick, colnames(pan))
  stopifnot(length(pick) > 0L)
  for (v in pick) pan[[v]] <- suppressWarnings(as.numeric(pan[[v]]))
  df_full <- pan[, pick, drop = FALSE]
  df_full <- df_full[rowSums(!is.na(df_full)) > 0L, , drop = FALSE]
  pan <- pan[rownames(df_full), , drop = FALSE]

  # Taxonomic tree from MSW93 columns
  tax <- data.frame(
    Order   = pan$MSW93_Order,
    Family  = pan$MSW93_Family,
    Genus   = pan$MSW93_Genus,
    Species = pan$Species_Key,
    row.names = pan$Species_Key,
    stringsAsFactors = FALSE
  )
  tree <- .build_tax_tree(tax, seed = PIGAUTO_CI_CONFIG$seed)
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(df_full)))
  df_full <- df_full[tree$tip.label, , drop = FALSE]

  .run_bench(df_full, tree, dataset = DATASET, out_dir = OUT_DIR,
             cfg = PIGAUTO_CI_CONFIG, log_transform = TRUE)
}, error = function(e) {
  .write_skip_marker(DATASET, OUT_DIR, conditionMessage(e))
})
