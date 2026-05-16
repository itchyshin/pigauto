#!/usr/bin/env Rscript
# script/gha/run_bench_leptraits.R
#
# CI wrapper for the LepTraits 1.0 butterfly bench.
# Shirey et al. 2022, Sci. Data 9: 382. Figshare 10.6084/m9.figshare.16974106.
# Tree built taxonomically from Family/Genus/Species (single Order, Lepidoptera).
#
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md.
#
# NOTE: LepTraits download URL is best-effort; if figshare schema drifts
# the wrapper falls through to the skip-marker.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  here <- "."
  source(file.path(here, "script", "gha", "_ci_config.R"))
  devtools::load_all(here, quiet = TRUE)
})

DATASET <- "leptraits"
OUT_DIR <- file.path("script", "gha", "results", DATASET)

tryCatch({
  # LepTraits consensus.csv from figshare 16974106
  csv_path <- .download_to_cache(
    "https://ndownloader.figshare.com/files/30884420",
    "leptraits_consensus.csv"
  )
  lt <- utils::read.csv(csv_path, stringsAsFactors = FALSE)

  trait_pick <- intersect(c("WS_L", "FW_L", "FlightDuration",
                              "NumberOfHostplantFamilies"),
                          colnames(lt))
  stopifnot(length(trait_pick) > 0L)
  tax_cols <- intersect(c("Family", "Genus", "Species"), colnames(lt))
  stopifnot(all(c("Family", "Genus", "Species") %in% tax_cols))

  ok <- !is.na(lt$Family) & nzchar(lt$Family) &
        !is.na(lt$Genus) & nzchar(lt$Genus) &
        !is.na(lt$Species) & nzchar(lt$Species)
  lt <- lt[ok, , drop = FALSE]
  lt$Species_Key <- gsub(" ", "_",
                          paste(lt$Genus, lt$Species, sep = "_"))
  lt <- lt[!duplicated(lt$Species_Key), , drop = FALSE]
  rownames(lt) <- lt$Species_Key

  for (v in trait_pick) lt[[v]] <- suppressWarnings(as.numeric(lt[[v]]))
  df_traits <- lt[, trait_pick, drop = FALSE]
  df_traits <- df_traits[rowSums(!is.na(df_traits)) > 0L, , drop = FALSE]
  lt <- lt[rownames(df_traits), , drop = FALSE]

  # LepTraits is all Lepidoptera — fake an "Order" col so .build_tax_tree works.
  tax <- data.frame(
    Order   = "Lepidoptera",
    Family  = lt$Family, Genus = lt$Genus,
    Species = lt$Species_Key,
    row.names = lt$Species_Key, stringsAsFactors = FALSE
  )
  tree <- .build_tax_tree(tax, seed = PIGAUTO_CI_CONFIG$seed)
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(df_traits)))
  df_traits <- df_traits[tree$tip.label, , drop = FALSE]

  .run_bench(df_traits, tree, dataset = DATASET, out_dir = OUT_DIR,
             cfg = PIGAUTO_CI_CONFIG, log_transform = TRUE)
}, error = function(e) {
  .write_skip_marker(DATASET, OUT_DIR, conditionMessage(e))
})
