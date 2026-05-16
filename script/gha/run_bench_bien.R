#!/usr/bin/env Rscript
# script/gha/run_bench_bien.R
#
# CI wrapper for the BIEN plants bench.
# BIEN R package + V.PhyloMaker2 (GitHub-only) — both heavy installs.
# On first runs we expect this to skip-marker; once the workflow proves
# stable on the lighter datasets, the maintainer can wire up a proper
# BIEN install step (BIEN, V.PhyloMaker2) in the workflow yaml.
#
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  here <- "."
  source(file.path(here, "script", "gha", "_ci_config.R"))
  devtools::load_all(here, quiet = TRUE)
})

DATASET <- "bien"
OUT_DIR <- file.path("script", "gha", "results", DATASET)

tryCatch({
  if (!requireNamespace("BIEN", quietly = TRUE) ||
      !requireNamespace("V.PhyloMaker2", quietly = TRUE)) {
    stop("BIEN and V.PhyloMaker2 are required but not installed. ",
         "Wire up the workflow's setup-r-dependencies step to include them ",
         "once the lighter-dataset jobs are stable.")
  }

  trait_picks <- c(
    "leaf area",
    "leaf area per leaf dry mass",
    "stem wood density",
    "seed mass",
    "maximum whole plant height"
  )
  bien_long <- BIEN::BIEN_trait_traits_per_species(trait = trait_picks)
  if (is.null(bien_long) || nrow(bien_long) == 0L) {
    stop("BIEN_trait_traits_per_species returned no rows")
  }
  bien_long$species_key <- gsub(" ", "_", bien_long$scrubbed_species_binomial)
  bien_long$value <- suppressWarnings(as.numeric(bien_long$trait_value))
  bien_long <- bien_long[!is.na(bien_long$value), , drop = FALSE]

  trait_map <- c(
    "leaf area"                            = "leaf_area",
    "leaf area per leaf dry mass"          = "sla",
    "stem wood density"                    = "wood_density",
    "seed mass"                            = "seed_mass",
    "maximum whole plant height"           = "height_m"
  )
  agg <- stats::aggregate(value ~ species_key + trait_name,
                            data = bien_long, FUN = function(x)
                              stats::median(x, na.rm = TRUE))
  agg$short <- trait_map[agg$trait_name]
  wide <- reshape(agg[, c("species_key","short","value")],
                   timevar = "short", idvar = "species_key", direction = "wide")
  colnames(wide) <- sub("^value\\.", "", colnames(wide))
  wide <- wide[rowSums(!is.na(wide[, -1, drop = FALSE])) > 0L, , drop = FALSE]
  rownames(wide) <- wide$species_key
  df_traits <- wide[, setdiff(colnames(wide), "species_key"), drop = FALSE]

  # V.PhyloMaker2 lookup
  species_list <- data.frame(
    species = gsub("_", " ", rownames(df_traits)),
    genus   = vapply(strsplit(gsub("_", " ", rownames(df_traits)), " ", fixed = TRUE),
                     `[`, character(1), 1L),
    family  = NA_character_,
    stringsAsFactors = FALSE
  )
  out <- V.PhyloMaker2::phylo.maker(species_list, scenarios = "S3")
  tree <- out$scenario.3
  tree$tip.label <- gsub(" ", "_", tree$tip.label)
  keep_sp <- intersect(rownames(df_traits), tree$tip.label)
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep_sp))
  df_traits <- df_traits[tree$tip.label, , drop = FALSE]

  .run_bench(df_traits, tree, dataset = DATASET, out_dir = OUT_DIR,
             cfg = PIGAUTO_CI_CONFIG, log_transform = TRUE)
}, error = function(e) {
  .write_skip_marker(DATASET, OUT_DIR, conditionMessage(e))
})
