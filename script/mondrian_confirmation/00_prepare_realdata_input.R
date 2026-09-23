#!/usr/bin/env Rscript
# Prepare a canonical observed-data input for the Stage-C Mondrian pilot.
#
# Usage:
#   Rscript 00_prepare_realdata_input.R pantheria output.rds cache_dir [n_subset]
#   Rscript 00_prepare_realdata_input.R fishbase   output.rds cache_dir [n_subset]
#
# The output deliberately contains native missing values, but the runner masks
# only originally observed values and retains that mask as its ground-truth
# receipt.  `cache_dir` must already contain the downloaded source data; this
# script never silently substitutes a synthetic tree or data set.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 4L) {
  stop("expected: dataset output.rds cache_dir [n_subset]", call. = FALSE)
}
dataset <- match.arg(args[[1L]], c("pantheria", "fishbase"))
out_file <- args[[2L]]
cache_dir <- args[[3L]]
n_subset <- if (length(args) == 4L) as.integer(args[[4L]]) else NA_integer_
if (!is.na(n_subset) && (!is.finite(n_subset) || n_subset < 30L)) {
  stop("n_subset must be at least 30", call. = FALSE)
}

suppressPackageStartupMessages(library(ape))

align_input <- function(data, tree, label) {
  if (is.null(rownames(data)) || anyDuplicated(rownames(data))) {
    stop(label, " data require unique row names", call. = FALSE)
  }
  overlap <- intersect(tree$tip.label, rownames(data))
  if (length(overlap) < 30L) stop(label, " has fewer than 30 matched tips", call. = FALSE)
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, overlap))
  data <- data[tree$tip.label, , drop = FALSE]
  if (!identical(rownames(data), tree$tip.label)) stop("alignment failed", call. = FALSE)
  list(data = data, tree = tree)
}

prepare_pantheria <- function() {
  pan_file <- file.path(cache_dir, "pantheria.txt")
  tree_file <- file.path(cache_dir, "mammal_tree.tre")
  if (!file.exists(pan_file) || !file.exists(tree_file)) {
    stop("PanTHERIA cache missing pantheria.txt or mammal_tree.tre; ",
         "run script/fetch_pantheria_and_tree.R first", call. = FALSE)
  }
  pan <- utils::read.table(pan_file, header = TRUE, sep = "\t", na.strings = "-999",
                           stringsAsFactors = FALSE, quote = "", comment.char = "")
  pan$species_key <- gsub("[^A-Za-z0-9_]", "",
                          paste(pan$MSW93_Genus, pan$MSW93_Species, sep = "_"))
  tree <- ape::read.tree(tree_file)
  tree$tip.label <- gsub("[^A-Za-z0-9_]", "", tree$tip.label)
  spec <- list(
    body_mass_g = list(src = "X5.1_AdultBodyMass_g", kind = "log"),
    head_body_length_mm = list(src = "X13.1_AdultHeadBodyLen_mm", kind = "log"),
    gestation_d = list(src = "X9.1_GestationLen_d", kind = "log"),
    litter_size = list(src = "X15.1_LitterSize", kind = "count"),
    max_longevity_m = list(src = "X17.1_MaxLongevity_m", kind = "log"),
    diet_breadth = list(src = "X6.1_DietBreadth", kind = "ordinal"),
    habitat_breadth = list(src = "X12.1_HabitatBreadth", kind = "ordinal"),
    terrestriality = list(src = "X12.2_Terrestriality", kind = "factor")
  )
  spec <- spec[vapply(spec, function(x) x$src %in% names(pan), logical(1))]
  out <- data.frame(row.names = pan$species_key)
  for (nm in names(spec)) {
    x <- suppressWarnings(as.numeric(pan[[spec[[nm]]$src]]))
    kind <- spec[[nm]]$kind
    if (kind == "log") {
      x[!is.finite(x) | x <= 0] <- NA_real_
      out[[nm]] <- log(x)
    } else if (kind == "count") {
      x[!is.finite(x)] <- NA_real_
      out[[nm]] <- as.integer(round(x))
    } else if (kind == "ordinal") {
      x[!is.finite(x)] <- NA_real_
      out[[nm]] <- factor(as.character(x), levels = as.character(sort(unique(stats::na.omit(x)))), ordered = TRUE)
    } else {
      x[!is.finite(x)] <- NA_real_
      out[[nm]] <- factor(as.character(x), levels = as.character(sort(unique(stats::na.omit(x)))))
    }
  }
  out <- out[rowSums(!is.na(out)) > 0L, , drop = FALSE]
  align_input(out, tree, "PanTHERIA")
}

prepare_fishbase <- function() {
  if (!requireNamespace("rfishbase", quietly = TRUE) || !requireNamespace("fishtree", quietly = TRUE)) {
    stop("FishBase preparation requires rfishbase and fishtree", call. = FALSE)
  }
  tree <- fishtree::fishtree_phylogeny()
  tree$tip.label <- gsub("_", " ", tree$tip.label)
  taxa <- rfishbase::load_taxa()[, c("SpecCode", "Species"), drop = FALSE]
  species <- rfishbase::species()
  keep <- intersect(c("SpecCode", "Length", "Weight", "BodyShapeI", "DepthRangeDeep", "Vulnerability"), names(species))
  species <- merge(taxa, species[, keep, drop = FALSE], by = "SpecCode", all.x = FALSE)
  ecology <- rfishbase::ecology()
  troph <- data.frame(SpecCode = integer(), Troph = numeric())
  if ("SpecCode" %in% names(ecology)) {
    if (all(c("DietTroph", "FoodTroph") %in% names(ecology))) {
      troph <- transform(ecology[, c("SpecCode", "DietTroph", "FoodTroph"), drop = FALSE],
                         Troph = ifelse(is.finite(DietTroph), DietTroph, FoodTroph))[, c("SpecCode", "Troph")]
    } else if ("DietTroph" %in% names(ecology)) {
      troph <- ecology[, c("SpecCode", "DietTroph"), drop = FALSE]; names(troph)[2L] <- "Troph"
    } else if ("FoodTroph" %in% names(ecology)) {
      troph <- ecology[, c("SpecCode", "FoodTroph"), drop = FALSE]; names(troph)[2L] <- "Troph"
    }
    troph <- troph[is.finite(troph$Troph) & !duplicated(troph$SpecCode), , drop = FALSE]
  }
  out <- merge(species, troph, by = "SpecCode", all.x = TRUE)
  out <- out[out$Species %in% tree$tip.label, , drop = FALSE]
  rownames(out) <- out$Species
  out <- out[, intersect(c("Length", "Weight", "DepthRangeDeep", "Vulnerability", "Troph", "BodyShapeI"), names(out)), drop = FALSE]
  for (nm in setdiff(names(out), "BodyShapeI")) out[[nm]] <- suppressWarnings(as.numeric(out[[nm]]))
  if ("BodyShapeI" %in% names(out)) {
    shape <- tolower(trimws(as.character(out$BodyShapeI)))
    shape[shape == ""] <- NA_character_
    keep_shape <- names(table(shape, useNA = "no"))[table(shape, useNA = "no") >= 50L]
    shape[!shape %in% keep_shape] <- NA_character_
    out$BodyShapeI <- droplevels(factor(shape))
  }
  out <- out[rowSums(!is.na(out)) > 0L, , drop = FALSE]
  align_input(out, tree, "FishBase")
}

prepared <- if (identical(dataset, "pantheria")) prepare_pantheria() else prepare_fishbase()
if (!is.na(n_subset) && n_subset < nrow(prepared$data)) {
  set.seed(20260818L)
  keep <- sort(sample(rownames(prepared$data), n_subset))
  prepared$tree <- ape::drop.tip(prepared$tree, setdiff(prepared$tree$tip.label, keep))
  prepared$data <- prepared$data[prepared$tree$tip.label, , drop = FALSE]
}
if (any(vapply(prepared$data, function(x) sum(!is.na(x)) < 20L, logical(1)))) {
  stop("prepared data contain a trait with fewer than 20 observed values", call. = FALSE)
}
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(c(prepared, list(dataset = dataset, prepared_at = format(Sys.time(), tz = "UTC", usetz = TRUE))), out_file)
cat(sprintf("Prepared %s: %d tips x %d traits -> %s\n", dataset, nrow(prepared$data), ncol(prepared$data), out_file))
print(vapply(prepared$data, function(x) sum(!is.na(x)), integer(1)))
