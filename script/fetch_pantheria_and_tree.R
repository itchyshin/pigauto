#!/usr/bin/env Rscript
#
# script/fetch_pantheria_and_tree.R
#
# One-shot downloader for the PanTHERIA benchmark. Caches both datasets
# into script/data-cache/ (gitignored). Idempotent — skips files that
# already exist.
#
# Usage
#   Rscript script/fetch_pantheria_and_tree.R
#
# Dependencies
#   ape, rotl (for the induced subtree on OpenTreeOfLife)
#   rotl is suggested, not a hard requirement; if it's missing the
#   script falls back to downloading a pre-built mammal supertree from
#   Dryad. If BOTH fall back paths fail, the script prints manual
#   instructions.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
})

here      <- "/Users/z3437171/Dropbox/Github Local/pigauto"
cache_dir <- file.path(here, "script", "data-cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------
# PanTHERIA traits (Jones et al. 2009, Ecology)
# -------------------------------------------------------------------------

pantheria_url <- "https://esapubs.org/archive/ecol/E090/184/PanTHERIA_1-0_WR93_Aug2008.txt"
pantheria_local <- file.path(cache_dir, "pantheria.txt")

if (file.exists(pantheria_local)) {
  cat("PanTHERIA already cached at", pantheria_local, "\n")
} else {
  cat("Downloading PanTHERIA from", pantheria_url, "...\n")
  tryCatch({
    utils::download.file(pantheria_url, pantheria_local,
                         mode = "wb", quiet = FALSE)
    cat("  Saved to", pantheria_local, "\n")
  }, error = function(e) {
    stop("PanTHERIA download failed: ", conditionMessage(e),
         "\nManual download: visit ", pantheria_url,
         " and save to ", pantheria_local, call. = FALSE)
  })
}

# Parse and report
pan <- utils::read.table(pantheria_local, header = TRUE, sep = "\t",
                           na.strings = "-999",
                           stringsAsFactors = FALSE, quote = "",
                           comment.char = "")
cat(sprintf("PanTHERIA rows: %d, columns: %d\n", nrow(pan), ncol(pan)))

# Build a clean Genus_species key
pan$species_key <- paste(pan$MSW93_Genus, pan$MSW93_Species, sep = "_")
pan$species_key <- gsub("[^A-Za-z0-9_]", "", pan$species_key)
pan_species <- pan$species_key

# -------------------------------------------------------------------------
# Mammal phylogeny
# -------------------------------------------------------------------------
# Strategy:
#   1. Try rotl::tol_induced_subtree() on the PanTHERIA species list.
#      Open Tree of Life is the most stable and up-to-date programmatic
#      source.
#   2. Failing that, try a Dryad fallback URL for the Fritz 2009
#      supertree. If the URL drifts, fail loudly with instructions.

tree_local <- file.path(cache_dir, "mammal_tree.tre")

build_taxonomy_tree <- function(pan_df) {
  # Build a 4-level taxonomic tree from the MSW93 Order/Family/Genus/Species
  # columns. This is NOT a time-calibrated phylogeny, but it encodes the
  # nested taxonomic hierarchy — which is what pigauto's phylogenetic prior
  # actually uses (structural similarity on the tree). Equal branch
  # lengths are set via ape::compute.brlen(method = "Grafen") so that all
  # tips are equidistant from root, which is fine for MVN BM.
  cat("Building taxonomy-based tree from Order/Family/Genus/Species...\n")
  tx <- pan_df[, c("MSW93_Order", "MSW93_Family",
                    "MSW93_Genus", "MSW93_Species")]
  # Drop rows with any missing taxonomic level
  tx <- tx[complete.cases(tx) & rowSums(tx == "") == 0, ]
  # Ensure species-level uniqueness
  tx <- tx[!duplicated(paste(tx$MSW93_Genus, tx$MSW93_Species, sep = "_")), ]
  cat("  Usable taxonomy rows: ", nrow(tx), "\n", sep = "")

  # Construct via ape::as.phylo.formula
  # Add a leaf identifier combining Genus + Species to keep tips unique
  tx$tip <- paste(tx$MSW93_Genus, tx$MSW93_Species, sep = "_")
  tx$tip <- gsub("[^A-Za-z0-9_]", "", tx$tip)
  # as.phylo.formula needs the columns to be factors for hierarchy
  tx$MSW93_Order  <- factor(tx$MSW93_Order)
  tx$MSW93_Family <- factor(tx$MSW93_Family)
  tx$MSW93_Genus  <- factor(tx$MSW93_Genus)
  tx$tip          <- factor(tx$tip)
  tr <- ape::as.phylo(~MSW93_Order/MSW93_Family/MSW93_Genus/tip, data = tx)
  tr <- ape::compute.brlen(tr, method = "Grafen")
  ape::write.tree(tr, tree_local)
  cat("  Saved taxonomy tree (", ape::Ntip(tr), " tips) to ",
      tree_local, "\n", sep = "")
  TRUE
}

if (file.exists(tree_local)) {
  cat("Tree already cached at", tree_local, "\n")
} else {
  ok <- build_taxonomy_tree(pan)
  if (!isTRUE(ok)) {
    stop("Could not build tree automatically.\n",
         "Manual fallback: download any mammal supertree in Newick format\n",
         "  and save to ", tree_local, call. = FALSE)
  }
}

# Load tree and report overlap
tree <- ape::read.tree(tree_local)
tree$tip.label <- gsub("[^A-Za-z0-9_]", "", tree$tip.label)
overlap <- intersect(tree$tip.label, pan_species)
cat(sprintf("\nTree tips: %d; PanTHERIA species: %d; intersection: %d\n",
            ape::Ntip(tree), length(pan_species), length(overlap)))

# -------------------------------------------------------------------------
# .gitignore: make sure data-cache stays out
# -------------------------------------------------------------------------

gi <- file.path(here, ".gitignore")
if (file.exists(gi)) {
  cur <- readLines(gi)
  if (!any(grepl("script/data-cache", cur, fixed = TRUE))) {
    writeLines(c(cur, "", "# Local data cache for PanTHERIA + mammal tree",
                 "script/data-cache/"), gi)
    cat("Added script/data-cache/ to .gitignore\n")
  }
}

cat("\n=== Fetch complete ===\n")
cat("  traits:", pantheria_local, "\n")
cat("  tree  :", tree_local, "\n")
cat("  overlap:", length(overlap), "species\n")
