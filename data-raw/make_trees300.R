## data-raw/make_trees300.R
## Creates trees300 (multiPhylo): 50 posterior trees for the avonet300 species.
## Run from the package root: Rscript data-raw/make_trees300.R
##
## Requires: megatrees (>= 0.1.4), ape, pigauto (for tree300 tip names).
## Source: megatrees::tree_bird_n100 — 100 randomly selected posterior trees
## from the BirdTree Hackett backbone (Jetz et al. 2012).

library(ape)
library(megatrees)

# ---- Get the 300 species names ------------------------------------------------
# Load our bundled MCC tree to get the canonical species set
tree300 <- get(load(here::here("data", "tree300.rda")))
our_spp <- tree300$tip.label
stopifnot(length(our_spp) == 300L)

# ---- Verify all 300 species exist in the megatree ----------------------------
mega_tips <- tree_bird_n100[[1]]$tip.label
n_match   <- sum(our_spp %in% mega_tips)
cat("Species matching megatree:", n_match, "/ 300\n")
stopifnot(n_match == 300L)

# ---- Sample 50 trees and prune -----------------------------------------------
set.seed(42)
idx <- sample(100, 50)

trees300 <- lapply(idx, function(i) {
  ape::keep.tip(tree_bird_n100[[i]], our_spp)
})
class(trees300) <- "multiPhylo"

cat("Number of trees:", length(trees300), "\n")
cat("All have 300 tips:",
    all(sapply(trees300, function(tr) length(tr$tip.label)) == 300), "\n")

# Quick sanity: edge-length variation across trees
el_vec <- sapply(trees300, function(tr) sum(tr$edge.length))
cat("Total edge-length range:", round(range(el_vec)), "\n")
cat("Total edge-length SD:   ", round(sd(el_vec), 1), "\n")

# ---- Save --------------------------------------------------------------------
save(trees300, file = here::here("data", "trees300.rda"), compress = "xz")
f_kb <- round(file.size(here::here("data", "trees300.rda")) / 1024)
cat("Saved data/trees300.rda (", f_kb, "KB)\n")
