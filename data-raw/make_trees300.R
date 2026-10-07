## data-raw/make_trees300.R
## Creates trees300 (multiPhylo): 50 posterior trees for the avonet300 species.
## Run from the package root: Rscript data-raw/make_trees300.R
##
## Requires: megatrees (>= 1.0.0), ape, pigauto (for tree300 tip names).
## Source: megatrees::get_tree_bird_n100(), which contains 50 Ericson-backbone
## and 50 Hackett-backbone BirdTree posterior trees (Jetz et al. 2012). This
## script samples 50 from the combined set, so trees300 contains both
## backbones. Original backbone citations:
## - Ericson PG et al. 2006. Diversification of Neoaves: integration of
##   molecular sequence data and fossils. Biology Letters 2:543-547.
##   doi:10.1098/rsbl.2006.0523.
## - Hackett SJ et al. 2008. A Phylogenomic Study of Birds Reveals Their
##   Evolutionary History. Science 320:1763-1768. doi:10.1126/science.1157704.
## BirdTree posterior-tree source: Jetz W et al. 2012. The global diversity
## of birds in space and time. Nature 491:444-448. doi:10.1038/nature11631.
## The megatrees package records an MIT license; redistribution rights for
## the underlying BirdTree data in pigauto remain open.

library(ape)
library(megatrees)

# ---- Get the 300 species names ------------------------------------------------
# Load our bundled MCC tree to get the canonical species set
tree300 <- get(load(here::here("data", "tree300.rda")))
our_spp <- tree300$tip.label
stopifnot(length(our_spp) == 300L)

# ---- Verify all 300 species exist in the megatree ----------------------------
tree_bird_n100 <- megatrees::get_tree_bird_n100()
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
