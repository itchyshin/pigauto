#!/usr/bin/env Rscript
#
# data-raw/make_ctmax_sim.R
#
# Generate a simulated multi-observation-per-species dataset for
# demonstrating pigauto's observation-level covariate imputation.
#
# The dataset mimics a critical thermal maximum (CTmax) study where
# each species has multiple measurements taken at different
# acclimation temperatures. The data-generating process is:
#
#   CTmax_ij = mu + phylo_i + beta * acclim_temp_ij + epsilon_ij
#
# where phylo_i ~ BM(tree), acclim_temp_ij varies within species,
# and beta is the acclimation response ratio (ARR).
#
# Output: data/ctmax_sim.rda
#         (uses tree300 as the underlying phylogeny)

library(ape)

set.seed(2026)

# Load the tree already bundled with pigauto
load("data/tree300.rda")
tree <- tree300
species <- tree$tip.label
n_species <- length(species)

# ---- Parameters -----------------------------------------------------------
mu    <- 38.0          # global mean CTmax (degrees C)
sigma <- 3.0           # BM process SD
beta  <- 0.10          # ARR: CTmax increases by 0.10 per degree acclimation
sigma_eps <- 1.5       # residual (within-species) SD

# ---- Phylogenetic species effects (BM) -----------------------------------
# Scale tree to get desired phylogenetic variance
tree_sc <- tree
tree_sc$edge.length <- tree_sc$edge.length / max(vcv(tree_sc)) * sigma^2
phylo_vals <- rTraitCont(tree_sc, model = "BM", sigma = 1.0)

# Ensure phylo_vals are aligned with species
phylo_vals <- phylo_vals[species]

# ---- Generate multi-obs data ---------------------------------------------
# Variable number of observations per species: Poisson(4) + 1, capped at 15
n_obs_per <- pmin(rpois(n_species, 4) + 1L, 15L)
total_obs <- sum(n_obs_per)

# Build the data frame
sp_col    <- rep(species, times = n_obs_per)
phylo_col <- rep(phylo_vals, times = n_obs_per)

# Acclimation temperature: 10-35 degrees, with within-species variation
acclim_temp <- numeric(total_obs)
idx <- 1L
for (i in seq_len(n_species)) {
  n_i <- n_obs_per[i]
  # Species-level baseline acclimation drawn from uniform [15, 30]
  base_temp <- runif(1, 15, 30)
  # Within-species variation: ±10 degrees around the baseline
  acclim_temp[idx:(idx + n_i - 1L)] <- base_temp + runif(n_i, -10, 10)
  idx <- idx + n_i
}

# CTmax: mu + phylo + beta * acclim_temp + epsilon
epsilon <- rnorm(total_obs, 0, sigma_eps)
ctmax   <- mu + phylo_col + beta * acclim_temp + epsilon

# ---- Introduce realistic missingness -------------------------------------
# 30% of species have NO observations (completely unobserved)
# 15% of remaining observations are missing at random
unobserved_sp <- sample(species, round(0.3 * n_species))
sp_mask <- sp_col %in% unobserved_sp

ctmax_miss <- ctmax
ctmax_miss[sp_mask] <- NA

# Additional MCAR missingness in observed species
obs_idx <- which(!sp_mask)
n_mcar  <- round(0.15 * length(obs_idx))
mcar_idx <- sample(obs_idx, n_mcar)
ctmax_miss[mcar_idx] <- NA

# ---- Assemble data frame -------------------------------------------------
ctmax_sim <- data.frame(
  species     = sp_col,
  acclim_temp = round(acclim_temp, 1),
  CTmax       = round(ctmax_miss, 2),
  stringsAsFactors = FALSE
)

# Store the truth for benchmarking (not in the main dataset)
attr(ctmax_sim, "ctmax_truth") <- round(ctmax, 2)

# ---- Summary -------------------------------------------------------------
cat(sprintf("Species:        %d\n", n_species))
cat(sprintf("Total obs:      %d\n", total_obs))
cat(sprintf("Obs/species:    %.1f (range %d-%d)\n",
            mean(n_obs_per), min(n_obs_per), max(n_obs_per)))
cat(sprintf("Unobserved spp: %d (%.0f%%)\n",
            length(unobserved_sp), 100 * length(unobserved_sp) / n_species))
cat(sprintf("Total NA cells: %d / %d (%.1f%%)\n",
            sum(is.na(ctmax_sim$CTmax)), nrow(ctmax_sim),
            100 * sum(is.na(ctmax_sim$CTmax)) / nrow(ctmax_sim)))
cat(sprintf("DGP: CTmax = %.1f + phylo + %.2f * acclim_temp + eps (sd=%.1f)\n",
            mu, beta, sigma_eps))

# ---- Save -----------------------------------------------------------------
save(ctmax_sim, file = "data/ctmax_sim.rda", compress = "xz")
cat("Wrote data/ctmax_sim.rda\n")
