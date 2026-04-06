# comparing Rphylopars and impute_phylo

library(torch)
library(ape)
library(MASS)     # mvrnorm
library(Matrix)   # nearPD
library(Metrics)  # rmse
library(ggplot2)
library(here)
library(Rphylopars)

#library(Rphylopars)
library(Rphylopars)
library(phytools) # for simulating pure-birth phylogenies
set.seed(21) # Set the seed for reproducible results
trait_cov <- matrix(c(4,2,2.2,2,3,1.5,2.2,1.5,2.5),nrow = 3,ncol = 3)
tree <- pbtree(n = 1000)
sim_data <- simtraits(v = trait_cov,tree = tree, nmissing = 50)
p_BM <- phylopars(trait_data = sim_data$trait_data,tree = sim_data$tree)
# 
# sim_data$original_X
# p_BM$anc_recon[1:50,] # Data with imputed species means
# p_BM$anc_var[1:20,] # Variances for each estimate
# p_BM$anc_recon[1:20,] - sqrt(p_BM$anc_var[1:20,])*1.96 # Lower 95% CI
# p_BM$anc_recon[1:20,] + sqrt(p_BM$anc_var[1:20,])*1.96 # Upper 95% CI
# 
# # pigauot

# 1. Find the indices of the missing entries in your data
miss_idx <- which(is.na(sim_data$trait_data), arr.ind = TRUE)
miss_idx[, 2] <- miss_idx[, 2] - 1
# 2. Extract the “true” values (before masking) at those positions
#    sim_data$original_X is the complete matrix returned by simtraits()
true_missing <- sim_data$original_X[,,1]
true_missing <- true_missing[ miss_idx ]
# 3. Extract the phylopars-imputed means for those same species×trait cells
anc <- p_BM$anc_recon[1:1000, ]
imputed_means <- as.matrix(anc)[ miss_idx ]


env_data <- data.frame(rep(0, dim(as.matrix(sim_data$trait_data))[1]))
star_tree <- tree
star_tree$edge.length[] <- max(branching.times(tree)) # make all tips equidistant




res <- impute_phylo(
  trait_data   = as.data.frame(sim_data$trait_data[,2:4]),
  phylo_tree   = star_tree,
  env_data     = env_data,
  species_id   = star_tree$tip.label,
  latent_dim   = 128,
  epochs       = 6000,
  n_samples    = 1000,
  lr           = 1e-2,
  patience     = 6000,
  ckpt_path    = "checkpoints",
  n_hidden_layers = 1
)

res$imputed_missing$imputed
true_missing
imputed_means

# correlaiton between these values
cor(res$imputed_missing$imputed, true_missing)
cor(true_missing, imputed_means)


cat("RMSE on masked entries:", rmse(true_missing,res$imputed_missing$imputed), "\n")
cat("RMSE on masked entries:", rmse(true_missing, imputed_means), "\n")
