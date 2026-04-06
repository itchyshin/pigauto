###############################################################################
## Simulate phylogenetically-correlated continuous *and* categorical traits
## plus environmental covariates, then impose 25% MCAR missingness
###############################################################################

library(ape)        # rtree()
library(MASS)       # mvrnorm()
library(Matrix)     # nearPD()

set.seed(123)

## 1. SETTINGS
n_sp   <- 300          # number of species
p_cts  <-  5           # continuous traits
p_bin  <-  3           # binary traits
p_ord  <-  1           # ordinal traits
p_mult <-  1           # multinomial traits (3-level)
p_env  <-  6           # environmental predictors

## 2. PHYLOGENY
tree   <- rtree(n_sp)
tree <- compute.brlen(tree) # ensure all branches are positive
V_phy  <- cov2cor(vcv(tree))

#library(ape)
# 
# make_star_phylogeny <- function(n_tips, tip_length = 1) {
#   # Create a phylo object for a star tree (all branches from the root)
#   tree <- list()
#   tree$edge <- cbind(rep(n_tips + 1, n_tips), 1:n_tips)
#   tree$Nnode <- 1
#   tree$tip.label <- paste0("sp", seq_len(n_tips))
#   tree$edge.length <- rep(tip_length, n_tips)
#   class(tree) <- "phylo"
#   tree
# }
# 
# # Example:
# tree2 <- make_star_phylogeny(n_tips = 300, tip_length = 1)

# Plot to confirm it's a star:
#plot(star_tree, show.tip.label = FALSE, main="Star Phylogeny")

# Example:
#star_tree <- make_star_phylogeny(n_tips = 300, tip_length = 1)

## 3. ENVIRONMENT
Sigma_env <- diag(p_env)
env_data  <- mvrnorm(n_sp, rep(0, p_env), Sigma_env)
colnames(env_data) <- paste0("env", 1:p_env)

## 4A. CONTINUOUS TRAITS
lambda <- rep(0.9, 5) # runif(p_cts, 0.2, 0.8)
Sigma_list <- lapply(lambda, function(lam) lam * V_phy + (1 - lam) * diag(n_sp))
Z_phy <- sapply(Sigma_list, function(Sig) mvrnorm(n = 1, mu = rep(0, n_sp), Sigma = Sig))
B_cts <- matrix(runif(p_cts * p_env, -1, 1), p_cts, p_env)
set.seed(42)
cts_traits <- Z_phy +
  as.matrix(env_data) %*% t(B_cts*0) + #  t(B_cts +
  matrix(rnorm(n_sp * p_cts, 0, 0.1), n_sp, p_cts)
colnames(cts_traits) <- paste0("cnt", 1:p_cts)

## 4B. BINARY TRAITS
lambda_bin <- rep(0.99, p_bin) #runif(p_bin, 0.2, 0.8)
Sigma_bin_list <- lapply(lambda_bin, function(lam) lam * V_phy + (1 - lam) * diag(n_sp))
Z_phy_bin <- sapply(Sigma_bin_list, function(Sig) mvrnorm(n = 1, mu = rep(0, n_sp), Sigma = Sig))
B_bin <- matrix(runif(p_bin * p_env, -1, 1), p_bin, p_env)
linpred_bin <- Z_phy_bin + as.matrix(env_data) %*% t(B_bin*2)
prob_bin <- plogis(linpred_bin)
bin_traits <- matrix(rbinom(n = n_sp * p_bin, size = 1, prob = as.vector(prob_bin)), nrow = n_sp, ncol = p_bin)
colnames(bin_traits) <- paste0("bin", 1:p_bin)

## 4C. ORDINAL TRAITS
lambda_ord <- runif(p_ord, 0.2, 0.8)
Sigma_ord_list <- lapply(lambda_ord, function(lam) lam * V_phy + (1 - lam) * diag(n_sp))
Z_phy_ord <- sapply(Sigma_ord_list, function(Sig) mvrnorm(n = 1, mu = rep(0, n_sp), Sigma = Sig))
B_ord <- matrix(runif(p_ord * p_env, -1, 1), p_ord, p_env)
set.seed(42)
resid_ord <- matrix(rnorm(n_sp * p_ord, 0, 0.1), n_sp, p_ord)
linpred_ord <- Z_phy_ord + as.matrix(env_data) %*% t(B_ord) + resid_ord
ord_traits <- apply(linpred_ord, 2, function(latent) {
  thr <- quantile(latent, probs = c(1/3, 2/3))
  cut(
    latent,
    breaks = c(-Inf, thr, Inf),
    labels = c("low", "mid", "high"),
    ordered_result = TRUE
  )
})

ord_traits <- factor(ord_traits, levels = c("low", "mid", "high"), labels = c("low", "mid", "high"), ordered = TRUE)
ord_traits <- as.data.frame(ord_traits)
colnames(ord_traits) <- paste0("ord", seq_len(p_ord))

## 4D. MULTINOMIAL TRAIT (3-level)
n_mult   <- p_mult
n_levels <- 3
phylo_cor <- 0.5
lambda_mult <- c(0.5, 0.7)
V1 <- lambda_mult[1] * V_phy + (1 - lambda_mult[1]) * diag(n_sp)
V2 <- lambda_mult[2] * V_phy + (1 - lambda_mult[2]) * diag(n_sp)
Cov12 <- phylo_cor * sqrt(lambda_mult[1] * lambda_mult[2]) * V_phy
top    <- cbind(V1, Cov12)
bottom <- cbind(t(Cov12), V2)
Sigma_mult <- rbind(top, bottom)
Z_mult_latent <- matrix(mvrnorm(1, mu = rep(0, 2 * n_sp), Sigma = Sigma_mult), nrow = n_sp, ncol = 2)
B_mult <- matrix(runif((n_levels - 1) * p_env, -1, 1), n_levels - 1, p_env)
linpred_mult <- Z_mult_latent + as.matrix(env_data) %*% t(B_mult) + matrix(rnorm(n_sp * (n_levels - 1), 0, 0.1), n_sp, n_levels - 1)
softmax <- function(x) exp(x) / sum(exp(x))
logit_mat <- cbind(0, linpred_mult)
probs <- t(apply(logit_mat, 1, softmax))
mult_trait <- apply(probs, 1, function(p) sample(1:3, 1, prob = p))
mult_traits <- factor(mult_trait, levels = 1:3, labels = c("A", "B", "C"))
mult_traits <- as.data.frame(mult_traits)
colnames(mult_traits) <- paste0("mult", seq_len(n_mult))

## 5. COMBINE TRAITS INTO DATA FRAME
traits_df <- data.frame(
  cts_traits,
  bin_traits,
  ord_traits,
  mult_traits
)

## 6. IMPOSING 25% MCAR MISSINGNESS — SUPPORTS ALL TYPES
miss_mat <- matrix(runif(n_sp * ncol(traits_df)) < 0.50, n_sp, ncol(traits_df))
traits_df_miss <- traits_df # copy

for (j in seq_along(traits_df)) {
  idx <- which(miss_mat[,j])
  if (is.numeric(traits_df[[j]]) || is.integer(traits_df[[j]])) {
    traits_df_miss[idx, j] <- NA
  } else if (is.factor(traits_df[[j]])) {
    # factors, incl. ordered and multinomial
    traits_df_miss[[j]][idx] <- NA
  }
}

## 7. READY FOR VGAE PIPELINE OR EXPLORATION
env_df    <- as.data.frame(env_data)

str(traits_df_miss[1:5, ])
str(env_df[1:5, ])

# Now traits_df_miss contains 25% MCAR missing data for all trait types!

#' Classify trait columns by variable type
#'
#' @param df A data.frame containing trait columns.
#' @param binary_threshold Integer: max unique values for binary (default: 2)
#' @return A list with elements: continuous, binary, ordinal, multinomial
#' @examples
#' classify_traits(traits_df_miss)
classify_traits <- function(df, binary_threshold = 2) {
  res <- list(continuous = character(0),
              binary     = character(0),
              ordinal    = character(0),
              multinomial= character(0))
  
  for (col in names(df)) {
    x <- df[[col]]
    # Ignore completely missing columns
    if (all(is.na(x))) next
    
    # Ordinal: ordered factor
    if (is.ordered(x)) {
      res$ordinal <- c(res$ordinal, col)
    } 
    # Multinomial: unordered factor (not binary)
    else if (is.factor(x)) {
      n_lev <- nlevels(x)
      if (n_lev > binary_threshold) {
        res$multinomial <- c(res$multinomial, col)
      } else {
        res$binary <- c(res$binary, col) # fallback: factor with 2 levels = binary
      }
    }
    # Binary: integer/numeric, only 2 unique non-missing values (e.g., 0,1)
    else if (is.integer(x) || is.numeric(x)) {
      uniq <- unique(x[!is.na(x)])
      if (length(uniq) == binary_threshold &&
          all(uniq %in% c(0, 1))) {
        res$binary <- c(res$binary, col)
      } else {
        res$continuous <- c(res$continuous, col)
      }
    } 
    # Fallback: treat as continuous
    else {
      res$continuous <- c(res$continuous, col)
    }
  }
  res
}

classify_traits(traits_df_miss)
