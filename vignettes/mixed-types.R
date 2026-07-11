## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = torch::torch_is_installed()
)


## ----simulate-----------------------------------------------------------------
library(pigauto)
library(ape)

set.seed(42)
n <- 60
tree <- rtree(n)

traits <- data.frame(
  row.names = tree$tip.label,
  mass      = exp(rnorm(n, 3, 0.5)),
  clutch    = as.integer(rpois(n, 3) + 1L),
  migr      = factor(sample(c("no", "yes"), n, replace = TRUE)),
  diet      = factor(sample(c("herb", "carn", "omni"), n, replace = TRUE)),
  threat    = ordered(sample(c("LC", "VU", "EN"), n, replace = TRUE),
                      levels = c("LC", "VU", "EN"))
)


## ----preprocess---------------------------------------------------------------
pd <- preprocess_traits(traits, tree)
print(pd)


## ----trait_map----------------------------------------------------------------
str(pd$trait_map, max.level = 1)
pd$trait_map$diet


## ----splits-------------------------------------------------------------------
spl <- make_missing_splits(pd$X_scaled, missing_frac = 0.20,
                           seed = 1, trait_map = pd$trait_map)
cat("Val cells (latent):", length(spl$val_idx), "\n")
cat("Test cells (latent):", length(spl$test_idx), "\n")


## ----baseline-----------------------------------------------------------------
bl <- fit_baseline(pd, tree, splits = spl)
dim(bl$mu)


## ----train, message = FALSE---------------------------------------------------
fit <- fit_pigauto(
  pd, tree,
  splits = spl,
  epochs = 200L,
  eval_every = 50L,
  patience = 5L,
  verbose = FALSE,
  seed = 1
)
print(fit)


## ----predict------------------------------------------------------------------
pred <- predict(fit, return_se = TRUE)
head(pred$imputed)


## ----probs--------------------------------------------------------------------
# Binary: probability of "yes"
head(pred$probabilities$migr)

# Categorical: probability of each diet class
head(pred$probabilities$diet)


## ----se-----------------------------------------------------------------------
head(pred$se)


## ----evaluate-----------------------------------------------------------------
ev <- evaluate_imputation(pred, pd$X_scaled, spl)
print(ev)


## ----mi_workflow, eval=FALSE--------------------------------------------------
# # Precompute nonlinear auxiliaries explicitly when required.
# analysis_data$z_sq <- analysis_data$z^2
# mi <- multi_impute_analysis(
#   data = analysis_data, formula = y ~ x + z,
#   missing = "x", model = "glm", m = 50L,
#   auxiliary = "z_sq", seed = 1L
# )
#
# # Fit a downstream model to each
# fits <- with_imputations(mi, function(d) {
#   glm(y ~ x + z, data = d, family = binomial)
# })
#
# # Pool with Rubin's rules
# pool_mi(fits)
