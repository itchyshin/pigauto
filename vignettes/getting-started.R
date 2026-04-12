## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 6,
  fig.height = 4
)


## ----install, eval=FALSE------------------------------------------------------
# install.packages("pigauto")
# 
# # torch backend (required; ~1 GB first-time download)
# torch::install_torch()


## ----load---------------------------------------------------------------------
library(pigauto)
library(ape)


## ----quickstart, eval=FALSE---------------------------------------------------
# result <- impute(traits, tree)
# 
# result$completed    # data.frame: observed values preserved, NAs filled
# result$imputed_mask # logical matrix: TRUE where a cell was imputed
# result$prediction$se # per-cell uncertainty (original units)


## ----data---------------------------------------------------------------------
data(avonet300, tree300, package = "pigauto")

head(avonet300)
ape::Ntip(tree300)


## ----preprocess---------------------------------------------------------------
traits <- avonet300[, -1]              # drop Species_Key column
rownames(traits) <- avonet300$Species_Key

pd <- preprocess_traits(traits, tree300, log_transform = TRUE)
print(pd)


## ----splits-------------------------------------------------------------------
splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.25, seed = 42)

cat("Total held-out cells:", length(splits$val_idx) + length(splits$test_idx), "\n")
cat("  Validation:", length(splits$val_idx), "\n")
cat("  Test:      ", length(splits$test_idx), "\n")


## ----baseline, eval=TRUE------------------------------------------------------
baseline <- fit_baseline(pd, tree300, splits = splits, model = "BM")

cat("Baseline object contains:\n")
cat("  mu matrix:", nrow(baseline$mu), "x", ncol(baseline$mu), "\n")
cat("  se matrix:", nrow(baseline$se), "x", ncol(baseline$se), "\n")


## ----graph, eval=TRUE---------------------------------------------------------
graph <- build_phylo_graph(
  tree300,
  k_eigen    = 8,
  sigma_mult = 0.5
)

cat("Graph: n =", graph$n, "species\n")
cat("Adjacency: [", nrow(graph$adj), "x", ncol(graph$adj), "]\n")
cat("Spectral coords: [", nrow(graph$coords), "x", ncol(graph$coords), "]\n")
cat("Kernel bandwidth sigma:", round(graph$sigma, 3), "\n")


## ----train, eval=FALSE--------------------------------------------------------
# fit <- fit_pigauto(
#   data            = pd,
#   tree            = tree300,
#   splits          = splits,
#   graph           = graph,
#   baseline        = baseline,
#   hidden_dim      = 64,
#   k_eigen         = 8,
#   dropout         = 0.10,
#   lr              = 3e-3,
#   epochs          = 2000,
#   corruption_rate = 0.55,
#   lambda_shrink   = 0.03,
#   eval_every      = 100,
#   patience        = 10,
#   verbose         = TRUE,
#   seed            = 1
# )


## ----print-fit, eval=FALSE----------------------------------------------------
# print(fit)


## ----predict, eval=FALSE------------------------------------------------------
# pred <- predict(fit, return_se = TRUE)
# 
# # pred$imputed: 300 x 4 matrix in original units
# # pred$se:      300 x 4 uncertainty matrix (original units)
# head(pred$imputed)


## ----evaluate, eval=FALSE-----------------------------------------------------
# # BM baseline RMSE on test cells
# eval_bm <- evaluate_imputation(baseline$mu, pd$X_scaled, splits)
# eval_bm[eval_bm$split == "test", c("trait", "n", "rmse", "pearson_r")]


## ----compare, eval=FALSE------------------------------------------------------
# # GNN test RMSE stored in fit object
# data.frame(
#   trait    = fit$trait_names,
#   bm_rmse  = eval_bm$rmse[eval_bm$split == "test"],
#   gnn_rmse = fit$test_rmse
# )


## ----plot-history, eval=FALSE-------------------------------------------------
# plot(fit, type = "history")


## ----plot-uncertainty, eval=FALSE---------------------------------------------
# plot_uncertainty(pred, trait_name = "Mass")


## ----gpu, eval=FALSE----------------------------------------------------------
# torch::cuda_is_available()       # NVIDIA GPU
# torch::backends_mps_is_available() # Apple Silicon GPU


## ----cache, eval=FALSE--------------------------------------------------------
# graph <- build_phylo_graph(tree300, k_eigen = 8, cache_path = "tree300_graph.rds")


## ----own-data, eval=FALSE-----------------------------------------------------
# tree   <- ape::read.tree("my_phylogeny.nwk")
# traits <- read_traits("my_traits.csv", species_col = "species")
# pd     <- preprocess_traits(traits, tree)
# # ...proceed as above

