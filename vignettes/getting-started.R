## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 6,
  fig.height = 4
)


## ----canonical-journey, eval=FALSE--------------------------------------------
# library(pigauto)
#
# traits <- read_traits("traits.csv")
# tree <- read_tree("tree.nwk")
# check_pigauto(traits, tree)
# result <- impute(traits, tree)
# completed <- completed_data(result)
# pigauto_report(result)


## ----declarations, eval=FALSE-------------------------------------------------
# # Ordered ecological states
# traits$threat <- ordered(traits$threat, levels = c("LC", "NT", "VU", "EN", "CR"))
#
# # Integer-valued continuous measurement; integer otherwise means count
# result <- impute(traits, tree, trait_types = c(length_mm = "continuous"))
#
# # Numeric 0-1 proportion and zero-inflated count require explicit declarations
# result <- impute(traits, tree, trait_types = c(cover = "proportion", parasites = "zi_count"))
#
# # A composition is complete by row and its observed components sum to one
# result <- impute(traits, tree,
#   multi_proportion_groups = list(diet = c("diet_insect", "diet_fruit", "diet_seed")))


## ----install, eval=FALSE------------------------------------------------------
# # CRAN release
# install.packages("pigauto")
#
# # Development version
# pak::pak("itchyshin/pigauto")
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
#
# # Conformal prediction intervals (nominal held-out diagnostic)
# pred$conformal_lower[["Mass"]]    # lower bound per species (original units)
# pred$conformal_upper[["Mass"]]    # upper bound per species (original units)
# pred$conformal_coverage           # empirical coverage on val set (target ≈ 0.95)


## ----evaluate, eval=FALSE-----------------------------------------------------
# evaluate_imputation(pred$imputed_latent, pd$X_scaled, splits)
# fit$r_cal


## ----plot-history, eval=FALSE-------------------------------------------------
# plot(fit, type = "history")


## ----plot-uncertainty, eval=FALSE---------------------------------------------
# plot(pred, type = "intervals", trait = "Mass")


## ----mi-workflow, eval=FALSE--------------------------------------------------
# analysis_data$z_sq <- analysis_data$z^2
#
# mi <- multi_impute_analysis(
#   data = analysis_data,
#   formula = y ~ x + z,
#   missing = "x",
#   model = "lm",
#   m = 50L,
#   auxiliary = "z_sq",
#   seed = 1L
# )
#
# fits <- with_imputations(mi, function(d) lm(y ~ x + z, data = d))
# pool_mi(fits)


## ----active-impute, eval=FALSE------------------------------------------------
# res <- impute(traits, tree)
#
# # Top-10 individual (species, trait) cells, ordered by expected
# # total variance reduction:
# suggest_next_observation(res, top_n = 10, by = "cell")
#
# # Or aggregated by species (sum of variance reductions across the
# # species' currently-missing continuous-family traits):
# suggest_next_observation(res, top_n = 10, by = "species")


## ----advanced-controls, eval=FALSE--------------------------------------------
# result <- impute(
#   traits, tree,
#   joint_solver = "rphylopars",
#   em_iterations = 2L,
#   joint_refine_iter = 1L,
#   conformal_split_val = TRUE
# )


## ----exact-route, eval=FALSE--------------------------------------------------
# result <- impute(traits, tree, predict_method = "exact")


## ----phylo-signal-gate, eval=FALSE--------------------------------------------
# result <- impute(traits, tree,
#                  phylo_signal_gate     = TRUE,   # default
#                  phylo_signal_threshold = 0.2,   # default; min lambda to keep BM/GNN
#                  phylo_signal_method   = "lambda")


## ----phylo-signal-gate-off, eval=FALSE----------------------------------------
# result <- impute(traits, tree, phylo_signal_gate = FALSE)


## ----covariate-helpers, eval=FALSE--------------------------------------------
# # Step 1: GBIF occurrence centroids (median lat/lon across cleaned points)
# gbif <- pull_gbif_centroids(species          = rownames(my_traits),
#                             cache_dir        = "cache/gbif",
#                             occurrence_limit = 500L)
# #> data.frame with columns species, centroid_lat, centroid_lon, n_occurrences
#
# # Step 2: WorldClim bioclim summaries (median + IQR of all 19 bio variables
# # across each species' GBIF occurrences)
# clim <- pull_worldclim_per_species(species             = rownames(my_traits),
#                                     gbif_cache_dir      = "cache/gbif",
#                                     worldclim_cache_dir = "cache/worldclim",
#                                     resolution          = "10m")
# #> data.frame with bio1..bio19 medians + IQRs + n_extracted per species
#
# # Step 3: hand the whole climate frame to impute() as covariates
# result <- impute(my_traits, my_tree, covariates = clim)


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
