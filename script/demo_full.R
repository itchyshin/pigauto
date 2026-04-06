# ============================================================================
# pigauto full demo
# ============================================================================
# Run this script interactively (section by section) to:
#   1. Install and load the package
#   2. Browse help pages
#   3. Run the BM baseline pipeline
#   4. Run the TabPFN baseline (tabular foundation model, no phylogeny)
#   5. Train the pigauto GNN on top of each baseline
#   6. Compare all methods
#
# Requirements:
#   - R >= 4.1, torch, ape, ggplot2, Rphylopars
#   - For TabPFN: reticulate + Python with tabimpute (see Section 4)
# ============================================================================


# -- 1. Install & load -------------------------------------------------------

# Install pigauto from local source (run once):
# devtools::install("/Users/z3437171/Dropbox/Github Local/pigauto")

library(pigauto)
library(ape)
library(ggplot2)

cat("torch available:", torch::torch_is_installed(), "\n")
cat("CUDA available: ", torch::cuda_is_available(), "\n")
cat("MPS available:  ", torch::backends_mps_is_available(), "\n")


# -- 2. Browse help pages ----------------------------------------------------

# Core functions
?preprocess_traits
?make_missing_splits
?fit_baseline
?fit_baseline_tabpfn
?build_phylo_graph
?fit_pigauto
?evaluate_imputation
?setup_tabpfn

# Bundled data
?avonet300
?tree300


# -- 3. Data preparation -----------------------------------------------------

data(avonet300, tree300, package = "pigauto")

str(avonet300)
str(tree300)

# Prepare trait matrix (species as rownames)
traits <- avonet300
rownames(traits) <- traits$Species_Key
traits$Species_Key <- NULL

head(traits)

# Preprocess: align to tree, log-transform, z-score
pd <- preprocess_traits(traits, tree300, log_transform = TRUE)
print(pd)

# Create train/val/test splits (25% of cells held out)
set.seed(42)
splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.25, seed = 42)

cat("Training cells: ", sum(splits$mask), "\n")
cat("Validation cells:", length(splits$val_idx), "\n")
cat("Test cells:      ", length(splits$test_idx), "\n")


# -- 4. Baseline 1: Brownian Motion (Rphylopars) ----------------------------

cat("\n--- Fitting BM baseline (uses phylogenetic tree) ---\n")
bl_bm <- fit_baseline(pd, tree300, splits = splits)

eval_bm <- evaluate_imputation(bl_bm$mu, pd$X_scaled, splits, bl_bm$se)
cat("\nBM baseline results:\n")
print(eval_bm)


# -- 5. Baseline 2: TabPFN (no phylogeny) -----------------------------------
#
# First-time setup (run once -- creates a Python virtualenv):
#   setup_tabpfn()
#
# If you already have tabimpute installed in another environment:
#   reticulate::use_virtualenv("your-env-name")
#   or: reticulate::use_condaenv("your-conda-env")

tabpfn_ok <- requireNamespace("reticulate", quietly = TRUE) &&
  tryCatch({
    reticulate::import("tabimpute.interface")
    TRUE
  }, error = function(e) FALSE)

if (tabpfn_ok) {
  cat("\n--- Fitting TabPFN baseline (NO phylogenetic tree) ---\n")
  bl_tab <- fit_baseline_tabpfn(pd, splits = splits, envname = NULL)

  eval_tab <- evaluate_imputation(bl_tab$mu, pd$X_scaled, splits)
  cat("\nTabPFN baseline results:\n")
  print(eval_tab)
} else {
  cat("\n[SKIPPED] TabPFN not available.\n")
  cat("To install: pigauto::setup_tabpfn()\n")
  bl_tab <- NULL
}


# -- 6. Build phylogenetic graph ---------------------------------------------

cat("\n--- Building phylogenetic graph ---\n")
graph <- build_phylo_graph(tree300, k_eigen = 8)
cat("Adjacency:       ", dim(graph$adj)[1], "x", dim(graph$adj)[2], "\n")
cat("Spectral coords: ", dim(graph$coords)[1], "x", dim(graph$coords)[2], "\n")


# -- 7. pigauto: BM + GNN ---------------------------------------------------

cat("\n--- Training pigauto (BM baseline + GNN) ---\n")
fit_bm_gnn <- fit_pigauto(
  data    = pd,
  tree    = tree300,
  splits  = splits,
  graph   = graph,
  baseline = bl_bm,
  epochs  = 2000L,
  verbose = TRUE,
  seed    = 1
)
print(fit_bm_gnn)

# val/test RMSE from the fit object (computed during training)
cat("\npigauto (BM+GNN) val  RMSE:", fit_bm_gnn$val_rmse, "\n")
cat("pigauto (BM+GNN) test RMSE:", fit_bm_gnn$test_rmse, "\n")


# -- 8. pigauto: TabPFN + GNN (if available) --------------------------------

if (!is.null(bl_tab)) {
  cat("\n--- Training pigauto (TabPFN baseline + GNN) ---\n")
  fit_tab_gnn <- fit_pigauto(
    data     = pd,
    tree     = tree300,
    splits   = splits,
    graph    = graph,
    baseline = bl_tab,
    epochs   = 2000L,
    verbose  = TRUE,
    seed     = 1
  )
  print(fit_tab_gnn)

  cat("\npigauto (TabPFN+GNN) val  RMSE:", fit_tab_gnn$val_rmse, "\n")
  cat("pigauto (TabPFN+GNN) test RMSE:", fit_tab_gnn$test_rmse, "\n")
}


# -- 9. Summary comparison --------------------------------------------------

cat("\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("COMPARISON SUMMARY (test split RMSE, z-score scale)\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")

results <- data.frame(
  method    = "BM (Rphylopars)",
  test_rmse = mean(eval_bm$rmse[eval_bm$split == "test"]),
  uses_tree = TRUE,
  uses_gnn  = FALSE
)

if (!is.null(bl_tab)) {
  results <- rbind(results, data.frame(
    method    = "TabPFN (no tree)",
    test_rmse = mean(eval_tab$rmse[eval_tab$split == "test"]),
    uses_tree = FALSE,
    uses_gnn  = FALSE
  ))
}

results <- rbind(results, data.frame(
  method    = "pigauto (BM + GNN)",
  test_rmse = fit_bm_gnn$test_rmse,
  uses_tree = TRUE,
  uses_gnn  = TRUE
))

if (!is.null(bl_tab)) {
  results <- rbind(results, data.frame(
    method    = "pigauto (TabPFN + GNN)",
    test_rmse = fit_tab_gnn$test_rmse,
    uses_tree = TRUE,
    uses_gnn  = TRUE
  ))
}

print(results, row.names = FALSE)

cat("\nKey questions answered:\n")
cat("  1. Does phylogeny help?    Compare BM vs TabPFN\n")
cat("  2. Does the GNN help?      Compare BM vs pigauto(BM+GNN)\n")
cat("  3. Does tree + FM combine? Compare all four\n")


# -- 10. Visualise -----------------------------------------------------------

# Training history (BM + GNN)
plot(fit_bm_gnn, type = "history")

# Scatter: observed vs predicted (val split)
plot(fit_bm_gnn, type = "scatter", split = "test")

# Uncertainty ribbons (requires predict output with SE)
pred_bm_gnn <- predict(fit_bm_gnn, return_se = TRUE)
for (tr in pd$trait_names) {
  p <- plot_uncertainty(pred_bm_gnn, trait_name = tr)
  print(p)
}

cat("\nDone. See results above.\n")
