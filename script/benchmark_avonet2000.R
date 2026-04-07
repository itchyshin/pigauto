# ============================================================================
# Benchmark: pigauto vs Rphylopars on AVONET 2000-species dataset
# ============================================================================
#
# Real-data benchmark using 2000 bird species, 8 continuous morphometric
# traits, and a pre-existing held-out truth set (~1594 masked cells).
#
# Data files (from BACE/dev/testing_data/):
#   - avonet_2000_masked.csv : 2000 species x 8 traits with introduced NAs
#   - avonet_2000_truth.csv  : held-out (species_tip, trait, true_value)
#   - Hackett_tree_2000.tre  : matching Hackett phylogeny (2000 tips)
#
# Methods compared:
#   1. BM (Rphylopars via pigauto)  -- phylogenetic BM baseline only
#   2. pigauto (BM + GNN)           -- residual GNN on top of BM baseline
#   3. Rphylopars (standalone)      -- direct Rphylopars on log-scale data
#   4. TabPFN (optional)            -- tabular foundation model (no phylogeny)
#
# All evaluation is on the original (un-logged, un-scaled) scale, matching
# the truth CSV values.
#
# Requirements: pigauto, ape, Rphylopars, ggplot2, torch
# ============================================================================


# -- 0. Configuration and paths -----------------------------------------------

# Path to the BACE testing data directory (relative to pigauto root)
data_dir <- file.path(
  dirname(getwd()),  # parent of pigauto

  "BACE", "dev", "testing_data"
)

# If running from a different location, set data_dir manually:
# data_dir <- "/Users/z3437171/Dropbox/Github Local/BACE/dev/testing_data"

# Output directory (inside pigauto/script/)
out_dir <- file.path(getwd(), "script")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("==========================================================\n")
cat("pigauto vs Rphylopars benchmark: AVONET 2000 species\n")
cat("==========================================================\n")
cat("Data directory:", data_dir, "\n")
cat("Output directory:", out_dir, "\n\n")


# -- 1. Load packages and data ------------------------------------------------

library(pigauto)
library(ape)
library(Rphylopars)
library(ggplot2)

cat("torch available:", torch::torch_is_installed(), "\n")
cat("CUDA available: ", torch::cuda_is_available(), "\n")
cat("MPS available:  ", torch::backends_mps_is_available(), "\n")

# Check TabPFN availability
TABPFN_OK <- tryCatch({
  requireNamespace("reticulate", quietly = TRUE) &&
    { reticulate::import("tabimpute.interface"); TRUE }
}, error = function(e) FALSE)
cat("TabPFN available:", TABPFN_OK, "\n\n")

# 1a. Masked trait data
masked_path <- file.path(data_dir, "avonet_2000_masked.csv")
stopifnot(file.exists(masked_path))
masked_df <- read.csv(masked_path, stringsAsFactors = FALSE)
cat("Masked data loaded:", nrow(masked_df), "species x",
    ncol(masked_df) - 1, "traits\n")

# 1b. Truth (held-out values)
truth_path <- file.path(data_dir, "avonet_2000_truth.csv")
stopifnot(file.exists(truth_path))
truth <- read.csv(truth_path, stringsAsFactors = FALSE)
cat("Truth data loaded:", nrow(truth), "held-out entries\n")
cat("  Traits in truth:", paste(sort(unique(truth$trait)), collapse = ", "), "\n")

# 1c. Phylogenetic tree
tree_path <- file.path(data_dir, "Hackett_tree_2000.tre")
stopifnot(file.exists(tree_path))
tree <- read_tree(tree_path)

# Compute branch lengths if absent or all-zero
if (is.null(tree$edge.length) || all(tree$edge.length == 0)) {
  cat("Computing Grafen branch lengths...\n")
  tree <- ape::compute.brlen(tree, method = "Grafen")
} else {
  cat("Tree already has branch lengths.\n")
}

cat("Tree:", length(tree$tip.label), "tips\n\n")

# Verify species alignment
stopifnot(all(masked_df$species %in% tree$tip.label))
stopifnot(length(tree$tip.label) == nrow(masked_df))


# -- 2. Prepare for pigauto ---------------------------------------------------

cat("--- Preparing data for pigauto ---\n")

# Set species as rownames, drop species column
traits_df <- masked_df
rownames(traits_df) <- traits_df$species
traits_df$species <- NULL

# Ensure all columns are numeric
for (nm in names(traits_df)) {
  traits_df[[nm]] <- as.numeric(traits_df[[nm]])
}

trait_names <- names(traits_df)
cat("Trait columns:", paste(trait_names, collapse = ", "), "\n")

# Preprocess: align to tree, auto-detect log columns, z-score
# log_transform = TRUE will auto-log columns with all-positive observed values.
# For AVONET: Mass, Wing.Length, Beak.Length_Culmen, Tarsus.Length,
#   Tail.Length, Range.Size will be logged (all positive).
# Centroid.Latitude and Centroid.Longitude will NOT be logged (negative values).
pd <- preprocess_traits(traits_df, tree, log_transform = TRUE)
print(pd)

# Report which traits were log-transformed
for (tm in pd$trait_map) {
  if (isTRUE(tm$log_transform)) {
    cat("  Log-transformed:", tm$name, "\n")
  } else {
    cat("  Not logged:     ", tm$name, "\n")
  }
}
cat("\n")


# -- 3. Build custom splits from pre-existing NA pattern + truth CSV -----------

cat("--- Building val/test splits from truth CSV ---\n")

# The avonet_2000_masked.csv already has NAs from a prior masking step.
# The truth CSV tells us which (species, trait) combinations have known values.
# We split those truth entries into val (25%) and test (75%).
# Other NAs (not in truth) are truly missing -- no evaluation is possible.

n_sp <- nrow(pd$X_scaled)
p_lat <- ncol(pd$X_scaled)

# Map each truth entry to its linear index in X_scaled
# X_scaled rows are in tree$tip.label order (from preprocess_traits)
# X_scaled columns correspond to latent columns via trait_map
species_order <- pd$species_names  # == tree$tip.label

# Build a lookup: trait name -> latent column index (1 col per continuous trait)
trait_to_latent_col <- setNames(integer(length(pd$trait_map)), names(pd$trait_map))
for (tm in pd$trait_map) {
  trait_to_latent_col[tm$name] <- tm$latent_cols[1]
}

# For each truth row, find its position in X_scaled
truth_linear_idx <- integer(nrow(truth))
truth_matched <- logical(nrow(truth))

for (i in seq_len(nrow(truth))) {
  sp <- truth$species_tip[i]
  tr <- truth$trait[i]

  row_idx <- match(sp, species_order)
  col_idx <- trait_to_latent_col[tr]

  if (!is.na(row_idx) && !is.na(col_idx)) {
    # Linear index: column-major (R default)
    truth_linear_idx[i] <- row_idx + (col_idx - 1L) * n_sp
    truth_matched[i] <- TRUE
  }
}

n_matched <- sum(truth_matched)
cat("Truth entries matched to X_scaled:", n_matched, "of", nrow(truth), "\n")

if (n_matched < nrow(truth)) {
  n_unmatched <- nrow(truth) - n_matched
  cat("  WARNING:", n_unmatched, "truth entries could not be matched.\n")
  cat("  Unmatched species:",
      head(truth$species_tip[!truth_matched], 5), "...\n")
}

# Keep only matched entries
truth_idx <- truth_linear_idx[truth_matched]
truth_sub <- truth[truth_matched, ]

# Split into val (25%) and test (75%)
set.seed(42)
n_truth <- length(truth_idx)
val_size <- floor(0.25 * n_truth)
perm <- sample(n_truth)
val_idx  <- truth_idx[perm[seq_len(val_size)]]
test_idx <- truth_idx[perm[(val_size + 1L):n_truth]]

# Also record which truth rows are val vs test (for later evaluation)
truth_sub$split <- NA_character_
truth_sub$split[perm[seq_len(val_size)]] <- "val"
truth_sub$split[perm[(val_size + 1L):n_truth]] <- "test"

# Build the mask: TRUE where values are observed (not NA) and not held out
mask <- !is.na(pd$X_scaled)
mask[c(val_idx, test_idx)] <- FALSE

splits <- list(
  val_idx  = val_idx,
  test_idx = test_idx,
  n        = n_sp,
  p        = p_lat,
  n_traits = length(pd$trait_map),
  mask     = mask
)

cat("Validation cells:", length(val_idx), "\n")
cat("Test cells:      ", length(test_idx), "\n")
cat("Observed (train):", sum(mask), "of", n_sp * p_lat, "total cells\n\n")


# -- 4. Run pigauto pipeline ---------------------------------------------------

cat("==========================================================\n")
cat("SECTION 4: pigauto pipeline\n")
cat("==========================================================\n\n")

# 4a. BM baseline
cat("--- Fitting BM baseline via Rphylopars ---\n")
bl <- fit_baseline(pd, tree, splits = splits)
cat("BM baseline fitted. mu:", dim(bl$mu)[1], "x", dim(bl$mu)[2], "\n\n")

# 4b. Build phylogenetic graph
cat("--- Building phylogenetic graph (2000 tips) ---\n")
graph <- build_phylo_graph(tree, k_eigen = 8)
cat("Graph: adj", dim(graph$adj)[1], "x", dim(graph$adj)[2],
    "| coords", dim(graph$coords)[1], "x", dim(graph$coords)[2], "\n\n")

# 4c. Train GNN
cat("--- Training pigauto GNN (epochs = 2000) ---\n")
fit <- tryCatch({
  fit_pigauto(
    data    = pd,
    tree    = tree,
    splits  = splits,
    graph   = graph,
    baseline = bl,
    epochs  = 2000L,
    verbose = TRUE,
    seed    = 42
  )
}, error = function(e) {
  cat("\n*** GNN training failed ***\n")
  cat("Error:", conditionMessage(e), "\n")
  cat("Continuing with BM baseline only.\n\n")
  NULL
})

# 4d. Predict with MC dropout (if training succeeded)
pred_pigauto <- NULL
if (!is.null(fit)) {
  cat("\n--- Generating pigauto predictions (MC dropout, 5 imputations) ---\n")
  pred_pigauto <- predict(fit, return_se = TRUE, n_imputations = 5L)
  cat("Predictions generated.\n\n")
}


# -- 5. Run standalone Rphylopars (comparison baseline) ------------------------

cat("==========================================================\n")
cat("SECTION 5: Standalone Rphylopars\n")
cat("==========================================================\n\n")

# Rphylopars operates on log-scale data, matching what pigauto does internally.
# We log-transform the same traits that pigauto auto-detected.

# Identify which traits were log-transformed
log_set <- character(0)
for (tm in pd$trait_map) {
  if (isTRUE(tm$log_transform)) log_set <- c(log_set, tm$name)
}

cat("Traits to log-transform for Rphylopars:", paste(log_set, collapse = ", "), "\n")

# Prepare data frame: species column + traits (already masked)
rp_df <- masked_df
rp_df <- rp_df[match(tree$tip.label, rp_df$species), ]  # align to tree order

# Log-transform the same traits
for (nm in log_set) {
  rp_df[[nm]] <- log(rp_df[[nm]])
}

# Rphylopars requires a "species" column as the first column
rp_input <- data.frame(
  species = rp_df$species,
  rp_df[, trait_names, drop = FALSE],
  stringsAsFactors = FALSE
)

cat("Running Rphylopars (BM, no phenotypic error)...\n")
rp_time <- system.time({
  rp_fit <- tryCatch(
    suppressWarnings(
      Rphylopars::phylopars(
        trait_data  = rp_input,
        tree        = tree,
        model       = "BM",
        pheno_error = FALSE
      )
    ),
    error = function(e) {
      cat("  Rphylopars failed:", conditionMessage(e), "\n")
      cat("  Retrying with branch-length jitter...\n")
      tree_j <- tree
      tree_j$edge.length <- tree_j$edge.length + 1e-6
      suppressWarnings(
        Rphylopars::phylopars(
          trait_data  = rp_input,
          tree        = tree_j,
          model       = "BM",
          pheno_error = FALSE
        )
      )
    }
  )
})
cat("Rphylopars completed in", round(rp_time["elapsed"], 1), "seconds.\n")

# Extract tip reconstructions (log scale)
tip_labels <- rp_fit$tree$tip.label
rp_recon_log <- rp_fit$anc_recon[tip_labels, trait_names, drop = FALSE]
rp_var_log   <- rp_fit$anc_var[tip_labels, trait_names, drop = FALSE]
rp_se_log    <- sqrt(pmax(rp_var_log, 0))

# Back-transform to original scale
rp_recon_orig <- rp_recon_log
for (nm in log_set) {
  rp_recon_orig[, nm] <- exp(rp_recon_log[, nm])
}
# Non-logged traits remain as-is (they were not transformed)

rp_recon_df <- as.data.frame(rp_recon_orig)
rp_recon_df$species <- rownames(rp_recon_df)
cat("Rphylopars predictions extracted for", nrow(rp_recon_df), "species.\n\n")


# -- 6. Back-transform pigauto predictions and merge with truth ----------------

cat("==========================================================\n")
cat("SECTION 6: Merge predictions with truth\n")
cat("==========================================================\n\n")

# Prepare comparison data frame
# Start with truth_sub (matched truth entries with val/test split assignment)
compare_df <- truth_sub

# 6a. BM baseline predictions (back-transform from latent scale)
bm_pred_original <- numeric(nrow(compare_df))
for (i in seq_len(nrow(compare_df))) {
  sp <- compare_df$species_tip[i]
  tr <- compare_df$trait[i]
  tm <- pd$trait_map[[tr]]
  row_idx <- match(sp, species_order)
  col_idx <- tm$latent_cols[1]

  # Reverse z-score
  val_log <- bl$mu[row_idx, col_idx] * tm$sd + tm$mean

  # Reverse log if applicable
  if (isTRUE(tm$log_transform)) {
    bm_pred_original[i] <- exp(val_log)
  } else {
    bm_pred_original[i] <- val_log
  }
}
compare_df$bm_pred <- bm_pred_original

# 6b. pigauto predictions (predict() returns original-scale values)
if (!is.null(pred_pigauto)) {
  pigauto_pred_original <- numeric(nrow(compare_df))
  pigauto_se_original   <- numeric(nrow(compare_df))

  for (i in seq_len(nrow(compare_df))) {
    sp <- compare_df$species_tip[i]
    tr <- compare_df$trait[i]
    row_idx <- match(sp, pred_pigauto$species_names)
    pigauto_pred_original[i] <- pred_pigauto$imputed[[tr]][row_idx]
    pigauto_se_original[i]   <- pred_pigauto$se[row_idx, tr]
  }
  compare_df$pigauto_pred <- pigauto_pred_original
  compare_df$pigauto_se   <- pigauto_se_original
}

# 6c. Standalone Rphylopars predictions
rp_pred_original <- numeric(nrow(compare_df))
for (i in seq_len(nrow(compare_df))) {
  sp <- compare_df$species_tip[i]
  tr <- compare_df$trait[i]
  row_idx <- match(sp, rownames(rp_recon_orig))
  rp_pred_original[i] <- rp_recon_orig[row_idx, tr]
}
compare_df$rphylopars_pred <- rp_pred_original

# 6d. TabPFN predictions (no phylogeny, if available)
bl_tab <- NULL
if (TABPFN_OK) {
  cat("--- Fitting TabPFN baseline (2000 species) ---\n")
  bl_tab <- tryCatch({
    fit_baseline_tabpfn(pd, splits = splits)
  }, error = function(e) {
    cat("  TabPFN failed (may exceed memory at 2000 species):",
        conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(bl_tab)) {
    tabpfn_pred_original <- numeric(nrow(compare_df))
    for (i in seq_len(nrow(compare_df))) {
      sp <- compare_df$species_tip[i]
      tr <- compare_df$trait[i]
      tm <- pd$trait_map[[tr]]
      row_idx <- match(sp, species_order)
      col_idx <- tm$latent_cols[1]
      val_log <- bl_tab$mu[row_idx, col_idx] * tm$sd + tm$mean
      if (isTRUE(tm$log_transform)) {
        tabpfn_pred_original[i] <- exp(val_log)
      } else {
        tabpfn_pred_original[i] <- val_log
      }
    }
    compare_df$tabpfn_pred <- tabpfn_pred_original
    cat("TabPFN predictions extracted.\n\n")
  }
}

cat("Comparison data frame:", nrow(compare_df), "rows\n")
cat("Columns:", paste(names(compare_df), collapse = ", "), "\n\n")


# -- 7. Compute comparison metrics per trait -----------------------------------

cat("==========================================================\n")
cat("SECTION 7: Compute metrics\n")
cat("==========================================================\n\n")

compute_metrics <- function(truth_vals, pred_vals, method_name, trait_name,
                            n_obs, se_vals = NULL) {
  ok <- is.finite(truth_vals) & is.finite(pred_vals)
  n <- sum(ok)
  if (n == 0) {
    return(data.frame(
      trait = trait_name, method = method_name, n = 0L,
      rmse = NA_real_, mae = NA_real_, pearson_r = NA_real_,
      coverage_95 = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  t <- truth_vals[ok]
  p <- pred_vals[ok]

  rmse_val <- sqrt(mean((t - p)^2))
  mae_val  <- mean(abs(t - p))
  r_val    <- if (n > 1) cor(t, p) else NA_real_

  cov95 <- NA_real_

  if (!is.null(se_vals)) {
    s <- se_vals[ok]
    lower <- p - 1.96 * s
    upper <- p + 1.96 * s
    cov95 <- mean(t >= lower & t <= upper)
  }

  data.frame(
    trait       = trait_name,
    method      = method_name,
    n           = n,
    rmse        = round(rmse_val, 4),
    mae         = round(mae_val, 4),
    pearson_r   = round(r_val, 4),
    coverage_95 = if (is.na(cov95)) NA_character_ else paste0(round(100 * cov95, 1), "%"),
    stringsAsFactors = FALSE
  )
}

# Compute per-trait metrics for test split only
test_df <- compare_df[compare_df$split == "test", ]

results_list <- list()
for (tr in sort(unique(test_df$trait))) {
  sub <- test_df[test_df$trait == tr, ]

  # BM (pigauto's internal baseline)
  results_list[[length(results_list) + 1L]] <- compute_metrics(
    sub$true_value, sub$bm_pred, "BM (pigauto)", tr, nrow(sub)
  )

  # pigauto BM+GNN
  if ("pigauto_pred" %in% names(sub)) {
    results_list[[length(results_list) + 1L]] <- compute_metrics(
      sub$true_value, sub$pigauto_pred, "pigauto (BM+GNN)", tr, nrow(sub),
      se_vals = sub$pigauto_se
    )
  }

  # Standalone Rphylopars
  results_list[[length(results_list) + 1L]] <- compute_metrics(
    sub$true_value, sub$rphylopars_pred, "Rphylopars", tr, nrow(sub)
  )

  # TabPFN (if available)
  if ("tabpfn_pred" %in% names(sub)) {
    results_list[[length(results_list) + 1L]] <- compute_metrics(
      sub$true_value, sub$tabpfn_pred, "TabPFN", tr, nrow(sub)
    )
  }
}

results_table <- do.call(rbind, results_list)
rownames(results_table) <- NULL

cat("Per-trait comparison (TEST split, original scale):\n")
cat("----------------------------------------------------------\n")
print(results_table, row.names = FALSE, right = FALSE)
cat("\n")

# Also compute overall (pooled across traits) metrics
cat("Overall metrics (pooled across all traits, TEST split):\n")
cat("----------------------------------------------------------\n")

overall_list <- list()
all_methods <- c("BM (pigauto)", "pigauto (BM+GNN)", "Rphylopars")
if ("tabpfn_pred" %in% names(compare_df)) {
  all_methods <- c(all_methods, "TabPFN")
}
for (method in all_methods) {
  method_rows <- results_table[results_table$method == method, ]
  if (nrow(method_rows) == 0) next

  # Weighted mean by n
  ns <- method_rows$n
  total_n <- sum(ns)
  w_rmse <- sum(ns * method_rows$rmse, na.rm = TRUE) / total_n
  w_mae  <- sum(ns * method_rows$mae, na.rm = TRUE) / total_n

  # Mean correlation (Fisher-z transform)
  rs <- as.numeric(method_rows$pearson_r)
  rs <- rs[is.finite(rs)]
  if (length(rs) > 0) {
    z <- atanh(pmin(pmax(rs, -0.999), 0.999))
    mean_r <- tanh(mean(z))
  } else {
    mean_r <- NA_real_
  }

  overall_list[[length(overall_list) + 1L]] <- data.frame(
    method    = method,
    n_total   = total_n,
    mean_rmse = round(w_rmse, 4),
    mean_mae  = round(w_mae, 4),
    mean_r    = round(mean_r, 4),
    stringsAsFactors = FALSE
  )
}

overall_table <- do.call(rbind, overall_list)
print(overall_table, row.names = FALSE, right = FALSE)
cat("\n")


# -- 8. Save results and generate figures --------------------------------------

cat("==========================================================\n")
cat("SECTION 8: Save results and figures\n")
cat("==========================================================\n\n")

# 8a. Save comparison table to CSV
csv_path <- file.path(out_dir, "benchmark_avonet2000_results.csv")
write.csv(results_table, csv_path, row.names = FALSE)
cat("Results saved to:", csv_path, "\n")

csv_overall_path <- file.path(out_dir, "benchmark_avonet2000_overall.csv")
write.csv(overall_table, csv_overall_path, row.names = FALSE)
cat("Overall metrics saved to:", csv_overall_path, "\n\n")

# 8b. Per-trait scatter plots: predicted vs truth, coloured by method

# Reshape compare_df to long format for plotting
scatter_long <- data.frame(stringsAsFactors = FALSE)

# BM
bm_part <- data.frame(
  species    = test_df$species_tip,
  trait      = test_df$trait,
  true_value = test_df$true_value,
  predicted  = test_df$bm_pred,
  method     = "BM (pigauto)",
  stringsAsFactors = FALSE
)
scatter_long <- rbind(scatter_long, bm_part)

# pigauto
if ("pigauto_pred" %in% names(test_df)) {
  pig_part <- data.frame(
    species    = test_df$species_tip,
    trait      = test_df$trait,
    true_value = test_df$true_value,
    predicted  = test_df$pigauto_pred,
    method     = "pigauto (BM+GNN)",
    stringsAsFactors = FALSE
  )
  scatter_long <- rbind(scatter_long, pig_part)
}

# Rphylopars
rp_part <- data.frame(
  species    = test_df$species_tip,
  trait      = test_df$trait,
  true_value = test_df$true_value,
  predicted  = test_df$rphylopars_pred,
  method     = "Rphylopars",
  stringsAsFactors = FALSE
)
scatter_long <- rbind(scatter_long, rp_part)

# TabPFN
if ("tabpfn_pred" %in% names(test_df)) {
  tab_part <- data.frame(
    species    = test_df$species_tip,
    trait      = test_df$trait,
    true_value = test_df$true_value,
    predicted  = test_df$tabpfn_pred,
    method     = "TabPFN",
    stringsAsFactors = FALSE
  )
  scatter_long <- rbind(scatter_long, tab_part)
}

# Remove rows with non-finite values
scatter_long <- scatter_long[is.finite(scatter_long$true_value) &
                               is.finite(scatter_long$predicted), ]

method_colours <- c(
  "BM (pigauto)"     = "#e41a1c",
  "pigauto (BM+GNN)" = "#377eb8",
  "Rphylopars"       = "#4daf4a",
  "TabPFN"           = "#984ea3"
)

p_scatter <- ggplot(scatter_long,
                    aes(x = true_value, y = predicted, colour = method)) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ trait, scales = "free") +
  scale_colour_manual(values = method_colours) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom") +
  labs(
    title = "Predicted vs Truth (test split, original scale)",
    x = "True value",
    y = "Predicted value",
    colour = "Method"
  )

scatter_path <- file.path(out_dir, "benchmark_avonet2000_scatter.pdf")
ggsave(scatter_path, p_scatter, width = 12, height = 10)
cat("Scatter plot saved to:", scatter_path, "\n")

# 8c. Uncertainty ribbon plot (pigauto only, per trait)
if (!is.null(pred_pigauto)) {
  # Build a truth lookup for overlay
  truth_lookup <- setNames(truth_sub$true_value,
                           paste(truth_sub$species_tip, truth_sub$trait, sep = "::"))

  ribbon_plots <- list()
  for (tr in trait_names) {
    pred_vals <- pred_pigauto$imputed[[tr]]
    se_vals   <- pred_pigauto$se[, tr]
    sp_names  <- pred_pigauto$species_names

    # Sort by predicted value
    idx <- order(pred_vals)

    df_rib <- data.frame(
      i     = seq_along(pred_vals),
      pred  = pred_vals[idx],
      lower = (pred_vals - 1.96 * se_vals)[idx],
      upper = (pred_vals + 1.96 * se_vals)[idx]
    )

    # Overlay truth where available
    truth_vals <- numeric(length(sp_names))
    truth_vals[] <- NA_real_
    for (j in seq_along(sp_names)) {
      key <- paste(sp_names[idx[j]], tr, sep = "::")
      if (key %in% names(truth_lookup)) {
        truth_vals[j] <- truth_lookup[key]
      }
    }
    df_rib$truth <- truth_vals

    p_rib <- ggplot(df_rib, aes(x = i)) +
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  fill = "#377eb8", alpha = 0.25) +
      geom_line(aes(y = pred), colour = "#377eb8", linewidth = 0.6) +
      geom_point(data = df_rib[is.finite(df_rib$truth), ],
                 aes(y = truth), colour = "#e41a1c", size = 1, alpha = 0.7) +
      theme_minimal(base_size = 11) +
      labs(
        title = paste("pigauto uncertainty:", tr),
        x = "Species (sorted by prediction)",
        y = tr
      )

    ribbon_plots[[tr]] <- p_rib
  }

  ribbon_path <- file.path(out_dir, "benchmark_avonet2000_uncertainty.pdf")
  pdf(ribbon_path, width = 10, height = 6)
  for (p in ribbon_plots) print(p)
  dev.off()
  cat("Uncertainty ribbon plots saved to:", ribbon_path, "\n")
}

# 8d. Bar chart: RMSE by trait x method
results_table$rmse_num <- as.numeric(results_table$rmse)
p_bar <- ggplot(results_table,
                aes(x = trait, y = rmse_num, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = method_colours) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  labs(
    title = "RMSE by Trait and Method (test split, original scale)",
    x = "Trait",
    y = "RMSE",
    fill = "Method"
  )

bar_path <- file.path(out_dir, "benchmark_avonet2000_rmse_bar.pdf")
ggsave(bar_path, p_bar, width = 10, height = 6)
cat("RMSE bar chart saved to:", bar_path, "\n")

# 8e. pigauto training history (if available)
if (!is.null(fit)) {
  p_hist <- plot(fit, type = "history")
  hist_path <- file.path(out_dir, "benchmark_avonet2000_training_history.pdf")
  ggsave(hist_path, p_hist, width = 8, height = 5)
  cat("Training history saved to:", hist_path, "\n")
}

cat("\n==========================================================\n")
cat("Benchmark complete.\n")
cat("==========================================================\n")
