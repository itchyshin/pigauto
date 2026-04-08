# ============================================================================
# Benchmark: pigauto mixed-type imputation on AVONET 300 species
# ============================================================================
#
# Evaluates mixed-type imputation (continuous + categorical + ordinal) using
# the bundled avonet300 dataset. Uses pigauto's own val/test masking to
# create held-out cells, then compares BM baseline vs pigauto (BM + GNN).
#
# Traits:
#   - Continuous (4): Mass, Beak.Length_Culmen, Tarsus.Length, Wing.Length
#   - Categorical (2): Trophic.Level (4 levels), Primary.Lifestyle (5 levels)
#   - Ordinal (1): Migration (3 levels)
#
# Methods:
#   1. BM baseline (Rphylopars for continuous, frequency for categorical)
#   2. pigauto (BM + GNN) with mixed-type loss
#
# Run from pigauto root:
#   source("script/benchmark_avonet_mixed.R")
# ============================================================================

library(pigauto)
library(ape)
library(ggplot2)

cat("==================================================================\n")
cat("pigauto mixed-type benchmark: AVONET 300 species\n")
cat("==================================================================\n\n")

# ---- Configuration ----------------------------------------------------------
N_REPS        <- 5L       # random-split replicates
MISSING_FRAC  <- 0.25
EPOCHS        <- 2000L
EVAL_EVERY    <- 100L
K_EIGEN       <- 8L
SEED_BASE     <- 100L

out_dir <- file.path(getwd(), "script")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# ---- Load data --------------------------------------------------------------
data(avonet300, tree300)

df <- avonet300
rownames(df) <- df$Species_Key
df$Species_Key <- NULL

cat("Dataset: 300 species x 7 traits\n")
cat("Trait types:\n")
for (col in names(df)) {
  cat(sprintf("  %-20s %s\n", col, paste(class(df[[col]]), collapse = "/")))
}
cat("\n")


# ---- Run replicates ---------------------------------------------------------
all_results <- list()

for (rep_id in seq_len(N_REPS)) {
  seed <- SEED_BASE + rep_id
  set.seed(seed)
  cat(sprintf("\n--- Replicate %d/%d (seed %d) ---\n", rep_id, N_REPS, seed))

  # 1. Preprocess
  pd <- preprocess_traits(df, tree300, log_transform = TRUE)
  cat(sprintf("  Latent matrix: %d x %d (%d latent cols)\n",
              nrow(pd$X_scaled), ncol(pd$X_scaled), pd$p_latent))

  # 2. Split
  splits <- make_missing_splits(pd$X_scaled, missing_frac = MISSING_FRAC,
                                seed = seed, trait_map = pd$trait_map)
  cat(sprintf("  Val cells: %d, Test cells: %d\n",
              length(splits$val_idx), length(splits$test_idx)))

  truth <- pd$X_scaled

  # 3. BM baseline
  cat("  Fitting BM baseline...\n")
  bl <- fit_baseline(pd, tree300, splits = splits)
  eval_bl <- tryCatch(
    evaluate_imputation(bl$mu, truth, splits, trait_map = pd$trait_map),
    error = function(e) { cat("  [ERROR]", conditionMessage(e), "\n"); NULL }
  )

  # 4. pigauto (BM + GNN)
  cat(sprintf("  Training pigauto (epochs = %d)...\n", EPOCHS))
  graph <- build_phylo_graph(tree300, k_eigen = K_EIGEN)
  fit <- tryCatch(
    fit_pigauto(pd, tree300, splits = splits, graph = graph,
                baseline = bl, epochs = EPOCHS, eval_every = EVAL_EVERY,
                verbose = FALSE, seed = seed),
    error = function(e) { cat("  [ERROR]", conditionMessage(e), "\n"); NULL }
  )

  eval_gnn <- NULL
  if (!is.null(fit)) {
    pred <- predict(fit, return_se = TRUE)
    eval_gnn <- tryCatch(
      evaluate_imputation(pred, truth, splits),
      error = function(e) { cat("  [ERROR]", conditionMessage(e), "\n"); NULL }
    )
  }

  # 5. Collect results
  collect <- function(eval_df, method_name) {
    if (is.null(eval_df)) return(NULL)
    test_df <- eval_df[eval_df$split == "test", , drop = FALSE]
    if (nrow(test_df) == 0) return(NULL)
    test_df$method <- method_name
    test_df$rep    <- rep_id
    test_df$seed   <- seed
    test_df
  }

  rep_results <- rbind(
    collect(eval_bl,  "BM_baseline"),
    collect(eval_gnn, "pigauto_GNN")
  )

  if (!is.null(rep_results) && nrow(rep_results) > 0) {
    all_results[[length(all_results) + 1]] <- rep_results
  }

  # Print rep summary
  if (!is.null(rep_results)) {
    for (method in c("BM_baseline", "pigauto_GNN")) {
      md <- rep_results[rep_results$method == method, ]
      if (nrow(md) == 0) next
      cat(sprintf("  %s:\n", method))
      for (i in seq_len(nrow(md))) {
        r <- md[i, ]
        if (r$type %in% c("continuous", "count")) {
          cat(sprintf("    %-20s RMSE=%.4f  r=%.4f\n",
                      r$trait, r$rmse, r$pearson_r))
        } else if (r$type %in% c("binary", "categorical")) {
          cat(sprintf("    %-20s Accuracy=%.3f\n", r$trait, r$accuracy))
        } else if (r$type == "ordinal") {
          cat(sprintf("    %-20s RMSE=%.4f  rho=%.4f\n",
                      r$trait, r$rmse, r$spearman_rho))
        }
      }
    }
  }
}


# ---- Aggregate results -------------------------------------------------------
results <- do.call(rbind, all_results)

if (nrow(results) > 0) {
  cat("\n\n")
  cat("==================================================================\n")
  cat("SUMMARY: Mean across", N_REPS, "replicates\n")
  cat("==================================================================\n")

  # Per-trait, per-method summary
  for (trait_nm in unique(results$trait)) {
    trait_df <- results[results$trait == trait_nm, ]
    trait_type <- trait_df$type[1]

    cat(sprintf("\n  %s (%s):\n", trait_nm, trait_type))

    for (method in c("BM_baseline", "pigauto_GNN")) {
      md <- trait_df[trait_df$method == method, ]
      if (nrow(md) == 0) next

      if (trait_type %in% c("continuous", "count")) {
        cat(sprintf("    %-15s RMSE=%.4f (sd=%.4f)  r=%.4f\n",
                    method,
                    mean(md$rmse, na.rm = TRUE),
                    stats::sd(md$rmse, na.rm = TRUE),
                    mean(md$pearson_r, na.rm = TRUE)))
      } else if (trait_type %in% c("binary", "categorical")) {
        cat(sprintf("    %-15s Accuracy=%.3f (sd=%.3f)\n",
                    method,
                    mean(md$accuracy, na.rm = TRUE),
                    stats::sd(md$accuracy, na.rm = TRUE)))
      } else if (trait_type == "ordinal") {
        cat(sprintf("    %-15s RMSE=%.4f  rho=%.4f\n",
                    method,
                    mean(md$rmse, na.rm = TRUE),
                    mean(md$spearman_rho, na.rm = TRUE)))
      }
    }

    # GNN improvement
    bl_df <- trait_df[trait_df$method == "BM_baseline", ]
    gn_df <- trait_df[trait_df$method == "pigauto_GNN", ]
    if (nrow(bl_df) > 0 && nrow(gn_df) > 0) {
      if (trait_type %in% c("continuous", "count", "ordinal")) {
        bl_rmse <- mean(bl_df$rmse, na.rm = TRUE)
        gn_rmse <- mean(gn_df$rmse, na.rm = TRUE)
        pct <- (bl_rmse - gn_rmse) / bl_rmse * 100
        cat(sprintf("    --> GNN improvement: %+.1f%% RMSE\n", pct))
      } else {
        bl_acc <- mean(bl_df$accuracy, na.rm = TRUE)
        gn_acc <- mean(gn_df$accuracy, na.rm = TRUE)
        cat(sprintf("    --> GNN accuracy change: %+.1f pp\n",
                    (gn_acc - bl_acc) * 100))
      }
    }
  }

  # Save results
  csv_path <- file.path(out_dir, "benchmark_avonet_mixed_results.csv")
  write.csv(results, csv_path, row.names = FALSE)
  cat(sprintf("\nResults saved to %s\n", csv_path))
}

cat("\nDone.\n")
