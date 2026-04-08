# ============================================================================
# Final benchmark: pigauto on AVONET 300 dataset
# ============================================================================
#
# Compares 3 methods across 5 replicates (seeds 101-105):
#   1. BM_baseline  -- phylogenetic BM imputation only
#   2. pigauto_NEW  -- BM + GNN with default settings (epochs=1000)
#   3. pigauto_LONG -- BM + GNN with extended training (epochs=3000)
#
# Traits (AVONET 300 species):
#   Continuous (4): Mass, Beak.Length_Culmen, Tarsus.Length, Wing.Length
#   Categorical (2): Trophic.Level, Primary.Lifestyle
#   Ordinal (1): Migration
#
# Run from pigauto root:
#   Rscript script/benchmark_final.R
# ============================================================================

cat("============================================================\n")
cat("pigauto Final Benchmark: AVONET 300 species\n")
cat("============================================================\n\n")

# ---- Load package ----------------------------------------------------------
devtools::load_all(".")
library(ape)

cat("torch available:", torch::torch_is_installed(), "\n\n")

# ---- Configuration ---------------------------------------------------------
N_REPS        <- 5L
SEEDS         <- 101:105
MISSING_FRAC  <- 0.25
K_EIGEN       <- "auto"
EVAL_EVERY    <- 50L

# Method-specific settings
EPOCHS_NEW    <- 1000L
PATIENCE_NEW  <- 10L

EPOCHS_LONG   <- 3000L
PATIENCE_LONG <- 15L

OUTPUT_CSV    <- "script/benchmark_final.csv"

# ---- Load data -------------------------------------------------------------
data(avonet300, tree300)

df <- avonet300
rownames(df) <- df$Species_Key
df$Species_Key <- NULL

cat("Dataset: 300 species x 7 traits\n")
cat("Trait types:\n")
for (col in names(df)) {
  cat(sprintf("  %-25s %s\n", col, paste(class(df[[col]]), collapse = "/")))
}
cat("\n")

# ---- Run replicates --------------------------------------------------------
all_results <- list()
t_start_all <- Sys.time()

for (rep_id in seq_len(N_REPS)) {
  seed <- SEEDS[rep_id]
  set.seed(seed)
  t_start_rep <- Sys.time()

  cat(sprintf("\n%s\n", strrep("=", 60)))
  cat(sprintf("Replicate %d/%d (seed %d)\n", rep_id, N_REPS, seed))
  cat(sprintf("%s\n", strrep("=", 60)))

  # ---- 4a. Preprocess traits ------------------------------------------------
  cat("  [1/6] Preprocessing traits...\n")
  pd <- preprocess_traits(df, tree300, log_transform = TRUE)
  cat(sprintf("         Latent matrix: %d x %d (%d latent cols)\n",
              nrow(pd$X_scaled), ncol(pd$X_scaled), pd$p_latent))

  # ---- 4b. Create splits ----------------------------------------------------
  cat("  [2/6] Creating train/val/test splits...\n")
  spl <- make_missing_splits(pd$X_scaled, missing_frac = MISSING_FRAC,
                             seed = seed, trait_map = pd$trait_map)
  cat(sprintf("         Val cells: %d, Test cells: %d\n",
              length(spl$val_idx), length(spl$test_idx)))

  truth <- pd$X_scaled

  # ---- 4c. Build graph (reused across methods) ------------------------------
  cat("  [3/6] Building phylogenetic graph...\n")
  graph <- build_phylo_graph(tree300, k_eigen = K_EIGEN)
  cat(sprintf("         k_eigen: %d (auto-selected)\n", ncol(graph$coords)))

  # ---- 4d. Fit baseline (reused across methods) -----------------------------
  cat("  [4/6] Fitting BM baseline...\n")
  bl <- tryCatch(
    fit_baseline(pd, tree300, splits = spl),
    error = function(e) {
      cat("         [ERROR] fit_baseline failed:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(bl)) {
    cat("  Skipping replicate (baseline failed).\n")
    next
  }

  # Evaluate baseline on test set
  eval_bl <- tryCatch(
    evaluate_imputation(bl$mu, truth, spl, trait_map = pd$trait_map),
    error = function(e) {
      cat("         [ERROR] evaluate baseline:", conditionMessage(e), "\n")
      NULL
    }
  )

  # ---- 4e. Fit pigauto_NEW (epochs=1000, patience=10) -----------------------
  cat(sprintf("  [5/6] Training pigauto_NEW (epochs=%d, patience=%d)...\n",
              EPOCHS_NEW, PATIENCE_NEW))

  fit_new <- tryCatch({
    fit_pigauto(
      data       = pd,
      tree       = tree300,
      splits     = spl,
      graph      = graph,
      baseline   = bl,
      use_attention = TRUE,
      k_eigen    = K_EIGEN,
      epochs     = EPOCHS_NEW,
      eval_every = EVAL_EVERY,
      patience   = PATIENCE_NEW,
      verbose    = FALSE,
      seed       = seed
    )
  }, error = function(e) {
    cat("         [ERROR] fit_pigauto (NEW):", conditionMessage(e), "\n")
    NULL
  })

  eval_new <- NULL
  pred_new <- NULL
  if (!is.null(fit_new)) {
    cat(sprintf("         Stopped at epoch %d (val_loss=%.4f)\n",
                max(fit_new$history$epoch), fit_new$val_rmse))
    pred_new <- tryCatch(
      predict(fit_new, return_se = TRUE),
      error = function(e) {
        cat("         [ERROR] predict (NEW):", conditionMessage(e), "\n")
        NULL
      }
    )
    if (!is.null(pred_new)) {
      eval_new <- tryCatch(
        evaluate_imputation(pred_new, truth, spl),
        error = function(e) {
          cat("         [ERROR] evaluate (NEW):", conditionMessage(e), "\n")
          NULL
        }
      )
    }
  }

  # ---- 4f. Fit pigauto_LONG (epochs=3000, patience=15) ----------------------
  cat(sprintf("  [6/6] Training pigauto_LONG (epochs=%d, patience=%d)...\n",
              EPOCHS_LONG, PATIENCE_LONG))

  fit_long <- tryCatch({
    fit_pigauto(
      data       = pd,
      tree       = tree300,
      splits     = spl,
      graph      = graph,
      baseline   = bl,
      use_attention = TRUE,
      k_eigen    = K_EIGEN,
      epochs     = EPOCHS_LONG,
      eval_every = EVAL_EVERY,
      patience   = PATIENCE_LONG,
      verbose    = FALSE,
      seed       = seed
    )
  }, error = function(e) {
    cat("         [ERROR] fit_pigauto (LONG):", conditionMessage(e), "\n")
    NULL
  })

  eval_long <- NULL
  pred_long <- NULL
  if (!is.null(fit_long)) {
    cat(sprintf("         Stopped at epoch %d (val_loss=%.4f)\n",
                max(fit_long$history$epoch), fit_long$val_rmse))
    pred_long <- tryCatch(
      predict(fit_long, return_se = TRUE),
      error = function(e) {
        cat("         [ERROR] predict (LONG):", conditionMessage(e), "\n")
        NULL
      }
    )
    if (!is.null(pred_long)) {
      eval_long <- tryCatch(
        evaluate_imputation(pred_long, truth, spl),
        error = function(e) {
          cat("         [ERROR] evaluate (LONG):", conditionMessage(e), "\n")
          NULL
        }
      )
    }
  }

  # ---- 4g. Collect results for this replicate --------------------------------
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
    collect(eval_bl,   "BM_baseline"),
    collect(eval_new,  "pigauto_NEW"),
    collect(eval_long, "pigauto_LONG")
  )

  if (!is.null(rep_results) && nrow(rep_results) > 0) {
    all_results[[length(all_results) + 1]] <- rep_results
  }

  # Print replicate summary
  t_rep <- round(difftime(Sys.time(), t_start_rep, units = "mins"), 1)
  cat(sprintf("\n  Replicate %d complete (%.1f min)\n", rep_id, t_rep))

  if (!is.null(rep_results)) {
    for (method in c("BM_baseline", "pigauto_NEW", "pigauto_LONG")) {
      md <- rep_results[rep_results$method == method, ]
      if (nrow(md) == 0) next
      cat(sprintf("    %s:\n", method))
      for (i in seq_len(nrow(md))) {
        r <- md[i, ]
        if (r$type %in% c("continuous", "count")) {
          cat(sprintf("      %-25s RMSE=%.4f  r=%.4f\n",
                      r$trait, r$rmse, r$pearson_r))
        } else if (r$type %in% c("binary", "categorical")) {
          cat(sprintf("      %-25s Accuracy=%.3f\n", r$trait, r$accuracy))
        } else if (r$type == "ordinal") {
          cat(sprintf("      %-25s RMSE=%.4f  rho=%.4f\n",
                      r$trait, r$rmse, r$spearman_rho))
        }
      }
    }
  }

  # Save intermediate progress
  if (length(all_results) > 0) {
    interim <- do.call(rbind, all_results)
    write.csv(interim, OUTPUT_CSV, row.names = FALSE)
    cat(sprintf("  [saved %d rows to %s]\n", nrow(interim), OUTPUT_CSV))
  }
}


# ============================================================================
# Aggregate results
# ============================================================================
results <- do.call(rbind, all_results)
rownames(results) <- NULL

t_total <- round(difftime(Sys.time(), t_start_all, units = "mins"), 1)
cat(sprintf("\n\nTotal runtime: %.1f minutes\n", t_total))

if (nrow(results) == 0) {
  cat("No results collected. Exiting.\n")
  quit(save = "no", status = 1)
}

write.csv(results, OUTPUT_CSV, row.names = FALSE)
cat(sprintf("Results saved to %s (%d rows)\n\n", OUTPUT_CSV, nrow(results)))


# ============================================================================
# Summary table
# ============================================================================
cat(strrep("=", 70), "\n")
cat("SUMMARY: Mean across", N_REPS, "replicates (test set)\n")
cat(strrep("=", 70), "\n\n")

# Per-trait, per-method summary
methods_order <- c("BM_baseline", "pigauto_NEW", "pigauto_LONG")

for (trait_nm in unique(results$trait)) {
  trait_df <- results[results$trait == trait_nm, ]
  trait_type <- trait_df$type[1]

  cat(sprintf("  %s (%s):\n", trait_nm, trait_type))

  for (method in methods_order) {
    md <- trait_df[trait_df$method == method, ]
    if (nrow(md) == 0) next

    if (trait_type %in% c("continuous", "count")) {
      cat(sprintf("    %-18s RMSE=%.4f (sd=%.4f)  r=%.4f  coverage=%.3f\n",
                  method,
                  mean(md$rmse, na.rm = TRUE),
                  stats::sd(md$rmse, na.rm = TRUE),
                  mean(md$pearson_r, na.rm = TRUE),
                  mean(md$coverage_95, na.rm = TRUE)))
    } else if (trait_type %in% c("binary", "categorical")) {
      cat(sprintf("    %-18s Accuracy=%.3f (sd=%.3f)  Brier=%.4f\n",
                  method,
                  mean(md$accuracy, na.rm = TRUE),
                  stats::sd(md$accuracy, na.rm = TRUE),
                  mean(md$brier, na.rm = TRUE)))
    } else if (trait_type == "ordinal") {
      cat(sprintf("    %-18s RMSE=%.4f (sd=%.4f)  rho=%.4f\n",
                  method,
                  mean(md$rmse, na.rm = TRUE),
                  stats::sd(md$rmse, na.rm = TRUE),
                  mean(md$spearman_rho, na.rm = TRUE)))
    }
  }

  # Improvement of LONG over baseline
  bl_df <- trait_df[trait_df$method == "BM_baseline", ]
  lo_df <- trait_df[trait_df$method == "pigauto_LONG", ]
  ne_df <- trait_df[trait_df$method == "pigauto_NEW", ]
  if (nrow(bl_df) > 0 && nrow(lo_df) > 0) {
    if (trait_type %in% c("continuous", "count", "ordinal")) {
      bl_rmse <- mean(bl_df$rmse, na.rm = TRUE)
      lo_rmse <- mean(lo_df$rmse, na.rm = TRUE)
      ne_rmse <- mean(ne_df$rmse, na.rm = TRUE)
      pct_new  <- (bl_rmse - ne_rmse) / bl_rmse * 100
      pct_long <- (bl_rmse - lo_rmse) / bl_rmse * 100
      cat(sprintf("    --> NEW  improvement: %+.1f%% RMSE\n", pct_new))
      cat(sprintf("    --> LONG improvement: %+.1f%% RMSE\n", pct_long))
    } else {
      bl_acc <- mean(bl_df$accuracy, na.rm = TRUE)
      lo_acc <- mean(lo_df$accuracy, na.rm = TRUE)
      ne_acc <- mean(ne_df$accuracy, na.rm = TRUE)
      cat(sprintf("    --> NEW  change: %+.1f pp accuracy\n",
                  (ne_acc - bl_acc) * 100))
      cat(sprintf("    --> LONG change: %+.1f pp accuracy\n",
                  (lo_acc - bl_acc) * 100))
    }
  }
  cat("\n")
}

# Overall summary by method
cat(strrep("-", 70), "\n")
cat("OVERALL SUMMARY:\n")
cat(strrep("-", 70), "\n")

for (method in methods_order) {
  md <- results[results$method == method, ]
  cont <- md[md$type %in% c("continuous", "count"), ]
  disc <- md[md$type %in% c("binary", "categorical"), ]
  ordi <- md[md$type == "ordinal", ]

  cat(sprintf("\n  %s:\n", method))
  if (nrow(cont) > 0) {
    cat(sprintf("    Continuous/count mean RMSE: %.4f\n",
                mean(cont$rmse, na.rm = TRUE)))
    cat(sprintf("    Continuous/count mean r:    %.4f\n",
                mean(cont$pearson_r, na.rm = TRUE)))
  }
  if (nrow(disc) > 0) {
    cat(sprintf("    Discrete mean accuracy:    %.3f\n",
                mean(disc$accuracy, na.rm = TRUE)))
  }
  if (nrow(ordi) > 0) {
    cat(sprintf("    Ordinal mean RMSE:         %.4f\n",
                mean(ordi$rmse, na.rm = TRUE)))
    cat(sprintf("    Ordinal mean rho:          %.4f\n",
                mean(ordi$spearman_rho, na.rm = TRUE)))
  }
}


# ============================================================================
# Generate HTML report using the last replicate
# ============================================================================
cat("\n\n")
cat(strrep("=", 70), "\n")
cat("GENERATING HTML REPORT\n")
cat(strrep("=", 70), "\n")

report_fit <- if (!is.null(fit_long)) fit_long else fit_new
if (!is.null(report_fit)) {
  tryCatch({
    pigauto_report(
      fit         = report_fit,
      data        = pd,
      splits      = spl,
      output_path = "script/benchmark_final_report.html",
      title       = "pigauto Final Benchmark Report",
      open        = FALSE
    )
    cat("Report saved to: script/benchmark_final_report.html\n")
  }, error = function(e) {
    cat("  [ERROR] pigauto_report failed:", conditionMessage(e), "\n")
  })
} else {
  cat("  No fit object available for report generation.\n")
}


cat("\n")
cat(strrep("=", 70), "\n")
cat("BENCHMARK COMPLETE\n")
cat(strrep("=", 70), "\n")
