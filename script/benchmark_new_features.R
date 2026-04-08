# ============================================================================
# Benchmark: OLD pigauto (no attention/calibration) vs NEW pigauto
#            (attention + calibration + conformal + phylo label propagation)
# ============================================================================
#
# Compares two pigauto configurations on the AVONET 300 mixed-type dataset:
#   OLD: use_attention = FALSE  (no gate calibration, no conformal)
#   NEW: use_attention = TRUE   (with gate calibration + conformal intervals)
#
# Evaluation:
#   - Continuous: RMSE, Pearson r
#   - Categorical: Accuracy
#   - Ordinal: RMSE, Spearman rho
#   - Conformal coverage (NEW only): fraction of test values within
#     conformal_lower/conformal_upper
# ============================================================================

devtools::load_all()
library(ape)

cat("==================================================================\n")
cat("Benchmark: OLD vs NEW pigauto on AVONET 300\n")
cat("==================================================================\n\n")

# ---- Configuration ----------------------------------------------------------
N_REPS        <- 3L
SEEDS         <- c(101L, 102L, 103L)
MISSING_FRAC  <- 0.25
EPOCHS        <- 500L
EVAL_EVERY    <- 50L
PATIENCE      <- 10L
K_EIGEN       <- 8L

out_dir <- file.path(getwd(), "script")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Load data --------------------------------------------------------------
data(avonet300, tree300)

df <- avonet300
rownames(df) <- df$Species_Key
df$Species_Key <- NULL

cat("Dataset: 300 species x", ncol(df), "traits\n")
cat("Trait types:\n")
for (col in names(df)) {
  cat(sprintf("  %-20s %s\n", col, paste(class(df[[col]]), collapse = "/")))
}
cat("\n")

# ---- Run replicates ---------------------------------------------------------
all_results <- list()

for (rep_id in seq_len(N_REPS)) {
  seed <- SEEDS[rep_id]
  set.seed(seed)
  cat(sprintf("\n=== Replicate %d/%d (seed %d) ===\n", rep_id, N_REPS, seed))

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
    error = function(e) { cat("  [BM ERROR]", conditionMessage(e), "\n"); NULL }
  )

  # 4. Build graph (shared between OLD and NEW)
  graph <- build_phylo_graph(tree300, k_eigen = K_EIGEN)

  # 5. OLD pigauto: use_attention = FALSE
  cat(sprintf("  Training OLD pigauto (attention=FALSE, epochs=%d)...\n", EPOCHS))
  fit_old <- tryCatch(
    fit_pigauto(pd, tree300, splits = splits, graph = graph,
                baseline = bl, epochs = EPOCHS, eval_every = EVAL_EVERY,
                patience = PATIENCE, use_attention = FALSE,
                verbose = FALSE, seed = seed),
    error = function(e) { cat("  [OLD ERROR]", conditionMessage(e), "\n"); NULL }
  )

  eval_old <- NULL
  conformal_coverage_old <- NULL
  if (!is.null(fit_old)) {
    pred_old <- predict(fit_old, return_se = TRUE)
    eval_old <- tryCatch(
      evaluate_imputation(pred_old, truth, splits),
      error = function(e) { cat("  [OLD EVAL ERROR]", conditionMessage(e), "\n"); NULL }
    )
  }

  # 6. NEW pigauto: use_attention = TRUE (gets calibration + conformal)
  cat(sprintf("  Training NEW pigauto (attention=TRUE, epochs=%d)...\n", EPOCHS))
  fit_new <- tryCatch(
    fit_pigauto(pd, tree300, splits = splits, graph = graph,
                baseline = bl, epochs = EPOCHS, eval_every = EVAL_EVERY,
                patience = PATIENCE, use_attention = TRUE,
                verbose = FALSE, seed = seed),
    error = function(e) { cat("  [NEW ERROR]", conditionMessage(e), "\n"); NULL }
  )

  eval_new <- NULL
  conformal_coverage_new <- NULL
  if (!is.null(fit_new)) {
    pred_new <- predict(fit_new, return_se = TRUE)
    eval_new <- tryCatch(
      evaluate_imputation(pred_new, truth, splits),
      error = function(e) { cat("  [NEW EVAL ERROR]", conditionMessage(e), "\n"); NULL }
    )

    # Conformal coverage: fraction of test values within conformal intervals
    if (!is.null(pred_new$conformal_lower) && !is.null(pred_new$conformal_upper)) {
      conformal_coverage_new <- list()
      n <- nrow(truth)

      # Build test mask matrix
      test_mat <- matrix(FALSE, n, ncol(truth))
      test_mat[splits$test_idx] <- TRUE

      for (tm in pd$trait_map) {
        nm <- tm$name
        if (!(tm$type %in% c("continuous", "count", "ordinal"))) next
        lc <- tm$latent_cols

        # Test rows for this trait
        test_rows <- which(test_mat[, lc[1]])
        if (length(test_rows) == 0) next

        # Get truth in original scale for coverage check
        truth_latent <- truth[test_rows, lc[1]]
        ok <- is.finite(truth_latent)
        if (sum(ok) == 0) next

        # Back-transform truth to original scale for comparison with conformal bounds
        if (tm$type == "continuous") {
          truth_orig <- truth_latent[ok] * tm$sd + tm$mean
          if (isTRUE(tm$log_transform)) truth_orig <- exp(truth_orig)
        } else if (tm$type == "count") {
          truth_orig <- expm1(truth_latent[ok] * tm$sd + tm$mean)
        } else if (tm$type == "ordinal") {
          truth_orig <- round(truth_latent[ok] * tm$sd + tm$mean)
        }

        lower <- pred_new$conformal_lower[test_rows[ok], nm]
        upper <- pred_new$conformal_upper[test_rows[ok], nm]

        coverage <- mean(truth_orig >= lower & truth_orig <= upper, na.rm = TRUE)
        conformal_coverage_new[[nm]] <- coverage
      }
    }
  }

  # 7. Collect results
  collect <- function(eval_df, method_name, coverage = NULL) {
    if (is.null(eval_df)) return(NULL)
    test_df <- eval_df[eval_df$split == "test", , drop = FALSE]
    if (nrow(test_df) == 0) return(NULL)
    test_df$method <- method_name
    test_df$rep    <- rep_id
    test_df$seed   <- seed

    # Add conformal coverage column
    test_df$conformal_coverage <- NA_real_
    if (!is.null(coverage)) {
      for (i in seq_len(nrow(test_df))) {
        nm <- test_df$trait[i]
        if (nm %in% names(coverage)) {
          test_df$conformal_coverage[i] <- coverage[[nm]]
        }
      }
    }
    test_df
  }

  rep_results <- rbind(
    collect(eval_bl,  "BM_baseline"),
    collect(eval_old, "pigauto_OLD"),
    collect(eval_new, "pigauto_NEW", conformal_coverage_new)
  )

  if (!is.null(rep_results) && nrow(rep_results) > 0) {
    all_results[[length(all_results) + 1]] <- rep_results
  }

  # Print rep summary
  if (!is.null(rep_results)) {
    for (method in c("BM_baseline", "pigauto_OLD", "pigauto_NEW")) {
      md <- rep_results[rep_results$method == method, ]
      if (nrow(md) == 0) next
      cat(sprintf("  %s:\n", method))
      for (i in seq_len(nrow(md))) {
        r <- md[i, ]
        if (r$type %in% c("continuous", "count")) {
          msg <- sprintf("    %-20s RMSE=%.4f  r=%.4f", r$trait, r$rmse, r$pearson_r)
          if (!is.na(r$conformal_coverage)) {
            msg <- paste0(msg, sprintf("  conformal_cov=%.3f", r$conformal_coverage))
          }
          cat(msg, "\n")
        } else if (r$type %in% c("binary", "categorical")) {
          cat(sprintf("    %-20s Accuracy=%.3f\n", r$trait, r$accuracy))
        } else if (r$type == "ordinal") {
          msg <- sprintf("    %-20s RMSE=%.4f  rho=%.4f", r$trait, r$rmse, r$spearman_rho)
          if (!is.na(r$conformal_coverage)) {
            msg <- paste0(msg, sprintf("  conformal_cov=%.3f", r$conformal_coverage))
          }
          cat(msg, "\n")
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

    for (method in c("BM_baseline", "pigauto_OLD", "pigauto_NEW")) {
      md <- trait_df[trait_df$method == method, ]
      if (nrow(md) == 0) next

      if (trait_type %in% c("continuous", "count")) {
        msg <- sprintf("    %-15s RMSE=%.4f (sd=%.4f)  r=%.4f",
                    method,
                    mean(md$rmse, na.rm = TRUE),
                    stats::sd(md$rmse, na.rm = TRUE),
                    mean(md$pearson_r, na.rm = TRUE))
        if (any(!is.na(md$conformal_coverage))) {
          msg <- paste0(msg, sprintf("  conf_cov=%.3f",
                                     mean(md$conformal_coverage, na.rm = TRUE)))
        }
        cat(msg, "\n")
      } else if (trait_type %in% c("binary", "categorical")) {
        cat(sprintf("    %-15s Accuracy=%.3f (sd=%.3f)\n",
                    method,
                    mean(md$accuracy, na.rm = TRUE),
                    stats::sd(md$accuracy, na.rm = TRUE)))
      } else if (trait_type == "ordinal") {
        msg <- sprintf("    %-15s RMSE=%.4f  rho=%.4f",
                    method,
                    mean(md$rmse, na.rm = TRUE),
                    mean(md$spearman_rho, na.rm = TRUE))
        if (any(!is.na(md$conformal_coverage))) {
          msg <- paste0(msg, sprintf("  conf_cov=%.3f",
                                     mean(md$conformal_coverage, na.rm = TRUE)))
        }
        cat(msg, "\n")
      }
    }

    # OLD vs NEW comparison
    old_df <- trait_df[trait_df$method == "pigauto_OLD", ]
    new_df <- trait_df[trait_df$method == "pigauto_NEW", ]
    if (nrow(old_df) > 0 && nrow(new_df) > 0) {
      if (trait_type %in% c("continuous", "count", "ordinal")) {
        old_rmse <- mean(old_df$rmse, na.rm = TRUE)
        new_rmse <- mean(new_df$rmse, na.rm = TRUE)
        pct <- (old_rmse - new_rmse) / old_rmse * 100
        cat(sprintf("    --> NEW vs OLD: %+.1f%% RMSE\n", pct))
      } else {
        old_acc <- mean(old_df$accuracy, na.rm = TRUE)
        new_acc <- mean(new_df$accuracy, na.rm = TRUE)
        cat(sprintf("    --> NEW vs OLD: %+.1f pp accuracy\n",
                    (new_acc - old_acc) * 100))
      }
    }

    # BM vs NEW comparison
    bl_df <- trait_df[trait_df$method == "BM_baseline", ]
    if (nrow(bl_df) > 0 && nrow(new_df) > 0) {
      if (trait_type %in% c("continuous", "count", "ordinal")) {
        bl_rmse <- mean(bl_df$rmse, na.rm = TRUE)
        new_rmse <- mean(new_df$rmse, na.rm = TRUE)
        pct <- (bl_rmse - new_rmse) / bl_rmse * 100
        cat(sprintf("    --> NEW vs BM:  %+.1f%% RMSE\n", pct))
      } else {
        bl_acc <- mean(bl_df$accuracy, na.rm = TRUE)
        new_acc <- mean(new_df$accuracy, na.rm = TRUE)
        cat(sprintf("    --> NEW vs BM:  %+.1f pp accuracy\n",
                    (new_acc - bl_acc) * 100))
      }
    }
  }

  # Conformal coverage summary (NEW only)
  new_results <- results[results$method == "pigauto_NEW", ]
  conf_traits <- new_results[!is.na(new_results$conformal_coverage), ]
  if (nrow(conf_traits) > 0) {
    cat("\n\n  Conformal 95% coverage (NEW pigauto, target = 0.95):\n")
    for (trait_nm in unique(conf_traits$trait)) {
      td <- conf_traits[conf_traits$trait == trait_nm, ]
      cat(sprintf("    %-20s coverage=%.3f (sd=%.3f)\n",
                  trait_nm,
                  mean(td$conformal_coverage, na.rm = TRUE),
                  stats::sd(td$conformal_coverage, na.rm = TRUE)))
    }
  }

  # Save results
  csv_path <- file.path(out_dir, "benchmark_new_features.csv")
  write.csv(results, csv_path, row.names = FALSE)
  cat(sprintf("\nResults saved to %s\n", csv_path))
}

cat("\nDone.\n")
