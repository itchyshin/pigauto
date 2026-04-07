# ============================================================================
# Simulation benchmark for pigauto
# ============================================================================
#
# Generates phylogenetic data with known truth via BACE's sim_bace() and
# pigauto's simulate_non_bm(), masks values with pigauto's own masking,
# then compares up to 4 imputation methods across 12 scenarios:
#
#   1. BM_baseline   — Rphylopars BM imputation (phylogenetic)
#   2. pigauto_GNN   — BM baseline + graph autoencoder correction
#   3. TabPFN        — Tabular foundation model (no phylogeny)
#   4. TabPFN_GNN    — TabPFN baseline + graph autoencoder correction
#
# Scenarios 1–8: BM-generated data (BACE sim_bace), all 5 trait types.
# Scenarios 9–12: non-BM data (OU, regime shift, nonlinear).
#
# Requirements:
#   - pigauto (loaded or installed)
#   - BACE   (devtools::load_all from local clone)
#   - torch, ape, Rphylopars
#   - ggplot2 (optional, for figures)
#
# Run from the pigauto package root:
#   source("script/benchmark_simulation.R")
# ============================================================================


# ---- Section 0: Configuration -----------------------------------------------

N_REPS        <- 3L      # replicates per scenario (increase to 10+ for real study)
MISSING_FRAC  <- 0.25    # fraction of cells held out for evaluation
EPOCHS        <- 500L    # GNN training epochs (increase to 2000 for real study)
EVAL_EVERY    <- 50L     # validation evaluation interval
K_EIGEN       <- 8L      # spectral node features
OUTPUT_CSV    <- "script/benchmark_results.csv"
USE_TABPFN    <- TRUE    # include TabPFN methods (skipped if unavailable)


# ---- Section 1: Load packages ------------------------------------------------

library(pigauto)
library(ape)

cat("Loading BACE via devtools::load_all...\n")
devtools::load_all("BACE/")

cat("torch available:", torch::torch_is_installed(), "\n")

# Check TabPFN availability
TABPFN_OK <- FALSE
if (USE_TABPFN) {
  TABPFN_OK <- tryCatch({
    requireNamespace("reticulate", quietly = TRUE) &&
      { reticulate::import("tabimpute.interface"); TRUE }
  }, error = function(e) FALSE)
}
cat("TabPFN available:", TABPFN_OK, "\n")


# ---- Section 2: adapt_sim_for_pigauto() -------------------------------------
#
# sim_bace returns a data.frame with a `species` column and potentially
# multiple rows per species (n_cases >= n_species).  pigauto needs:
#   - one row per species (species-level means for duplicated species)
#   - species as rownames, no species column
#   - proper R types: factor (binary/categorical), ordered (ordinal),
#     integer (count), numeric (continuous)
#
# BACE type encoding in complete_data:
#   "gaussian"      -> numeric
#   "poisson"       -> integer
#   "binary"        -> ordered factor with integer levels (0, 1)
#   "thresholdK"    -> ordered factor with integer levels (1, ..., K)
#   "multinomialK"  -> (unordered) factor with LETTERS levels

adapt_sim_for_pigauto <- function(sim_result, scenario) {
  cd   <- sim_result$complete_data
  tree <- sim_result$tree

  # All variable types: response first, then predictors
  all_types <- c(scenario$response_type, scenario$predictor_types)
  var_names <- sim_result$params$var_names  # e.g. c("y", "x1", "x2", ...)
  trait_cols <- var_names  # columns to keep as traits

  # -- Aggregate to one row per species (species-level mean / mode) -----------
  spp <- cd$species
  unique_spp <- sort(unique(spp))

  agg <- data.frame(row.names = unique_spp)

  for (k in seq_along(trait_cols)) {
    col_name <- trait_cols[k]
    col_type <- all_types[k]
    vals     <- cd[[col_name]]

    if (col_type == "gaussian") {
      # Mean per species
      sp_vals <- tapply(as.numeric(vals), spp, mean)
      agg[[col_name]] <- as.numeric(sp_vals[unique_spp])

    } else if (col_type == "poisson") {
      # Round of mean per species (keep as integer)
      sp_vals <- tapply(as.numeric(vals), spp, mean)
      agg[[col_name]] <- as.integer(round(sp_vals[unique_spp]))

    } else if (col_type == "binary") {
      # BACE outputs ordered factor with levels "0", "1"
      # pigauto expects unordered factor with 2 levels
      num_vals <- as.integer(as.character(vals))
      sp_vals  <- tapply(num_vals, spp, function(x) {
        tb <- table(x)
        as.integer(names(tb)[which.max(tb)])
      })
      int_agg <- sp_vals[unique_spp]
      agg[[col_name]] <- factor(int_agg, levels = c(0L, 1L))

    } else if (grepl("^multinomial", col_type)) {
      # BACE outputs unordered factor with LETTERS levels -- mode per species
      sp_vals <- tapply(as.character(vals), spp, function(x) {
        tb <- table(x)
        names(tb)[which.max(tb)]
      })
      levs <- levels(vals)
      if (is.null(levs)) levs <- sort(unique(as.character(vals)))
      agg[[col_name]] <- factor(sp_vals[unique_spp], levels = levs)

    } else if (grepl("^threshold", col_type)) {
      # BACE outputs ordered factor with integer levels (1, ..., K)
      # pigauto expects ordered factor
      num_vals <- as.integer(as.character(vals))
      sp_vals  <- tapply(num_vals, spp, function(x) {
        tb <- table(x)
        as.integer(names(tb)[which.max(tb)])
      })
      int_agg <- sp_vals[unique_spp]
      n_cats  <- as.integer(gsub("threshold", "", col_type))
      levs    <- as.character(seq_len(n_cats))
      agg[[col_name]] <- ordered(as.character(int_agg), levels = levs)
    }
  }

  # -- Prune tree to species present in data -----------------------------------
  keep <- intersect(tree$tip.label, unique_spp)
  if (length(keep) < length(unique_spp)) {
    # Some species in data may not be in tree; drop them
    agg <- agg[rownames(agg) %in% keep, , drop = FALSE]
  }
  drop_tips <- setdiff(tree$tip.label, keep)
  if (length(drop_tips) > 0) {
    tree <- ape::drop.tip(tree, drop_tips)
  }

  list(traits = agg, tree = tree)
}


# ---- Section 3: run_one_replicate() ------------------------------------------

run_one_replicate <- function(scenario, rep_id) {
  sc_id <- scenario$id
  seed  <- rep_id * 1000L + sc_id

  cat(sprintf("\n--- Scenario %s (%s), rep %d, seed %d ---\n",
              sc_id, scenario$description, rep_id, seed))

  # -- 3a. Simulate data --------------------------------------------------------
  sim_type <- scenario$sim_type %||% "BM"

  adapted <- NULL
  if (sim_type == "BM") {
    # BACE sim_bace path (scenarios 1-8)
    n_vars   <- 1L + length(scenario$predictor_types)
    ps_vec   <- rep(scenario$phylo_signal, length.out = n_vars)
    miss_vec <- rep(0, n_vars)

    set.seed(seed)
    sim <- tryCatch(
      sim_bace(
        response_type  = scenario$response_type,
        predictor_types = scenario$predictor_types,
        phylo_signal   = ps_vec,
        n_cases        = scenario$n_species,
        n_species      = scenario$n_species,
        missingness    = miss_vec
      ),
      error = function(e) {
        cat("  [ERROR] sim_bace failed:", conditionMessage(e), "\n")
        NULL
      }
    )
    if (is.null(sim)) return(NULL)

    adapted <- tryCatch(
      adapt_sim_for_pigauto(sim, scenario),
      error = function(e) {
        cat("  [ERROR] adapt_sim_for_pigauto failed:", conditionMessage(e), "\n")
        NULL
      }
    )
  } else {
    # Non-BM simulation path (scenarios 9+)
    set.seed(seed)
    tree <- ape::rtree(scenario$n_species)
    traits <- tryCatch(
      simulate_non_bm(
        tree       = tree,
        n_traits   = scenario$n_traits %||% 4L,
        scenario   = sim_type,
        alpha      = scenario$alpha %||% 2.0,
        sigma      = scenario$sigma %||% 1.0,
        shift_magnitude = scenario$shift_magnitude %||% 2.0,
        seed       = seed
      ),
      error = function(e) {
        cat("  [ERROR] simulate_non_bm failed:", conditionMessage(e), "\n")
        NULL
      }
    )
    if (is.null(traits)) return(NULL)
    adapted <- list(traits = traits, tree = tree)
  }
  if (is.null(adapted)) return(NULL)

  cat(sprintf("  Adapted data: %d species, %d traits\n",
              nrow(adapted$traits), ncol(adapted$traits)))

  # -- 3c. Preprocess -----------------------------------------------------------
  pd <- tryCatch(
    preprocess_traits(adapted$traits, adapted$tree, log_transform = FALSE),
    error = function(e) {
      cat("  [ERROR] preprocess_traits failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(pd)) return(NULL)

  cat(sprintf("  Latent matrix: %d x %d\n", nrow(pd$X_scaled), ncol(pd$X_scaled)))

  # -- 3d. Create splits --------------------------------------------------------
  missing_frac <- scenario$missing_frac %||% MISSING_FRAC
  splits <- make_missing_splits(pd$X_scaled, missing_frac = missing_frac,
                                seed = seed, trait_map = pd$trait_map)

  cat(sprintf("  Val cells: %d, Test cells: %d\n",
              length(splits$val_idx), length(splits$test_idx)))

  # Store truth for evaluation
  truth <- pd$X_scaled

  # -- 3e. Method 1: BM baseline ------------------------------------------------
  cat("  Fitting BM baseline...\n")
  bl <- tryCatch(
    fit_baseline(pd, adapted$tree, splits = splits),
    error = function(e) {
      cat("  [ERROR] fit_baseline failed:", conditionMessage(e), "\n")
      NULL
    }
  )

  eval_bl <- NULL
  if (!is.null(bl)) {
    eval_bl <- tryCatch(
      evaluate_imputation(bl$mu, truth, splits, trait_map = pd$trait_map),
      error = function(e) {
        cat("  [ERROR] evaluate baseline failed:", conditionMessage(e), "\n")
        NULL
      }
    )
  }

  # -- 3f. Method 2: pigauto (BM + GNN) ----------------------------------------
  eval_gnn <- NULL
  tryCatch({
    cat("  Building phylo graph...\n")
    graph <- build_phylo_graph(adapted$tree, k_eigen = K_EIGEN)

    cat(sprintf("  Training pigauto (epochs = %d)...\n", EPOCHS))
    fit <- fit_pigauto(
      data     = pd,
      tree     = adapted$tree,
      splits   = splits,
      graph    = graph,
      baseline = bl,
      epochs   = EPOCHS,
      eval_every = EVAL_EVERY,
      verbose  = FALSE,
      seed     = seed
    )

    cat("  Predicting...\n")
    pred <- predict(fit, return_se = TRUE)

    eval_gnn <- evaluate_imputation(pred, truth, splits)
  }, error = function(e) {
    cat("  [ERROR] pigauto GNN failed:", conditionMessage(e), "\n")
  })

  # -- 3g. Method 3: TabPFN (no phylogeny) --------------------------------------
  eval_tab <- NULL
  if (TABPFN_OK) {
    tryCatch({
      cat("  Fitting TabPFN baseline...\n")
      bl_tab <- fit_baseline_tabpfn(pd, splits = splits)
      eval_tab <- evaluate_imputation(bl_tab$mu, truth, splits,
                                      trait_map = pd$trait_map)
    }, error = function(e) {
      cat("  [ERROR] TabPFN baseline failed:", conditionMessage(e), "\n")
    })
  }

  # -- 3h. Method 4: pigauto (TabPFN + GNN) ------------------------------------
  eval_tab_gnn <- NULL
  if (TABPFN_OK && !is.null(eval_tab)) {
    tryCatch({
      # Build graph if not already done
      if (!exists("graph", inherits = FALSE)) {
        graph <- build_phylo_graph(adapted$tree, k_eigen = K_EIGEN)
      }
      cat(sprintf("  Training pigauto with TabPFN baseline (epochs = %d)...\n",
                  EPOCHS))
      fit_tg <- fit_pigauto(
        data     = pd,
        tree     = adapted$tree,
        splits   = splits,
        graph    = graph,
        baseline = bl_tab,
        epochs   = EPOCHS,
        eval_every = EVAL_EVERY,
        verbose  = FALSE,
        seed     = seed
      )
      cat("  Predicting (TabPFN+GNN)...\n")
      pred_tg <- predict(fit_tg, return_se = TRUE)
      eval_tab_gnn <- evaluate_imputation(pred_tg, truth, splits)
    }, error = function(e) {
      cat("  [ERROR] pigauto TabPFN+GNN failed:", conditionMessage(e), "\n")
    })
  }

  # -- 3i. Collect results -------------------------------------------------------
  collect_results <- function(eval_df, method_name) {
    if (is.null(eval_df)) return(NULL)
    # Only keep test split
    df <- eval_df[eval_df$split == "test", , drop = FALSE]
    if (nrow(df) == 0) return(NULL)

    # Melt metric columns to long format
    metric_cols <- c("rmse", "pearson_r", "coverage_95", "mae",
                     "spearman_rho", "accuracy", "brier")
    rows <- list()
    for (mc in metric_cols) {
      vals <- df[[mc]]
      ok   <- !is.na(vals)
      if (!any(ok)) next
      rows[[length(rows) + 1L]] <- data.frame(
        scenario    = sc_id,
        description = scenario$description,
        rep         = rep_id,
        method      = method_name,
        trait       = df$trait[ok],
        type        = df$type[ok],
        metric      = mc,
        value       = vals[ok],
        stringsAsFactors = FALSE
      )
    }
    if (length(rows) == 0) return(NULL)
    do.call(rbind, rows)
  }

  rbind(
    collect_results(eval_bl,      "BM_baseline"),
    collect_results(eval_gnn,     "pigauto_GNN"),
    collect_results(eval_tab,     "TabPFN"),
    collect_results(eval_tab_gnn, "TabPFN_GNN")
  )
}


# ---- Section 4: Define scenarios ---------------------------------------------

scenarios <- list(
  list(
    id              = 1L,
    description     = "All continuous, low signal",
    response_type   = "gaussian",
    predictor_types = c("gaussian", "gaussian", "gaussian"),
    phylo_signal    = 0.3,
    n_species       = 150L
  ),
  list(
    id              = 2L,
    description     = "All continuous, high signal",
    response_type   = "gaussian",
    predictor_types = c("gaussian", "gaussian", "gaussian"),
    phylo_signal    = 0.7,
    n_species       = 150L
  ),
  list(
    id              = 3L,
    description     = "Mixed (cont+binary+count), low signal",
    response_type   = "gaussian",
    predictor_types = c("gaussian", "binary", "poisson"),
    phylo_signal    = 0.3,
    n_species       = 150L
  ),
  list(
    id              = 4L,
    description     = "Mixed (cont+binary+count), high signal",
    response_type   = "gaussian",
    predictor_types = c("gaussian", "binary", "poisson"),
    phylo_signal    = 0.7,
    n_species       = 150L
  ),
  list(
    id              = 5L,
    description     = "Mixed + categorical",
    response_type   = "gaussian",
    predictor_types = c("gaussian", "multinomial4", "binary"),
    phylo_signal    = 0.3,
    n_species       = 150L
  ),
  list(
    id              = 6L,
    description     = "High missingness (50%)",
    response_type   = "gaussian",
    predictor_types = c("gaussian", "gaussian", "gaussian"),
    phylo_signal    = 0.3,
    n_species       = 150L,
    missing_frac    = 0.50
  ),
  list(
    id              = 7L,
    description     = "Near-zero phylogenetic signal",
    response_type   = "gaussian",
    predictor_types = c("gaussian", "gaussian", "gaussian"),
    phylo_signal    = 0.05,
    n_species       = 150L
  ),
  list(
    id              = 8L,
    description     = "All 5 types",
    response_type   = "poisson",
    predictor_types = c("gaussian", "binary", "multinomial3", "threshold3"),
    phylo_signal    = 0.4,
    n_species       = 150L
  ),

  # ---- Non-BM scenarios (9-12): where the GNN should add value ----------------

  list(
    id              = 9L,
    description     = "OU, moderate pull (alpha=2)",
    sim_type        = "OU",
    alpha           = 2.0,
    sigma           = 1.0,
    n_species       = 150L,
    n_traits        = 4L
  ),
  list(
    id              = 10L,
    description     = "OU, strong pull (alpha=5)",
    sim_type        = "OU",
    alpha           = 5.0,
    sigma           = 1.0,
    n_species       = 150L,
    n_traits        = 4L
  ),
  list(
    id              = 11L,
    description     = "Regime shift (bimodal)",
    sim_type        = "regime_shift",
    shift_magnitude = 2.0,
    sigma           = 1.0,
    n_species       = 150L,
    n_traits        = 4L
  ),
  list(
    id              = 12L,
    description     = "Non-linear correlations",
    sim_type        = "nonlinear",
    sigma           = 1.0,
    n_species       = 150L,
    n_traits        = 4L
  )
)


# ---- Section 5: Run all scenarios x replicates --------------------------------

cat("\n")
cat(strrep("=", 70), "\n")
cat("PIGAUTO SIMULATION BENCHMARK\n")
cat(sprintf("Scenarios: %d | Replicates: %d | Epochs: %d | Missing: %.0f%%\n",
            length(scenarios), N_REPS, EPOCHS, MISSING_FRAC * 100))
cat(strrep("=", 70), "\n")

all_results <- list()
counter     <- 0L

for (sc in scenarios) {
  for (rep_id in seq_len(N_REPS)) {
    counter <- counter + 1L
    cat(sprintf("\n[%d/%d] ", counter, length(scenarios) * N_REPS))

    res <- tryCatch(
      run_one_replicate(sc, rep_id),
      error = function(e) {
        cat("  [FATAL] Replicate failed:", conditionMessage(e), "\n")
        NULL
      }
    )

    if (!is.null(res) && nrow(res) > 0) {
      all_results[[length(all_results) + 1L]] <- res
    }

    # Save intermediate progress
    if (length(all_results) > 0) {
      interim <- do.call(rbind, all_results)
      write.csv(interim, OUTPUT_CSV, row.names = FALSE)
    }
  }
}

# Final combined data.frame
if (length(all_results) > 0) {
  results_df <- do.call(rbind, all_results)
  rownames(results_df) <- NULL
  write.csv(results_df, OUTPUT_CSV, row.names = FALSE)
  cat(sprintf("\nResults saved to %s (%d rows)\n", OUTPUT_CSV, nrow(results_df)))
} else {
  results_df <- data.frame()
  cat("\nNo results collected.\n")
}


# ---- Section 6: Summary tables -----------------------------------------------

if (nrow(results_df) > 0) {

  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("SUMMARY: Mean +/- SD per (scenario, method, type, metric)\n")
  cat(strrep("=", 70), "\n")

  summary_df <- aggregate(
    value ~ scenario + description + method + type + metric,
    data = results_df,
    FUN  = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
  )
  # aggregate with c() creates a matrix column; unpack it
  summary_df <- cbind(
    summary_df[, c("scenario", "description", "method", "type", "metric")],
    mean = summary_df$value[, "mean"],
    sd   = summary_df$value[, "sd"]
  )
  summary_df$sd[is.na(summary_df$sd)] <- 0

  summary_df <- summary_df[order(summary_df$scenario, summary_df$type,
                                 summary_df$metric, summary_df$method), ]
  print(summary_df, row.names = FALSE)

  # -- Compute improvement: each method vs BM baseline --------------------------
  cat("\n")
  cat(strrep("-", 70), "\n")
  cat("IMPROVEMENT: % change vs BM baseline (all methods)\n")
  cat("  Positive = method better; for RMSE/MAE/Brier lower is better\n")
  cat(strrep("-", 70), "\n")

  lower_better  <- c("rmse", "mae", "brier")
  higher_better <- c("pearson_r", "accuracy", "spearman_rho")

  bl_rows   <- summary_df[summary_df$method == "BM_baseline", ]
  merge_key <- c("scenario", "type", "metric")

  # Compare each non-baseline method to BM
  other_methods <- setdiff(unique(summary_df$method), "BM_baseline")
  comp_list <- list()
  for (meth in other_methods) {
    m_rows <- summary_df[summary_df$method == meth, ]
    cc <- merge(
      bl_rows[, c(merge_key, "description", "mean")],
      m_rows[, c(merge_key, "mean")],
      by = merge_key, suffixes = c("_bl", "_method")
    )
    cc$method <- meth
    cc$improvement_pct <- ifelse(
      cc$metric %in% lower_better,
      (cc$mean_bl - cc$mean_method) / abs(cc$mean_bl) * 100,
      (cc$mean_method - cc$mean_bl) / pmax(abs(cc$mean_bl), 1e-8) * 100
    )
    comp_list[[length(comp_list) + 1L]] <- cc
  }
  comp <- do.call(rbind, comp_list)
  comp <- comp[order(comp$scenario, comp$type, comp$metric, comp$method), ]

  cat("\n")
  print(comp[, c("scenario", "description", "method", "type", "metric",
                  "mean_bl", "mean_method", "improvement_pct")],
        row.names = FALSE, digits = 3)

} else {
  cat("\nNo results to summarise.\n")
}


# ---- Section 7: Figures (if ggplot2 available) --------------------------------

if (nrow(results_df) > 0 && requireNamespace("ggplot2", quietly = TRUE)) {

  library(ggplot2)

  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("GENERATING FIGURES\n")
  cat(strrep("=", 70), "\n")

  # -- 7a. RMSE by method x scenario for continuous traits -----------------------
  rmse_cont <- results_df[results_df$metric == "rmse" &
                             results_df$type == "continuous", ]

  if (nrow(rmse_cont) > 0) {
    rmse_agg <- aggregate(value ~ scenario + description + method,
                          data = rmse_cont, FUN = mean, na.rm = TRUE)

    p1 <- ggplot(rmse_agg,
                 aes(x = factor(scenario), y = value, fill = method)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      labs(x = "Scenario", y = "Test RMSE (latent scale)",
           title = "Continuous traits: RMSE by method and scenario",
           fill = "Method") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    print(p1)

    ggsave("script/benchmark_rmse_continuous.pdf", p1,
           width = 8, height = 5)
    cat("  Saved: script/benchmark_rmse_continuous.pdf\n")
  }

  # -- 7b. Accuracy by method x scenario for binary/categorical traits -----------
  acc_disc <- results_df[results_df$metric == "accuracy" &
                            results_df$type %in% c("binary", "categorical"), ]

  if (nrow(acc_disc) > 0) {
    acc_agg <- aggregate(value ~ scenario + description + method + type,
                         data = acc_disc, FUN = mean, na.rm = TRUE)

    p2 <- ggplot(acc_agg,
                 aes(x = factor(scenario), y = value, fill = method)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      facet_wrap(~ type) +
      labs(x = "Scenario", y = "Accuracy",
           title = "Discrete traits: Accuracy by method and scenario",
           fill = "Method") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    print(p2)

    ggsave("script/benchmark_accuracy_discrete.pdf", p2,
           width = 8, height = 5)
    cat("  Saved: script/benchmark_accuracy_discrete.pdf\n")
  }

  # -- 7c. Improvement plot (% improvement over BM baseline, all methods) --------
  if (exists("comp") && nrow(comp) > 0) {
    comp_plot <- comp[comp$metric %in% c(lower_better, higher_better), ]

    if (nrow(comp_plot) > 0) {
      comp_plot$label <- paste0("S", comp_plot$scenario, ": ", comp_plot$type)

      p3 <- ggplot(comp_plot,
                   aes(x = reorder(label, improvement_pct),
                       y = improvement_pct, fill = metric)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        facet_wrap(~ method) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
        coord_flip() +
        labs(x = NULL, y = "% improvement over BM baseline",
             title = "Relative improvement by scenario, type, and method",
             fill = "Metric") +
        theme_minimal(base_size = 11) +
        theme(legend.position = "bottom")
      print(p3)

      ggsave("script/benchmark_improvement.pdf", p3,
             width = 10, height = 8)
      cat("  Saved: script/benchmark_improvement.pdf\n")
    }
  }

} else {
  cat("\nSkipping figures (ggplot2 not available).\n")
}


cat("\n")
cat(strrep("=", 70), "\n")
cat("BENCHMARK COMPLETE\n")
cat(strrep("=", 70), "\n")
