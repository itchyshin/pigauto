#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# pigauto vs TabPFN benchmark
#
# Compares four methods on the standard simulation battery:
#   1. baseline           - phylogenetic label propagation (fit_baseline)
#   2. pigauto            - GNN on the phylogenetic baseline (current default)
#   3. tabpfn             - TabPFN (tabimpute) alone, no phylogeny
#   4. pigauto_tabpfn     - GNN on top of TabPFN as baseline
#
# The goal is to answer: does a tabular foundation model plus phylogenetic
# information beat either alone?  Does it beat the simple phylogenetic
# label-propagation baseline?
#
# Architecture note: R torch (for the GNN) and Python torch (bundled with
# tabimpute) are ABI-incompatible in the same process.  TabPFN is therefore
# invoked as a subprocess via script/run_tabpfn.py, which keeps the two
# libtorch namespaces completely isolated.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
})

# ---- Config ---------------------------------------------------------------
N_SPECIES  <- 150L
SCENARIOS  <- c("BM", "OU", "nonlinear", "regime_shift", "mixed")
MISS_FRAC  <- 0.25
N_REPS     <- 2L
EPOCHS     <- 600L

TABPFN_PY   <- "/Users/z3437171/.virtualenvs/r-tabpfn/bin/python"
TABPFN_SCR  <- "script/run_tabpfn.py"

METHODS <- c("baseline", "pigauto", "tabpfn", "pigauto_tabpfn")

# ---- TabPFN via subprocess ------------------------------------------------
# Mirrors fit_baseline_tabpfn() but isolates Python torch in a subprocess.
tabpfn_baseline <- function(data, splits) {
  X <- data$X_scaled
  # Mask val + test cells so TabPFN only sees the training cells
  if (!is.null(splits)) {
    X[splits$val_idx]  <- NA
    X[splits$test_idx] <- NA
  }

  tmp_in  <- tempfile(fileext = ".csv")
  tmp_out <- tempfile(fileext = ".csv")
  on.exit({ unlink(tmp_in); unlink(tmp_out) }, add = TRUE)

  write.csv(as.data.frame(X), tmp_in, row.names = FALSE)

  status <- suppressWarnings(
    system2(TABPFN_PY, c(TABPFN_SCR, tmp_in, tmp_out),
            stdout = TRUE, stderr = TRUE)
  )
  if (!file.exists(tmp_out)) {
    stop("TabPFN subprocess failed:\n", paste(status, collapse = "\n"))
  }

  mu <- as.matrix(read.csv(tmp_out))
  dimnames(mu) <- list(data$species_names, data$trait_names)

  # Approximate SE as per-trait residual SD on training cells
  obs_train <- !is.na(data$X_scaled)
  if (!is.null(splits)) {
    obs_train[splits$val_idx]  <- FALSE
    obs_train[splits$test_idx] <- FALSE
  }
  se <- matrix(NA_real_, nrow(mu), ncol(mu), dimnames = dimnames(mu))
  for (j in seq_len(ncol(mu))) {
    keep <- obs_train[, j] & is.finite(data$X_scaled[, j]) & is.finite(mu[, j])
    if (sum(keep) > 1L) {
      se[, j] <- stats::sd(data$X_scaled[keep, j] - mu[keep, j])
    }
  }

  list(mu = mu, se = se)
}

# ---- Evaluate a raw prediction matrix in the same format as evaluate() ----
eval_pred_matrix <- function(pred_mat, data, splits, method_name) {
  trait_map <- data$trait_map
  truth     <- data$X_scaled
  n         <- nrow(truth)
  test_idx  <- splits$test_idx
  row_i     <- ((test_idx - 1L) %% n) + 1L
  col_j     <- ceiling(test_idx / n)

  eval_test_cells(
    pred_latent  = pred_mat,
    pred_obj     = NULL,
    truth_latent = truth,
    row_i        = row_i,
    col_j        = col_j,
    trait_map    = trait_map,
    n            = n,
    method       = method_name,
    conformal_scores = NULL
  )
}

# ---- Main loop ------------------------------------------------------------
all_scenarios <- c("BM", "OU", "regime_shift", "nonlinear", "mixed")
all_results   <- list()
idx <- 0L

cat("\n=== pigauto vs TabPFN benchmark ===\n")
cat(sprintf("n_species=%d, reps=%d, missing=%.2f, epochs=%d\n\n",
            N_SPECIES, N_REPS, MISS_FRAC, EPOCHS))

for (scen in SCENARIOS) {
  cat(sprintf("=== Scenario: %s ===\n", scen))

  for (rep in seq_len(N_REPS)) {
    rep_seed <- rep * 100L + match(scen, all_scenarios)
    cat(sprintf("  Rep %d/%d (seed %d)\n", rep, N_REPS, rep_seed))

    set.seed(rep_seed)
    tree <- ape::rtree(N_SPECIES)
    df <- if (scen == "mixed") {
      simulate_mixed_traits(tree, seed = rep_seed)
    } else if (scen == "BM") {
      simulate_bm_traits(tree, 3L, seed = rep_seed)
    } else {
      simulate_non_bm(tree, n_traits = 3L, scenario = scen, seed = rep_seed)
    }

    pd  <- preprocess_traits(df, tree, log_transform = FALSE)
    spl <- make_missing_splits(pd$X_scaled, missing_frac = MISS_FRAC,
                               seed = rep_seed, trait_map = pd$trait_map)

    # (1) Phylogenetic baseline
    bl_phy <- fit_baseline(pd, tree, splits = spl)

    # (3) TabPFN baseline (subprocess — before loading R torch would be safer
    #     but we can call a subprocess even after torch is loaded)
    bl_tab <- tryCatch(
      tabpfn_baseline(pd, spl),
      error = function(e) {
        cat("    tabpfn error: ", conditionMessage(e), "\n")
        NULL
      }
    )

    # (2) pigauto on phylogenetic baseline
    fit_phy <- fit_pigauto(pd, tree, splits = spl, baseline = bl_phy,
                           epochs = EPOCHS, verbose = FALSE, seed = rep_seed)
    ev_phy  <- evaluate(fit_phy, data = pd, splits = spl)   # baseline + pigauto

    ev_tab <- NULL
    ev_pig_tab <- NULL
    if (!is.null(bl_tab)) {
      ev_tab <- eval_pred_matrix(bl_tab$mu, pd, spl, method_name = "tabpfn")

      # (4) pigauto on TabPFN baseline (stacked)
      fit_tab <- tryCatch(
        fit_pigauto(pd, tree, splits = spl, baseline = bl_tab,
                    epochs = EPOCHS, verbose = FALSE, seed = rep_seed),
        error = function(e) {
          cat("    pigauto_tabpfn error: ", conditionMessage(e), "\n")
          NULL
        }
      )
      if (!is.null(fit_tab)) {
        ev_stack <- evaluate(fit_tab, data = pd, splits = spl)
        # Drop the "baseline" rows (which would duplicate tabpfn rows) and
        # rename pigauto -> pigauto_tabpfn.
        ev_pig_tab <- ev_stack[ev_stack$method == "pigauto", ]
        ev_pig_tab$method <- "pigauto_tabpfn"
      }
    }

    rep_df <- rbind(ev_phy, ev_tab, ev_pig_tab)
    rep_df$scenario <- scen
    rep_df$rep      <- rep

    idx <- idx + 1L
    all_results[[idx]] <- rep_df
  }
}

results <- do.call(rbind, all_results)
rownames(results) <- NULL

saveRDS(list(results = results,
             scenarios = SCENARIOS,
             n_species = N_SPECIES,
             n_reps = N_REPS,
             methods = METHODS),
        "script/bench_tabpfn.rds")

# ---- Report ----------------------------------------------------------------
cat("\n\n=== Summary ===\n")
cat(strrep("=", 72), "\n")

for (scen in SCENARIOS) {
  cat("\nScenario:", scen, "\n")
  cat(strrep("-", 60), "\n")
  sub <- results[results$scenario == scen, ]

  cont_rmse <- sub[sub$metric == "rmse" & sub$type == "continuous", ]
  if (nrow(cont_rmse) > 0L) {
    agg <- aggregate(value ~ method, data = cont_rmse, FUN = mean)
    agg <- agg[match(METHODS, agg$method), ]
    agg <- agg[!is.na(agg$method), ]
    for (i in seq_len(nrow(agg))) {
      cat(sprintf("  %-18s  RMSE: %.4f\n", agg$method[i], agg$value[i]))
    }
  }

  acc <- sub[sub$metric == "accuracy" &
               sub$type %in% c("binary", "categorical"), ]
  if (nrow(acc) > 0L) {
    agg <- aggregate(value ~ method, data = acc, FUN = mean)
    agg <- agg[match(METHODS, agg$method), ]
    agg <- agg[!is.na(agg$method), ]
    for (i in seq_len(nrow(agg))) {
      cat(sprintf("  %-18s  Accuracy: %.1f%%\n",
                  agg$method[i], 100 * agg$value[i]))
    }
  }
}

cat("\n", strrep("=", 72), "\n", sep = "")
cat("Relative RMSE improvement vs phylogenetic baseline\n")
cat(strrep("-", 60), "\n")
cat(sprintf("  %-15s", "scenario"))
for (m in setdiff(METHODS, "baseline")) cat(sprintf(" %16s", m))
cat("\n")

for (scen in SCENARIOS) {
  sub <- results[results$scenario == scen &
                   results$metric == "rmse" &
                   results$type == "continuous", ]
  if (nrow(sub) == 0L) next
  bl_mean <- mean(sub$value[sub$method == "baseline"], na.rm = TRUE)
  cat(sprintf("  %-15s", scen))
  for (m in setdiff(METHODS, "baseline")) {
    m_mean <- mean(sub$value[sub$method == m], na.rm = TRUE)
    if (!is.finite(m_mean)) {
      cat(sprintf(" %16s", "-"))
    } else {
      pct <- 100 * (bl_mean - m_mean) / bl_mean
      cat(sprintf(" %+14.1f%%", pct))
    }
  }
  cat("\n")
}

cat("\nAccuracy change (percentage points) vs phylogenetic baseline\n")
cat(strrep("-", 60), "\n")
cat(sprintf("  %-15s", "scenario"))
for (m in setdiff(METHODS, "baseline")) cat(sprintf(" %16s", m))
cat("\n")

for (scen in SCENARIOS) {
  sub <- results[results$scenario == scen &
                   results$metric == "accuracy" &
                   results$type %in% c("binary", "categorical"), ]
  if (nrow(sub) == 0L) next
  bl_mean <- mean(sub$value[sub$method == "baseline"], na.rm = TRUE)
  cat(sprintf("  %-15s", scen))
  for (m in setdiff(METHODS, "baseline")) {
    m_mean <- mean(sub$value[sub$method == m], na.rm = TRUE)
    if (!is.finite(m_mean)) {
      cat(sprintf(" %16s", "-"))
    } else {
      pp <- 100 * (m_mean - bl_mean)
      cat(sprintf(" %+14.1fpp", pp))
    }
  }
  cat("\n")
}

cat("\nSaved: script/bench_tabpfn.rds\n")
