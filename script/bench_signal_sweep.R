#!/usr/bin/env Rscript
#
# script/bench_signal_sweep.R
#
# Phase 8 MVP: Pagel's lambda signal-strength sweep × mixed-type traits.
# Tests whether pigauto discriminates over the baseline across the
# phylogenetic-signal range, and where each method starts winning.
#
# Design
#   - Tree:        ape::rcoal(300) — bundled tree300 scale.
#   - Traits:      2 continuous + 1 binary + 1 categorical K=3.
#   - Pagel's λ:   {0.1, 0.3, 0.5, 0.7, 0.9, 1.0}. V_λ = λ·V + (1-λ)·I.
#   - missing:     MCAR at 30% per trait.
#   - Methods:     mean_baseline, pigauto_LP, pigauto_default, pigauto_em5.
#   - Replicates:  3 per cell.
#
# Output
#   script/bench_signal_sweep.rds
#   script/bench_signal_sweep.md
#
# Run:
#   Rscript script/bench_signal_sweep.R

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  library(MASS)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto",
    quiet = TRUE
  )
})

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "bench_signal_sweep.rds")
out_md  <- file.path(here, "script", "bench_signal_sweep.md")

# -------------------------------------------------------------------------
# Config
# -------------------------------------------------------------------------

CONFIG <- list(
  lambdas      = c(0.1, 0.3, 0.5, 0.7, 0.9, 1.0),
  n_species    = 300L,
  miss_frac    = 0.30,
  n_reps       = 3L,
  methods      = c("mean_baseline", "pigauto_LP",
                   "pigauto_default", "pigauto_em5")
)

# -------------------------------------------------------------------------
# Data generator: Pagel's λ on mixed-type traits
# -------------------------------------------------------------------------

# Build an n x n covariance matrix at Pagel's λ on the tree:
#   V_λ = λ · V_BM + (1 − λ) · diag(diag(V_BM))
# Positive-definite for any λ in [0, 1].
lambda_vcv <- function(tree, lambda) {
  V <- ape::vcv(tree)   # BM variance-covariance
  D <- diag(diag(V))
  lambda * V + (1 - lambda) * D
}

# Sample n species × p traits from MVN(0, V_λ) stacked iid across traits.
# Returns an n × p matrix with tip.label rownames.
sim_mvn_traits <- function(tree, V_lambda, p, seed) {
  set.seed(seed)
  n <- length(tree$tip.label)
  L <- t(chol(V_lambda + diag(1e-8, n)))   # ridge for numerical safety
  out <- matrix(0, n, p)
  for (j in seq_len(p)) out[, j] <- as.numeric(L %*% stats::rnorm(n))
  rownames(out) <- tree$tip.label
  out
}

# Simulate 2 continuous + 1 binary + 1 categorical K=3 at Pagel's λ.
# Returns a data.frame with:
#   c1, c2  (numeric, z-scored latent)
#   b1      (factor "0"/"1", threshold at 0)
#   cat1    (factor A/B/C, argmax of 3 latent continua)
sim_lambda_mixed <- function(tree, lambda, seed) {
  V_lambda <- lambda_vcv(tree, lambda)

  # 2 continuous + 1 binary latent + 3 categorical latent = 6 columns
  latent <- sim_mvn_traits(tree, V_lambda, p = 6, seed = seed)
  c1 <- latent[, 1]
  c2 <- latent[, 2]
  b1 <- factor(ifelse(latent[, 3] > 0, "1", "0"), levels = c("0", "1"))
  cat_scores <- latent[, 4:6]
  colnames(cat_scores) <- c("A", "B", "C")
  cat_labels <- apply(cat_scores, 1L, function(r) names(r)[which.max(r)])
  cat1 <- factor(cat_labels, levels = c("A", "B", "C"))

  df <- data.frame(c1 = c1, c2 = c2, b1 = b1, cat1 = cat1,
                    row.names = tree$tip.label)
  df
}

# -------------------------------------------------------------------------
# Missingness injection (MCAR per trait)
# -------------------------------------------------------------------------

inject_mcar <- function(df, miss_frac, seed) {
  set.seed(seed)
  n <- nrow(df)
  mask <- lapply(names(df), function(v) {
    sample.int(n, ceiling(n * miss_frac))
  })
  names(mask) <- names(df)
  for (v in names(df)) df[mask[[v]], v] <- NA
  list(data = df, mask = mask)
}

# -------------------------------------------------------------------------
# Method implementations
# -------------------------------------------------------------------------

# mean_baseline: column mean for continuous, modal class for discrete.
method_mean <- function(df_miss, tree, truth, mask) {
  out <- df_miss
  for (v in names(out)) {
    if (is.factor(out[[v]])) {
      mode <- names(sort(table(out[[v]], useNA = "no"), decreasing = TRUE))[1]
      out[[v]][mask[[v]]] <- factor(mode, levels = levels(out[[v]]))
    } else {
      out[[v]][mask[[v]]] <- mean(out[[v]], na.rm = TRUE)
    }
  }
  out
}

# pigauto_LP: force the legacy label-propagation path (no Level-C).
# We achieve this by temporarily unloading Rphylopars via a stub; the
# cleanest portable way is to call impute() with a dummy data where
# joint_mvn_available() naturally returns FALSE.  Since we can't easily
# stub that in-script, instead we just call impute() with a config that
# sends binary/categorical through LP by setting em_iterations = 0L
# AND stripping categoricals from the joint path — BUT in v0.9.1+ the
# joint path always fires if Rphylopars is available. Alternative:
# compare against the single-trait per-column BM path by using
# simulate_benchmark's "LP-only" flag.  For the MVP we approximate
# pigauto_LP as mean_baseline on discrete + per-column BM on continuous
# — that's the "pre-Level-C" regime.  Not exact; see spec Phase 8.1 for
# the full treatment.
method_pigauto_lp_approx <- function(df_miss, tree, truth, mask) {
  out <- df_miss
  # Continuous cols: per-column BM via pigauto's impute() with single col at a time
  for (v in names(out)) {
    if (!is.factor(out[[v]])) {
      one <- out[, v, drop = FALSE]
      res <- tryCatch(
        pigauto::impute(one, tree, epochs = 50L, verbose = FALSE,
                         n_imputations = 1L, seed = 1L,
                         missing_frac = 0),
        error = function(e) NULL
      )
      if (!is.null(res)) {
        out[[v]] <- res$completed[[v]]
      }
    }
  }
  # Discrete: mode (legacy LP-on-single-trait is effectively this)
  for (v in names(out)) {
    if (is.factor(out[[v]])) {
      mode <- names(sort(table(out[[v]], useNA = "no"), decreasing = TRUE))[1]
      out[[v]][mask[[v]]] <- factor(mode, levels = levels(out[[v]]))
    }
  }
  out
}

# pigauto_default: em_iterations = 0L.
# n_imputations=20 gives MC-dropout draws for coverage computation; the
# fitted res is attached as an attribute so eval_cell can read conformal
# intervals + imputed_datasets.
method_pigauto_default <- function(df_miss, tree, truth, mask, seed) {
  res <- pigauto::impute(df_miss, tree, epochs = 300L, verbose = FALSE,
                           n_imputations = 20L, seed = seed,
                           missing_frac = 0.25, em_iterations = 0L)
  out <- res$completed
  attr(out, "res_obj") <- res
  out
}

# pigauto_em5: em_iterations = 5L (Phase 6 diagonal).
method_pigauto_em5 <- function(df_miss, tree, truth, mask, seed) {
  res <- pigauto::impute(df_miss, tree, epochs = 300L, verbose = FALSE,
                           n_imputations = 20L, seed = seed,
                           missing_frac = 0.25, em_iterations = 5L)
  out <- res$completed
  attr(out, "res_obj") <- res
  out
}

# -------------------------------------------------------------------------
# Evaluation — RMSE / Pearson r / accuracy + both coverage types
# -------------------------------------------------------------------------

eval_cell <- function(truth, completed, mask) {
  res_obj <- attr(completed, "res_obj")
  lo      <- if (!is.null(res_obj)) res_obj$prediction$conformal_lower else NULL
  hi      <- if (!is.null(res_obj)) res_obj$prediction$conformal_upper else NULL
  mi_list <- if (!is.null(res_obj)) res_obj$prediction$imputed_datasets else NULL

  rows <- list()
  for (v in names(truth)) {
    idx <- mask[[v]]
    if (length(idx) == 0L) next
    t_v <- truth[[v]][idx]
    c_v <- completed[[v]][idx]
    if (is.factor(t_v)) {
      acc <- mean(as.character(c_v) == as.character(t_v), na.rm = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(
        trait = v, type = "discrete",
        metric = "accuracy", value = acc
      )
    } else {
      rmse <- sqrt(mean((t_v - c_v)^2, na.rm = TRUE))
      pear <- suppressWarnings(stats::cor(t_v, c_v, use = "complete.obs"))
      rows[[length(rows) + 1L]] <- data.frame(
        trait = v, type = "continuous",
        metric = "rmse", value = rmse
      )
      rows[[length(rows) + 1L]] <- data.frame(
        trait = v, type = "continuous",
        metric = "pearson_r", value = pear
      )
      t_num <- as.numeric(t_v)
      if (!is.null(lo) && !is.null(hi) && v %in% colnames(lo)) {
        lo_v <- lo[idx, v]; hi_v <- hi[idx, v]
        valid <- is.finite(lo_v) & is.finite(hi_v) & is.finite(t_num)
        if (any(valid)) {
          hits <- t_num[valid] >= lo_v[valid] & t_num[valid] <= hi_v[valid]
          rows[[length(rows) + 1L]] <- data.frame(
            trait = v, type = "continuous",
            metric = "coverage95_conformal", value = mean(hits))
        }
      }
      if (!is.null(mi_list) && length(mi_list) > 1L && v %in% names(mi_list[[1]])) {
        draws_mat <- vapply(mi_list, function(d) as.numeric(d[idx, v]),
                             numeric(length(idx)))
        if (!is.matrix(draws_mat))
          draws_mat <- matrix(draws_mat, ncol = length(mi_list))
        if (nrow(draws_mat) == length(idx) && ncol(draws_mat) > 1L) {
          q_lo <- apply(draws_mat, 1L, stats::quantile, probs = 0.025, na.rm = TRUE)
          q_hi <- apply(draws_mat, 1L, stats::quantile, probs = 0.975, na.rm = TRUE)
          valid <- is.finite(q_lo) & is.finite(q_hi) & is.finite(t_num)
          if (any(valid)) {
            hits <- t_num[valid] >= q_lo[valid] & t_num[valid] <= q_hi[valid]
            rows[[length(rows) + 1L]] <- data.frame(
              trait = v, type = "continuous",
              metric = "coverage95_mcdropout", value = mean(hits))
          }
        }
      }
    }
  }
  do.call(rbind, rows)
}

# -------------------------------------------------------------------------
# Main loop
# -------------------------------------------------------------------------

cat(sprintf("Phase 8 signal sweep: %d λ × %d reps × %d methods = %d cells\n",
            length(CONFIG$lambdas), CONFIG$n_reps, length(CONFIG$methods),
            length(CONFIG$lambdas) * CONFIG$n_reps * length(CONFIG$methods)))

results <- list()
script_start <- proc.time()[["elapsed"]]
for (rep_i in seq_len(CONFIG$n_reps)) {
  set.seed(rep_i)
  tree <- ape::rcoal(CONFIG$n_species,
                      tip.label = paste0("sp", seq_len(CONFIG$n_species)))
  for (lam in CONFIG$lambdas) {
    sim_seed <- 1000L * rep_i + as.integer(lam * 10)
    df_truth <- sim_lambda_mixed(tree, lam, seed = sim_seed)
    miss <- inject_mcar(df_truth, CONFIG$miss_frac, seed = sim_seed + 1)
    df_miss <- miss$data
    mask    <- miss$mask

    for (meth in CONFIG$methods) {
      t0 <- proc.time()[["elapsed"]]
      completed <- tryCatch(
        switch(meth,
          mean_baseline    = method_mean(df_miss, tree, df_truth, mask),
          pigauto_LP       = method_pigauto_lp_approx(df_miss, tree, df_truth, mask),
          pigauto_default  = method_pigauto_default(df_miss, tree, df_truth, mask, sim_seed + 2),
          pigauto_em5      = method_pigauto_em5(df_miss, tree, df_truth, mask, sim_seed + 2)
        ),
        error = function(e) { message("  ", meth, " @ λ=", lam, ": ", conditionMessage(e)); NULL }
      )
      wall <- proc.time()[["elapsed"]] - t0

      if (is.null(completed)) next

      ev <- eval_cell(df_truth, completed, mask)
      ev$method <- meth
      ev$lambda <- lam
      ev$rep    <- rep_i
      ev$wall_s <- wall
      results[[length(results) + 1L]] <- ev

      cat(sprintf("  rep=%d λ=%.1f %-18s %.1fs\n", rep_i, lam, meth, wall))
      saveRDS(list(results = do.call(rbind, results),
                    config = CONFIG,
                    script_wall = proc.time()[["elapsed"]] - script_start),
                out_rds)
    }
  }
}

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------

all_res <- do.call(rbind, results)
saveRDS(list(results = all_res, config = CONFIG,
              script_wall = proc.time()[["elapsed"]] - script_start),
         out_rds)

summary_tbl <- aggregate(value ~ method + lambda + trait + metric,
                          data = all_res, FUN = mean)

md <- c(
  "# Phase 8 MVP: signal-strength sweep (Pagel's λ)",
  "",
  sprintf("n_species = %d, miss_frac = %.2f, n_reps = %d",
          CONFIG$n_species, CONFIG$miss_frac, CONFIG$n_reps),
  sprintf("Total wall: %.1f min",
          (proc.time()[["elapsed"]] - script_start) / 60),
  "",
  "## Per-trait, per-λ means (averaged over reps)",
  "",
  "```",
  capture.output(print(summary_tbl, row.names = FALSE, max = 500L)),
  "```"
)
writeLines(md, out_md)
cat("\n=== DONE ===\n")
cat("  rds :", out_rds, "\n")
cat("  md  :", out_md, "\n")
