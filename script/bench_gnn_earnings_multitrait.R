#!/usr/bin/env Rscript
# script/bench_gnn_earnings_multitrait.R
#
# Day 2 of the GNN earnings experiment: tests pigauto's claimed
# architectural advantage of capturing CROSS-TRAIT structure that
# per-trait or linear-cross-trait baselines cannot.
#
# Day 1 (v1/v2) used univariate y = phylo + f(cov) + eps and
# found pigauto loses to phylolm-BLUP on linear/nonlinear and only
# wins modestly on interactive.  That sim does not exercise the
# GNN's main differentiator: cross-trait sharing.
#
# DGP for Day 2:
#
#   - Generate K = 4 correlated continuous traits via multivariate
#     BM on the tree, with cross-trait correlation matrix Sigma.
#   - Add a NONLINEAR cross-trait coupling on trait 4:
#       y4 = phylo_4 + nonlin_couple * sin(2*y1) * exp(0.3 * y2)
#     (the joint-MVN baseline can capture the LINEAR cross-trait
#     correlation but NOT the nonlinear sin*exp interaction.)
#   - 30% MCAR mask on each trait.
#
# Methods compared:
#   1. column_mean per trait
#   2. lm per trait with other traits as predictors  (linear cross-trait)
#   3. Rphylopars joint MVN BLUP                      (linear cross-trait + phylo)
#      (the smart baseline that captures BM cross-trait correlation)
#   4. pigauto with phylogeny + multi-trait imputation
#
# Pre-registered decision rule:
#   pigauto earns its keep on Day 2 if
#     mean(RMSE_pigauto, on trait 4)
#       < 0.90 * mean(RMSE_Rphylopars_BLUP, on trait 4)
#   AND pooled SE-bounded.  Trait 4 is the focal one because that's
#   where the nonlinear coupling lives -- linear methods CAN'T fit it
#   even with cross-trait info.
#
# Sweep:
#   nonlin_couple in {0, 0.5, 1.0}    (strength of nonlinear coupling)
#   tree:           tree300 only      (300 species, fast)
#   reps:           5
#
# This tests whether the GNN can do something Rphylopars truly cannot.
# If pigauto wins here, the architecture has a real advantage.

suppressPackageStartupMessages({
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH", unset = "")
  if (nzchar(pkg_path)) devtools::load_all(pkg_path, quiet = TRUE)
  else library(pigauto)
  library(ape)
  if (!requireNamespace("Rphylopars", quietly = TRUE))
    stop("Rphylopars required for Day 2 multitrait sim")
  library(Rphylopars)
})
options(warn = -1L)

SEED      <- 2026L
MISS_FRAC <- 0.30
EPOCHS    <- as.integer(Sys.getenv("PIGAUTO_BENCH_EPOCHS", "200"))
N_IMP     <- 10L
N_REPS    <- 5L
COUPLES   <- c(0.0, 0.5, 1.0)
N_TRAITS  <- 4L
TAG       <- Sys.getenv("PIGAUTO_BENCH_TAG", "multitrait_tree300")

cat_line <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# Load tree (default tree300)
tree_rda <- Sys.getenv("PIGAUTO_TREE_RDA", "")
tree_obj <- Sys.getenv("PIGAUTO_TREE_OBJ", "tree300")
if (nzchar(tree_rda)) { load(tree_rda); tree <- get(tree_obj) }
else { data(list = tree_obj); tree <- get(tree_obj) }
set.seed(SEED)
if (!ape::is.binary(tree)) tree <- ape::multi2di(tree, random = TRUE)
if (is.null(tree$edge.length))
  tree <- ape::compute.brlen(tree, method = "Grafen")
zero_e <- tree$edge.length <= 0
if (any(zero_e)) tree$edge.length[zero_e] <- 1e-6
n <- length(tree$tip.label)
cat_line(sprintf("multitrait sim: tree=%s, n=%d, K=%d traits, couples=%s, reps=%d",
                  TAG, n, N_TRAITS,
                  paste(COUPLES, collapse=","), N_REPS))

# DGP -----------------------------------------------------------------
sim_multivariate_bm <- function(tree, K, seed) {
  set.seed(seed)
  # Random K x K cross-trait correlation matrix Sigma
  A <- matrix(stats::rnorm(K * K), K, K)
  Sigma <- crossprod(A) / K + diag(K) * 0.5
  Sigma <- stats::cov2cor(Sigma)  # standardize to correlations
  # Multivariate BM via independent eigen-decomposition
  R <- stats::cov2cor(ape::vcv(tree))
  L_R <- t(chol(R + 1e-6 * diag(nrow(R))))
  L_S <- t(chol(Sigma + 1e-6 * diag(K)))
  Z <- matrix(stats::rnorm(n * K), n, K)
  Y <- L_R %*% Z %*% t(L_S)   # n x K  with row-cov R, col-cov Sigma
  colnames(Y) <- paste0("y", 1:K)
  rownames(Y) <- tree$tip.label
  list(Y = Y, Sigma = Sigma)
}

apply_nonlinear_coupling <- function(Y, couple) {
  if (couple == 0) return(Y)
  Y_new <- Y
  # y4 += couple * sin(2*y1) * exp(0.3*y2)  (nonlinear cross-trait)
  delta <- couple * sin(2 * Y[, 1]) * exp(0.3 * Y[, 2])
  Y_new[, 4] <- Y[, 4] + delta
  Y_new
}

# Methods -------------------------------------------------------------
fit_column_mean <- function(Y, mask) {
  out <- Y
  for (j in seq_len(ncol(Y))) {
    out[mask[, j], j] <- mean(Y[!mask[, j], j], na.rm = TRUE)
  }
  out
}

fit_lm_crosstrait <- function(Y, mask) {
  # For each trait j, regress on the OTHER traits (using observed cells where
  # other traits are observed too).
  out <- Y
  K <- ncol(Y)
  for (j in seq_len(K)) {
    obs_j <- !mask[, j]
    other <- setdiff(seq_len(K), j)
    df_o <- as.data.frame(Y[obs_j, , drop = FALSE])
    fit <- tryCatch(
      stats::lm(stats::as.formula(sprintf("%s ~ %s", colnames(Y)[j],
                                              paste(colnames(Y)[other], collapse=" + "))),
                  data = df_o),
      error = function(e) NULL)
    if (is.null(fit)) {
      out[mask[, j], j] <- mean(Y[obs_j, j], na.rm = TRUE)
    } else {
      df_h <- as.data.frame(Y)
      pred <- as.numeric(stats::predict(fit, newdata = df_h))
      out[mask[, j], j] <- pred[mask[, j]]
    }
  }
  out
}

fit_rphylopars_blup <- function(Y, mask, tree) {
  # Build Rphylopars-format data: species + traits, NA where masked
  Y_obs <- Y
  for (j in seq_len(ncol(Y))) Y_obs[mask[, j], j] <- NA
  df <- data.frame(species = rownames(Y), Y_obs, stringsAsFactors = FALSE)
  fit <- tryCatch(
    Rphylopars::phylopars(df, tree = tree, model = "BM",
                            phylo_correlated = TRUE,
                            pheno_correlated = TRUE,
                            REML = TRUE),
    error = function(e) NULL)
  if (is.null(fit)) {
    cat_line("  Rphylopars failed; falling back to column mean")
    return(fit_column_mean(Y, mask))
  }
  out <- fit$anc_recon[rownames(Y), colnames(Y)]
  # Replace observed cells with truth (Rphylopars predicts both observed and missing)
  for (j in seq_len(ncol(Y))) out[!mask[, j], j] <- Y[!mask[, j], j]
  as.matrix(out)
}

fit_pigauto_multitrait <- function(Y, mask, tree) {
  Y_obs <- Y
  for (j in seq_len(ncol(Y))) Y_obs[mask[, j], j] <- NA
  df <- as.data.frame(Y_obs)
  rownames(df) <- rownames(Y)
  res <- tryCatch(
    pigauto::impute(df, tree, epochs = EPOCHS, n_imputations = N_IMP,
                      verbose = FALSE, seed = SEED, safety_floor = TRUE),
    error = function(e) { cat_line("  pigauto error: ", conditionMessage(e)); NULL })
  if (is.null(res)) return(fit_column_mean(Y, mask))
  as.matrix(res$completed)
}

# Eval ---------------------------------------------------------------
rmse_per_trait <- function(pred, truth, mask) {
  K <- ncol(truth)
  out <- numeric(K)
  for (j in seq_len(K)) {
    idx <- which(mask[, j] & is.finite(pred[, j]) & is.finite(truth[, j]))
    if (length(idx) < 5L) { out[j] <- NA_real_; next }
    out[j] <- sqrt(mean((pred[idx, j] - truth[idx, j])^2))
  }
  names(out) <- colnames(truth)
  out
}

# Main sweep ---------------------------------------------------------
results <- list()
for (couple in COUPLES) for (rep_id in seq_len(N_REPS)) {
  rep_seed <- SEED + 1000L * rep_id + round(couple * 100)
  cat_line(sprintf("=== couple=%.1f rep=%d (seed=%d) ===",
                    couple, rep_id, rep_seed))

  sim <- sim_multivariate_bm(tree, N_TRAITS, rep_seed)
  Y_true <- apply_nonlinear_coupling(sim$Y, couple)

  # MCAR mask
  set.seed(rep_seed + 99L)
  mask <- matrix(FALSE, n, N_TRAITS, dimnames = list(NULL, colnames(Y_true)))
  for (j in seq_len(N_TRAITS)) {
    mask[sample.int(n, round(MISS_FRAC * n)), j] <- TRUE
  }

  preds <- list()
  walls <- numeric()

  t0 <- proc.time()[["elapsed"]]
  preds$column_mean <- fit_column_mean(Y_true, mask)
  walls["column_mean"] <- proc.time()[["elapsed"]] - t0

  t0 <- proc.time()[["elapsed"]]
  preds$lm_crosstrait <- fit_lm_crosstrait(Y_true, mask)
  walls["lm_crosstrait"] <- proc.time()[["elapsed"]] - t0

  t0 <- proc.time()[["elapsed"]]
  preds$rphylopars_blup <- fit_rphylopars_blup(Y_true, mask, tree)
  walls["rphylopars_blup"] <- proc.time()[["elapsed"]] - t0
  cat_line(sprintf("  rphylopars_blup %.1fs", walls["rphylopars_blup"]))

  t0 <- proc.time()[["elapsed"]]
  preds$pigauto_multitrait <- fit_pigauto_multitrait(Y_true, mask, tree)
  walls["pigauto_multitrait"] <- proc.time()[["elapsed"]] - t0
  cat_line(sprintf("  pigauto_multitrait %.1fs", walls["pigauto_multitrait"]))

  for (m in names(preds)) {
    r <- rmse_per_trait(preds[[m]], Y_true, mask)
    for (j in seq_along(r)) {
      results[[length(results) + 1L]] <- data.frame(
        couple = couple, rep_id = rep_id, method = m,
        trait_j = j, trait_name = names(r)[j],
        rmse = r[[j]], wall_s = walls[[m]],
        stringsAsFactors = FALSE)
    }
  }
}

results_df <- do.call(rbind, results)
out_rds <- sprintf("script/bench_gnn_earnings_%s.rds", TAG)
out_md  <- sprintf("script/bench_gnn_earnings_%s.md",  TAG)
saveRDS(list(results = results_df,
              meta = list(tag = TAG, n = n, K = N_TRAITS,
                            couples = COUPLES, n_reps = N_REPS,
                            epochs = EPOCHS, seed = SEED, miss_frac = MISS_FRAC,
                            decision_rule = "pigauto wins on Day 2 if RMSE_pigauto on trait 4 < 0.9 * RMSE_Rphylopars BLUP and gap > 1.96 SE")),
         out_rds, compress = "xz")

# Summary
agg <- aggregate(rmse ~ method + couple + trait_name,
                  data = results_df,
                  FUN = function(v) c(mean = mean(v, na.rm=TRUE),
                                       se   = stats::sd(v, na.rm=TRUE) / sqrt(sum(!is.na(v)))))
agg <- do.call(data.frame, agg)
colnames(agg)[ncol(agg) - 1L] <- "mean_rmse"
colnames(agg)[ncol(agg)]      <- "se_rmse"
agg <- agg[order(agg$trait_name, agg$couple, agg$method), ]

writeLines(c(
  sprintf("# Day 2: multitrait GNN earnings (%s, n=%d, K=%d)", TAG, n, N_TRAITS),
  "",
  "Per (couple, trait, method): mean RMSE +/- SE on held-out cells.",
  "Decisive trait: y4 (where the nonlinear sin*exp coupling lives).",
  "",
  "```",
  capture.output(print(agg, row.names = FALSE)),
  "```"), out_md)

cat_line("=== DONE ===")
cat_line("  rds: ", out_rds)
cat_line("  md : ", out_md)
