#!/usr/bin/env Rscript
# script/bench_gnn_earnings_v2.R
#
# Phase 4 of the GNN earnings experiment.  Addresses four flaws
# in v1 (script/bench_gnn_earnings.R) flagged by the simulation
# audit:
#
#  1. v1 used phylolm(model="BM") which assumes zero residual noise --
#     misspecified vs the DGP's eps term.  v2 uses model="lambda"
#     and uses the lambda-fitted covariance in the BLUP step.
#
#  2. v1 had a degenerate cell at alpha=0.4, beta=0.6 (zero noise
#     floor).  v2 caps alpha + beta <= 0.9 so eps variance >= 0.1.
#
#  3. v1 used n_reps=3 -- simulation SE on the same order as the
#     0.10 RMSE-ratio threshold.  v2 uses n_reps=5 and reports
#     mean ± SE per cell.
#
#  4. v1 generated covariates iid (rnorm), so covariates were
#     phylogenetically independent -- best case for pigauto.  v2
#     adds a phylo-structured covariate condition where cov1 is
#     drawn from BM on the same tree (worst case for pigauto:
#     covariate is phylo-redundant).  Real comparative biology data
#     usually sits between these two extremes.
#
# Pre-registered decision rule:
#   pigauto + cov "earns its keep" relative to phylolm-lambda BLUP
#   if mean(RMSE_pigauto_cov_sfT) < 0.90 * mean(RMSE_phylolm_lambda_blup)
#   AND the absolute gap exceeds 1.96 * pooled simulation SE,
#   on the cells with f in {nonlinear, interactive} and beta = 0.4.
#
# Sweep:
#   alpha    in {0.1, 0.4, 0.7}    (phylo signal: weak / mid / strong)
#   beta     in {0, 0.2, 0.4}      (cov signal: none / weak / mid)
#   f_type   in {linear, nonlinear, interactive}
#   cov_type in {iid, phylo_structured}
#   reps     = 5
#   skip the (alpha=0.7, beta=0.4) cell (alpha+beta > 0.9 violates noise floor)
#
# Methods (6, same as v1 but phylolm now uses lambda):
#   1. column_mean
#   2. lm
#   3. phylolm_lambda_blup  (linear cov + lambda-fitted BM mix BLUP)
#   4. pigauto_no_cov
#   5. pigauto_cov_sfT
#   6. pigauto_cov_sfF
#
# Tree controlled by env vars (same as v1).

suppressPackageStartupMessages({
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH", unset = "")
  if (nzchar(pkg_path)) devtools::load_all(pkg_path, quiet = TRUE)
  else library(pigauto)
  library(ape)
  library(phylolm)
})
options(warn = -1L)

SEED         <- 2026L
MISS_FRAC    <- 0.30
EPOCHS       <- as.integer(Sys.getenv("PIGAUTO_BENCH_EPOCHS", "150"))
N_IMP        <- 10L
N_TARGET     <- as.integer(Sys.getenv("PIGAUTO_BENCH_N", "0"))
TAG          <- Sys.getenv("PIGAUTO_BENCH_TAG", "v2_tree300")

ALPHAS       <- c(0.1, 0.4, 0.7)
BETAS        <- c(0.0, 0.2, 0.4)
F_TYPES      <- c("linear", "nonlinear", "interactive")
COV_TYPES    <- c("iid", "phylo_structured")
N_REPS       <- 5L

# valid (alpha, beta) cells: noise floor >= 0.1 always
valid_ab <- function(a, b) (1 - a - b) >= 0.1 - 1e-9

cat_line <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "),
                                ..., "\n", sep = "")

# ---------- load phylogeny (same as v1) ----------
tree_rda <- Sys.getenv("PIGAUTO_TREE_RDA", "")
tree_obj <- Sys.getenv("PIGAUTO_TREE_OBJ", "tree300")
if (nzchar(tree_rda)) {
  load(tree_rda)
  tree_full <- get(tree_obj)
} else {
  data(list = tree_obj)
  tree_full <- get(tree_obj)
}

set.seed(SEED)
if (!ape::is.binary(tree_full)) tree_full <- ape::multi2di(tree_full, random = TRUE)
if (!ape::is.rooted(tree_full)) tree_full <- ape::root.phylo(tree_full,
                                                                outgroup = 1L,
                                                                resolve.root = TRUE)
if (is.null(tree_full$edge.length))
  tree_full <- ape::compute.brlen(tree_full, method = "Grafen")
zero_e <- tree_full$edge.length <= 0
if (any(zero_e)) tree_full$edge.length[zero_e] <- 1e-6

if (N_TARGET > 0L && N_TARGET < length(tree_full$tip.label)) {
  set.seed(SEED + 7L)
  keep <- sample(tree_full$tip.label, N_TARGET)
  tree <- ape::keep.tip(tree_full, keep)
} else {
  tree <- tree_full
}
n <- length(tree$tip.label)
cat_line(sprintf("v2: tree=%s, n=%d, alphas=%s, betas=%s, f=%s, cov=%s, reps=%d",
                  TAG, n,
                  paste(ALPHAS, collapse=","),
                  paste(BETAS,  collapse=","),
                  paste(F_TYPES, collapse=","),
                  paste(COV_TYPES, collapse=","),
                  N_REPS))

# Precompute correlation matrix once
R_full <- stats::cov2cor(ape::vcv(tree))

# ---------- DGP helpers ----------
sim_phylo <- function(tree, seed) {
  set.seed(seed)
  v <- ape::rTraitCont(tree, model = "BM", sigma = 1, root.value = 0)
  v <- as.numeric(v[tree$tip.label])
  (v - mean(v)) / stats::sd(v)
}

sim_cov_iid <- function(n, seed) {
  set.seed(seed + 1L)
  m <- matrix(stats::rnorm(n * 5L), nrow = n)
  colnames(m) <- paste0("cov", 1:5)
  m
}

# Phylo-structured: cov1 is its own BM draw on the same tree (different
# seed, but same tree structure -> sister species share cov values).
# cov2-5 are iid noise.  Standardise each column.
sim_cov_phylo <- function(tree, seed) {
  n <- length(tree$tip.label)
  set.seed(seed + 31L)
  c1 <- ape::rTraitCont(tree, model = "BM", sigma = 1, root.value = 0)
  c1 <- as.numeric(c1[tree$tip.label])
  c1 <- (c1 - mean(c1)) / stats::sd(c1)
  set.seed(seed + 32L)
  rest <- matrix(stats::rnorm(n * 4L), nrow = n)
  m <- cbind(c1, rest)
  colnames(m) <- paste0("cov", 1:5)
  m
}

f_apply <- function(cov_mat, f_type) {
  c1 <- cov_mat[, 1]; c2 <- cov_mat[, 2]
  fx <- switch(f_type,
                linear      = c1,
                nonlinear   = sin(2 * c1) * exp(0.3 * c2),
                interactive = c1 * c2 + 0.5 * c1^2,
                stop("unknown f_type"))
  (fx - mean(fx)) / stats::sd(fx)
}

sim_y <- function(tree, alpha, beta, f_type, cov_type, seed) {
  n <- length(tree$tip.label)
  phylo <- sim_phylo(tree, seed)
  cov_m <- if (cov_type == "iid") sim_cov_iid(n, seed)
           else                    sim_cov_phylo(tree, seed)
  fx    <- f_apply(cov_m, f_type)
  set.seed(seed + 2L)
  eps   <- stats::rnorm(n)
  resid_w <- 1 - alpha - beta
  y <- sqrt(alpha) * phylo + sqrt(beta) * fx + sqrt(resid_w) * eps
  list(y = y, cov = cov_m)
}

# ---------- methods ----------
fit_lm <- function(y, cov_mat, mask) {
  obs <- !mask
  df_o <- data.frame(y = y[obs], cov_mat[obs, , drop = FALSE])
  fit <- stats::lm(y ~ ., data = df_o)
  df_h <- data.frame(cov_mat)
  as.numeric(stats::predict(fit, newdata = df_h))
}

# Phylolm with model="lambda" + BLUP correction using the lambda-fitted
# mixed covariance V_lambda = lambda * R + (1 - lambda) * I.
fit_phylolm_lambda_blup <- function(y, cov_mat, tree, R_full, mask) {
  obs <- !mask
  df_o <- data.frame(y = y[obs], cov_mat[obs, , drop = FALSE])
  rownames(df_o) <- tree$tip.label[obs]
  tree_o <- ape::keep.tip(tree, tree$tip.label[obs])

  fit <- tryCatch(
    phylolm::phylolm(y ~ cov1 + cov2 + cov3 + cov4 + cov5,
                       data = df_o, phy = tree_o, model = "lambda",
                       lower.bound = 0.0, upper.bound = 1.0),
    error = function(e) NULL)
  if (is.null(fit))
    return(rep(mean(y[obs], na.rm = TRUE), length(y)))

  lam <- as.numeric(fit$optpar)        # phylolm stores lambda in $optpar
  if (!is.finite(lam) || lam < 0 || lam > 1) lam <- 1.0   # fall back to BM

  # Fixed-effect prediction for ALL tips
  X_all <- cbind(1, cov_mat)
  beta_hat <- as.numeric(coef(fit))
  yhat_fixed <- as.numeric(X_all %*% beta_hat)

  # Build lambda-mixed covariance matrix (using correlation R, since
  # we'll multiply by sigma2 implicitly through the residual scale)
  o_idx <- which(obs); h_idx <- which(mask); n_o <- length(o_idx)
  V_full <- lam * R_full + (1 - lam) * diag(nrow(R_full))
  V_oo <- V_full[o_idx, o_idx, drop = FALSE]
  V_ho <- V_full[h_idx, o_idx, drop = FALSE]

  resid_o <- y[o_idx] - yhat_fixed[o_idx]
  V_oo_inv <- tryCatch(
    solve(V_oo + 1e-6 * diag(n_o)),
    error = function(e) MASS::ginv(V_oo))

  blup_h <- as.numeric(V_ho %*% (V_oo_inv %*% resid_o))
  pred <- yhat_fixed
  pred[h_idx] <- yhat_fixed[h_idx] + blup_h
  pred
}

fit_pigauto <- function(y, cov_mat, tree, mask, with_cov, sf) {
  df <- data.frame(y = y)
  rownames(df) <- tree$tip.label
  df$y[mask] <- NA_real_
  cov_df <- if (with_cov) {
    cdf <- as.data.frame(cov_mat)
    rownames(cdf) <- tree$tip.label
    cdf
  } else NULL
  res <- tryCatch(
    pigauto::impute(df, tree, covariates = cov_df,
                      epochs = EPOCHS, n_imputations = N_IMP,
                      verbose = FALSE, seed = SEED,
                      safety_floor = sf),
    error = function(e) NULL)
  if (is.null(res)) return(rep(NA_real_, length(y)))
  as.numeric(res$completed$y)
}

# ---------- evaluation ----------
rmse_pearson <- function(pred, truth, mask) {
  ok <- mask & is.finite(pred) & is.finite(truth)
  if (sum(ok) < 5L) return(c(rmse = NA_real_, r = NA_real_))
  d <- pred[ok] - truth[ok]
  r <- if (stats::sd(pred[ok]) < 1e-10) NA_real_
       else suppressWarnings(stats::cor(pred[ok], truth[ok]))
  c(rmse = sqrt(mean(d^2)), r = r)
}

# ---------- enumerate cells ----------
cells <- list()
for (a in ALPHAS) for (b in BETAS) {
  if (!valid_ab(a, b)) next
  for (ft in F_TYPES) for (ct in COV_TYPES) for (rp in seq_len(N_REPS)) {
    cells[[length(cells) + 1L]] <- list(alpha = a, beta = b, f_type = ft,
                                          cov_type = ct, rep_id = rp)
  }
}
n_cells <- length(cells)
cat_line(sprintf("Total cells: %d (skipping invalid alpha+beta>0.9)", n_cells))

# ---------- main sweep ----------
results <- list()
for (i in seq_len(n_cells)) {
  cc <- cells[[i]]
  rep_seed <- SEED + 1000L * cc$rep_id +
                round(cc$alpha * 100) +
                round(cc$beta * 10) +
                match(cc$f_type, F_TYPES) +
                100L * (match(cc$cov_type, COV_TYPES) - 1L)
  cat_line(sprintf("=== %d/%d: a=%.1f b=%.1f f=%s cov=%s rep=%d (seed=%d) ===",
                    i, n_cells,
                    cc$alpha, cc$beta, cc$f_type, cc$cov_type, cc$rep_id, rep_seed))

  sim <- sim_y(tree, cc$alpha, cc$beta, cc$f_type, cc$cov_type, rep_seed)
  y_true <- sim$y
  set.seed(rep_seed + 99L)
  mask <- rep(FALSE, n)
  mask[sample.int(n, round(MISS_FRAC * n))] <- TRUE

  preds <- list()
  walls <- numeric()

  t0 <- proc.time()[["elapsed"]]
  preds$column_mean <- {
    out <- rep(NA_real_, n); out[mask] <- mean(y_true[!mask]); out
  }
  walls["column_mean"] <- proc.time()[["elapsed"]] - t0

  t0 <- proc.time()[["elapsed"]]
  preds$lm <- fit_lm(y_true, sim$cov, mask)
  walls["lm"] <- proc.time()[["elapsed"]] - t0

  t0 <- proc.time()[["elapsed"]]
  preds$phylolm_lambda_blup <- fit_phylolm_lambda_blup(y_true, sim$cov, tree,
                                                        R_full, mask)
  walls["phylolm_lambda_blup"] <- proc.time()[["elapsed"]] - t0

  t0 <- proc.time()[["elapsed"]]
  preds$pigauto_no_cov <- fit_pigauto(y_true, sim$cov, tree, mask,
                                        with_cov = FALSE, sf = TRUE)
  walls["pigauto_no_cov"] <- proc.time()[["elapsed"]] - t0
  cat_line(sprintf("  pigauto_no_cov %.1fs", walls["pigauto_no_cov"]))

  t0 <- proc.time()[["elapsed"]]
  preds$pigauto_cov_sfT <- fit_pigauto(y_true, sim$cov, tree, mask,
                                         with_cov = TRUE, sf = TRUE)
  walls["pigauto_cov_sfT"] <- proc.time()[["elapsed"]] - t0
  cat_line(sprintf("  pigauto_cov_sfT %.1fs", walls["pigauto_cov_sfT"]))

  t0 <- proc.time()[["elapsed"]]
  preds$pigauto_cov_sfF <- fit_pigauto(y_true, sim$cov, tree, mask,
                                         with_cov = TRUE, sf = FALSE)
  walls["pigauto_cov_sfF"] <- proc.time()[["elapsed"]] - t0
  cat_line(sprintf("  pigauto_cov_sfF %.1fs", walls["pigauto_cov_sfF"]))

  for (m in names(preds)) {
    rp <- rmse_pearson(preds[[m]], y_true, mask)
    results[[length(results) + 1L]] <- data.frame(
      cell_id = i, alpha = cc$alpha, beta = cc$beta,
      f_type = cc$f_type, cov_type = cc$cov_type,
      rep_id = cc$rep_id, method = m,
      rmse = rp[["rmse"]], pearson_r = rp[["r"]],
      n_held = sum(mask), wall_s = walls[[m]],
      stringsAsFactors = FALSE)
  }
}

results_df <- do.call(rbind, results)

out_rds <- sprintf("script/bench_gnn_earnings_%s.rds", TAG)
out_md  <- sprintf("script/bench_gnn_earnings_%s.md",  TAG)

saveRDS(list(results = results_df,
              meta = list(tag = TAG, n = n,
                            alphas = ALPHAS, betas = BETAS,
                            f_types = F_TYPES, cov_types = COV_TYPES,
                            n_reps = N_REPS, epochs = EPOCHS,
                            seed = SEED, miss_frac = MISS_FRAC,
                            decision_rule = "pigauto wins if RMSE_pigauto < 0.9 * RMSE_phylolm_lambda AND gap > 1.96 * pooled SE on f in {nonlinear, interactive} & beta=0.4 cells")),
         out_rds, compress = "xz")

# Mean ± SE per (alpha, beta, f, cov, method)
agg <- do.call(rbind, by(results_df,
  list(results_df$alpha, results_df$beta, results_df$f_type,
        results_df$cov_type, results_df$method),
  function(d) {
    data.frame(alpha = d$alpha[1], beta = d$beta[1],
                f_type = d$f_type[1], cov_type = d$cov_type[1],
                method = d$method[1],
                mean_rmse = mean(d$rmse, na.rm = TRUE),
                se_rmse   = stats::sd(d$rmse, na.rm = TRUE) / sqrt(sum(!is.na(d$rmse))),
                mean_r    = mean(d$pearson_r, na.rm = TRUE),
                n_reps_ok = sum(!is.na(d$rmse)))
  }))
agg <- agg[order(agg$cov_type, agg$alpha, agg$beta, agg$f_type, agg$method), ]

md <- c(
  sprintf("# GNN earnings sim v2 (%s, n=%d)", TAG, n),
  "",
  "## Pre-registered decision rule",
  "",
  "pigauto + cov earns its keep if:",
  "  mean(RMSE_pigauto_cov_sfT) < 0.90 * mean(RMSE_phylolm_lambda_blup)",
  "  AND gap > 1.96 * pooled SE",
  "  on f in {nonlinear, interactive} and beta = 0.4",
  "",
  sprintf("## Sweep: alpha=%s, beta=%s, f=%s, cov=%s, n_reps=%d (skip alpha+beta>0.9)",
            paste(ALPHAS, collapse=","),
            paste(BETAS,  collapse=","),
            paste(F_TYPES, collapse=","),
            paste(COV_TYPES, collapse=","), N_REPS),
  "",
  "## Mean RMSE ± SE per cell, all methods",
  "",
  "```",
  capture.output(print(agg, row.names = FALSE)),
  "```")
writeLines(md, out_md)

cat_line("=== DONE ===")
cat_line("  rds: ", out_rds)
cat_line("  md : ", out_md)
