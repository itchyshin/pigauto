#!/usr/bin/env Rscript
# script/bench_gnn_earnings.R
#
# The "GNN earnings" experiment.  Tests whether pigauto's neural-network
# architecture provides value beyond what a smart linear baseline
# (phylolm BLUP under BM) already captures.
#
# Semi-synthetic DGP on a real phylogeny:
#
#   y = sqrt(alpha) * phylo + sqrt(beta) * f(cov) + sqrt(1-alpha-beta) * eps
#
# where:
#   - phylo ~ scale( BM(tree, sigma=1, root=0) )
#   - cov   = 5 standard normal covariates (cov1..cov5)
#   - f(.) ∈ {linear, nonlinear, interactive}:
#       linear      : f = cov1
#       nonlinear   : f = sin(2*cov1) * exp(0.3 * cov2)
#       interactive : f = cov1 * cov2 + 0.5 * cov1^2
#   - eps  ~ N(0, 1)
# Components scaled so total Var(y) ~= 1; alpha + beta + (1-alpha-beta) = 1.
#
# Sweep:
#   alpha ∈ {0.4}  (mid phylogenetic signal)  -- fixed for budget
#   beta  ∈ {0, 0.2, 0.4, 0.6}
#   f_type ∈ {linear, nonlinear, interactive}
#   reps   = 3
#
# 6 methods compared:
#   1. column_mean  (floor)
#   2. lm           (covariates only, no phylogeny)
#   3. phylolm_blup (linear cov + BM BLUP -- the smart linear baseline)
#   4. pigauto_no_cov   (phylogeny only via GNN)
#   5. pigauto_cov_sfT  (full architecture, safety_floor = TRUE)
#   6. pigauto_cov_sfF  (full architecture, safety_floor = FALSE)
#
# Mask: 30% MCAR of the simulated trait per cell.
#
# Tree input controlled by env vars:
#   PIGAUTO_TREE_RDA       path to .rda containing one tree (default: tree300)
#   PIGAUTO_TREE_OBJ       object name in the .rda
#   PIGAUTO_BENCH_N        subsample size (default: full tree)
#   PIGAUTO_BENCH_EPOCHS   pigauto epochs (default 150)
#   PIGAUTO_BENCH_TAG      output tag for filenames
#
# Output:
#   script/bench_gnn_earnings_<tag>.rds
#   script/bench_gnn_earnings_<tag>.md

suppressPackageStartupMessages({
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH", unset = "")
  if (nzchar(pkg_path)) devtools::load_all(pkg_path, quiet = TRUE)
  else library(pigauto)
  library(ape)
  library(phylolm)
})
options(warn = -1L)

# ---------- config ----------
SEED         <- 2026L
MISS_FRAC    <- 0.30
EPOCHS       <- as.integer(Sys.getenv("PIGAUTO_BENCH_EPOCHS", "150"))
N_IMP        <- 10L
N_TARGET     <- as.integer(Sys.getenv("PIGAUTO_BENCH_N", "0"))   # 0 = full
TAG          <- Sys.getenv("PIGAUTO_BENCH_TAG", "tree300")

ALPHA        <- 0.4
BETAS        <- c(0.0, 0.2, 0.4, 0.6)
F_TYPES      <- c("linear", "nonlinear", "interactive")
N_REPS       <- 3L

cat_line <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "),
                                ..., "\n", sep = "")

# ---------- load phylogeny ----------
tree_rda <- Sys.getenv("PIGAUTO_TREE_RDA", "")
tree_obj <- Sys.getenv("PIGAUTO_TREE_OBJ", "tree300")

if (nzchar(tree_rda)) {
  load(tree_rda)
  tree_full <- get(tree_obj)
} else {
  data(list = tree_obj)
  tree_full <- get(tree_obj)
}
if (!inherits(tree_full, "phylo"))
  stop("Loaded object is not a phylo")

# Sanitise: ensure binary, rooted, positive branch lengths
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
cat_line(sprintf("Loaded tree: tag=%s, n=%d tips, binary=%s, rooted=%s",
                  TAG, n, ape::is.binary(tree), ape::is.rooted(tree)))

# ---------- DGP helpers ----------
sim_phylo <- function(tree, seed) {
  set.seed(seed)
  v <- ape::rTraitCont(tree, model = "BM", sigma = 1, root.value = 0)
  v <- as.numeric(v[tree$tip.label])
  (v - mean(v)) / stats::sd(v)
}
sim_cov <- function(n, seed) {
  set.seed(seed + 1L)
  matrix(stats::rnorm(n * 5L), nrow = n,
          dimnames = list(NULL, paste0("cov", 1:5)))
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
sim_y <- function(tree, alpha, beta, f_type, seed) {
  n <- length(tree$tip.label)
  phylo <- sim_phylo(tree, seed)
  cov_m <- sim_cov(n, seed)
  fx    <- f_apply(cov_m, f_type)
  set.seed(seed + 2L)
  eps   <- stats::rnorm(n)
  resid_w <- max(0, 1 - alpha - beta)
  y <- sqrt(alpha) * phylo + sqrt(beta) * fx + sqrt(resid_w) * eps
  list(y = y, phylo = phylo, cov = cov_m, f = fx, eps = eps)
}

# ---------- methods ----------

#' Fit each method on observed cells, predict held-out cells.
#' Returns list with named numeric vectors of predictions, length n,
#' filled only at held-out positions.

fit_column_mean <- function(y_obs, n) {
  rep(mean(y_obs, na.rm = TRUE), n)
}

fit_lm <- function(y, cov_mat, mask) {
  obs <- !mask
  df_o <- data.frame(y = y[obs], cov_mat[obs, , drop = FALSE])
  fit <- stats::lm(y ~ ., data = df_o)
  df_h <- data.frame(cov_mat)
  pred <- as.numeric(stats::predict(fit, newdata = df_h))
  pred
}

fit_phylolm_blup <- function(y, cov_mat, tree, mask) {
  obs <- !mask
  df_o <- data.frame(y = y[obs], cov_mat[obs, , drop = FALSE])
  rownames(df_o) <- tree$tip.label[obs]
  tree_o <- ape::keep.tip(tree, tree$tip.label[obs])
  fit <- tryCatch(
    phylolm::phylolm(y ~ cov1 + cov2 + cov3 + cov4 + cov5,
                       data = df_o, phy = tree_o, model = "BM"),
    error = function(e) NULL)
  if (is.null(fit)) return(rep(mean(y[obs]), length(y)))

  # Fixed-effect prediction for ALL tips
  X_all <- cbind(1, cov_mat)
  beta_hat <- as.numeric(coef(fit))
  yhat_fixed <- as.numeric(X_all %*% beta_hat)

  # BLUP correction: y_h = yhat_fixed_h + R_ho R_oo^-1 (y_o - yhat_fixed_o)
  R <- stats::cov2cor(ape::vcv(tree))
  o_idx <- which(obs); h_idx <- which(mask)
  resid_o <- y[o_idx] - yhat_fixed[o_idx]
  R_ho <- R[h_idx, o_idx, drop = FALSE]
  R_oo <- R[o_idx, o_idx, drop = FALSE]
  R_oo_inv <- tryCatch(
    solve(R_oo + 1e-6 * diag(length(o_idx))),
    error = function(e) MASS::ginv(R_oo))
  blup_h <- as.numeric(R_ho %*% (R_oo_inv %*% resid_o))
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

# ---------- main sweep ----------
results <- list()
row_id  <- 1L
total_cells <- length(BETAS) * length(F_TYPES) * N_REPS
cell_id <- 0L

for (beta in BETAS) for (f_type in F_TYPES) for (rep_id in seq_len(N_REPS)) {
  cell_id <- cell_id + 1L
  rep_seed <- SEED + 1000L * rep_id + round(beta * 100) +
                match(f_type, F_TYPES)
  cat_line(sprintf("=== cell %d/%d: beta=%.1f f=%s rep=%d (seed=%d) ===",
                    cell_id, total_cells, beta, f_type, rep_id, rep_seed))

  sim <- sim_y(tree, ALPHA, beta, f_type, rep_seed)
  y_true <- sim$y

  # 30% MCAR
  set.seed(rep_seed + 99L)
  mask <- rep(FALSE, n)
  mask[sample.int(n, round(MISS_FRAC * n))] <- TRUE

  preds <- list()
  walls <- numeric()

  t0 <- proc.time()[["elapsed"]]
  preds$column_mean <- {
    out <- rep(NA_real_, n)
    out[mask] <- mean(y_true[!mask])
    out
  }
  walls["column_mean"] <- proc.time()[["elapsed"]] - t0

  t0 <- proc.time()[["elapsed"]]
  preds$lm <- fit_lm(y_true, sim$cov, mask)
  walls["lm"] <- proc.time()[["elapsed"]] - t0

  t0 <- proc.time()[["elapsed"]]
  preds$phylolm_blup <- fit_phylolm_blup(y_true, sim$cov, tree, mask)
  walls["phylolm_blup"] <- proc.time()[["elapsed"]] - t0

  t0 <- proc.time()[["elapsed"]]
  preds$pigauto_no_cov <- fit_pigauto(y_true, sim$cov, tree, mask,
                                        with_cov = FALSE, sf = TRUE)
  walls["pigauto_no_cov"] <- proc.time()[["elapsed"]] - t0
  cat_line(sprintf("  pigauto_no_cov done in %.1fs", walls["pigauto_no_cov"]))

  t0 <- proc.time()[["elapsed"]]
  preds$pigauto_cov_sfT <- fit_pigauto(y_true, sim$cov, tree, mask,
                                         with_cov = TRUE, sf = TRUE)
  walls["pigauto_cov_sfT"] <- proc.time()[["elapsed"]] - t0
  cat_line(sprintf("  pigauto_cov_sfT done in %.1fs", walls["pigauto_cov_sfT"]))

  t0 <- proc.time()[["elapsed"]]
  preds$pigauto_cov_sfF <- fit_pigauto(y_true, sim$cov, tree, mask,
                                         with_cov = TRUE, sf = FALSE)
  walls["pigauto_cov_sfF"] <- proc.time()[["elapsed"]] - t0
  cat_line(sprintf("  pigauto_cov_sfF done in %.1fs", walls["pigauto_cov_sfF"]))

  for (m in names(preds)) {
    rp <- rmse_pearson(preds[[m]], y_true, mask)
    results[[length(results) + 1L]] <- data.frame(
      cell_id = cell_id, beta = beta, f_type = f_type, rep_id = rep_id,
      method  = m,
      rmse    = rp[["rmse"]], pearson_r = rp[["r"]],
      n_held  = sum(mask), wall_s = walls[[m]],
      stringsAsFactors = FALSE)
    row_id <- row_id + 1L
  }
}

results_df <- do.call(rbind, results)

out_rds <- sprintf("script/bench_gnn_earnings_%s.rds", TAG)
out_md  <- sprintf("script/bench_gnn_earnings_%s.md",  TAG)

saveRDS(list(results = results_df,
              meta = list(tag = TAG, n = n, alpha = ALPHA,
                            betas = BETAS, f_types = F_TYPES,
                            n_reps = N_REPS, epochs = EPOCHS,
                            seed = SEED, miss_frac = MISS_FRAC)),
         out_rds, compress = "xz")

# Markdown summary
agg <- aggregate(cbind(rmse, pearson_r) ~ method + beta + f_type,
                  data = results_df, FUN = function(x) mean(x, na.rm = TRUE))
agg <- agg[order(agg$f_type, agg$beta, agg$method), ]

md <- c(
  sprintf("# GNN earnings sim (%s, n=%d)", TAG, n),
  "",
  sprintf("DGP: y = sqrt(%.1f)*phylo + sqrt(beta)*f(cov) + sqrt(residual)*eps",
            ALPHA),
  sprintf("Sweep: beta ∈ {%s}, f ∈ {%s}, %d reps, 30%% MCAR.",
            paste(BETAS, collapse=", "),
            paste(F_TYPES, collapse=", "), N_REPS),
  "",
  "## Mean across reps (RMSE, Pearson r)",
  "",
  "```",
  capture.output(print(agg, row.names = FALSE)),
  "```")
writeLines(md, out_md)

cat_line("=== DONE ===")
cat_line("  rds: ", out_rds)
cat_line("  md : ", out_md)
