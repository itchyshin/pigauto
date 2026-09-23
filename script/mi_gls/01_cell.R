#!/usr/bin/env Rscript
# script/mi_gls/01_cell.R
#
# One (regime, rep) job for the MI-GLS-attenuation campaign
# (arc/mi-gls-attenuation). Usage:
#   OPENBLAS_NUM_THREADS=1 Rscript script/mi_gls/01_cell.R <regime_id> <rep> <outdir>
#
# Simulates bivariate BM (x, y), rho = 0.7, on an n-tip tree at the regime's
# Pagel lambda (simulation only -- downstream imputation and BOTH downstream
# models always use the lambda = 1 (original) tree; see
# script/mi_gls/regimes.R). Applies the regime's missingness mechanism, then
# compares methods x downstream models. Writes one .rds:
#   <outdir>/regime_<regime_id>_rep_<rep>.rds
#
# Methods:
#   complete   -- reference: fit on the (unobserved-in-practice) truth
#   single     -- impute() default, one completed dataset, naive model SE
#   mi_conf_pc -- multi_impute(conformal, predict_method = "per_column")
#   mi_conf_ex -- multi_impute(conformal, predict_method = "exact")
#   mi_dropout -- multi_impute(mc_dropout)
#   draw_cond  -- draw_conditional_bm() (R/draws_conditional.R prototype)
#   oracle     -- true-model proper MI (only when regime$run_oracle, i.e.
#                 lambda == 1 & missing == "x_only"; see
#                 script/mondrian_confirmation/13_mi_gls_attenuation_diag.R)
#
# Downstream: OLS lm(y ~ x); PGLS gls(y ~ x, corBrownian(tree)). Pooled via
# pigauto::pool_mi() for multi-draw methods; single/complete report the raw
# model SE/df (no pooling -- there is only one dataset). True slope is 0.7
# for BOTH downstream models: with unit-variance bivariate BM (rho = 0.7),
# each (x_i, y_i) pair is marginally Bivariate Normal regardless of the
# phylogenetic correlation across rows, so the population regression
# coefficient is rho for OLS as well as PGLS -- only the standard errors /
# efficiency differ, not the estimand.
#
# epochs default 500 (the full campaign); override via env var
# MI_GLS_EPOCHS for a cheap smoke run (e.g. 150). m (draws) defaults to 20,
# override via MI_GLS_M.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("expected: regime_id rep outdir", call. = FALSE)
regime_id <- as.integer(args[[1L]])
rep_i     <- as.integer(args[[2L]])
outdir    <- args[[3L]]
if (!is.finite(regime_id) || !is.finite(rep_i)) {
  stop("regime_id and rep must be integers", call. = FALSE)
}

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressMessages({ devtools::load_all(quiet = TRUE); library(ape); library(nlme) })
if (requireNamespace("torch", quietly = TRUE)) {
  try(torch::torch_set_num_threads(1L), silent = TRUE)
  try(torch::torch_set_num_interop_threads(1L), silent = TRUE)
}

source(file.path("script", "mi_gls", "regimes.R"))
reg <- regimes[regimes$regime_id == regime_id, ]
if (nrow(reg) != 1L) stop("unknown regime_id: ", regime_id, call. = FALSE)

epochs    <- as.integer(Sys.getenv("MI_GLS_EPOCHS", "500"))
m         <- as.integer(Sys.getenv("MI_GLS_M", "20"))
true_beta <- reg$rho

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out_f <- file.path(outdir, sprintf("regime_%d_rep_%d.rds", regime_id, rep_i))
if (file.exists(out_f)) { cat("SKIP", out_f, "\n"); quit(save = "no") }

seed <- 20260923L + regime_id * 10000L + rep_i
set.seed(seed)

# ---- DGP --------------------------------------------------------------------
n <- reg$n
tree <- ape::rtree(n)
sim_tree <- if (reg$lambda == 1) tree else transform_tree_pagel(tree, reg$lambda)
V_sim <- ape::vcv(sim_tree); V_sim <- V_sim / max(V_sim)
Lc  <- chol(V_sim + 1e-8 * diag(n))
Sig <- matrix(reg$rho, 2, 2); diag(Sig) <- 1
Z <- t(Lc) %*% matrix(stats::rnorm(n * 2), n, 2) %*% chol(Sig)
truth <- data.frame(row.names = tree$tip.label, x = Z[, 1], y = Z[, 2])

# ---- missingness --------------------------------------------------------------
mar_phylo_mask <- function(tree, m_miss = 0.3) {
  n_tip <- ape::Ntip(tree)
  nodes <- (n_tip + 2L):(n_tip + tree$Nnode)
  sizes <- vapply(nodes, function(nd)
    length(ape::extract.clade(tree, nd)$tip.label), integer(1))
  cand <- nodes[sizes >= floor(0.15 * n_tip) & sizes <= ceiling(0.35 * n_tip)]
  picked <- if (length(cand) >= 2L) sample(cand, 2L) else
    nodes[order(abs(sizes - 0.25 * n_tip))[1:2]]
  in_clade <- rep(FALSE, n_tip); names(in_clade) <- tree$tip.label
  for (nd in picked) in_clade[ape::extract.clade(tree, nd)$tip.label] <- TRUE
  pvec <- ifelse(in_clade[tree$tip.label], 7, 1)
  pvec <- pvec * (m_miss * n_tip) / sum(pvec)
  pvec <- pmin(pvec, 0.95)
  mask <- stats::runif(n_tip) < pvec
  if (sum(!mask) < 20L) {
    keep <- sample(which(mask), sum(mask) - (n_tip - 20L))
    mask[keep] <- FALSE
  }
  mask
}
mcar_mask <- function(n_tip, m_miss = 0.3) stats::runif(n_tip) < m_miss

make_mask <- function() {
  if (identical(reg$mechanism, "MCAR")) mcar_mask(n) else mar_phylo_mask(tree)
}

mask_x <- make_mask()
mask_y <- if (identical(reg$missing, "both")) make_mask() else rep(FALSE, n)

df <- truth
df$x[mask_x] <- NA
df$y[mask_y] <- NA

# ---- downstream fit + pooling helpers ----------------------------------------
fit_downstream <- function(dat) {
  dat$species <- rownames(dat)
  ols <- tryCatch(stats::lm(y ~ x, data = dat), error = function(e) NULL)
  pgls <- tryCatch(
    nlme::gls(y ~ x, correlation = ape::corBrownian(phy = tree, form = ~species),
              data = dat, method = "ML"),
    error = function(e) NULL)
  list(ols = ols, pgls = pgls)
}

covers <- function(est, se, dfr, truth_val = true_beta, conf = 0.95) {
  if (!is.finite(est) || !is.finite(se) || !is.finite(dfr) || dfr <= 0) return(NA)
  tcrit <- stats::qt(1 - (1 - conf) / 2, dfr)
  truth_val >= est - tcrit * se & truth_val <= est + tcrit * se
}

na_row <- function(kind) {
  data.frame(downstream = kind, estimate = NA_real_, se = NA_real_,
            df = NA_real_, covered = NA)
}

single_row <- function(model_fit, kind) {
  if (is.null(model_fit)) return(na_row(kind))
  if (kind == "ols") {
    co  <- summary(model_fit)$coefficients
    est <- unname(co["x", "Estimate"]); se <- unname(co["x", "Std. Error"])
    dfr <- model_fit$df.residual
  } else {
    tt  <- summary(model_fit)$tTable
    est <- unname(tt["x", "Value"]); se <- unname(tt["x", "Std.Error"])
    dfr <- model_fit$dims$N - model_fit$dims$p
  }
  data.frame(downstream = kind, estimate = est, se = se, df = dfr,
            covered = covers(est, se, dfr))
}

pool_row <- function(fits_list, kind) {
  fits <- lapply(fits_list, `[[`, kind)
  ok <- !vapply(fits, is.null, logical(1))
  if (sum(ok) < 2L) return(na_row(kind))
  pooled <- tryCatch(pigauto::pool_mi(fits[ok]), error = function(e) e)
  if (inherits(pooled, "error")) return(na_row(kind))
  row <- pooled[pooled$term == "x", , drop = FALSE]
  data.frame(downstream = kind, estimate = row$estimate, se = row$std.error,
            df = row$df, covered = covers(row$estimate, row$std.error, row$df))
}

single_row_or_pool <- function(fits_list, kind) {
  if (length(fits_list) == 1L) single_row(fits_list[[1L]][[kind]], kind)
  else pool_row(fits_list, kind)
}

method_result <- function(method_name, datasets, wall_s) {
  fits_list <- lapply(datasets, fit_downstream)
  rbind(
    cbind(method = method_name, wall_s = wall_s,
         single_row_or_pool(fits_list, "ols")),
    cbind(method = method_name, wall_s = wall_s,
         single_row_or_pool(fits_list, "pgls"))
  )
}

error_result <- function(method_name, wall_s) {
  cbind(method = method_name, wall_s = wall_s,
       rbind(na_row("ols"), na_row("pgls")))
}

results <- list()

# ---- complete-data reference --------------------------------------------------
t0 <- proc.time()[["elapsed"]]
complete_result <- method_result("complete", list(truth), NA_real_)
complete_result$wall_s <- proc.time()[["elapsed"]] - t0
results$complete <- complete_result

# ---- single imputation ----------------------------------------------------------
t0 <- proc.time()[["elapsed"]]
single_res <- tryCatch(
  impute(df, tree, epochs = epochs, verbose = FALSE, seed = seed),
  error = function(e) e)
wall <- proc.time()[["elapsed"]] - t0
results$single <- if (!inherits(single_res, "error")) {
  method_result("single", list(single_res$completed), wall)
} else {
  error_result("single", wall)
}

# ---- multi_impute conformal, predict_method per_column / exact ---------------
run_mi_conformal <- function(predict_method) {
  method_name <- paste0("mi_conf_", predict_method)
  t0 <- proc.time()[["elapsed"]]
  mi <- tryCatch(
    multi_impute(df, tree, m = m, draws_method = "conformal",
                conformal_method = "split", predict_method = predict_method,
                epochs = epochs, verbose = FALSE, seed = seed, gnn = TRUE),
    error = function(e) e)
  wall <- proc.time()[["elapsed"]] - t0
  if (inherits(mi, "error")) return(error_result(method_name, wall))
  method_result(method_name, mi$datasets, wall)
}
results$mi_conf_pc <- run_mi_conformal("per_column")
results$mi_conf_ex <- run_mi_conformal("exact")

# ---- multi_impute mc_dropout -----------------------------------------------------
t0 <- proc.time()[["elapsed"]]
mi_drop <- tryCatch(
  multi_impute(df, tree, m = m, draws_method = "mc_dropout",
              epochs = epochs, verbose = FALSE, seed = seed, gnn = TRUE),
  error = function(e) e)
wall <- proc.time()[["elapsed"]] - t0
results$mi_dropout <- if (inherits(mi_drop, "error")) {
  error_result("mi_dropout", wall)
} else {
  method_result("mi_dropout", mi_drop$datasets, wall)
}

# ---- draw_conditional_bm() (R/draws_conditional.R prototype) -----------------
t0 <- proc.time()[["elapsed"]]
draw_res <- tryCatch({
  pd <- preprocess_traits(df, tree, log_transform = FALSE)
  draw_conditional_bm(list(data = pd, tree = tree), m = m, seed = seed)
}, error = function(e) e)
wall <- proc.time()[["elapsed"]] - t0
if (inherits(draw_res, "error")) {
  results$draw_cond <- error_result("draw_cond", wall)
  draw_cond_error <- conditionMessage(draw_res)
} else {
  results$draw_cond <- method_result("draw_cond", draw_res$datasets, wall)
  draw_cond_error <- NA_character_
}

# ---- oracle proper MI (true model; lambda == 1 & missing == "x_only" only) --
if (isTRUE(reg$run_oracle)) {
  t0 <- proc.time()[["elapsed"]]
  im <- which(mask_x)
  io <- c(setdiff(seq_len(n), im), n + seq_len(n))
  Vt <- ape::vcv(tree); Vt <- Vt / max(Vt)     # oracle always uses the TRUE (lambda = 1) tree
  C <- kronecker(Sig, Vt)
  obs_vec <- c(truth$x[-im], truth$y)
  mu_c <- C[im, io, drop = FALSE] %*% solve(C[io, io, drop = FALSE], obs_vec)
  V_c  <- C[im, im, drop = FALSE] -
    C[im, io, drop = FALSE] %*% solve(C[io, io, drop = FALSE], C[io, im, drop = FALSE])
  Lc2 <- t(chol(V_c + 1e-9 * diag(length(im))))
  oracle_datasets <- lapply(seq_len(m), function(k) {
    d <- truth
    d$x[im] <- as.numeric(mu_c + Lc2 %*% stats::rnorm(length(im)))
    d
  })
  wall <- proc.time()[["elapsed"]] - t0
  results$oracle <- method_result("oracle", oracle_datasets, wall)
}

# ---- assemble + save -----------------------------------------------------------
res_df <- do.call(rbind, results)
rownames(res_df) <- NULL

out <- list(
  regime_id = regime_id, rep = rep_i, seed = seed, regime = reg,
  epochs = epochs, m = m, true_beta = true_beta,
  n_missing_x = sum(mask_x), n_missing_y = sum(mask_y),
  draw_cond_error = draw_cond_error,
  results = res_df
)
saveRDS(out, out_f)

cat(sprintf("OK regime=%d rep=%d n=%d lambda=%.1f mech=%s missing=%s\n",
           regime_id, rep_i, n, reg$lambda, reg$mechanism, reg$missing))
print(res_df)
