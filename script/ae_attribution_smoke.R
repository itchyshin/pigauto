#!/usr/bin/env Rscript
# script/ae_attribution_smoke.R
#
# Item #3 (2026-04-29): "Is the GNN delta actually doing work, or is
# pigauto's win all from the analytical pipeline?"
#
# Method: gate-zeroing.  Fix 4 (commit fb4461e) unblocked the
# manual-override path, and the strict val-floor fix (commit d651c6c)
# guarantees that gate-zeroed pigauto recovers the pure-baseline
# performance.  We compare:
#
#   1. pigauto_full     -- calibrated gate, full pipeline (BM + GNN)
#   2. pigauto_gnn_off  -- gate forced to 0 (BM baseline only, no GNN)
#   3. phylolm_blup     -- standalone phylolm BLUP (no pigauto pipeline)
#   4. column_mean      -- trivial floor
#
# If pigauto_full > pigauto_gnn_off significantly:  GNN earning its keep.
# If pigauto_full == pigauto_gnn_off:                GNN contributes 0.
# If pigauto_gnn_off == phylolm_blup:                pipeline = phylolm.
# If pigauto_gnn_off > phylolm_blup:                 pipeline > phylolm.
#
# Single-obs DGP with n=200 species, 30% MCAR, lambda in {0.1, 0.3, 0.5}
# crossed with f-type in {linear, nonlinear, interactive} and beta=0.5.
# 3 reps per cell.  9 cells * 3 reps = 27 fits each method.

suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
  devtools::load_all(quiet = TRUE)
})

set.seed(202604291L)

N_SPECIES <- 200L
N_REPS    <- 3L
MISS_FRAC <- 0.30
SIGMA_RES <- 0.5
EPOCHS    <- 80L

LAMBDAS <- c(0.1, 0.3, 0.5)
F_TYPES <- c("linear", "nonlinear", "interactive")

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) cat(sprintf("[%5.0fs] ", proc.time()[["elapsed"]] - t0),
                                ..., "\n", sep = "")

# DGP
sim_one <- function(n_species, lambda, f_type, beta, seed) {
  set.seed(seed)
  tree <- ape::rtree(n_species)
  sp <- tree$tip.label
  phylo_named <- ape::rTraitCont(tree, model = "BM",
                                   sigma = sqrt(lambda), root.value = 0)
  phylo <- as.numeric(phylo_named[sp])  # reorder by sp first, THEN strip names
  names(phylo) <- sp
  Z <- matrix(rnorm(n_species * 5L), n_species, 5L)
  colnames(Z) <- paste0("z", 1:5)
  contrib <- switch(f_type,
    linear      = as.numeric(Z %*% c(0.5, 0.3, 0.2, 0.1, -0.1)),
    nonlinear   = sin(Z[, 1]) * exp(0.3 * Z[, 2]) + 0.5 * tanh(Z[, 3]) +
                    0.3 * sin(2 * Z[, 4]) - 0.2 * cos(Z[, 5]),
    interactive = Z[, 1] * Z[, 2] + 0.5 * Z[, 3]^2 + 0.3 * Z[, 1] * Z[, 4] -
                    0.2 * Z[, 2] * Z[, 5]
  )
  contrib <- as.numeric(scale(contrib))
  y <- phylo + beta * contrib + rnorm(n_species, 0, SIGMA_RES)

  truth <- y
  k <- floor(MISS_FRAC * n_species)
  mask_idx <- sample(seq_len(n_species), k)
  y_obs <- y; y_obs[mask_idx] <- NA

  list(tree = tree, sp = sp, truth = truth, y_obs = y_obs, mask = mask_idx,
        Z = as.data.frame(Z))
}

score <- function(pred, truth, mask) {
  e <- pred[mask] - truth[mask]
  c(rmse = sqrt(mean(e^2, na.rm = TRUE)),
    r    = stats::cor(pred[mask], truth[mask], use = "complete.obs"))
}

# Methods
m_column_mean <- function(d) {
  rep(mean(d$y_obs, na.rm = TRUE), length(d$y_obs))
}

m_phylolm_blup <- function(d) {
  obs_idx <- which(!is.na(d$y_obs))
  fit <- phylolm::phylolm(
    y ~ z1 + z2 + z3 + z4 + z5,
    data = data.frame(y = d$y_obs[obs_idx], d$Z[obs_idx, ]),
    phy = ape::keep.tip(d$tree, d$sp[obs_idx]),
    model = "lambda"
  )
  beta_hat <- coef(fit)
  Xm <- cbind(1, as.matrix(d$Z))
  Xm_pred <- Xm %*% beta_hat[match(c("(Intercept)",
                                       paste0("z", 1:5)), names(beta_hat))]

  # Phylogenetic prediction = X_m * beta + R_mo R_oo^-1 (y_o - X_o * beta)
  R <- ape::vcv(d$tree, model = "Brownian")
  R <- R / max(R)   # normalize
  ord <- match(d$sp, rownames(R))
  R <- R[ord, ord]
  m_idx <- d$mask
  o_idx <- setdiff(seq_along(d$sp), m_idx)
  R_mo <- R[m_idx, o_idx, drop = FALSE]
  R_oo <- R[o_idx, o_idx, drop = FALSE] + diag(1e-6, length(o_idx))
  resid_o <- d$y_obs[o_idx] - Xm_pred[o_idx]
  pred <- as.numeric(Xm_pred)
  pred[m_idx] <- as.numeric(Xm_pred[m_idx] +
                              R_mo %*% solve(R_oo, resid_o))
  pred
}

fit_pigauto_pair <- function(d, seed) {
  ddf <- data.frame(row.names = d$sp,
                    y = d$y_obs)
  cov_df <- d$Z
  res <- pigauto::impute(
    traits = ddf, tree = d$tree,
    covariates = cov_df, missing_frac = 0.20,
    epochs = EPOCHS, verbose = FALSE, seed = seed,
    safety_floor = TRUE
  )
  fit <- res$fit
  pred_full <- res$completed$y

  # Override the GNN gate to 0 but KEEP the calibrated BM / MEAN split.
  # This isolates the GNN delta's contribution while preserving whatever
  # the strict val-floor picked for the baseline-vs-mean blend.  If the
  # full pipeline chose r_BM = 1 (pure BM), gnn_off matches it exactly;
  # if it chose r_MEAN = 1 (pure mean), gnn_off matches that.
  p_lat <- as.integer(fit$model_config$input_dim)
  fit_off <- fit
  rb <- fit$r_cal_bm
  rg <- fit$r_cal_gnn
  rm_ <- fit$r_cal_mean
  if (is.null(rb) || length(rb) != p_lat) rb <- rep(1, p_lat)
  if (is.null(rm_) || length(rm_) != p_lat) rm_ <- rep(0, p_lat)
  if (is.null(rg) || length(rg) != p_lat) rg <- rep(0, p_lat)

  # Renormalise: r_bm + r_mean takes the GNN's weight.  Allocate the
  # GNN weight proportionally between BM and MEAN so the new corner
  # stays on the calibrated (BM, MEAN) ridge.
  s <- rb + rm_
  rb_new <- ifelse(s > 0, rb / s, 1)
  rm_new <- ifelse(s > 0, rm_ / s, 0)
  fit_off$calibrated_gates <- rep(0, p_lat)
  fit_off$r_cal_gnn        <- rep(0, p_lat)
  fit_off$r_cal_bm         <- rb_new
  fit_off$r_cal_mean       <- rm_new
  pred_obj <- predict(fit_off, return_se = FALSE)
  pred_internal <- pred_obj$imputed$y

  # Predict returns internal (tree-tip) order.  Convert back via input_row_order.
  iro <- res$data$input_row_order
  pred_off <- numeric(length(d$y_obs))
  inv_map <- match(seq_along(iro), iro)
  for (i in seq_along(d$y_obs)) {
    if (!is.na(d$y_obs[i])) {
      pred_off[i] <- d$y_obs[i]
    } else {
      pred_off[i] <- pred_internal[inv_map[i]]
    }
  }
  list(pred_full = pred_full, pred_gnn_off = pred_off,
        gates = list(r_cal_gnn = fit$r_cal_gnn,
                     r_cal_bm  = fit$r_cal_bm,
                     r_cal_mean = fit$r_cal_mean))
}

# Main loop
cells <- expand.grid(lambda = LAMBDAS, f_type = F_TYPES,
                       rep_id = seq_len(N_REPS),
                       KEEP.OUT.ATTRS = FALSE,
                       stringsAsFactors = FALSE)

results <- list()
for (i in seq_len(nrow(cells))) {
  row <- cells[i, ]
  seed <- 31000L + i
  log_line(sprintf("Cell %d/%d: f=%s, lam=%.2f, rep=%d",
                    i, nrow(cells), row$f_type, row$lambda, row$rep_id))
  d <- tryCatch(sim_one(N_SPECIES, row$lambda, row$f_type, 0.5, seed),
                 error = function(e) { log_line("  DGP ERR: ", e$message); NULL })
  if (is.null(d)) next

  push_result <- function(method, pred) {
    s <- tryCatch(score(pred, d$truth, d$mask),
                   error = function(e) c(rmse = NA, r = NA))
    results[[length(results) + 1L]] <- data.frame(
      lambda = row$lambda, f_type = row$f_type, rep = row$rep_id,
      method = method, rmse = s["rmse"], r = s["r"]
    )
  }

  # Baselines
  pred_cm <- tryCatch(m_column_mean(d), error = function(e) rep(NA, N_SPECIES))
  push_result("column_mean", pred_cm)
  pred_pl <- tryCatch(m_phylolm_blup(d), error = function(e) {
    log_line("  phylolm ERR: ", e$message); rep(NA, N_SPECIES) })
  push_result("phylolm_blup", pred_pl)

  # Pigauto pair
  ab <- tryCatch(fit_pigauto_pair(d, seed),
                  error = function(e) { log_line("  pigauto ERR: ", e$message); NULL })
  if (!is.null(ab)) {
    push_result("pigauto_full",    ab$pred_full)
    push_result("pigauto_gnn_off", ab$pred_gnn_off)
    rg_str <- if (is.numeric(ab$gates$r_cal_gnn)) {
      paste(round(ab$gates$r_cal_gnn, 3), collapse = ",")
    } else "NULL"
    log_line(sprintf("  full=%.3f gnn_off=%.3f phylolm=%.3f r_cal_gnn=%s",
                      score(ab$pred_full, d$truth, d$mask)["rmse"],
                      score(ab$pred_gnn_off, d$truth, d$mask)["rmse"],
                      score(pred_pl, d$truth, d$mask)["rmse"],
                      rg_str))
  }
}

df <- do.call(rbind, results)
saveRDS(df, "script/ae_attribution_smoke.rds")

# Summary
cat("\n========= AE-attribution summary =========\n")
agg <- aggregate(rmse ~ f_type + lambda + method, df, mean)
agg <- agg[order(agg$f_type, agg$lambda, agg$method), ]
print(agg)

# Per-cell pigauto_full vs pigauto_gnn_off RMSE delta
cat("\n========= GNN contribution (pigauto_full - pigauto_gnn_off) =========\n")
w <- reshape(df[df$method %in% c("pigauto_full", "pigauto_gnn_off"),
                  c("lambda","f_type","rep","method","rmse")],
              idvar = c("lambda","f_type","rep"), timevar = "method",
              direction = "wide")
w$gnn_delta <- w$rmse.pigauto_full - w$rmse.pigauto_gnn_off
w$gnn_pct   <- 100 * w$gnn_delta / w$rmse.pigauto_gnn_off
print(w)

cat("\n========= Mean GNN delta by cell =========\n")
ag <- aggregate(cbind(gnn_delta, gnn_pct) ~ f_type + lambda, w, mean)
print(ag)
cat("\nNegative gnn_delta = GNN reduces RMSE = GNN earning its keep.\n")
cat("Positive gnn_delta = GNN HURTS = strict val-floor should have caught it.\n")
log_line("DONE")
