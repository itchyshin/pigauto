#!/usr/bin/env Rscript
# script/bench_sim_bace_nonlinear.R
#
# GNN-ability sweep: BACE-simulated DGP + post-hoc nonlinear/interactive
# response injection x pigauto (Fix A-H + row-fix) vs analytical baselines.
#
# WHY THIS BENCH EXISTS
#
# bench_sim_bace_pigauto.R uses purely linear DGPs (y = X*beta + species + phylo
# + obs noise).  On linear DGPs, phylolm-lambda is mathematically optimal and
# pigauto's GNN should close its gate to MATCH (not beat) phylolm.  The
# safety_floor and gate-closure are pigauto's value there: robustness, not
# lift.
#
# To demonstrate the GNN's ABILITY, we need DGPs where the BM/GLS baseline
# mathematically can't capture the signal.  This script injects nonlinear /
# interactive structure on top of BACE-sim's linear gaussian output:
#
#   f_type = "linear"      -> y = (BACE-sim linear y) + 0  (control)
#   f_type = "nonlinear"   -> y = (BACE-sim base y) + beta * scale(sin(2*x1) * exp(0.3*x2))
#   f_type = "interactive" -> y = (BACE-sim base y) + beta * scale(x1*x2 + 0.5*x1^2)
#
# In linear cells, pigauto should match phylolm-lambda (gate-closed parity).
# In nonlinear / interactive cells, pigauto should clearly BEAT phylolm-lambda
# AND beat a nonlinear OLS baseline (lm with poly + interactions, no phylo).
#
# METHODS (5)
#   1. column_mean             -- loss floor
#   2. lm                      -- lm(y ~ x1 + x2 + x3) (linear, no phylo)
#   3. lm_nonlinear            -- lm(y ~ poly(x1, 2) + poly(x2, 2) + x3 + x1:x2)
#                                  (smart practitioner with poly + interaction features, no phylo)
#   4. phylolm_lambda_blup     -- GLS regression + lambda-fitted BM BLUP
#                                  (linear features + phylo, the "linear smart" baseline)
#   5. pigauto_cov_sfT         -- pigauto::impute() with covariates, safety_floor=TRUE
#                                  (nonlinear via GNN + phylo, what we are testing)
#
# DECISIVE PREDICTIONS
#
#   linear cells:        pigauto matches phylolm_lambda within 10 %; both clearly
#                         beat lm and lm_nonlinear.
#   nonlinear cells:     pigauto BEATS phylolm_lambda by >= 10 %; pigauto matches
#                         or beats lm_nonlinear; lm clearly worst.
#   interactive cells:   pigauto matches lm_nonlinear within 10 % (since that
#                         baseline has the bilinear feature) but BEATS
#                         phylolm_lambda by >= 10 %.
#
# Tier (env var PIGAUTO_TIER): smoke / medium / full.
#
# Run:
#   PIGAUTO_TIER=smoke   PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_sim_bace_nonlinear.R
#   PIGAUTO_TIER=medium  PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_sim_bace_nonlinear.R

suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH",
                          unset = "/Users/z3437171/Dropbox/Github Local/pigauto")
  devtools::load_all(pkg_path, quiet = TRUE)
  source(file.path(pkg_path, "script", "sim_bace_dgp.R"))
})
options(warn = -1L)

TIER  <- toupper(Sys.getenv("PIGAUTO_TIER", "smoke"))
HERE  <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(HERE, "script",
                      sprintf("bench_sim_bace_nonlinear_%s.rds", tolower(TIER)))
out_log <- file.path(HERE, "script",
                      sprintf("bench_sim_bace_nonlinear_%s.log", tolower(TIER)))

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) cat(sprintf("[%6.0fs] ", proc.time()[["elapsed"]] - t0),
                                ..., "\n", sep = "")

# ---------------- Sweep grid (per tier) -----------------------------------

if (TIER == "SMOKE") {
  # 24 cells: tests basic feasibility, ~30-40 min wall.
  N_SPECIES <- c(100L)
  MULTI_OBS <- c(1L, 4L)
  PHYLO_SIG <- c(0.4)
  BETA_STR  <- c(0.0, 0.5)
  F_TYPES   <- c("linear", "nonlinear", "interactive")
  N_REPS    <- 2L
} else if (TIER == "MEDIUM") {
  # 216 cells: ~12 hr wall.
  N_SPECIES <- c(100L, 500L)
  MULTI_OBS <- c(1L, 4L)
  PHYLO_SIG <- c(0.2, 0.6)
  BETA_STR  <- c(0.0, 0.3, 0.6)
  F_TYPES   <- c("linear", "nonlinear", "interactive")
  N_REPS    <- 3L
} else if (TIER == "FULL") {
  # 1296 cells: not recommended without HPC.
  N_SPECIES <- c(100L, 500L, 1500L)
  MULTI_OBS <- c(1L, 4L)
  PHYLO_SIG <- c(0.2, 0.5, 0.8)
  BETA_STR  <- c(0.0, 0.3, 0.6)
  F_TYPES   <- c("linear", "nonlinear", "interactive")
  N_REPS    <- 3L
} else {
  stop("PIGAUTO_TIER must be one of: smoke, medium, full")
}

MISS_FRAC <- 0.30
EPOCHS    <- as.integer(Sys.getenv("PIGAUTO_BENCH_EPOCHS", "150"))

log_line("=========================================================")
log_line("Tier: ", TIER, "   (NONLINEAR DGP)")
log_line("Sweep: n_species=", paste(N_SPECIES, collapse=","),
          " | multi_obs=", paste(MULTI_OBS, collapse=","),
          " | phylo_signal=", paste(PHYLO_SIG, collapse=","),
          " | beta=", paste(BETA_STR, collapse=","),
          " | f_types=", paste(F_TYPES, collapse=","),
          " | reps=", N_REPS)
log_line("Output: ", out_rds)

# ---------------- Helper: nonlinear DGP wrapper ----------------------------

#' Generate a sim_bace_dgp() cell, then post-hoc inject a nonlinear /
#' interactive contribution into the response.
sim_nonlinear_dgp <- function(n_species, multi_obs_ratio, phylo_signal,
                                beta_strength, f_type, seed) {
  # Pull base from BACE with beta_resp = 0 so the response has only
  # phylo + species RE + obs noise (no linear cov effect).  We add the
  # nonlinear contribution on top.
  base <- sim_bace_dgp(
    n_species          = n_species,
    multi_obs_ratio    = multi_obs_ratio,
    phylo_signal       = phylo_signal,
    beta_resp_strength = 0,
    response_type      = "gaussian",
    n_predictors       = 3L,
    miss_frac          = MISS_FRAC,
    seed               = seed
  )

  # Build numeric design matrix from predictors (factors -> 0/1 codes).
  pred_cols <- base$predictor_names
  X_full <- vapply(pred_cols, function(p) {
    v <- base$df_complete[[p]]
    if (is.factor(v)) as.numeric(v) - 1 else as.numeric(v)
  }, numeric(nrow(base$df_complete)))
  X_full <- as.matrix(X_full)

  # Compute the nonlinear contribution.
  contrib <- switch(f_type,
    linear      = X_full[, 1],
    nonlinear   = sin(2 * X_full[, 1]) * exp(0.3 * X_full[, 2]),
    interactive = X_full[, 1] * X_full[, 2] + 0.5 * X_full[, 1]^2,
    stop("Unknown f_type: ", f_type)
  )
  contrib <- as.numeric(scale(contrib))   # mean 0, sd 1

  # Inject into response (both df_complete and df_observed; re-apply mask)
  delta <- beta_strength * contrib
  base$df_complete$y <- base$df_complete$y + delta
  base$df_observed$y <- base$df_observed$y + delta
  base$df_observed$y[base$mask] <- NA

  base$meta$f_type        <- f_type
  base$meta$beta_strength <- beta_strength
  base
}

# ---------------- Method implementations ----------------------------------

score <- function(pred, truth, mask) {
  ok <- mask & is.finite(truth) & is.finite(pred)
  if (sum(ok) < 5L) return(c(rmse = NA_real_, r = NA_real_))
  p <- pred[ok]; t <- truth[ok]
  rmse <- sqrt(mean((p - t)^2))
  r <- if (stats::sd(p) > 1e-10) suppressWarnings(stats::cor(p, t)) else NA_real_
  c(rmse = rmse, r = r)
}

method_column_mean <- function(d) {
  cm <- mean(d$df_observed$y, na.rm = TRUE)
  rep(cm, nrow(d$df_observed))
}

method_lm <- function(d) {
  ddf <- d$df_observed
  pred_n <- d$predictor_names
  fit_df <- ddf[!is.na(ddf$y), c("y", pred_n), drop = FALSE]
  fmla <- stats::as.formula(paste("y ~", paste(pred_n, collapse = " + ")))
  fit <- tryCatch(stats::lm(fmla, data = fit_df), error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, nrow(ddf)))
  predict(fit, newdata = ddf)
}

method_lm_nonlinear <- function(d) {
  ddf <- d$df_observed
  pred_n <- d$predictor_names
  # Build nonlinear formula: poly(x1, 2) + poly(x2, 2) + (other linear) + x1:x2
  # x3 may be binary -> just add linearly.
  if (length(pred_n) < 2L) return(method_lm(d))
  # Skip factor predictors from poly() to avoid errors
  fit_df <- ddf[!is.na(ddf$y), c("y", pred_n), drop = FALSE]
  poly_terms <- c()
  linear_terms <- c()
  for (p in pred_n) {
    v <- fit_df[[p]]
    if (is.factor(v) || length(unique(v)) < 4L) {
      linear_terms <- c(linear_terms, p)
    } else {
      poly_terms <- c(poly_terms, sprintf("poly(%s, 2, raw = TRUE)", p))
    }
  }
  rhs <- paste(c(poly_terms, linear_terms), collapse = " + ")
  if (length(poly_terms) >= 2L) {
    # Add bilinear interaction between the first two continuous predictors
    inter_pair <- pred_n[!sapply(fit_df[, pred_n], function(v)
      is.factor(v) || length(unique(v)) < 4L)][1:2]
    if (all(!is.na(inter_pair))) {
      rhs <- paste0(rhs, " + ", inter_pair[1], ":", inter_pair[2])
    }
  }
  fmla <- stats::as.formula(paste("y ~", rhs))
  fit <- tryCatch(stats::lm(fmla, data = fit_df), error = function(e) NULL)
  if (is.null(fit)) return(method_lm(d))
  tryCatch(predict(fit, newdata = ddf), error = function(e) rep(NA_real_, nrow(ddf)))
}

method_phylolm_blup <- function(d) {
  # Aggregate to species level, fit phylolm with model="lambda", BLUP held-out.
  ddf <- d$df_observed
  pred_n <- d$predictor_names
  multi_obs <- d$meta$multi_obs_ratio > 1L

  if (multi_obs) {
    sp_means_y <- tapply(ddf$y, ddf$species,
                          function(v) { v <- v[!is.na(v)]; if (length(v)) mean(v) else NA })
    sp_means_X <- as.data.frame(lapply(pred_n, function(p) {
      v <- ddf[[p]]; if (is.factor(v)) v <- as.numeric(v) - 1
      tapply(v, ddf$species, function(x) mean(x, na.rm = TRUE))
    }))
    names(sp_means_X) <- pred_n
    sp_df <- data.frame(species = names(sp_means_y),
                         y = as.numeric(sp_means_y),
                         sp_means_X, stringsAsFactors = FALSE)
  } else {
    X_num <- as.data.frame(lapply(pred_n, function(p) {
      v <- ddf[[p]]; if (is.factor(v)) as.numeric(v) - 1 else as.numeric(v)
    }))
    names(X_num) <- pred_n
    sp_df <- data.frame(species = ddf$species, y = ddf$y, X_num,
                         stringsAsFactors = FALSE)
  }

  # Subset to species in tree, observed y
  sp_df <- sp_df[sp_df$species %in% d$tree$tip.label, , drop = FALSE]
  sp_df_obs <- sp_df[!is.na(sp_df$y), , drop = FALSE]
  if (nrow(sp_df_obs) < length(pred_n) + 5L) {
    return(rep(NA_real_, nrow(ddf)))
  }
  rownames(sp_df_obs) <- sp_df_obs$species

  # Subset tree to observed species + held-out species
  tree_obs <- ape::keep.tip(d$tree,
                              intersect(sp_df_obs$species, d$tree$tip.label))

  fmla <- stats::as.formula(paste("y ~", paste(pred_n, collapse = " + ")))
  fit <- tryCatch(
    phylolm::phylolm(fmla, data = sp_df_obs, phy = tree_obs, model = "lambda"),
    error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, nrow(ddf)))

  # Predict species-level (point estimate via fitted intercept + beta * X + BLUP residual)
  beta_hat <- coef(fit)
  Xmat <- model.matrix(fmla, sp_df)  # all species
  fixed <- as.numeric(Xmat %*% beta_hat)

  # Residuals on observed species, BLUP for missing
  R <- ape::vcv(d$tree)
  R <- stats::cov2cor(R)
  R <- R[sp_df$species, sp_df$species]
  obs_idx <- which(!is.na(sp_df$y))
  miss_idx <- which(is.na(sp_df$y))
  if (length(miss_idx) == 0L) {
    sp_pred <- sp_df$y
  } else {
    e_obs <- sp_df$y[obs_idx] - fixed[obs_idx]
    R_oo <- R[obs_idx, obs_idx, drop = FALSE]
    R_mo <- R[miss_idx, obs_idx, drop = FALSE]
    blup <- as.numeric(R_mo %*% solve(R_oo + diag(1e-6, nrow(R_oo)), e_obs))
    sp_pred <- numeric(nrow(sp_df))
    sp_pred[obs_idx] <- fixed[obs_idx] + e_obs   # = sp_df$y[obs_idx]
    sp_pred[miss_idx] <- fixed[miss_idx] + blup
  }
  names(sp_pred) <- sp_df$species

  # Map back to obs level for multi-obs
  if (multi_obs) sp_pred[as.character(ddf$species)] else sp_pred[as.character(ddf$species)]
}

method_pigauto_cov_sfT <- function(d) {
  ddf <- d$df_observed
  pred_n <- d$predictor_names

  # Encode binary factor predictor as numeric
  cov_df <- as.data.frame(lapply(pred_n, function(p) {
    v <- ddf[[p]]; if (is.factor(v)) as.numeric(v) - 1 else as.numeric(v)
  }))
  names(cov_df) <- pred_n

  multi_obs <- d$meta$multi_obs_ratio > 1L

  res <- if (multi_obs) {
    pigauto::impute(traits = ddf[, c("species", "y"), drop = FALSE],
                     tree = d$tree, species_col = "species",
                     covariates = cov_df, missing_frac = 0.0,
                     epochs = EPOCHS, verbose = FALSE,
                     seed = sample.int(1e7L, 1L))
  } else {
    # single-obs: rownames must be species
    df_single <- ddf[, c("y", pred_n), drop = FALSE]
    rownames(df_single) <- ddf$species
    pigauto::impute(traits = df_single[, "y", drop = FALSE],
                     tree = d$tree,
                     covariates = cov_df, missing_frac = 0.0,
                     epochs = EPOCHS, verbose = FALSE,
                     seed = sample.int(1e7L, 1L))
  }
  if (multi_obs) res$completed$y else {
    # single-obs result has rownames = tree$tip.label.  Map to ddf order.
    res$completed[as.character(ddf$species), "y"]
  }
}

# ---------------- Main loop -----------------------------------------------

cells <- expand.grid(
  n_species  = N_SPECIES,
  multi_obs  = MULTI_OBS,
  phylo_sig  = PHYLO_SIG,
  beta_str   = BETA_STR,
  f_type     = F_TYPES,
  rep_id     = seq_len(N_REPS),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_cells <- nrow(cells)
log_line("Total cells: ", total_cells)

results <- list()
for (i in seq_len(total_cells)) {
  row <- cells[i, ]
  seed <- 7000L + i
  log_line(sprintf("Cell %d/%d: f=%s, n=%d, multi=%d, ps=%.1f, beta=%.1f, rep=%d (seed=%d)",
                    i, total_cells, row$f_type, row$n_species, row$multi_obs,
                    row$phylo_sig, row$beta_str, row$rep_id, seed))
  d <- tryCatch(
    sim_nonlinear_dgp(n_species = row$n_species,
                        multi_obs_ratio = row$multi_obs,
                        phylo_signal = row$phylo_sig,
                        beta_strength = row$beta_str,
                        f_type = row$f_type, seed = seed),
    error = function(e) { log_line("  DGP ERROR: ", e$message); NULL })
  if (is.null(d)) next

  truth <- d$df_complete$y

  out <- list()
  for (m_name in c("column_mean", "lm", "lm_nonlinear",
                    "phylolm_lambda_blup", "pigauto_cov_sfT")) {
    fn <- get(paste0("method_", m_name))
    pred <- tryCatch(fn(d),
                      error = function(e) {
                        log_line("  ", m_name, " ERROR: ", e$message)
                        rep(NA_real_, nrow(d$df_observed))
                      })
    s <- score(pred, truth, d$mask)
    out[[m_name]] <- list(rmse = s["rmse"], r = s["r"])
  }

  for (m in names(out)) {
    results[[length(results) + 1L]] <- data.frame(
      n_species = row$n_species, multi_obs = row$multi_obs,
      phylo_signal = row$phylo_sig, beta_strength = row$beta_str,
      f_type = row$f_type, rep = row$rep_id,
      method = m,
      rmse = out[[m]]$rmse, r = out[[m]]$r,
      stringsAsFactors = FALSE
    )
  }

  # Incremental save
  if (i %% 4L == 0L || i == total_cells) {
    df <- do.call(rbind, results)
    saveRDS(df, out_rds)
    log_line(sprintf("  [checkpoint] saved %d rows", nrow(df)))
  }
  invisible(gc(verbose = FALSE))
}

df <- do.call(rbind, results)
saveRDS(df, out_rds)
log_line(sprintf("=== DONE: %d cells, %d rows ===", total_cells, nrow(df)))
log_line(sprintf("  rds: %s", out_rds))
