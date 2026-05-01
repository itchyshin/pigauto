#!/usr/bin/env Rscript
# script/bench_sim_bace_pigauto.R
#
# Comprehensive sweep: BACE-simulated DGP (sim_bace_dgp.R) x pigauto (Fix A-H)
# vs analytical baselines (column_mean, lm, phylolm-lambda BLUP).
#
# Methods compared per cell:
#   1. column_mean         -- global mean (loss floor)
#   2. lm                  -- lm(y ~ X), no phylogeny
#   3. phylolm_lambda_blup -- GLS regression + lambda-fitted BM BLUP (gaussian only)
#   4. pigauto_cov_sfT     -- pigauto::impute() with covariates, safety_floor=TRUE
#
# Tier (env var PIGAUTO_TIER): smoke / medium / full.
#
# Output (saved incrementally):
#   script/bench_sim_bace_pigauto_<TIER>.rds  -- tidy results data.frame
#   script/bench_sim_bace_pigauto_<TIER>.log  -- progress log
#
# Run:
#   PIGAUTO_TIER=smoke   PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_sim_bace_pigauto.R
#   PIGAUTO_TIER=medium  PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_sim_bace_pigauto.R
#   PIGAUTO_TIER=full    PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_sim_bace_pigauto.R

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
                      sprintf("bench_sim_bace_pigauto_%s.rds", tolower(TIER)))
out_log <- file.path(HERE, "script",
                      sprintf("bench_sim_bace_pigauto_%s.log", tolower(TIER)))

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) cat(sprintf("[%6.0fs] ", proc.time()[["elapsed"]] - t0),
                                ..., "\n", sep = "", file = stdout(), append = TRUE)

# ---------------- Sweep grid (per tier) -----------------------------------

if (TIER == "SMOKE") {
  N_SPECIES <- c(100L)
  MULTI_OBS <- c(1L, 4L)
  PHYLO_SIG <- c(0.2, 0.6)
  BETA_STR  <- c(0.0, 0.5)
  RESP_TYPE <- c("gaussian")
  N_PRED    <- c(3L)
  N_REPS    <- 2L
} else if (TIER == "MEDIUM") {
  N_SPECIES <- c(100L, 500L)
  MULTI_OBS <- c(1L, 4L)
  PHYLO_SIG <- c(0.2, 0.6)
  BETA_STR  <- c(0.0, 0.3, 0.6)
  RESP_TYPE <- c("gaussian", "binary")
  N_PRED    <- c(3L)
  N_REPS    <- 3L
} else if (TIER == "FULL") {
  N_SPECIES <- c(100L, 500L, 1500L)
  MULTI_OBS <- c(1L, 4L)
  PHYLO_SIG <- c(0.2, 0.5, 0.8)
  BETA_STR  <- c(0.0, 0.3, 0.6)
  RESP_TYPE <- c("gaussian", "binary", "threshold4")
  N_PRED    <- c(3L, 5L)
  N_REPS    <- 3L
} else {
  stop("PIGAUTO_TIER must be one of: smoke, medium, full")
}

MISS_FRAC <- 0.30
EPOCHS    <- as.integer(Sys.getenv("PIGAUTO_BENCH_EPOCHS", "150"))

log_line("=========================================================")
log_line("Tier: ", TIER)
log_line("Sweep: n_species=", paste(N_SPECIES, collapse=","),
          " | multi_obs=", paste(MULTI_OBS, collapse=","),
          " | phylo_signal=", paste(PHYLO_SIG, collapse=","),
          " | beta=", paste(BETA_STR, collapse=","),
          " | response=", paste(RESP_TYPE, collapse=","),
          " | n_pred=", paste(N_PRED, collapse=","),
          " | reps=", N_REPS)
log_line("Output: ", out_rds)

# ---------------- Method implementations -----------------------------------

# Compute RMSE / accuracy on held-out cells
score <- function(pred, truth, mask, response_type) {
  ok <- mask & is.finite(truth)
  if (response_type == "gaussian") {
    if (sum(ok) < 5L) return(c(rmse = NA_real_, acc = NA_real_, r = NA_real_))
    p <- as.numeric(pred[ok]); t <- as.numeric(truth[ok])
    rmse <- sqrt(mean((p - t)^2, na.rm = TRUE))
    r <- if (stats::sd(p, na.rm = TRUE) > 1e-10)
           suppressWarnings(stats::cor(p, t, use = "complete.obs"))
         else NA_real_
    c(rmse = rmse, acc = NA_real_, r = r)
  } else {
    # binary / threshold (categorical)
    p <- pred[ok]; t <- as.integer(truth[ok])
    # If pred is probability/continuous, threshold at 0.5 for binary,
    # round to nearest integer for threshold
    if (response_type == "binary") {
      p_class <- as.integer(as.numeric(p) >= 0.5)
    } else {
      p_class <- as.integer(round(as.numeric(p)))
    }
    acc <- mean(p_class == t, na.rm = TRUE)
    c(rmse = NA_real_, acc = acc, r = NA_real_)
  }
}

# Method 1: column mean
fit_column_mean <- function(sim) {
  obs <- !sim$mask
  m <- if (sim$response_type == "gaussian") {
    mean(sim$df_observed[[sim$response_name]][obs], na.rm = TRUE)
  } else {
    # mode for discrete
    tab <- table(sim$df_observed[[sim$response_name]][obs])
    as.numeric(names(tab)[which.max(tab)])
  }
  pred <- rep(m, nrow(sim$df_observed))
  pred
}

# Method 2: linear model (no phylogeny) — gaussian only really
fit_lm <- function(sim) {
  obs <- !sim$mask
  df_o <- sim$df_observed[obs, c(sim$response_name, sim$predictor_names),
                          drop = FALSE]
  if (sim$response_type == "gaussian") {
    fit <- tryCatch(stats::lm(stats::as.formula(paste0(sim$response_name,
                                  " ~ ", paste(sim$predictor_names, collapse=" + "))),
                                data = df_o),
                     error = function(e) NULL)
    if (is.null(fit)) return(rep(mean(df_o[[sim$response_name]], na.rm=TRUE),
                                  nrow(sim$df_observed)))
    df_h <- sim$df_observed[, sim$predictor_names, drop = FALSE]
    as.numeric(stats::predict(fit, newdata = df_h))
  } else {
    # logistic / multinomial via glm — keep simple: classify majority
    tab <- table(df_o[[sim$response_name]])
    rep(as.numeric(names(tab)[which.max(tab)]), nrow(sim$df_observed))
  }
}

# Method 3: phylolm-lambda BLUP — gaussian only.  When multi_obs > 1 the
# tree has fewer tips than rows; aggregate to species means then
# kriging-predict at species level, then broadcast back to obs.
fit_phylolm_lambda_blup <- function(sim) {
  if (sim$response_type != "gaussian") return(rep(NA_real_, nrow(sim$df_observed)))

  df_full <- sim$df_observed   # has species + y + x1...xK
  resp <- sim$response_name
  preds <- sim$predictor_names

  # Aggregate to species level (mean of obs per species).  Multi-obs
  # collapses to species means for the kriging step.  The hold-out
  # mask is at obs level; a species is "held-out" if ANY obs of it is
  # held-out (we'll predict at species level then broadcast).
  by_sp <- split(df_full, df_full$species)
  spp <- names(by_sp)

  y_sp <- vapply(by_sp, function(d) {
    v <- d[[resp]]
    v <- v[!is.na(v)]
    if (length(v) == 0L) NA_real_ else mean(v)
  }, numeric(1))
  X_sp <- t(vapply(by_sp, function(d) {
    sapply(preds, function(p) mean(suppressWarnings(as.numeric(d[[p]])), na.rm = TRUE))
  }, numeric(length(preds))))
  rownames(X_sp) <- spp

  obs_sp <- !is.na(y_sp)
  spp_obs <- spp[obs_sp]
  if (length(spp_obs) < length(preds) + 5L) {
    return(rep(mean(y_sp, na.rm = TRUE), nrow(sim$df_observed)))
  }

  tree <- ape::keep.tip(sim$tree, spp)
  tree_o <- ape::keep.tip(tree, spp_obs)
  df_o <- data.frame(y = y_sp[spp_obs], X_sp[spp_obs, , drop = FALSE])
  rownames(df_o) <- spp_obs

  fit <- tryCatch(
    phylolm::phylolm(stats::as.formula(paste0("y ~ ", paste(preds, collapse=" + "))),
                       data = df_o, phy = tree_o, model = "lambda",
                       lower.bound = 0.0, upper.bound = 1.0),
    error = function(e) NULL)
  if (is.null(fit)) return(rep(mean(y_sp, na.rm = TRUE), nrow(sim$df_observed)))

  lam <- as.numeric(fit$optpar)
  if (!is.finite(lam) || lam < 0 || lam > 1) lam <- 1.0

  X_all <- cbind(1, X_sp)
  colnames(X_all)[1] <- "(Intercept)"
  beta_hat <- as.numeric(coef(fit))
  yhat_fixed <- as.numeric(X_all %*% beta_hat)

  R_full <- stats::cov2cor(ape::vcv(tree))
  o_idx <- match(spp_obs, spp); h_idx <- which(!obs_sp)
  if (length(h_idx) == 0L) {
    pred_sp <- yhat_fixed
  } else {
    V_full <- lam * R_full + (1 - lam) * diag(nrow(R_full))
    V_oo <- V_full[o_idx, o_idx]; V_ho <- V_full[h_idx, o_idx]
    V_oo_inv <- tryCatch(
      solve(V_oo + 1e-6 * diag(length(o_idx))),
      error = function(e) MASS::ginv(V_oo))
    resid_o <- y_sp[o_idx] - yhat_fixed[o_idx]
    blup_h <- as.numeric(V_ho %*% (V_oo_inv %*% resid_o))
    pred_sp <- yhat_fixed
    pred_sp[h_idx] <- yhat_fixed[h_idx] + blup_h
  }
  names(pred_sp) <- spp

  # Broadcast species-level predictions back to obs level
  as.numeric(pred_sp[as.character(df_full$species)])
}

# Method 4: pigauto with covariates (Fix A-H)
fit_pigauto_cov <- function(sim, epochs = EPOCHS) {
  df_full <- sim$df_observed   # has species + y + x1...xK
  resp <- sim$response_name
  preds <- sim$predictor_names
  multi_obs <- sim$meta$multi_obs_ratio > 1L

  # Build pigauto traits + covariates frames
  if (multi_obs) {
    traits_df <- df_full[, c("species", resp), drop = FALSE]
    cov_df    <- df_full[, preds, drop = FALSE]
    species_arg <- "species"
  } else {
    traits_df <- df_full[, resp, drop = FALSE]
    rownames(traits_df) <- df_full$species
    cov_df <- df_full[, preds, drop = FALSE]
    rownames(cov_df) <- df_full$species
    species_arg <- NULL
  }

  # Trait type: pigauto auto-detects continuous numeric, but binary needs
  # factor and threshold/multinomial need ordered/factor
  if (sim$response_type == "binary") {
    traits_df[[resp]] <- factor(traits_df[[resp]], levels = c(0, 1))
  } else if (grepl("^threshold", sim$response_type)) {
    traits_df[[resp]] <- ordered(traits_df[[resp]])
  }

  res <- tryCatch(
    pigauto::impute(
      traits      = traits_df,
      tree        = sim$tree,
      covariates  = cov_df,
      species_col = species_arg,
      epochs      = epochs,
      n_imputations = 5L,
      verbose     = FALSE,
      seed        = 1L,
      safety_floor = TRUE
    ),
    error = function(e) {
      log_line("    pigauto error: ", conditionMessage(e)); NULL })
  if (is.null(res)) return(rep(NA_real_, nrow(df_full)))

  # Extract response-trait predictions; pigauto's $completed has same shape
  # as the input traits_df.  Missing cells are now imputed.
  pred <- res$completed[[resp]]
  if (is.factor(pred)) pred <- as.numeric(as.character(pred))
  pred
}

# ---------------- Main sweep loop ------------------------------------------

results <- list()
cell_id <- 0L
total_cells <- length(N_SPECIES) * length(MULTI_OBS) * length(PHYLO_SIG) *
                length(BETA_STR) * length(RESP_TYPE) * length(N_PRED) * N_REPS
log_line("Total cells: ", total_cells)

for (rt in RESP_TYPE)
for (np in N_PRED)
for (ns in N_SPECIES)
for (mo in MULTI_OBS)
for (ps in PHYLO_SIG)
for (bs in BETA_STR)
for (rep_id in seq_len(N_REPS)) {
  cell_id <- cell_id + 1L
  cell_seed <- 1000L + 17L * cell_id
  log_line(sprintf("Cell %d/%d: response=%s, n_pred=%d, n=%d, multi=%d, ps=%.1f, beta=%.1f, rep=%d (seed=%d)",
                    cell_id, total_cells, rt, np, ns, mo, ps, bs, rep_id, cell_seed))

  # Generate the data
  sim <- tryCatch(
    sim_bace_dgp(n_species = ns, multi_obs_ratio = mo,
                  phylo_signal = ps, beta_resp_strength = bs,
                  response_type = rt, n_predictors = np,
                  miss_frac = MISS_FRAC, seed = cell_seed),
    error = function(e) { log_line("  sim_bace_dgp error: ", conditionMessage(e)); NULL })
  if (is.null(sim)) next

  truth <- sim$df_complete[[sim$response_name]]
  cell_meta <- list(cell_id = cell_id, response_type = rt, n_predictors = np,
                     n_species = ns, multi_obs = mo, phylo_signal = ps,
                     beta_strength = bs, rep_id = rep_id, seed = cell_seed,
                     n_held = sum(sim$mask), n_obs = nrow(sim$df_observed))

  # Method 1: column_mean
  t1 <- proc.time()[["elapsed"]]
  pred_cm <- fit_column_mean(sim)
  s_cm <- score(pred_cm, truth, sim$mask, rt)
  results[[length(results) + 1L]] <- c(cell_meta, list(method = "column_mean",
    rmse = s_cm[["rmse"]], acc = s_cm[["acc"]], r = s_cm[["r"]],
    wall_s = proc.time()[["elapsed"]] - t1))

  # Method 2: lm
  t1 <- proc.time()[["elapsed"]]
  pred_lm <- fit_lm(sim)
  s_lm <- score(pred_lm, truth, sim$mask, rt)
  results[[length(results) + 1L]] <- c(cell_meta, list(method = "lm",
    rmse = s_lm[["rmse"]], acc = s_lm[["acc"]], r = s_lm[["r"]],
    wall_s = proc.time()[["elapsed"]] - t1))

  # Method 3: phylolm-lambda BLUP (gaussian only)
  if (rt == "gaussian") {
    t1 <- proc.time()[["elapsed"]]
    pred_phy <- fit_phylolm_lambda_blup(sim)
    s_phy <- score(pred_phy, truth, sim$mask, rt)
    results[[length(results) + 1L]] <- c(cell_meta, list(method = "phylolm_lambda_blup",
      rmse = s_phy[["rmse"]], acc = s_phy[["acc"]], r = s_phy[["r"]],
      wall_s = proc.time()[["elapsed"]] - t1))
  }

  # Method 4: pigauto with covariates (Fix A-H)
  t1 <- proc.time()[["elapsed"]]
  pred_pig <- fit_pigauto_cov(sim)
  s_pig <- score(pred_pig, truth, sim$mask, rt)
  results[[length(results) + 1L]] <- c(cell_meta, list(method = "pigauto_cov_sfT",
    rmse = s_pig[["rmse"]], acc = s_pig[["acc"]], r = s_pig[["r"]],
    wall_s = proc.time()[["elapsed"]] - t1))

  # Save incremental results every 4 cells (in case bench is killed)
  if (cell_id %% 4L == 0L) {
    saveRDS(do.call(rbind, lapply(results, function(r) {
      d <- as.data.frame(r, stringsAsFactors = FALSE); d
    })), out_rds, compress = "xz")
    log_line(sprintf("  [checkpoint] saved %d rows", length(results)))
  }
}

# Final save
results_df <- do.call(rbind, lapply(results, function(r) {
  d <- as.data.frame(r, stringsAsFactors = FALSE); d
}))
saveRDS(results_df, out_rds, compress = "xz")
log_line(sprintf("=== DONE: %d cells, %d rows ===", cell_id, nrow(results_df)))
log_line("  rds: ", out_rds)
