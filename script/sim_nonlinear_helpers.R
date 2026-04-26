# script/sim_nonlinear_helpers.R
#
# Shared helpers for the nonlinear-DGP and transformer-ablation benches.
# Sourced by:
#   - script/bench_sim_bace_nonlinear.R
#   - script/bench_transformer_ablation.R
#
# Defines:
#   sim_nonlinear_dgp()        -- BACE-sim base + post-hoc nonlinear injection
#   score()                    -- RMSE + Pearson r on held-out cells
#   method_column_mean()       -- loss-floor baseline
#   method_lm()                -- linear OLS, no phylo
#   method_lm_nonlinear()      -- poly(2) + bilinear interaction OLS, no phylo
#   method_phylolm_blup()      -- phylolm-lambda BLUP, the linear-smart baseline
#
# Pure helpers; no top-level work.  Both calling scripts must already
# load pigauto + ape + phylolm + source script/sim_bace_dgp.R.
# Both must set EPOCHS and MISS_FRAC (or fall back to defaults below).

if (!exists("EPOCHS"))    EPOCHS    <- 150L
if (!exists("MISS_FRAC")) MISS_FRAC <- 0.30

#' Generate a sim_bace_dgp() cell, then post-hoc inject a nonlinear /
#' interactive contribution into the response.
#'
#' @param n_species      integer
#' @param multi_obs_ratio integer, 1 = single-obs, > 1 = multi-obs
#' @param phylo_signal   numeric in [0, 1)
#' @param beta_strength  numeric. SDs of nonlinear contribution to add.
#' @param f_type         "linear" | "nonlinear" | "interactive"
#' @param seed           integer
#' @return list with tree, df_complete, df_observed, mask, predictor_names,
#'         response_name, meta (incl. f_type and beta_strength)
sim_nonlinear_dgp <- function(n_species, multi_obs_ratio, phylo_signal,
                                beta_strength, f_type, seed) {
  base <- sim_bace_dgp(
    n_species          = n_species,
    multi_obs_ratio    = multi_obs_ratio,
    phylo_signal       = phylo_signal,
    beta_resp_strength = 0,    # zero linear y-on-x effect
    response_type      = "gaussian",
    n_predictors       = 3L,
    miss_frac          = MISS_FRAC,
    seed               = seed
  )
  pred_cols <- base$predictor_names
  X_full <- vapply(pred_cols, function(p) {
    v <- base$df_complete[[p]]
    if (is.factor(v)) as.numeric(v) - 1 else as.numeric(v)
  }, numeric(nrow(base$df_complete)))
  X_full <- as.matrix(X_full)

  contrib <- switch(f_type,
    linear      = X_full[, 1],
    nonlinear   = sin(2 * X_full[, 1]) * exp(0.3 * X_full[, 2]),
    interactive = X_full[, 1] * X_full[, 2] + 0.5 * X_full[, 1]^2,
    stop("Unknown f_type: ", f_type)
  )
  contrib <- as.numeric(scale(contrib))
  delta <- beta_strength * contrib

  base$df_complete$y <- base$df_complete$y + delta
  base$df_observed$y <- base$df_observed$y + delta
  base$df_observed$y[base$mask] <- NA

  base$meta$f_type        <- f_type
  base$meta$beta_strength <- beta_strength
  base
}

#' RMSE + Pearson r on held-out cells.
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
  if (length(pred_n) < 2L) return(method_lm(d))
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

  sp_df <- sp_df[sp_df$species %in% d$tree$tip.label, , drop = FALSE]
  sp_df_obs <- sp_df[!is.na(sp_df$y), , drop = FALSE]
  if (nrow(sp_df_obs) < length(pred_n) + 5L) {
    return(rep(NA_real_, nrow(ddf)))
  }
  rownames(sp_df_obs) <- sp_df_obs$species
  tree_obs <- ape::keep.tip(d$tree,
                              intersect(sp_df_obs$species, d$tree$tip.label))

  fmla <- stats::as.formula(paste("y ~", paste(pred_n, collapse = " + ")))
  fit <- tryCatch(
    phylolm::phylolm(fmla, data = sp_df_obs, phy = tree_obs, model = "lambda"),
    error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, nrow(ddf)))

  beta_hat <- coef(fit)
  # Build the prediction design matrix using a NO-LHS formula so model.matrix
  # doesn't drop rows where y is NA (it would otherwise leave `fixed` with the
  # wrong length and the subsequent miss_idx indexing returns NAs).
  fmla_rhs <- stats::as.formula(paste("~", paste(pred_n, collapse = " + ")))
  Xmat <- model.matrix(fmla_rhs, sp_df)
  fixed <- as.numeric(Xmat %*% beta_hat)

  R <- ape::vcv(d$tree); R <- stats::cov2cor(R); R <- R[sp_df$species, sp_df$species]
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
    sp_pred[obs_idx] <- sp_df$y[obs_idx]
    sp_pred[miss_idx] <- fixed[miss_idx] + blup
  }
  names(sp_pred) <- sp_df$species
  sp_pred[as.character(ddf$species)]
}
