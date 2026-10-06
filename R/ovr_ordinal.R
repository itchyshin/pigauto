# PROTOTYPE (feat/discrete-lambda-ordinal, opt-in, default off): cumulative
# decomposition of an ordinal trait. K ordered classes -> K-1 synthetic
# binary traits "class >= k", each fitted through the same threshold-joint
# machinery the one-vs-rest categorical path uses (so each gets its own
# lambda under options(pigauto.discrete_lambda = "estimate")). Enabled with
# options(pigauto.ordinal_method = "cumulative").

#' Build a pigauto_data whose ordinal column is replaced by "class >= k"
#'
#' The ordinal latent column is z-scored integer class; the integer class is
#' recovered as round(z * sd + mean). The synthetic binary column is 0/1 with
#' NA preserved, written in place of the ordinal column.
#' @keywords internal
#' @noRd
build_cumulative_pd <- function(data, tm_ord, k) {
  col <- tm_ord$latent_cols
  cls <- round(data$X_scaled[, col] * tm_ord$sd + tm_ord$mean)
  pd_k <- data
  pd_k$X_scaled[, col] <- ifelse(is.na(cls), NA_real_, as.numeric(cls >= k))
  new_tm <- lapply(data$trait_map, function(tm) {
    if (identical(tm$latent_cols, tm_ord$latent_cols) &&
        identical(tm$type, "ordinal")) {
      list(name = paste0(tm_ord$name %||% "ord", "_ge_", k),
           type = "binary", n_latent = 1L, latent_cols = col,
           levels = c("no", "yes"), mean = 0, sd = 1)
    } else tm
  })
  pd_k$trait_map <- new_tm
  attr(pd_k, "synthetic_bin_col") <- col
  pd_k
}

#' Fit K-1 cumulative binary threshold-joint baselines for an ordinal trait
#'
#' @param tm_ord trait_map entry of the ordinal trait (type "ordinal").
#' @return numeric matrix (n_species x (K-1)); column k is P(class >= k)
#'   for k = 1..K-1 (classes coded 0..K-1); NA where a fit failed. NULL if
#'   every fit failed.
#' @keywords internal
#' @noRd
fit_cumulative_ordinal_fits <- function(data, tree, tm_ord,
                                         splits = NULL, graph = NULL,
                                         soft_aggregate = FALSE,
                                         joint_solver = "inhouse",
                                         predict_method = "exact",
                                         joint_refine_iter = 0L,
                                         predict_method_explicit = NULL) {
  if (is.null(predict_method_explicit)) {
    predict_method_explicit <- !missing(predict_method)
  }
  stopifnot(joint_mvn_available())
  K <- length(tm_ord$levels)
  if (K < 2L) return(NULL)
  spp <- if (!is.null(data$species_names)) data$species_names else rownames(data$X_scaled)
  probs <- matrix(NA_real_, nrow = length(spp), ncol = K - 1L,
                  dimnames = list(spp, paste0("ge_", seq_len(K - 1L))))
  for (k in seq_len(K - 1L)) {
    pd_k <- build_cumulative_pd(data, tm_ord, k)
    jt <- tryCatch(
      fit_joint_threshold_baseline(pd_k, tree, splits = splits, graph = graph,
                                    soft_aggregate = soft_aggregate,
                                    joint_solver = joint_solver,
                                    predict_method = predict_method,
                                    joint_refine_iter = joint_refine_iter,
                                    predict_method_explicit = predict_method_explicit),
      error = function(e) NULL
    )
    if (is.null(jt)) next
    bin_idx <- which(jt$liab_types == "binary" &
                       jt$liab_cols == attr(pd_k, "synthetic_bin_col"))
    if (length(bin_idx) == 0L) next
    dec <- decode_binary_liability(mu_liab = jt$mu_liab[, bin_idx],
                                    se_liab = jt$se_liab[, bin_idx])
    probs[, k] <- dec$p
  }
  if (all(is.na(probs))) return(NULL)
  probs
}

#' Turn cumulative probabilities P(Y >= k) into class probabilities
#'
#' Failed (all-NA) columns are filled by linear interpolation in k between
#' neighbours (1 at k = 0, 0 at k = K). Monotone non-increasing in k via
#' cumulative minimum; P(Y = k) = P(Y>=k) - P(Y>=k+1), clipped at 0 and
#' renormalised.
#' @param cum_probs n x (K-1) matrix of P(Y >= k), k = 1..K-1.
#' @return list(probs = n x K matrix, class = integer 0..K-1 most probable)
#' @keywords internal
#' @noRd
decode_cumulative_ordinal <- function(cum_probs) {
  n <- nrow(cum_probs); K <- ncol(cum_probs) + 1L
  probs <- matrix(NA_real_, n, K)
  for (i in seq_len(n)) {
    g <- c(1, cum_probs[i, ], 0)          # P(Y >= 0..K)
    ok <- !is.na(g)
    if (!all(ok)) g <- stats::approx(which(ok), g[ok], seq_along(g))$y
    g <- pmin(1, pmax(0, g))
    g <- cummin(g)
    p <- pmax(g[-length(g)] - g[-1], 0)
    s <- sum(p)
    probs[i, ] <- if (s > 0) p / s else rep(1 / K, K)
  }
  list(probs = probs, class = max.col(probs, ties.method = "first") - 1L)
}
