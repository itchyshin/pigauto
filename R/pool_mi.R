#' Pool downstream model fits across multiple imputations (Rubin's rules)
#'
#' Combine regression coefficients from `M` model fits — one per imputed
#' dataset — into a single pooled table using Rubin's rules. The pooled
#' standard errors properly account for *both* within-imputation sampling
#' variance and between-imputation variance, so downstream inference
#' propagates the uncertainty introduced by imputation.
#'
#' @param fits A list of model fits of length `M >= 2`. Any model class
#'   implementing `coef()` and `vcov()` works (e.g. [`stats::lm`],
#'   [`nlme::gls`], [`lme4::lmer`], `glmmTMB::glmmTMB`, `phylolm::phylolm`,
#'   `phylolm::phyloglm`). The output of [with_imputations()] is accepted
#'   directly. `MCMCglmm` fits are rejected — see Details.
#' @param conf.level Confidence level for the pooled interval (default
#'   `0.95`).
#' @param coef_fun Function extracting a named numeric coefficient vector
#'   from one fit. Defaults to [stats::coef()]. Supply a custom extractor
#'   for models where `coef()` returns a list (e.g. `glmmTMB` with multiple
#'   components) — see Examples.
#' @param vcov_fun Function extracting the variance-covariance matrix of
#'   the fixed-effect coefficients. Defaults to [stats::vcov()]. Must
#'   return a square matrix whose row/column names match `coef_fun(fit)`.
#' @param df_fun Optional function returning the complete-data residual
#'   degrees of freedom `nu_com` from one fit. When supplied, pooled
#'   degrees of freedom use the Barnard & Rubin (1999) small-sample
#'   correction, which is less biased for short series. When `NULL` (the
#'   default) the classical Rubin (1987) formula is used.
#'
#' @return A data.frame with one row per coefficient and columns:
#'   \describe{
#'     \item{`term`}{Coefficient name.}
#'     \item{`estimate`}{Pooled point estimate (`mean` across fits).}
#'     \item{`std.error`}{Pooled standard error `sqrt(T)` where
#'       `T = W + (1 + 1/M) * B`.}
#'     \item{`df`}{Pooled degrees of freedom (Barnard-Rubin if `df_fun`
#'       supplied, else classical Rubin).}
#'     \item{`statistic`}{`estimate / std.error`.}
#'     \item{`p.value`}{Two-sided p-value from a t distribution on `df`.}
#'     \item{`conf.low`, `conf.high`}{Pooled `conf.level` interval.}
#'     \item{`fmi`}{Fraction of missing information.}
#'     \item{`riv`}{Relative increase in variance due to non-response.}
#'   }
#'
#' @details
#' Let \eqn{\hat\theta_i} be the coefficient vector from fit *i* and
#' \eqn{U_i = \mathrm{vcov}(\mathrm{fit}_i)}, for \eqn{i = 1, \ldots, M}.
#' Rubin's rules (Rubin 1987) give
#' \deqn{\bar\theta = M^{-1} \sum_i \hat\theta_i}
#' \deqn{W = M^{-1} \sum_i \mathrm{diag}(U_i)}
#' \deqn{B = (M-1)^{-1} \sum_i (\hat\theta_i - \bar\theta)^2}
#' \deqn{T = W + (1 + 1/M) B}
#' with pooled standard error \eqn{\sqrt{T}}. The relative increase in
#' variance is \eqn{r = (1 + 1/M) B / W}, the classical pooled df is
#' \eqn{\nu_{\text{old}} = (M-1)(1 + 1/r)^2}, and the fraction of missing
#' information is
#' \deqn{\mathrm{fmi} = (r + 2/(\nu + 3)) / (r + 1).}
#' When `df_fun` returns finite complete-data df `nu_com`, the
#' Barnard-Rubin (1999) correction combines
#' \eqn{\nu_{\text{obs}} = ((\nu_{\text{com}}+1)/(\nu_{\text{com}}+3))
#' \nu_{\text{com}} (1 - \lambda)} with `nu_old` via
#' \eqn{\nu_{\text{BR}} = 1/(1/\nu_{\text{old}} + 1/\nu_{\text{obs}})}.
#'
#' **MCMCglmm fits** are rejected because Rubin's rules are not the right
#' tool for posterior samples: variance decomposition does not generalise
#' cleanly to posterior distributions. For a Bayesian pigauto workflow
#' (pigauto as imputer, MCMCglmm as inference engine), concatenate the
#' posterior samples across imputations manually — stack `fit$Sol` and
#' `fit$VCV` row-wise with `do.call(rbind, ...)` and wrap the result in
#' `coda::as.mcmc()`. See section 7 of the `pigauto_workflow_mixed`
#' tutorial (`system.file("doc", "pigauto_workflow_mixed.html",
#' package = "pigauto")`) for a worked example. For an integrated
#' chained-equation MCMC workflow where imputation and inference happen
#' in a single engine, use the companion BACE package end-to-end
#' (`BACE::bace()` + `BACE::pool_posteriors()`).
#'
#' @references
#' Rubin DB (1987). *Multiple Imputation for Nonresponse in Surveys.*
#' Wiley.
#'
#' Barnard J, Rubin DB (1999). "Small-sample degrees of freedom with
#' multiple imputation." *Biometrika* 86(4): 948–955.
#'
#' Nakagawa S, Freckleton RP (2008). "Missing inaction: the dangers of
#' ignoring missing data." *Trends in Ecology & Evolution* 23(11):
#' 592–596.
#'
#' Nakagawa S, Freckleton RP (2011). "Model averaging, missing data and
#' multiple imputation: a case study for behavioural ecology."
#' *Behavioral Ecology and Sociobiology* 65(1): 103–116.
#'
#' @seealso [multi_impute()], [with_imputations()]
#'
#' @examples
#' \dontrun{
#' # Typical workflow
#' mi   <- multi_impute(traits, tree, m = 100)
#' fits <- with_imputations(mi, function(d) lm(y ~ x1 + x2, data = d))
#' pool_mi(fits)
#'
#' # glmmTMB with custom extractors (conditional component only)
#' pool_mi(
#'   fits,
#'   coef_fun = function(f) glmmTMB::fixef(f)$cond,
#'   vcov_fun = function(f) vcov(f)$cond
#' )
#' }
#'
#' @export
pool_mi <- function(fits,
                    conf.level = 0.95,
                    coef_fun   = stats::coef,
                    vcov_fun   = stats::vcov,
                    df_fun     = NULL) {

  # ---- Validate input ----
  if (!is.list(fits)) {
    stop("`fits` must be a list of model fits.", call. = FALSE)
  }

  # Drop and warn about captured errors (from with_imputations() with
  # error handling on). These are stored as try-error objects.
  is_err <- vapply(fits, function(x) inherits(x, "try-error") ||
                     inherits(x, "pigauto_mi_error"),
                   logical(1))
  if (any(is_err)) {
    warning(sprintf("Dropping %d fit%s that failed in with_imputations().",
                    sum(is_err), if (sum(is_err) == 1L) "" else "s"),
            call. = FALSE)
    fits <- fits[!is_err]
  }

  M <- length(fits)
  if (M < 2L) {
    stop("Need at least 2 fits to pool; got ", M, ".", call. = FALSE)
  }

  # MCMCglmm detection — Rubin's rules don't apply cleanly to posteriors.
  if (inherits(fits[[1]], "MCMCglmm")) {
    stop(
      "pool_mi() applies Rubin's rules and is not appropriate for ",
      "MCMCglmm fits. For a Bayesian pigauto workflow, concatenate ",
      "posterior samples across imputations directly: ",
      "`coda::as.mcmc(do.call(rbind, lapply(fits, function(f) f$Sol)))`. ",
      "See the 'pigauto_workflow_mixed' tutorial section 7 for a worked ",
      "example. For an integrated chained-equation MCMC workflow, use ",
      "BACE::bace() + BACE::pool_posteriors() end-to-end.",
      call. = FALSE
    )
  }

  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      conf.level <= 0 || conf.level >= 1) {
    stop("`conf.level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }

  # ---- Extract coefficients and variance-covariance ----
  coefs <- lapply(fits, coef_fun)

  # Sanity: all elements must be numeric vectors with consistent names.
  coef_classes <- vapply(coefs, function(x) is.numeric(x) && !is.null(names(x)),
                         logical(1))
  if (!all(coef_classes)) {
    stop("`coef_fun()` must return a named numeric vector for every fit. ",
         "If your model returns a list (e.g. glmmTMB), supply a custom ",
         "`coef_fun` — see the Examples section of ?pool_mi.",
         call. = FALSE)
  }

  nm_ref <- names(coefs[[1]])
  nm_ok  <- vapply(coefs, function(x) identical(names(x), nm_ref),
                   logical(1))
  if (!all(nm_ok)) {
    # Diagnostic: show which fits deviated.
    offenders <- which(!nm_ok)
    stop(
      "Coefficient names differ across fits. Rubin's rules require a ",
      "common set of terms. First offending fit: index ", offenders[1],
      ". Reference names: ", paste(nm_ref, collapse = ", "),
      ". Offender names: ", paste(names(coefs[[offenders[1]]]),
                                  collapse = ", "),
      call. = FALSE
    )
  }

  # M x p matrix of coefficients
  coef_mat <- do.call(rbind, lapply(coefs, function(x) unname(x[nm_ref])))

  # Per-fit diagonal of vcov (within-imputation variance of each coef)
  vars_mat <- matrix(NA_real_, nrow = M, ncol = length(nm_ref))
  for (i in seq_len(M)) {
    V <- vcov_fun(fits[[i]])
    if (is.null(V) || !is.matrix(V) || any(dim(V) < length(nm_ref))) {
      stop("`vcov_fun()` did not return a valid matrix for fit ", i, ".",
           call. = FALSE)
    }
    diag_V <- diag(V)
    # Match by name if possible, else fall back to position.
    if (!is.null(names(diag_V)) && all(nm_ref %in% names(diag_V))) {
      vars_mat[i, ] <- diag_V[nm_ref]
    } else if (length(diag_V) == length(nm_ref)) {
      vars_mat[i, ] <- diag_V
    } else {
      stop("`vcov_fun()` returned a matrix of the wrong size for fit ", i,
           " (", length(diag_V), " rows vs ", length(nm_ref),
           " coefficients).", call. = FALSE)
    }
  }

  # ---- Rubin's rules ----
  theta_bar <- colMeans(coef_mat)
  W <- colMeans(vars_mat)                         # within-imputation
  # Between-imputation variance: apply(var) uses (M-1) denominator.
  if (M > 1L) {
    B <- apply(coef_mat, 2, stats::var)
  } else {
    B <- rep(0, length(theta_bar))
  }
  total_var <- W + (1 + 1 / M) * B                # total variance T
  se_pool   <- sqrt(total_var)

  # Relative increase in variance. Guard against W == 0 (can happen when
  # a coefficient has zero within-imputation variance, e.g. a constrained
  # parameter); treat r as Inf in that case, which makes lambda == 1 and
  # fmi == 1.
  r <- ifelse(W > 0, (1 + 1 / M) * B / W, Inf)
  lambda <- ifelse(total_var > 0, (1 + 1 / M) * B / total_var, 1)

  # ---- Degrees of freedom ----
  # Classical Rubin (1987): v_old = (M-1)(1 + 1/r)^2
  #   With r -> Inf (all variance between-imputation) this collapses to
  #   v_old = M - 1, which is the right limit for pure MC noise.
  v_old <- ifelse(
    is.finite(r) & r > 0,
    (M - 1) * (1 + 1 / r)^2,
    M - 1
  )

  # Barnard-Rubin (1999) refinement if complete-data df is available.
  if (!is.null(df_fun)) {
    nu_com <- tryCatch(df_fun(fits[[1]]), error = function(e) NA_real_)
    if (is.numeric(nu_com) && length(nu_com) == 1L && is.finite(nu_com) &&
        nu_com > 0) {
      v_obs <- ((nu_com + 1) / (nu_com + 3)) * nu_com * (1 - lambda)
      v_bar <- 1 / (1 / v_old + 1 / v_obs)
    } else {
      v_bar <- v_old
    }
  } else {
    v_bar <- v_old
  }

  # ---- Fraction of missing information and inference ----
  fmi <- ifelse(
    is.finite(r),
    (r + 2 / (v_bar + 3)) / (r + 1),
    1
  )

  t_stat <- theta_bar / se_pool
  # Two-sided p-value. Guard against se_pool == 0 (constant coefficient).
  p_val <- ifelse(
    se_pool > 0,
    2 * stats::pt(-abs(t_stat), df = v_bar),
    NA_real_
  )
  alpha <- 1 - conf.level
  q <- stats::qt(1 - alpha / 2, df = v_bar)
  conf_lo <- theta_bar - q * se_pool
  conf_hi <- theta_bar + q * se_pool

  out <- data.frame(
    term      = nm_ref,
    estimate  = theta_bar,
    std.error = se_pool,
    df        = v_bar,
    statistic = t_stat,
    p.value   = p_val,
    conf.low  = conf_lo,
    conf.high = conf_hi,
    fmi       = fmi,
    riv       = r,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  attr(out, "m") <- M
  attr(out, "conf.level") <- conf.level
  class(out) <- c("pigauto_pooled", "data.frame")
  out
}


#' @export
print.pigauto_pooled <- function(x, digits = 4, ...) {
  M <- attr(x, "m")
  cl <- attr(x, "conf.level")
  cat(sprintf("Pooled estimates from %d multiply-imputed fits (Rubin's rules)\n",
              M))
  cat(sprintf("Confidence level: %.0f%%\n\n", 100 * cl))

  # Pretty-print with rounded numerics.
  df_show <- x
  num_cols <- vapply(df_show, is.numeric, logical(1))
  df_show[num_cols] <- lapply(df_show[num_cols], function(v) {
    formatC(v, digits = digits, format = "fg", flag = "#")
  })
  class(df_show) <- "data.frame"
  print(df_show, row.names = FALSE)

  # Interpretation hint on fraction of missing information.
  max_fmi <- suppressWarnings(max(x$fmi, na.rm = TRUE))
  if (is.finite(max_fmi) && max_fmi > 0.5) {
    cat(sprintf(
      "\nNote: max fmi = %.2f. Consider M >= %d imputations for stable SEs.\n",
      max_fmi, max(100L, ceiling(100 * max_fmi))
    ))
  }
  invisible(x)
}
