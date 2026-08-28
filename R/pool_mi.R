#' Pool downstream model fits across multiple imputations (Rubin's rules)
#'
#' Combine regression coefficients from `M` model fits -- one per imputed
#' dataset -- into a single pooled table using Rubin's rules. The pooled
#' standard errors combine within-imputation sampling variance and
#' between-imputation variance. Correct Rubin arithmetic does not make an
#' incompatible imputation model inferentially valid. For pigauto inference,
#' fits should come from the documented [multi_impute_analysis()] workflow.
#' The backend is experimental and supports fixed-effect coefficients only. Variance
#' components, correlations, random-effect predictions, BLUPs/conditional
#' modes, latent loadings, and other structured parameters are unsupported.
#'
#' @param fits A list of model fits of length `M >= 2`. Automatic fixed-effect
#'   adapters are provided for [`stats::lm`], [`stats::glm`], `nlme::gls`,
#'   `nlme::lme`, `lme4::merMod`, `glmmTMB::glmmTMB`, `drmTMB`, and
#'   `gllvmTMB_multi` fits. Other classes implementing `coef()` and `vcov()`
#'   also work. The output of [with_imputations()] is accepted directly.
#'   `MCMCglmm` fits are rejected -- see Details.
#'   Bare user-supplied fit lists are accepted with an explicit warning because
#'   their imputation provenance cannot be verified.
#' @param conf.level Confidence level for the pooled interval (default
#'   `0.95`).
#' @param coef_fun Optional function extracting a named numeric fixed-effect
#'   vector from one fit. `NULL` uses the automatic class adapter. Custom
#'   coefficient and covariance extractors can be supplied independently;
#'   callers must ensure custom extractors select fixed effects only.
#' @param vcov_fun Optional function extracting the fixed-effect covariance
#'   matrix. `NULL` uses the automatic class adapter. Base matrices and
#'   `Matrix` objects are accepted; the selected covariance block must be
#'   square, finite, symmetric, and have non-negative diagonal entries.
#' @param df_fun Optional function returning the complete-data residual
#'   degrees of freedom `nu_com` from one fit. When supplied, pooled
#'   degrees of freedom use the Barnard & Rubin (1999) small-sample
#'   correction, which is less biased for short series. When `NULL` (the
#'   default) the classical Rubin (1987) formula is used.
#' @param tidy_fun Optional function returning a data.frame with unique
#'   `term`, numeric `estimate`, and non-negative finite `std.error` columns.
#'   This is an alternative to `coef_fun` and `vcov_fun`, not a supplement;
#'   combining them is an error. Callers must ensure it selects fixed effects
#'   only.
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
#' With no between-imputation variation (`B = 0`), the classical limit is
#' `df = Inf`, `riv = 0`, and `fmi = 0`. If finite complete-data df are
#' supplied, Barnard--Rubin instead retains finite observed-data df and its
#' small-sample FMI adjustment. A completely deterministic quantity
#' (`B = W = 0`) always has zero FMI and infinite df.
#'
#' Only fixed-effect coefficients are pooled in version 0.10.0. Random-effect
#' variances and correlations, BLUPs/conditional modes, latent loadings, and
#' other structured parameters require parameter-specific transformations and
#' are not supported by the automatic `pool_mi()` adapters. Custom extractors
#' are an expert escape hatch and cannot be inspected by pigauto; using them to
#' select unsupported structured parameters is outside the documented scope.
#' The `glmmTMB` adapter selects conditional fixed effects only. The `drmTMB`
#' adapter includes named distributional fixed-effect blocks such as regression
#' coefficients for `mu` and `sigma`; those are fixed coefficients, not
#' random-effect variance components.
#'
#' **MCMCglmm fits** are rejected because Rubin's rules are not the right
#' tool for posterior samples: variance decomposition does not generalise
#' cleanly to posterior distributions. No MCMCglmm downstream workflow is
#' supported by the initial analysis-aware backend.
#'
#' @references
#' Rubin DB (1987). *Multiple Imputation for Nonresponse in Surveys.*
#' Wiley.
#'
#' Barnard J, Rubin DB (1999). "Small-sample degrees of freedom with
#' multiple imputation." *Biometrika* 86(4): 948-955.
#'
#' Nakagawa S, Freckleton RP (2008). "Missing inaction: the dangers of
#' ignoring missing data." *Trends in Ecology & Evolution* 23(11):
#' 592-596.
#'
#' Nakagawa S, Freckleton RP (2011). "Model averaging, missing data and
#' multiple imputation: a case study for behavioural ecology."
#' *Behavioral Ecology and Sociobiology* 65(1): 103-116.
#'
#' @seealso [multi_impute_analysis()], [with_imputations()]
#'
#' @examples
#' \donttest{
#' dat <- data.frame(y = stats::rnorm(30L), x = stats::rnorm(30L),
#'                   z = stats::rnorm(30L))
#' dat$x[seq(3L, 30L, by = 5L)] <- NA_real_
#' mi <- multi_impute_analysis(
#'   data = dat, formula = y ~ x + z, missing = "x",
#'   model = "lm", m = 2L
#' )
#' fits <- with_imputations(mi, function(d) lm(y ~ x + z, data = d))
#' pool_mi(fits)
#' }
#'
#' @export
pool_mi <- function(fits,
                    conf.level = 0.95,
                    coef_fun   = NULL,
                    vcov_fun   = NULL,
                    df_fun     = NULL,
                    tidy_fun   = NULL) {

  # ---- Validate input ----
  if (!is.list(fits)) {
    stop("`fits` must be a list of model fits.", call. = FALSE)
  }

  workflow <- attr(fits, "mi_workflow")
  unverified_bare <- FALSE
  if (identical(workflow, "pigauto_diagnostic_mi") ||
      inherits(fits, "pigauto_diagnostic_mi")) {
    stop("Prediction-diagnostic completions from `multi_impute()` cannot be ",
         "pooled. Use `multi_impute_analysis()` for the supported narrow ",
         "workflow.", call. = FALSE)
  } else if (identical(workflow, "pigauto_tree_sensitivity_diagnostic") ||
             inherits(fits, "pigauto_tree_sensitivity_diagnostic") ||
             inherits(fits, "pigauto_mi_trees")) {
    stop("Tree-sensitivity diagnostics from `multi_impute_trees()` cannot be ",
         "pooled. Tree-aware downstream pooling is unsupported; use ",
         "`multi_impute_analysis()` for the supported narrow workflow.",
         call. = FALSE)
  } else if (inherits(fits, "pigauto_mi")) {
    stop("Legacy `pigauto_mi` fit lists have no analysis-aware provenance. ",
         "Refit with `multi_impute_analysis()` and `with_imputations()` ",
         "before pooling.", call. = FALSE)
  } else if (inherits(fits, "pigauto_mi_fits")) {
    if (is.null(workflow)) {
      stop("Legacy `pigauto_mi_fits` objects have no analysis-aware ",
           "provenance. Refit with `multi_impute_analysis()` and ",
           "`with_imputations()` before pooling.", call. = FALSE)
    }
    if (!identical(workflow, "pigauto_analysis_mi_v1")) {
      stop("`fits` carries unknown multiple-imputation provenance and cannot be ",
           "pooled safely.", call. = FALSE)
    }
  } else if (identical(workflow, "pigauto_analysis_mi_v1")) {
    stop("`mi_workflow = \"pigauto_analysis_mi_v1\"` is inconsistent ",
         "provenance for a fit list not created by `with_imputations()`.",
         call. = FALSE)
  } else if (is.null(workflow)) {
    unverified_bare <- TRUE
  } else {
    stop("`fits` carries unknown multiple-imputation provenance and cannot be ",
         "pooled safely.", call. = FALSE)
  }

  tree_index <- attr(fits, "tree_index")
  n_trees <- attr(fits, "n_trees")
  m_per_tree <- attr(fits, "m_per_tree")

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
    if (!is.null(tree_index)) tree_index <- tree_index[!is_err]
  }

  M <- length(fits)
  if (M < 2L) {
    stop("Need at least 2 fits to pool; got ", M, ".", call. = FALSE)
  }

  # MCMCglmm detection -- Rubin's rules don't apply cleanly to posteriors.
  if (inherits(fits[[1]], "MCMCglmm")) {
    stop(
      "pool_mi() applies Rubin's rules and is not appropriate for ",
      "MCMCglmm fits. For a Bayesian pigauto workflow, concatenate ",
      "posterior samples across imputations directly: ",
      "`coda::as.mcmc(do.call(rbind, lapply(fits, function(f) f$Sol)))`. ",
      "For the frequentist Rubin's-rules path, see ",
      "`vignette('mixed-types', package = 'pigauto')`.",
      call. = FALSE
    )
  }

  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("`conf.level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }

  if (!is.null(tidy_fun) &&
      (!is.null(coef_fun) || !is.null(vcov_fun))) {
    stop("`tidy_fun` cannot be combined with `coef_fun` or `vcov_fun`.",
         call. = FALSE)
  }
  if (!is.null(coef_fun) && !is.function(coef_fun)) {
    stop("`coef_fun` must be a function or `NULL`.", call. = FALSE)
  }
  if (!is.null(vcov_fun) && !is.function(vcov_fun)) {
    stop("`vcov_fun` must be a function or `NULL`.", call. = FALSE)
  }
  if (!is.null(df_fun) && !is.function(df_fun)) {
    stop("`df_fun` must be a function or `NULL`.", call. = FALSE)
  }
  if (!is.null(tidy_fun) && !is.function(tidy_fun)) {
    stop("`tidy_fun` must be a function or `NULL`.", call. = FALSE)
  }

  # ---- Extract coefficients and variance-covariance ----
  if (!is.null(tidy_fun)) {
    tidy_tables <- lapply(seq_len(M), function(i) {
      .pool_mi_validate_tidy(tidy_fun(fits[[i]]), i)
    })
    coefs <- lapply(tidy_tables, function(x) {
      stats::setNames(x$estimate, x$term)
    })
  } else {
    coef_extractor <- if (is.null(coef_fun)) .pool_mi_auto_coef else coef_fun
    vcov_extractor <- if (is.null(vcov_fun)) .pool_mi_auto_vcov else vcov_fun
    coefs <- lapply(fits, coef_extractor)
  }

  # Sanity: all elements must be numeric vectors with consistent names.
  coef_classes <- vapply(coefs, function(x) {
    is.numeric(x) && length(x) > 0L && !is.null(names(x)) &&
      all(nzchar(names(x))) && !anyDuplicated(names(x)) && all(is.finite(x))
  }, logical(1))
  if (!all(coef_classes)) {
    stop("The coefficient extractor must return a finite, uniquely named ",
         "numeric vector for every fit.",
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
  if (!is.null(tidy_fun)) {
    vars_mat[,] <- do.call(rbind, lapply(tidy_tables, function(x) {
      stats::setNames(x$std.error^2, x$term)[nm_ref]
    }))
  } else {
    for (i in seq_len(M)) {
      V <- .pool_mi_validate_vcov(vcov_extractor(fits[[i]]), nm_ref, i)
      vars_mat[i, ] <- diag(V)
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

  # Relative increase in variance, including both W == 0 limits. A fully
  # deterministic quantity (W == B == 0) has no missing information; a
  # quantity with W == 0 and B > 0 has all of its uncertainty between draws.
  r <- rep(NA_real_, length(W))
  r[W > 0] <- (1 + 1 / M) * B[W > 0] / W[W > 0]
  r[W == 0 & B > 0] <- Inf
  r[W == 0 & B == 0] <- 0
  lambda <- ifelse(total_var > 0,
                   (1 + 1 / M) * B / total_var,
                   0)

  # ---- Degrees of freedom ----
  # Classical Rubin (1987): v_old = (M-1)(1 + 1/r)^2
  #   With r -> Inf (all variance between-imputation) this collapses to
  #   v_old = M - 1, which is the right limit for pure MC noise.
  v_old <- rep(Inf, length(r))
  positive_r <- is.finite(r) & r > 0
  v_old[positive_r] <- (M - 1) * (1 + 1 / r[positive_r])^2
  v_old[is.infinite(r)] <- M - 1

  # Barnard-Rubin (1999) refinement if complete-data df is available.
  if (!is.null(df_fun)) {
    nu_com <- df_fun(fits[[1]])
    if (!is.numeric(nu_com) || length(nu_com) != 1L ||
        !is.finite(nu_com) || nu_com <= 0) {
      stop("`df_fun()` must return one finite positive number.",
           call. = FALSE)
    }
    v_obs <- ((nu_com + 1) / (nu_com + 3)) * nu_com * (1 - lambda)
    v_bar <- v_old
    has_observed_df <- v_obs > 0 & total_var > 0
    v_bar[has_observed_df] <-
      1 / (1 / v_old[has_observed_df] + 1 / v_obs[has_observed_df])
  } else {
    v_bar <- v_old
  }

  # ---- Fraction of missing information and inference ----
  fmi <- rep(1, length(r))
  finite_r <- is.finite(r)
  fmi[finite_r] <-
    (r[finite_r] + 2 / (v_bar[finite_r] + 3)) / (r[finite_r] + 1)
  fmi[total_var == 0] <- 0

  t_stat <- ifelse(se_pool > 0, theta_bar / se_pool, NA_real_)
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
  if (!is.null(workflow)) attr(out, "mi_workflow") <- workflow
  if (!is.null(tree_index)) {
    attr(out, "tree_index") <- tree_index
    attr(out, "n_trees") <- n_trees
    attr(out, "m_per_tree") <- m_per_tree
  }
  class(out) <- c("pigauto_pooled", "data.frame")
  if (unverified_bare) {
    warning("pool_mi() received a bare fit list with unverified provenance; ",
            "Rubin arithmetic will run, but the imputation workflow was not ",
            "validated by pigauto.", call. = FALSE)
  }
  out
}


.pool_mi_auto_coef <- function(fit) {
  if (inherits(fit, "gllvmTMB_multi")) {
    td <- .pool_mi_gllvm_tidy(fit)
    return(stats::setNames(td$estimate, td$term))
  }

  if (inherits(fit, "drmTMB")) {
    blocks <- stats::coef(fit)
    if (!is.list(blocks) || is.null(names(blocks))) {
      stop("The drmTMB fixed-effect adapter could not find named ",
           "distributional coefficient blocks.", call. = FALSE)
    }
    values <- unlist(blocks, use.names = FALSE)
    block_names <- rep(names(blocks), lengths(blocks))
    term_names <- unlist(lapply(blocks, names), use.names = FALSE)
    names(values) <- paste0(block_names, ":", term_names)
    return(values)
  }

  if (inherits(fit, "glmmTMB")) {
    return(glmmTMB::fixef(fit)$cond)
  }
  if (inherits(fit, "merMod")) {
    return(lme4::fixef(fit))
  }
  if (inherits(fit, "lme")) {
    return(nlme::fixef(fit))
  }

  stats::coef(fit)
}


.pool_mi_auto_vcov <- function(fit) {
  if (inherits(fit, "gllvmTMB_multi")) {
    td <- .pool_mi_gllvm_tidy(fit)
    V <- diag(td$std.error^2, nrow = nrow(td))
    dimnames(V) <- list(td$term, td$term)
    return(V)
  }
  if (inherits(fit, "glmmTMB")) {
    return(stats::vcov(fit)$cond)
  }
  stats::vcov(fit)
}


.pool_mi_gllvm_tidy <- function(fit) {
  package <- sub("_multi$", "", class(fit)[[1L]])
  namespace <- tryCatch(
    asNamespace(package),
    error = function(e) NULL
  )
  method <- if (is.null(namespace)) NULL else utils::getS3method(
    "tidy", "gllvmTMB_multi", optional = TRUE, envir = namespace
  )
  if (is.null(method)) {
    stop("The gllvmTMB fixed-effect adapter requires the package that created ",
         "the fit to be installed and to provide tidy.gllvmTMB_multi().",
         call. = FALSE)
  }
  .pool_mi_validate_tidy(method(fit, effects = "fixed"), "gllvmTMB")
}


.pool_mi_validate_tidy <- function(x, fit_index) {
  label <- paste0("fit ", fit_index)
  required <- c("term", "estimate", "std.error")
  if (!is.data.frame(x) || !all(required %in% names(x)) || nrow(x) == 0L) {
    stop("`tidy_fun()` must return a non-empty data.frame with `term`, ",
         "`estimate`, and `std.error` columns for ", label, ".",
         call. = FALSE)
  }
  if (!is.character(x$term) || anyNA(x$term) || any(!nzchar(x$term)) ||
      anyDuplicated(x$term)) {
    stop("`tidy_fun()` returned missing, empty, or duplicate terms for ",
         label, ".", call. = FALSE)
  }
  if (!is.numeric(x$estimate) || any(!is.finite(x$estimate)) ||
      !is.numeric(x$std.error) || any(!is.finite(x$std.error)) ||
      any(x$std.error < 0)) {
    stop("`tidy_fun()` must return finite numeric estimates and finite ",
         "non-negative standard errors for ", label, ".", call. = FALSE)
  }
  x[, required, drop = FALSE]
}


.pool_mi_validate_vcov <- function(V, terms, fit_index) {
  if (!is.matrix(V) && !inherits(V, "Matrix")) {
    stop("The covariance extractor must return a base matrix or Matrix ",
         "object for fit ", fit_index, ".", call. = FALSE)
  }
  V <- as.matrix(V)
  if (!is.numeric(V) || nrow(V) != ncol(V)) {
    stop("The covariance matrix must be numeric and square for fit ",
         fit_index, ".", call. = FALSE)
  }

  rn <- rownames(V)
  cn <- colnames(V)
  if (xor(is.null(rn), is.null(cn))) {
    stop("The covariance matrix must have both row and column names, or ",
         "neither, for fit ", fit_index, ".", call. = FALSE)
  }
  if (!is.null(rn)) {
    if (anyDuplicated(rn) || anyDuplicated(cn) || !setequal(rn, cn) ||
        !all(terms %in% rn)) {
      stop("Covariance row/column names do not align with coefficient ",
           "terms for fit ", fit_index, ".", call. = FALSE)
    }
    V <- V[terms, terms, drop = FALSE]
  } else if (nrow(V) != length(terms)) {
    stop("An unnamed covariance matrix must have exactly one row per ",
         "coefficient for fit ", fit_index, ".", call. = FALSE)
  } else {
    dimnames(V) <- list(terms, terms)
  }

  if (any(!is.finite(V))) {
    stop("The selected covariance block contains non-finite values for fit ",
         fit_index, ".", call. = FALSE)
  }
  symmetry_error <- max(abs(V - t(V)))
  symmetry_scale <- max(1, max(abs(V)))
  if (symmetry_error > 1e-8 * symmetry_scale) {
    stop("The selected covariance block is not symmetric for fit ",
         fit_index, ".", call. = FALSE)
  }
  if (any(diag(V) < 0)) {
    stop("The covariance diagonal contains negative values for fit ",
         fit_index, ".", call. = FALSE)
  }
  V
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
