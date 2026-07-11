# Harness-only positive controls for the v0.10.0 MI validation.

.numeric_outcome <- function(y) {
  if (is.factor(y)) return(as.numeric(y) - 1)
  as.numeric(y)
}

.capture_warnings <- function(expr) {
  warnings <- character(0)
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings))
}

.validate_completed_datasets <- function(datasets, dgp, m, engine) {
  if (!is.list(datasets) || length(datasets) != as.integer(m)) {
    stop(engine, " returned ", length(datasets), " datasets; expected ", m,
         ".", call. = FALSE)
  }
  missing <- !is.finite(dgp$observed$x)
  observed <- !missing
  for (i in seq_along(datasets)) {
    data <- datasets[[i]]
    if (!is.data.frame(data) || nrow(data) != nrow(dgp$observed)) {
      stop(engine, " dataset ", i, " does not preserve row count.",
           call. = FALSE)
    }
    if (!("x" %in% names(data)) || any(!is.finite(data$x))) {
      stop(engine, " dataset ", i, " has missing or non-finite x.",
           call. = FALSE)
    }
    if (!isTRUE(all.equal(data$x[observed], dgp$observed$x[observed],
                          tolerance = 0))) {
      stop(engine, " dataset ", i, " changed observed x values.",
           call. = FALSE)
    }
    if (dgp$dgp == "lmer" &&
        !identical(as.character(data$species),
                   as.character(dgp$observed$species))) {
      stop(engine, " dataset ", i, " changed species row order.",
           call. = FALSE)
    }
  }
  invisible(datasets)
}

.draw_oracle_conditional_imputations <- function(dgp, m, seed) {
  m <- as.integer(m)
  if (length(m) != 1L || is.na(m) || m < 2L) {
    stop("m must be an integer >= 2.", call. = FALSE)
  }
  missing <- !is.finite(dgp$observed$x)
  if (!any(missing)) return(rep(list(dgp$observed), m))

  mu_prior <- dgp$x_prior_mean[missing]
  sd_prior <- rep(dgp$x_prior_sd, sum(missing))
  y <- .numeric_outcome(dgp$observed$y[missing])
  z <- dgp$observed$z[missing]
  species <- if (dgp$dgp == "lmer") dgp$observed$species[missing] else NULL

  draw_one <- function(i) {
    set.seed(as.integer(seed + i))
    if (dgp$dgp %in% c("lm", "lmer")) {
      beta_x <- 0.60
      beta_z <- -0.40
      random_intercept <- if (dgp$dgp == "lmer") {
        unname(dgp$random_intercept[species])
      } else {
        rep(0, length(y))
      }
      outcome_offset <- 1 + beta_z * z + random_intercept
      posterior_var <- 1 / (1 / sd_prior^2 + beta_x^2)
      posterior_mean <- posterior_var * (
        mu_prior / sd_prior^2 + beta_x * (y - outcome_offset)
      )
      x_draw <- stats::rnorm(
        length(y), mean = posterior_mean, sd = sqrt(posterior_var)
      )
    } else {
      x_draw <- numeric(length(y))
      for (j in seq_along(y)) {
        accepted <- FALSE
        while (!accepted) {
          proposal <- stats::rnorm(1L, mu_prior[j], sd_prior[j])
          probability <- stats::plogis(-0.20 + 0.80 * proposal - 0.50 * z[j])
          likelihood <- if (y[j] == 1) probability else 1 - probability
          accepted <- stats::runif(1L) <= likelihood
        }
        x_draw[j] <- proposal
      }
    }
    completed <- dgp$observed
    completed$x[missing] <- x_draw
    completed
  }

  lapply(seq_len(m), draw_one)
}

.draw_bayes_norm_lm_imputations <- function(dgp, m, seed) {
  if (dgp$dgp != "lm") {
    stop("The analytic normal-regression imputer is only defined for lm.",
         call. = FALSE)
  }
  m <- as.integer(m)
  if (length(m) != 1L || is.na(m) || m < 2L) {
    stop("m must be an integer >= 2.", call. = FALSE)
  }
  data <- dgp$observed
  observed <- is.finite(data$x)
  missing <- !observed
  design <- stats::model.matrix(~ y + z + I(z^2), data = data)
  q_observed <- design[observed, , drop = FALSE]
  q_missing <- design[missing, , drop = FALSE]
  x_observed <- data$x[observed]
  qr_fit <- qr(q_observed)
  if (qr_fit$rank != ncol(q_observed)) {
    stop("The analytic lm imputation design is rank deficient.", call. = FALSE)
  }
  beta_hat <- as.numeric(qr.coef(qr_fit, x_observed))
  residual <- x_observed - as.numeric(q_observed %*% beta_hat)
  residual_df <- nrow(q_observed) - ncol(q_observed)
  if (residual_df <= 0L) {
    stop("The analytic lm imputer has no residual degrees of freedom.",
         call. = FALSE)
  }
  rss <- sum(residual^2)
  if (!is.finite(rss) || rss <= 0) {
    stop("The analytic lm imputer requires positive finite residual variance.",
         call. = FALSE)
  }
  # The QR fit supplies rank detection and stable coefficients; Cholesky gives
  # the covariance in the original (unpivoted) coefficient order.
  xtx_inverse <- chol2inv(chol(crossprod(q_observed)))

  elapsed <- system.time({
    datasets <- lapply(seq_len(m), function(i) {
      set.seed(as.integer(seed + i))
      sigma2 <- rss / stats::rchisq(1L, df = residual_df)
      beta <- beta_hat + as.numeric(
        t(chol(sigma2 * xtx_inverse)) %*% stats::rnorm(length(beta_hat))
      )
      completed <- data
      completed$x[missing] <- as.numeric(q_missing %*% beta) +
        stats::rnorm(sum(missing), sd = sqrt(sigma2))
      completed
    })
  })[["elapsed"]]
  .validate_completed_datasets(datasets, dgp, m, "bayes_norm_lm")
  list(
    datasets = datasets,
    diagnostics = list(
      engine = "bayes_norm_lm", warnings = character(0),
      elapsed_seconds = unname(elapsed), residual_df = residual_df,
      imputation_formula = "x ~ y + z + I(z^2)"
    )
  )
}

.draw_smcfcs_imputations <- function(dgp, m, seed, numit = 20L,
                                     rjlimit = 10000L) {
  if (dgp$dgp != "glm") {
    stop("SMCFCS is used only for the logit glm validation cell.",
         call. = FALSE)
  }
  if (!requireNamespace("smcfcs", quietly = TRUE)) {
    stop("The standard single-level comparator requires package 'smcfcs'.",
         call. = FALSE)
  }
  data <- dgp$observed[, c("x", "y", "z"), drop = FALSE]
  data$z_sq <- data$z^2
  data$y <- .numeric_outcome(data$y)
  method <- c(x = "norm", y = "", z = "", z_sq = "")
  predictor_matrix <- matrix(
    0, nrow = ncol(data), ncol = ncol(data),
    dimnames = list(names(data), names(data))
  )
  predictor_matrix["x", c("z", "z_sq")] <- 1

  set.seed(as.integer(seed))
  elapsed <- system.time({
    captured <- .capture_warnings(smcfcs::smcfcs(
      originaldata = data,
      smtype = "logistic",
      smformula = "y ~ x + z",
      method = method,
      predictorMatrix = predictor_matrix,
      m = as.integer(m),
      numit = as.integer(numit),
      rjlimit = as.integer(rjlimit),
      noisy = FALSE
    ))
  })[["elapsed"]]
  fit <- captured$value
  datasets <- lapply(fit$impDatasets, function(completed) {
    completed$z_sq <- NULL
    completed$y <- dgp$observed$y
    completed
  })
  .validate_completed_datasets(datasets, dgp, m, "smcfcs")
  list(
    datasets = datasets,
    diagnostics = list(
      engine = "smcfcs", warnings = captured$warnings,
      elapsed_seconds = unname(elapsed), numit = as.integer(numit),
      rjlimit = as.integer(rjlimit),
      sm_coef_iter = fit$smCoefIter
    )
  )
}

.draw_jomo_lmer_imputations <- function(dgp, m, seed, nburn = 1000L,
                                         nbetween = 100L) {
  if (!requireNamespace("jomo", quietly = TRUE)) {
    stop("The standard multilevel comparator requires package 'jomo'.",
         call. = FALSE)
  }
  data <- dgp$observed[, c("y", "x", "z", "species"), drop = FALSE]
  data$z_sq <- data$z^2
  data$species <- factor(data$species)
  level <- rep(1L, ncol(data))

  set.seed(as.integer(seed))
  elapsed <- system.time({
    captured <- .capture_warnings(jomo::jomo.smc(
      formula = y ~ x + z + (1 | species),
      data = data,
      level = level,
      nburn = as.integer(nburn),
      nbetween = as.integer(nbetween),
      nimp = as.integer(m),
      output = 0,
      model = "lmer"
    ))
  })[["elapsed"]]
  long <- captured$value
  datasets <- lapply(seq_len(m), function(i) {
    completed <- long[long$Imputation == i, , drop = FALSE]
    completed <- completed[order(completed$id), c("y", "x", "z", "z_sq"),
                           drop = FALSE]
    completed$z_sq <- NULL
    completed$species <- as.character(dgp$observed$species)
    rownames(completed) <- NULL
    completed
  })
  .validate_completed_datasets(datasets, dgp, m, "jomo.smc")
  list(
    datasets = datasets,
    diagnostics = list(
      engine = "jomo", warnings = captured$warnings,
      elapsed_seconds = unname(elapsed), nburn = as.integer(nburn),
      nbetween = as.integer(nbetween)
    )
  )
}

.draw_standard_smc_imputations <- function(dgp, m, seed,
                                             smcfcs_numit = 20L,
                                             smcfcs_rjlimit = 10000L,
                                             jomo_nburn = 1000L,
                                             jomo_nbetween = 100L) {
  if (!exists("multi_impute_analysis", mode = "function", inherits = TRUE)) {
    stop("Load the current pigauto source before running package controls.",
         call. = FALSE)
  }
  data <- dgp$observed
  data$z_sq <- data$z^2
  formula <- switch(
    dgp$dgp,
    lm = y ~ x + z,
    glm = y ~ x + z,
    lmer = y ~ x + z + (1 | species)
  )
  control <- switch(
    dgp$dgp,
    lm = list(),
    glm = list(numit = as.integer(smcfcs_numit),
               rjlimit = as.integer(smcfcs_rjlimit)),
    lmer = list(nburn = as.integer(jomo_nburn),
                nbetween = as.integer(jomo_nbetween))
  )
  elapsed <- system.time({
    mi <- multi_impute_analysis(
      data = data, formula = formula, missing = "x", model = dgp$dgp,
      m = as.integer(m), auxiliary = "z_sq", seed = as.integer(seed),
      control = control
    )
  })[["elapsed"]]
  datasets <- lapply(mi$datasets, function(completed) {
    completed$z_sq <- NULL
    completed
  })
  .validate_completed_datasets(datasets, dgp, m, mi$engine)
  list(
    datasets = datasets,
    diagnostics = c(
      list(
        engine = mi$engine,
        warnings = mi$diagnostics$warnings,
        elapsed_seconds = unname(elapsed)
      ),
      control
    )
  )
}
