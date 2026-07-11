#' Analysis-aware multiple imputation for narrow regression models
#'
#' Generate multiple imputations conditional on a specified downstream
#' regression model. Validity requires MAR, a correctly specified conditional
#' imputation model, and adequate SMCFCS or jomo convergence. This experimental
#' interface is deliberately narrow: exactly one continuous covariate may
#' contain missing values, while the outcome, all other predictors, auxiliary
#' variables, and any grouping variable must be fully observed.
#'
#' @param data A data frame containing the analysis variables.
#' @param formula A two-sided additive model formula. For `model = "lmer"`,
#'   exactly one random-intercept term of the form `(1 | group)` is required.
#' @param missing A single character string naming the continuous numeric
#'   covariate to impute. It must occur as a main-effect predictor in
#'   `formula` and must be the only column in `data` containing missing values.
#' @param model One of `"lm"`, `"glm"`, or `"lmer"`. The `"glm"` route is
#'   restricted to binomial-logit regression. The `"lmer"` route is restricted
#'   to a Gaussian random-intercept model.
#' @param m Integer number of completed datasets. Must be at least two.
#' @param auxiliary Character vector naming fully observed numeric columns used
#'   only by the imputation model. They cannot duplicate formula variables.
#'   Derived terms must be created explicitly in `data` before calling this
#'   function.
#' @param seed Integer random seed.
#' @param control Named list of engine controls. The `"lm"` route accepts no
#'   controls. The `"glm"` route accepts `numit` (default 20) and `rjlimit`
#'   (default 100000). The `"lmer"` route accepts `nburn` (default 1000) and
#'   `nbetween` (default 100).
#'
#' @return An object with classes `"pigauto_analysis_mi"` and `"pigauto_mi"`.
#'   Its `datasets` component is a list of `m` completed data frames compatible
#'   with [with_imputations()] and [pool_mi()]. The object also records the
#'   analysis model, formula, engine, controls, engine version, seed, runtime
#'   provenance, and imputed-variable metadata. Any engine warning aborts the
#'   call instead of returning an inference-ready object.
#'
#' @section Experimental scope:
#' This function does not fit pigauto's phylogenetic baseline or graph neural
#' network and does not propagate posterior-tree uncertainty. The `"lm"` route
#' uses proper Bayesian Normal-regression draws. The `"glm"` route requires the
#' optional `smcfcs` package, and the `"lmer"` route requires the optional
#' `jomo` and `lme4` packages. Interactions, transformed terms, random slopes,
#' nested or crossed random effects, multiple incomplete variables, missing
#' outcomes, and pooling of variance components or BLUPs are unsupported.
#'
#' @examples
#' set.seed(11)
#' dat <- data.frame(
#'   y = stats::rnorm(80),
#'   x = stats::rnorm(80),
#'   z = stats::rnorm(80)
#' )
#' dat$x[seq(4, 80, by = 5)] <- NA_real_
#' mi <- multi_impute_analysis(
#'   dat, y ~ x + z, missing = "x", model = "lm", m = 5, seed = 12
#' )
#' fits <- with_imputations(
#'   mi, function(d) stats::lm(y ~ x + z, data = d), .progress = FALSE
#' )
#' pool_mi(fits)
#'
#' @seealso [multi_impute()], [with_imputations()], [pool_mi()]
#' @export
multi_impute_analysis <- function(data, formula, missing,
                                  model = c("lm", "glm", "lmer"),
                                  m = 50L, auxiliary = character(),
                                  seed = 1L, control = list()) {
  model <- match.arg(model)
  validated <- .analysis_mi_validate(
    data = data, formula = formula, missing = missing, model = model,
    m = m, auxiliary = auxiliary, seed = seed, control = control
  )

  captured <- switch(
    model,
    lm = .analysis_mi_capture_warnings(
      .analysis_mi_bayes_norm(
        validated$data, validated$fixed_formula, validated$missing,
        validated$auxiliary, validated$m, validated$seed
      )
    ),
    glm = {
      .analysis_mi_require("smcfcs", "model = \"glm\"")
      .analysis_mi_capture_warnings(
        .analysis_mi_smcfcs(
          validated$data, validated$fixed_formula, validated$missing,
          validated$auxiliary, validated$m, validated$seed,
          validated$control
        )
      )
    },
    lmer = {
      .analysis_mi_require("lme4", "model = \"lmer\"")
      .analysis_mi_require("jomo", "model = \"lmer\"")
      .analysis_mi_capture_warnings(
        .analysis_mi_jomo(
          validated$data, validated$formula, validated$missing,
          validated$auxiliary, validated$group, validated$m,
          validated$seed, validated$control
        )
      )
    }
  )
  if (model %in% c("glm", "lmer")) {
    .analysis_mi_abort_on_warnings(captured$warnings, captured$messages, model)
  }

  datasets <- .analysis_mi_validate_datasets(
    captured$value, validated$data, validated$missing, validated$m
  )
  engine <- switch(model, lm = "bayes_norm", glm = "smcfcs", lmer = "jomo")
  dependency_versions <- switch(
    model,
    lm = c(R = as.character(getRversion())),
    glm = c(smcfcs = as.character(utils::packageVersion("smcfcs"))),
    lmer = c(
      jomo = as.character(utils::packageVersion("jomo")),
      lme4 = as.character(utils::packageVersion("lme4"))
    )
  )
  engine_version <- if (engine == "bayes_norm") {
    unname(dependency_versions[["R"]])
  } else {
    unname(dependency_versions[[engine]])
  }

  structure(
    list(
      datasets = datasets,
      m = validated$m,
      draws_method = "analysis_aware",
      model = model,
      engine = engine,
      formula = formula,
      missing = validated$missing,
      auxiliary = validated$auxiliary,
      control = validated$control,
      seed = validated$seed,
      engine_version = engine_version,
      dependency_versions = dependency_versions,
      package_version = as.character(utils::packageVersion("pigauto")),
      r_version = R.version.string,
      platform = R.version$platform,
      rng_kind = RNGkind(),
      diagnostics = list(
        warnings = captured$warnings,
        messages = captured$messages,
        output = captured$output
      ),
      imputed_mask = is.na(validated$data[[validated$missing]]),
      data_original = validated$data,
      call = match.call()
    ),
    class = c("pigauto_analysis_mi", "pigauto_mi", "list")
  )
}

.analysis_mi_validate <- function(data, formula, missing, model, m, auxiliary,
                                  seed, control) {
  if (!is.data.frame(data) || nrow(data) < 2L) {
    stop("`data` must be a data frame with at least two rows.", call. = FALSE)
  }
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be a two-sided model formula.", call. = FALSE)
  }
  if (!is.character(missing) || length(missing) != 1L || is.na(missing) ||
      !nzchar(missing) || !(missing %in% names(data))) {
    stop("`missing` must name exactly one column in `data`.", call. = FALSE)
  }
  if (!is.double(data[[missing]]) || is.factor(data[[missing]])) {
    stop("`missing` must name a double-precision continuous covariate.",
         call. = FALSE)
  }
  missing_rows <- is.na(data[[missing]])
  if (!any(missing_rows)) {
    stop("`missing` must contain at least one missing value.", call. = FALSE)
  }
  if (any(!is.finite(data[[missing]][!missing_rows]))) {
    stop("Observed values of `missing` must be finite.", call. = FALSE)
  }
  other_na <- vapply(data[setdiff(names(data), missing)], anyNA, logical(1))
  if (any(other_na)) {
    stop(
      "Only `", missing, "` may contain missing values; incomplete columns: ",
      paste(names(other_na)[other_na], collapse = ", "), ".", call. = FALSE
    )
  }

  auxiliary <- as.character(auxiliary)
  if (anyNA(auxiliary) || any(!nzchar(auxiliary)) || anyDuplicated(auxiliary)) {
    stop("`auxiliary` must contain unique, non-missing column names.",
         call. = FALSE)
  }
  absent_aux <- setdiff(auxiliary, names(data))
  if (length(absent_aux)) {
    stop("Unknown `auxiliary` columns: ", paste(absent_aux, collapse = ", "),
         ".", call. = FALSE)
  }
  if (missing %in% auxiliary) {
    stop("`missing` cannot also be an auxiliary column.", call. = FALSE)
  }

  m <- .analysis_mi_integer(m, "m", lower = 2L)
  seed <- .analysis_mi_integer(seed, "seed", lower = 0L)
  control <- .analysis_mi_control(control, model)

  if (model == "lmer") {
    .analysis_mi_require("lme4", "validation of model = \"lmer\"")
    bars <- suppressWarnings(lme4::findbars(formula))
    if (length(bars) != 1L || deparse(bars[[1]][[2]]) != "1" ||
        !is.symbol(bars[[1]][[3]])) {
      stop("`model = \"lmer\"` requires exactly one `(1 | group)` term.",
           call. = FALSE)
    }
    group <- as.character(bars[[1]][[3]])
    fixed_formula <- suppressWarnings(lme4::nobars(formula))
  } else {
    if (length(.analysis_mi_findbars_safe(formula)) > 0L) {
      stop("Random-effect terms are supported only for `model = \"lmer\"`.",
           call. = FALSE)
    }
    group <- NULL
    fixed_formula <- formula
  }

  fixed <- .analysis_mi_validate_fixed_formula(fixed_formula, data, missing)
  response <- fixed$response
  predictors <- fixed$predictors
  overlap <- intersect(auxiliary, c(response, predictors, group))
  if (length(overlap)) {
    stop("`auxiliary` columns must not duplicate formula variables: ",
         paste(overlap, collapse = ", "), ".", call. = FALSE)
  }
  involved <- unique(c(response, predictors, auxiliary, group))
  absent <- setdiff(involved, names(data))
  if (length(absent)) {
    stop("Formula variables not found in `data`: ", paste(absent, collapse = ", "),
         ".", call. = FALSE)
  }
  if (any(vapply(data[involved], anyNA, logical(1)) & involved != missing)) {
    stop("Outcome, other predictors, auxiliaries, and group must be fully observed.",
         call. = FALSE)
  }
  nonfinite <- vapply(data[involved], function(x) {
    is.numeric(x) && any(!is.finite(x[!is.na(x)]))
  }, logical(1))
  if (any(nonfinite)) {
    stop("Model and auxiliary numeric columns must contain only finite values.",
         call. = FALSE)
  }
  fixed_inputs <- unique(c(setdiff(predictors, missing), auxiliary))
  unsupported <- fixed_inputs[!vapply(data[fixed_inputs], is.numeric,
                                       logical(1))]
  if (length(unsupported)) {
    stop("Fixed-effect predictors and auxiliaries must be numeric; unsupported: ",
         paste(unsupported, collapse = ", "), ".", call. = FALSE)
  }

  if (model %in% c("lm", "lmer")) {
    if (!is.numeric(data[[response]]) || any(!is.finite(data[[response]]))) {
      stop("Gaussian outcomes must be fully observed and finite numeric values.",
           call. = FALSE)
    }
  } else {
    .analysis_mi_binary_outcome(data[[response]])
  }
  if (!is.null(group)) {
    if (anyNA(data[[group]]) || length(unique(data[[group]])) < 2L) {
      stop("The random-intercept grouping variable must have at least two groups.",
           call. = FALSE)
    }
  }

  list(
    data = data, formula = formula, fixed_formula = fixed_formula,
    response = response, predictors = predictors, missing = missing,
    auxiliary = auxiliary, group = group, m = m, seed = seed, control = control
  )
}

.analysis_mi_findbars_safe <- function(formula) {
  if (!grepl("\\|", paste(deparse(formula), collapse = " "), fixed = FALSE)) {
    return(list())
  }
  if (!requireNamespace("lme4", quietly = TRUE)) return(list(quote(`|`(x, y))))
  suppressWarnings(lme4::findbars(formula))
}

.analysis_mi_validate_fixed_formula <- function(formula, data, missing) {
  if (!is.symbol(formula[[2]])) {
    stop("The formula outcome must be an untransformed column name.",
         call. = FALSE)
  }
  tt <- stats::terms(formula, data = data)
  labels <- attr(tt, "term.labels")
  orders <- attr(tt, "order")
  if (length(labels) == 0L || attr(tt, "intercept") != 1L ||
      any(orders != 1L) || any(labels == ".") ||
      any(!(labels %in% names(data)))) {
    stop("Only additive, untransformed main-effect predictors with an intercept are supported.",
         call. = FALSE)
  }
  response <- as.character(formula[[2]])
  if (!(missing %in% labels)) {
    stop("`missing` must appear as a main-effect predictor in `formula`.",
         call. = FALSE)
  }
  list(response = response, predictors = labels)
}

.analysis_mi_control <- function(control, model) {
  if (!is.list(control) || (length(control) && is.null(names(control))) ||
      anyDuplicated(names(control))) {
    stop("`control` must be a uniquely named list.", call. = FALSE)
  }
  defaults <- switch(
    model,
    lm = list(),
    glm = list(numit = 20L, rjlimit = 100000L),
    lmer = list(nburn = 1000L, nbetween = 100L)
  )
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unsupported `control` entries for model ", model, ": ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  out <- utils::modifyList(defaults, control)
  if (length(out)) {
    out <- Map(function(value, name) {
      .analysis_mi_integer(value, paste0("control$", name), lower = 1L)
    }, out, names(out))
  }
  out
}

.analysis_mi_integer <- function(x, name, lower) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < lower || x > .Machine$integer.max) {
    stop("`", name, "` must be a whole number >= ", lower, ".", call. = FALSE)
  }
  as.integer(x)
}

.analysis_mi_binary_outcome <- function(x) {
  if (is.factor(x)) {
    x <- droplevels(x)
    if (nlevels(x) != 2L) {
      stop("The `glm` outcome must have exactly two observed levels.",
           call. = FALSE)
    }
    return(as.integer(x) - 1L)
  }
  if (!(is.numeric(x) || is.logical(x))) {
    stop("The `glm` outcome must be coded 0/1 or be a two-level factor.",
         call. = FALSE)
  }
  values <- sort(unique(as.numeric(x)))
  if (length(values) != 2L || !identical(values, c(0, 1))) {
    stop("The `glm` outcome must be coded 0/1 or be a two-level factor.",
         call. = FALSE)
  }
  as.integer(x)
}

.analysis_mi_bayes_norm <- function(data, formula, missing, auxiliary, m, seed) {
  response <- as.character(formula[[2]])
  predictors <- attr(stats::terms(formula, data = data), "term.labels")
  imputation_predictors <- unique(c(response, setdiff(predictors, missing), auxiliary))
  imputation_formula <- stats::reformulate(imputation_predictors, response = missing)
  mf <- stats::model.frame(imputation_formula, data = data,
                           na.action = stats::na.pass, drop.unused.levels = FALSE)
  tt <- stats::delete.response(stats::terms(mf))
  X <- stats::model.matrix(tt, mf)
  x <- data[[missing]]
  observed <- !is.na(x)
  X_observed <- X[observed, , drop = FALSE]
  qr_x <- qr(X_observed)
  if (qr_x$rank != ncol(X_observed)) {
    stop("The Bayesian Normal imputation design matrix is rank deficient.",
         call. = FALSE)
  }
  condition_number <- kappa(X_observed, exact = TRUE)
  if (!is.finite(condition_number) || condition_number > 1e8) {
    stop("The Bayesian Normal imputation design matrix is ill-conditioned.",
         call. = FALSE)
  }
  fit <- stats::lm.fit(X_observed, x[observed])
  df <- sum(observed) - qr_x$rank
  sse <- sum(fit$residuals^2)
  if (df < 1L || !is.finite(sse) || sse <= 0) {
    stop("Bayesian Normal imputation needs positive residual degrees of freedom and variance.",
         call. = FALSE)
  }
  R <- chol(crossprod(X_observed))
  X_missing <- X[!observed, , drop = FALSE]

  set.seed(seed)
  lapply(seq_len(m), function(i) {
    sigma2 <- sse / stats::rchisq(1L, df = df)
    beta <- fit$coefficients + sqrt(sigma2) *
      backsolve(R, stats::rnorm(ncol(X_observed)))
    draw <- as.numeric(X_missing %*% beta) +
      stats::rnorm(sum(!observed), sd = sqrt(sigma2))
    completed <- data
    completed[[missing]][!observed] <- draw
    completed
  })
}

.analysis_mi_smcfcs <- function(data, formula, missing, auxiliary, m, seed,
                                control) {
  response <- as.character(formula[[2]])
  predictors <- attr(stats::terms(formula, data = data), "term.labels")
  columns <- unique(c(response, predictors, auxiliary))
  engine_data <- data[, columns, drop = FALSE]
  engine_data[[response]] <- .analysis_mi_binary_outcome(engine_data[[response]])
  method <- stats::setNames(rep("", ncol(engine_data)), names(engine_data))
  method[[missing]] <- "norm"
  predictor_matrix <- matrix(
    0, nrow = ncol(engine_data), ncol = ncol(engine_data),
    dimnames = list(names(engine_data), names(engine_data))
  )
  imputation_predictors <- unique(c(setdiff(predictors, missing), auxiliary))
  predictor_matrix[missing, imputation_predictors] <- 1

  set.seed(seed)
  fit <- smcfcs::smcfcs(
    originaldata = engine_data,
    smtype = "logistic",
    smformula = paste(deparse(formula), collapse = " "),
    method = method,
    predictorMatrix = predictor_matrix,
    m = m,
    numit = control$numit,
    rjlimit = control$rjlimit,
    noisy = FALSE
  )
  if (!is.list(fit$impDatasets) || length(fit$impDatasets) != m) {
    stop("smcfcs did not return the requested number of imputations.",
         call. = FALSE)
  }
  lapply(fit$impDatasets, function(engine_completed) {
    if (!is.data.frame(engine_completed) || nrow(engine_completed) != nrow(data)) {
      stop("smcfcs returned an imputation with incompatible rows.",
           call. = FALSE)
    }
    completed <- data
    rows <- is.na(data[[missing]])
    completed[[missing]][rows] <- as.numeric(engine_completed[[missing]][rows])
    completed
  })
}

.analysis_mi_jomo <- function(data, formula, missing, auxiliary, group, m,
                              seed, control) {
  columns <- unique(c(all.vars(formula), auxiliary))
  engine_data <- data[, columns, drop = FALSE]
  engine_data[[group]] <- factor(engine_data[[group]])
  level <- rep(1L, ncol(engine_data))

  set.seed(seed)
  long <- jomo::jomo.smc(
    formula = formula, data = engine_data, level = level,
    nburn = control$nburn, nbetween = control$nbetween, nimp = m,
    output = 0, model = "lmer"
  )
  required <- c("Imputation", "id", missing)
  if (!is.data.frame(long) || !all(required %in% names(long))) {
    stop("jomo returned an incompatible long-format result.", call. = FALSE)
  }
  lapply(seq_len(m), function(i) {
    one <- long[long$Imputation == i, , drop = FALSE]
    if (nrow(one) != nrow(data) || anyDuplicated(one$id) ||
        !setequal(one$id, seq_len(nrow(data)))) {
      stop("jomo returned an imputation with incompatible row identifiers.",
           call. = FALSE)
    }
    one <- one[order(one$id), , drop = FALSE]
    completed <- data
    rows <- is.na(data[[missing]])
    completed[[missing]][rows] <- as.numeric(one[[missing]][rows])
    completed
  })
}

.analysis_mi_capture_warnings <- function(expr) {
  warnings <- character()
  messages <- character()
  value <- NULL
  output <- utils::capture.output(
    value <- withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    type = "output"
  )
  list(
    value = value,
    warnings = warnings,
    messages = messages,
    output = output[nzchar(output)]
  )
}

.analysis_mi_abort_on_warnings <- function(warnings, messages, model) {
  if (!length(warnings)) return(invisible(TRUE))
  detail <- paste(unique(warnings), collapse = " | ")
  stop(
    "The analysis-aware ", model, " engine emitted a warning, so no ",
    "inference-ready object was returned: ", detail,
    if (length(messages)) paste0(" Engine messages: ",
                                 paste(unique(messages), collapse = " | "))
    else "",
    call. = FALSE
  )
}

.analysis_mi_validate_datasets <- function(datasets, original, missing, m) {
  if (!is.list(datasets) || length(datasets) != m) {
    stop("The imputation engine did not return exactly `m` datasets.",
         call. = FALSE)
  }
  missing_rows <- is.na(original[[missing]])
  for (i in seq_along(datasets)) {
    completed <- datasets[[i]]
    if (!is.data.frame(completed) || !identical(dim(completed), dim(original)) ||
        !identical(names(completed), names(original)) ||
        !identical(rownames(completed), rownames(original))) {
      stop("Imputed dataset ", i, " does not preserve the input shape and names.",
           call. = FALSE)
    }
    classes <- vapply(completed, function(x) paste(class(x), collapse = "/"),
                      character(1))
    original_classes <- vapply(original, function(x) paste(class(x), collapse = "/"),
                               character(1))
    if (!identical(classes, original_classes)) {
      stop("Imputed dataset ", i, " does not preserve input column classes.",
           call. = FALSE)
    }
    if (any(!is.finite(completed[[missing]]))) {
      stop("Imputed dataset ", i, " contains non-finite values in `missing`.",
           call. = FALSE)
    }
    if (!identical(completed[[missing]][!missing_rows],
                   original[[missing]][!missing_rows])) {
      stop("Imputed dataset ", i, " modified observed values.", call. = FALSE)
    }
    unchanged <- setdiff(names(original), missing)
    if (!identical(completed[unchanged], original[unchanged])) {
      stop("Imputed dataset ", i, " modified non-imputed columns.",
           call. = FALSE)
    }
  }
  datasets
}

.analysis_mi_require <- function(package, route) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "The analysis-aware ", route, " route requires optional package '",
      package, "'. Install it with install.packages(\"", package, "\").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @export
print.pigauto_analysis_mi <- function(x, ...) {
  warning_count <- length(x$diagnostics$warnings %||% character())
  cat("pigauto experimental analysis-aware multiple imputation\n")
  cat(sprintf("  M        : %d imputations\n", x$m))
  cat(sprintf("  Model    : %s\n", x$model))
  cat(sprintf("  Engine   : %s %s\n", x$engine, x$engine_version))
  cat(sprintf("  Formula  : %s\n", paste(deparse(x$formula), collapse = " ")))
  cat(sprintf("  Imputed  : %s (%d cells)\n", x$missing,
              sum(x$imputed_mask)))
  cat(sprintf("  Warnings : %d captured\n", warning_count))
  cat("\n  Fit downstream models:   with_imputations(mi, fit_fun)\n")
  cat("  Pool fixed effects with: pool_mi(fits)\n")
  invisible(x)
}
