#' Fit a downstream model on every imputed dataset
#'
#' Apply a user-supplied model-fitting function `.f` to each of the `M`
#' complete datasets stored in an analysis-aware MI object and return the list
#' of fits. This is the supported fitting entry for
#' [multi_impute_analysis()] only. Prediction-diagnostic and posterior-tree
#' sensitivity completions stop with an actionable error before any model is
#' fitted. The analysis-aware workflow is experimental and limited to
#' fixed-effect coefficients and their covariance matrices.
#'
#' @param mi A `pigauto_analysis_mi` object returned by
#'   [multi_impute_analysis()]. Legacy `pigauto_mi` objects, diagnostic
#'   completions, tree-sensitivity completions, and bare dataset lists are not
#'   supported fitting inputs.
#' @param .f A function of the form `function(dataset, ...)` that fits
#'   a model to one complete data.frame and returns a model object.
#'   [pool_mi()] supplies automatic fixed-effect adapters for its documented
#'   model classes. Other classes require `coef()` and `vcov()` methods that
#'   return compatible fixed-effect quantities, or explicit extractor
#'   functions supplied to `pool_mi()`; extractability does not imply that an
#'   unlisted analysis model has passed the analysis-aware validation gate.
#' @param ... Additional arguments passed to `.f` for every imputation.
#' @param .progress Logical. Show a text progress indicator (default
#'   `TRUE` in interactive sessions).
#' @param .on_error One of `"continue"` (default) or `"stop"`. When
#'   `"continue"`, errors from `.f` are captured per imputation and the
#'   loop proceeds; a warning at the end summarises failures. When
#'   `"stop"`, the first error aborts the entire run.
#'
#' @return A list of length `M` with class `"pigauto_mi_fits"`. Each
#'   element is either a model fit or, if `.f` errored on that
#'   imputation and `.on_error = "continue"`, an object of class
#'   `"pigauto_mi_error"` containing the captured condition. [pool_mi()]
#'   filters error elements automatically.
#'
#' @seealso [multi_impute_analysis()], [pool_mi()]
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
#' fits <- with_imputations(mi, function(d) stats::lm(y ~ x + z, data = d))
#' pool_mi(fits)
#' }
#'
#' @export
with_imputations <- function(mi, .f, ...,
                             .progress = interactive(),
                             .on_error = c("continue", "stop")) {

  .on_error <- match.arg(.on_error)

  if (!is.function(.f)) {
    stop("`.f` must be a function of one argument (a data.frame).",
         call. = FALSE)
  }

  workflow <- if (is.list(mi)) mi$mi_workflow else NULL
  if (identical(workflow, "pigauto_diagnostic_mi") ||
      inherits(mi, "pigauto_diagnostic_mi")) {
    stop("`multi_impute()` returns prediction-diagnostic completions, not ",
         "analysis-aware MI. Use `multi_impute_analysis()` before ",
         "`with_imputations()`.", call. = FALSE)
  }
  if (identical(workflow, "pigauto_tree_sensitivity_diagnostic") ||
      inherits(mi, "pigauto_tree_sensitivity_diagnostic") ||
      inherits(mi, "pigauto_mi_trees")) {
    stop("`multi_impute_trees()` returns tree-sensitivity diagnostics, not ",
         "analysis-aware MI. Tree-aware downstream pooling is unsupported; ",
         "use `multi_impute_analysis()` for the supported narrow workflow.",
         call. = FALSE)
  }
  if (inherits(mi, "pigauto_mi") && !inherits(mi, "pigauto_analysis_mi")) {
    stop("Legacy `pigauto_mi` objects have no analysis-aware provenance. ",
         "Create a new object with `multi_impute_analysis()` before ",
         "`with_imputations()`.", call. = FALSE)
  }
  if (!inherits(mi, "pigauto_analysis_mi") ||
      !identical(workflow, "pigauto_analysis_mi_v1")) {
    stop("`with_imputations()` requires an analysis-aware object from ",
         "`multi_impute_analysis()`.", call. = FALSE)
  }
  datasets <- mi$datasets

  M <- length(datasets)
  if (M < 2L) {
    stop("`mi` contains fewer than 2 datasets; nothing to pool later.",
         call. = FALSE)
  }

  dots <- list(...)

  fits <- vector("list", M)
  failures <- integer(0)

  for (i in seq_len(M)) {
    if (.progress) {
      msg <- sprintf("\r  fitting imputation %d/%d", i, M)
      message(msg, appendLF = FALSE)
    }

    dataset_i <- datasets[[i]]
    call_args <- c(list(dataset_i), dots)

    result <- tryCatch(
      do.call(.f, call_args),
      error = function(e) {
        structure(
          list(index = i, message = conditionMessage(e), call = conditionCall(e)),
          class = c("pigauto_mi_error", "condition")
        )
      }
    )

    if (inherits(result, "pigauto_mi_error")) {
      failures <- c(failures, i)
      if (.on_error == "stop") {
        if (.progress) message("")  # finish the progress line
        stop("`.f` failed on imputation ", i, ": ", result$message,
             call. = FALSE)
      }
    }

    fits[[i]] <- result
  }

  if (.progress) message("")  # newline after progress updates

  if (length(failures) > 0L) {
    warning(sprintf(
      "%d of %d fits failed (%s). They will be dropped by pool_mi().",
      length(failures), M,
      paste0("indices: ", paste(failures, collapse = ", "))
    ), call. = FALSE)
  }

  attr(fits, "n_fits")   <- M
  attr(fits, "n_failed") <- length(failures)
  attr(fits, "failed")   <- failures
  attr(fits, "mi_workflow") <- workflow
  class(fits) <- c("pigauto_mi_fits", "list")
  fits
}


#' @export
print.pigauto_mi_fits <- function(x, ...) {
  n <- attr(x, "n_fits") %||% length(x)
  nf <- attr(x, "n_failed") %||% 0L
  ok <- n - nf
  cat(sprintf("pigauto multi-imputation fits: %d/%d successful\n", ok, n))
  if (nf > 0L) {
    failed <- attr(x, "failed")
    cat(sprintf("  failed at indices: %s\n",
                paste(failed, collapse = ", ")))
  }
  if (ok > 0L) {
    # Show the class of the first successful fit.
    first_ok <- NULL
    for (i in seq_along(x)) {
      if (!inherits(x[[i]], "pigauto_mi_error")) {
        first_ok <- x[[i]]
        break
      }
    }
    if (!is.null(first_ok)) {
      cat(sprintf("  model class: %s\n",
                  paste(class(first_ok), collapse = " / ")))
    }
  }
  cat("\n  Pool with: pool_mi(fits)\n")
  invisible(x)
}


# Local null-coalesce so we don't require rlang at load time (rlang is
# already in Imports for .data, but this keeps the file self-contained).
`%||%` <- function(a, b) if (is.null(a)) b else a
