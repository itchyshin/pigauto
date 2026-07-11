#' Fit a downstream model on every imputed dataset
#'
#' Apply a user-supplied model-fitting function `.f` to each of the `M`
#' complete datasets stored in a `pigauto_mi` object and return the list
#' of fits. For inference, use an object returned by
#' [multi_impute_analysis()]; conformal-width, Brownian/MC-dropout, and PMM
#' prediction-diagnostic draws are unsupported. The analysis-aware workflow is
#' experimental and is
#' limited to fixed-effect coefficients and their covariance matrices.
#'
#' @param mi A `pigauto_mi` object returned by [multi_impute_analysis()]. Plain
#'   lists of data.frames are also accepted and treated as the `datasets`
#'   slot directly.
#' @param .f A function of the form `function(dataset, ...)` that fits
#'   a model to one complete data.frame and returns a model object.
#'   [pool_mi()] supplies automatic fixed-effect adapters for its documented
#'   model classes. Other classes require `coef()` and `vcov()` methods that
#'   return compatible fixed-effect quantities, or explicit extractor
#'   functions supplied to `pool_mi()`; extractability does not imply that an
#'   unlisted analysis model has passed the analysis-aware validation gate.
#'   When `mi` comes from [multi_impute_trees()], `.f` may also declare
#'   explicit `tree`, `tree_index`, or `imputation` arguments; these are
#'   filled with the posterior tree object, its index in `mi$trees`, and
#'   the stochastic-completion index. The dataset also carries matching
#'   attributes. This metadata support is for prediction-sensitivity
#'   diagnostics only; tree-aware downstream inference is unsupported.
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
#' \dontrun{
#' mi <- multi_impute_analysis(
#'   data = df, formula = y ~ x + z, missing = "x",
#'   model = "lm", m = 50L
#' )
#'
#' fits <- with_imputations(mi, function(d) stats::lm(y ~ x + z, data = d))
#'
#' pool_mi(fits)
#'
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

  # Accept either a pigauto_mi object or a bare list of data.frames.
  if (inherits(mi, "pigauto_mi")) {
    datasets <- mi$datasets
  } else if (is.list(mi) && all(vapply(mi, is.data.frame, logical(1)))) {
    datasets <- mi
  } else {
    stop("`mi` must be a `pigauto_mi` object (from multi_impute()) or a ",
         "list of data.frames.", call. = FALSE)
  }

  M <- length(datasets)
  if (M < 2L) {
    stop("`mi` contains fewer than 2 datasets; nothing to pool later.",
         call. = FALSE)
  }

  tree_index <- NULL
  trees <- NULL
  n_trees <- NULL
  m_per_tree <- NULL
  if (inherits(mi, "pigauto_mi_trees")) {
    tree_index <- mi$tree_index
    trees <- mi$trees
    n_trees <- mi$n_trees
    m_per_tree <- mi$m_per_tree
    if (is.null(tree_index) || length(tree_index) != M) {
      stop("`mi$tree_index` must have one entry per imputed dataset.",
           call. = FALSE)
    }
    if (anyNA(tree_index) || any(tree_index < 1L)) {
      stop("`mi$tree_index` must contain positive tree indices.",
           call. = FALSE)
    }
    if (!is.null(trees) && any(tree_index > length(trees))) {
      stop("`mi$trees` does not contain every tree referenced by `mi$tree_index`.",
           call. = FALSE)
    }
  }

  dots <- list(...)
  f_formals <- names(formals(.f))
  add_if_declared <- function(args, nm, value) {
    if (nm %in% f_formals && !(nm %in% names(args))) args[[nm]] <- value
    args
  }

  fits <- vector("list", M)
  failures <- integer(0)

  for (i in seq_len(M)) {
    if (.progress) {
      msg <- sprintf("\r  fitting imputation %d/%d", i, M)
      message(msg, appendLF = FALSE)
    }

    dataset_i <- datasets[[i]]
    call_args <- c(list(dataset_i), dots)
    if (!is.null(tree_index)) {
      t_idx <- tree_index[[i]]
      tree_i <- if (!is.null(trees)) trees[[t_idx]] else NULL
      attr(dataset_i, "tree_index") <- t_idx
      attr(dataset_i, "imputation") <- i
      if (!is.null(tree_i)) attr(dataset_i, "tree") <- tree_i
      call_args[[1L]] <- dataset_i
      call_args <- add_if_declared(call_args, "tree", tree_i)
      call_args <- add_if_declared(call_args, "tree_index", t_idx)
      call_args <- add_if_declared(call_args, "imputation", i)
    }

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
  if (!is.null(tree_index)) {
    attr(fits, "tree_index") <- tree_index
    attr(fits, "n_trees") <- n_trees
    attr(fits, "m_per_tree") <- m_per_tree
  }
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
