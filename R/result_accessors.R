#' Extract completed trait data
#'
#' @param result A `pigauto_result`.
#' @return Exactly the stored `result$completed` data.frame.
#' @examples
#' \dontrun{completed_data(result)}
#' @export
completed_data <- function(result) {
  if (!inherits(result, "pigauto_result")) stop("`result` must be a pigauto_result.", call. = FALSE)
  result$completed
}

.describe_pigauto_result <- function(result) {
  data.frame(
    role = c("completed", "prediction", "uncertainty", "conformal", "inference"),
    text = c("Completed data preserve observed cells and fill only modeled missing cells.",
             "Predictions are all-cell diagnostic outputs, not a second completed dataset.",
             "Uncertainty is type-dependent; discrete values are not Gaussian standard errors.",
             "Conformal bounds are nominal diagnostic bounds and do not make a broad certification claim.",
             "A pigauto_result does not authorize inference; use multi_impute_analysis() only in its documented narrow route."),
    stringsAsFactors = FALSE
  )
}

#' Summarise a pigauto result
#'
#' @param object A `pigauto_result`.
#' @param ... Unused.
#' @return A small `summary_pigauto_result` object containing input-check,
#'   completion, output-role, and evaluation summaries.
#' @export
summary.pigauto_result <- function(object, ...) {
  mask <- object$imputed_mask
  trait_cols <- if (!is.null(mask)) colnames(mask) else character()
  target_rows <- rep(TRUE, nrow(object$completed))
  if (inherits(object$check, "pigauto_check")) {
    data_only <- object$check$species$data_only$names
    ids <- if (is.null(object$check$species$species_col)) rownames(object$completed) else as.character(object$completed[[object$check$species$species_col]])
    target_rows <- !(ids %in% data_only)
  }
  completed_traits <- if (length(trait_cols)) object$completed[, trait_cols, drop = FALSE] else NULL
  successful_fills <- if (length(trait_cols)) mask & !is.na(as.matrix(completed_traits)) else NULL
  per_trait <- if (length(trait_cols)) stats::setNames(
    as.integer(colSums(successful_fills[target_rows, , drop = FALSE])), trait_cols
  ) else integer()
  unresolved <- if (length(trait_cols)) sum(is.na(completed_traits[target_rows, , drop = FALSE])) else NA_integer_
  target_cells <- if (length(trait_cols)) sum(target_rows) * length(trait_cols) else NA_integer_
  filled <- if (length(per_trait)) sum(successful_fills[target_rows, , drop = FALSE]) else 0L
  observed <- if (is.na(target_cells)) NA_integer_ else target_cells - filled - unresolved
  structure(list(check = object$check %||% "Preflight check unavailable",
                 counts = list(filled = filled, observed = observed, unresolved = unresolved),
                 per_trait = data.frame(trait = names(per_trait), filled = unname(per_trait), stringsAsFactors = FALSE),
                 output_roles = .describe_pigauto_result(object), evaluation = object$evaluation),
            class = "summary_pigauto_result")
}

#' Print a pigauto result summary
#'
#' @param x A `summary_pigauto_result` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.summary_pigauto_result <- function(x, ...) {
  cat("pigauto result summary\n")
  cat("  observed cells:", x$counts$observed, "\n")
  cat("  filled cells:", x$counts$filled, "\n")
  cat("  unresolved cells:", x$counts$unresolved, "\n")
  if (is.data.frame(x$output_roles)) cat("  roles: completed=filled data; prediction=diagnostic; inference=not authorized\n")
  invisible(x)
}
