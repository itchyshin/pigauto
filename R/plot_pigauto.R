#' Plot training history (ggplot2, deprecated)
#'
#' @description
#' Deprecated. Use the base-R S3 method instead:
#' \code{plot(fit, type = "history")}.
#'
#' Produces a ggplot2 validation loss training curve.
#'
#' @param x object of class \code{"pigauto_fit"}.
#' @param ... ignored.
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes geom_point geom_abline facet_wrap
#' @importFrom ggplot2 theme_minimal labs geom_line geom_ribbon
#' @importFrom rlang .data
#' @export
plot_history_gg <- function(x, ...) {
  if (!inherits(x, "pigauto_fit")) {
    stop("'x' must be a pigauto_fit object.")
  }
  if (nrow(x$history) == 0L) stop("No training history available.")

  val_col <- if ("val_loss" %in% names(x$history)) "val_loss" else "val_rmse"
  df <- x$history[is.finite(x$history[[val_col]]), ]

  ggplot2::ggplot(df, ggplot2::aes(x = .data$epoch, y = .data[[val_col]])) +
    ggplot2::geom_line(colour = "#1f78b4") +
    ggplot2::geom_point(size = 1.5, colour = "#1f78b4") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Training history",
      x = "Epoch", y = "Validation loss"
    )
}


#' Plot uncertainty ribbons for imputed trait values
#'
#' @description
#' Plots 95\% credible intervals around imputed values, sorted by predicted
#' value.  If \code{truth} is supplied, observed values are overlaid as
#' points.  Works for continuous, count, and ordinal traits.  For binary
#' traits, plots predicted probabilities with CI ribbon.
#'
#' @param pred_result list with \code{imputed} and \code{se} components
#'   (output of \code{predict.pigauto_fit}).
#' @param truth matrix or data.frame of true values (same scale as
#'   \code{pred_result$imputed}), or \code{NULL}.
#' @param trait_name character. Which trait to plot (must match column name).
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point
#' @importFrom ggplot2 theme_minimal labs geom_hline
#' @importFrom rlang .data
#' @export
plot_uncertainty <- function(pred_result, truth = NULL, trait_name) {
  # Determine trait type from trait_map if available
  tm_info <- NULL
  if (!is.null(pred_result$trait_map)) {
    for (tm in pred_result$trait_map) {
      if (tm$name == trait_name) { tm_info <- tm; break }
    }
  }

  # Get predicted values (numeric)
  if (is.data.frame(pred_result$imputed)) {
    pred_j <- pred_result$imputed[[trait_name]]
    if (is.null(pred_j)) {
      stop("'", trait_name, "' not found in pred_result$imputed.")
    }
    # For binary: use probabilities if available
    if (!is.null(tm_info) && tm_info$type == "binary" &&
        !is.null(pred_result$probabilities[[trait_name]])) {
      pred_j <- pred_result$probabilities[[trait_name]]
    } else {
      pred_j <- as.numeric(pred_j)
    }
  } else {
    if (!trait_name %in% colnames(pred_result$imputed)) {
      stop("'", trait_name, "' not found in pred_result$imputed.")
    }
    pred_j <- pred_result$imputed[, trait_name]
  }

  # SE
  se_j <- pred_result$se[, trait_name]
  if (is.null(se_j)) stop("No SE available for trait '", trait_name, "'.")

  idx <- order(pred_j)

  df <- data.frame(
    i     = seq_along(pred_j)[idx],
    pred  = pred_j[idx],
    lower = (pred_j - 1.96 * se_j)[idx],
    upper = (pred_j + 1.96 * se_j)[idx]
  )

  # Determine y-axis label
  y_label <- trait_name
  if (!is.null(tm_info)) {
    if (tm_info$type == "binary") y_label <- paste(trait_name, "(probability)")
    if (tm_info$type == "categorical") y_label <- paste(trait_name, "(entropy)")
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$i)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      fill = "#1f78b4", alpha = 0.25
    ) +
    ggplot2::geom_line(ggplot2::aes(y = .data$pred),
                       colour = "#1f78b4", linewidth = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Imputation uncertainty:", trait_name),
      x = "Species (sorted by prediction)", y = y_label
    )

  # Add reference line for binary
  if (!is.null(tm_info) && tm_info$type == "binary") {
    p <- p + ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
                                  colour = "grey40")
  }

  if (!is.null(truth)) {
    if (is.data.frame(truth)) {
      truth_j <- as.numeric(truth[[trait_name]])
    } else if (is.matrix(truth) && trait_name %in% colnames(truth)) {
      truth_j <- truth[, trait_name]
    } else {
      truth_j <- NULL
    }
    if (!is.null(truth_j)) {
      df$truth <- truth_j[idx]
      p <- p + ggplot2::geom_point(
        data = df[is.finite(df$truth), ],
        ggplot2::aes(y = .data$truth),
        colour = "grey30", size = 1.2, alpha = 0.6
      )
    }
  }

  p
}
