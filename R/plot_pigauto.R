#' Plot training history or observed-vs-predicted scatter
#'
#' @description
#' Produces either a validation RMSE training curve (\code{type = "history"})
#' or a scatter plot of BM baseline predictions against themselves for the
#' selected split (\code{type = "scatter"}).
#'
#' @param x object of class \code{"pigauto_fit"}.
#' @param type \code{"scatter"} (default) or \code{"history"}.
#' @param split which split to use for scatter: \code{"val"} (default) or
#'   \code{"test"}.
#' @param ... ignored.
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' plot(fit)
#' plot(fit, type = "history")
#' }
#' @importFrom ggplot2 ggplot aes geom_point geom_abline facet_wrap
#' @importFrom ggplot2 theme_minimal labs geom_line geom_ribbon
#' @importFrom rlang .data
#' @export
plot.pigauto_fit <- function(x, type = "scatter", split = "val", ...) {
  type  <- match.arg(type,  c("scatter", "history"))
  split <- match.arg(split, c("val", "test"))

  if (type == "history") {
    if (nrow(x$history) == 0L) stop("No training history available.")
    df <- x$history[is.finite(x$history$val_rmse), ]
    ggplot2::ggplot(df, ggplot2::aes(x = .data$epoch, y = .data$val_rmse)) +
      ggplot2::geom_line(colour = "#1f78b4") +
      ggplot2::geom_point(size = 1.5, colour = "#1f78b4") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Training history",
        x = "Epoch", y = "Validation RMSE"
      )
  } else {
    splits <- x$splits
    if (is.null(splits)) stop("No splits stored; cannot produce scatter plot.")
    idx <- if (split == "val") splits$val_idx else splits$test_idx
    if (length(idx) == 0L) stop("No ", split, " cells available.")

    n       <- length(x$species_names)
    row_j   <- ((idx - 1L) %% n) + 1L
    col_j   <- ceiling(idx / n)

    bm_vals <- mapply(
      \(r, c) x$norm$means[[c]] + x$norm$sds[[c]] * x$baseline$mu[r, c],
      row_j, col_j
    )
    trait_labels <- x$trait_names[col_j]

    df <- data.frame(
      bm    = bm_vals,
      trait = trait_labels,
      stringsAsFactors = FALSE
    )

    ggplot2::ggplot(df, ggplot2::aes(x = .data$bm, y = .data$bm)) +
      ggplot2::geom_point(alpha = 0.3, size = 1) +
      ggplot2::geom_abline(linetype = "dashed", colour = "grey40") +
      ggplot2::facet_wrap(~ .data$trait, scales = "free") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("BM baseline:", split, "split"),
        x = "BM prediction (original scale)", y = "BM prediction (original scale)"
      )
  }
}


#' Plot uncertainty ribbons for imputed trait values
#'
#' @description
#' Plots 95\% credible intervals around imputed values, sorted by predicted
#' value. If \code{truth} is supplied, observed values are overlaid as points.
#'
#' @param pred_result list with \code{imputed} and \code{se} matrices
#'   (output of \code{predict.pigauto_fit}).
#' @param truth matrix of true values (same scale as \code{pred_result}), or
#'   \code{NULL}.
#' @param trait_name character. Which trait to plot (must match column name).
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point
#' @importFrom ggplot2 theme_minimal labs
#' @importFrom rlang .data
#' @export
plot_uncertainty <- function(pred_result, truth = NULL, trait_name) {
  if (!trait_name %in% colnames(pred_result$imputed)) {
    stop("'", trait_name, "' not found in pred_result$imputed.")
  }

  pred_j <- pred_result$imputed[, trait_name]
  se_j   <- pred_result$se[, trait_name]
  idx    <- order(pred_j)

  df <- data.frame(
    i     = seq_along(pred_j)[idx],
    pred  = pred_j[idx],
    lower = (pred_j - 1.96 * se_j)[idx],
    upper = (pred_j + 1.96 * se_j)[idx]
  )

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
      x = "Species (sorted by prediction)", y = trait_name
    )

  if (!is.null(truth)) {
    truth_j <- truth[, trait_name]
    df$truth <- truth_j[idx]
    p <- p + ggplot2::geom_point(
      data = df[is.finite(df$truth), ],
      ggplot2::aes(y = .data$truth),
      colour = "grey30", size = 1.2, alpha = 0.6
    )
  }

  p
}
