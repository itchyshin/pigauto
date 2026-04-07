#' Plot training history or observed-vs-predicted scatter
#'
#' @description
#' Produces either a validation loss training curve (\code{type = "history"})
#' or a scatter plot of observed vs predicted values for the selected split
#' (\code{type = "scatter"}).
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

    # Support both old (val_rmse) and new (val_loss) column names
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
  } else {
    splits <- x$splits
    if (is.null(splits)) stop("No splits stored; cannot produce scatter plot.")
    idx <- if (split == "val") splits$val_idx else splits$test_idx
    if (length(idx) == 0L) stop("No ", split, " cells available.")

    # Only plot continuous/count/ordinal traits (numeric in latent space)
    trait_map <- x$trait_map
    n <- length(x$species_names)
    p <- if (!is.null(x$baseline$mu)) ncol(x$baseline$mu) else length(x$trait_names)

    row_j <- ((idx - 1L) %% n) + 1L
    col_j <- ceiling(idx / n)

    # Determine which latent columns are plottable (continuous-like)
    if (!is.null(trait_map)) {
      cont_lc <- integer(0)
      cont_names_map <- character(0)
      for (tm in trait_map) {
        if (tm$type %in% c("continuous", "count", "ordinal")) {
          cont_lc <- c(cont_lc, tm$latent_cols)
          cont_names_map <- c(cont_names_map, rep(tm$name, tm$n_latent))
        }
      }
      keep <- col_j %in% cont_lc
      if (sum(keep) == 0L) {
        stop("No continuous/count/ordinal traits to scatter-plot.")
      }
      row_j <- row_j[keep]
      col_j <- col_j[keep]
      trait_labels <- cont_names_map[match(col_j, cont_lc)]
    } else {
      trait_labels <- x$trait_names[col_j]
    }

    bm_vals <- mapply(function(r, c) x$baseline$mu[r, c], row_j, col_j)

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
        x = "BM prediction (latent scale)",
        y = "BM prediction (latent scale)"
      )
  }
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
