# Base R plotting functions for pigauto objects
#
# These functions use only base R graphics (no ggplot2 dependency).
# They supersede the ggplot2-based methods in plot_pigauto.R.
# To avoid an S3 method conflict, the old plot.pigauto_fit in
# plot_pigauto.R should be removed or renamed once these are adopted.


# -- Colour constants ---------------------------------------------------------

.pigauto_colours <- list(
  gnn      = "#2b8cbe",
  baseline = "#636363",
  worse    = "#e34a33",
  better   = "#31a354",
  # Trait type palette
  continuous  = "#31a354",
  count       = "#2b8cbe",
  ordinal     = "#e6ab02",
  binary      = "#e34a33",
  categorical = "#7570b3"
)


# -- 1. plot.pigauto_fit ------------------------------------------------------

#' Plot diagnostics for a fitted pigauto model
#'
#' S3 \code{plot} method for \code{pigauto_fit} objects, using base R
#' graphics.  Three plot types are available: training history,
#' calibrated gate values, and conformal prediction scores.
#'
#' @param x An object of class \code{"pigauto_fit"}.
#' @param type Character.
#'   \code{"history"} (default): 2x2 panel of training loss components
#'     (reconstruction, shrinkage, gate regularisation) and validation
#'     loss over epochs.
#'   \code{"gates"}: bar plot of calibrated gate values per trait,
#'     coloured by trait type (green = continuous, blue = count,
#'     orange = ordinal, red = binary, purple = categorical).
#'   \code{"conformal"}: bar plot of conformal prediction scores per
#'     trait with a reference line at the median score.
#' @param ... Additional arguments passed to base plotting functions.
#' @return Invisible \code{NULL}. Called for its side effect (plotting).
#' @examples
#' \dontrun{
#' fit <- fit_pigauto(data, tree)
#' plot(fit)
#' plot(fit, type = "history")
#' plot(fit, type = "gates")
#' plot(fit, type = "conformal")
#' }
#' @importFrom graphics par plot plot.new points lines abline barplot axis text
#'   legend mtext segments arrows rect polygon boxplot hist
#' @importFrom grDevices adjustcolor
#' @importFrom stats predict setNames
#' @export
plot.pigauto_fit <- function(x, type = "history", ...) {
  type <- match.arg(type, c("history", "gates", "conformal"))

  if (type == "history") {
    .plot_history(x, ...)
  } else if (type == "gates") {
    .plot_gates(x, ...)
  } else {
    .plot_conformal_fit(x, ...)
  }
  invisible(NULL)
}


# -- History panel (2x2) -----------------------------------------------------

.plot_history <- function(x, ...) {
  h <- x$history
  if (is.null(h) || nrow(h) == 0L) {
    stop("No training history available in this pigauto_fit object.")
  }

  op <- par(mfrow = c(2, 2), mar = c(4, 4.5, 2.5, 1), cex.lab = 1.1,
            cex.main = 1.1, family = "")
  on.exit(par(op), add = TRUE)

  col_line <- .pigauto_colours$gnn
  col_fill <- adjustcolor(col_line, alpha.f = 0.15)

  # Helper for each panel
  .panel <- function(y, ylab, main) {
    ok <- is.finite(y)
    if (sum(ok) == 0L) {
      plot.new()
      mtext(paste(main, "(no data)"), side = 3, line = 0.5, cex = 0.9)
      return(invisible(NULL))
    }
    ep <- h$epoch[ok]
    yy <- y[ok]
    plot(ep, yy, type = "n", xlab = "Epoch", ylab = ylab, main = main,
         las = 1, bty = "l", ...)
    # Light fill under line
    polygon(c(ep, rev(ep)),
            c(rep(par("usr")[3], length(ep)), rev(yy)),
            col = col_fill, border = NA)
    lines(ep, yy, col = col_line, lwd = 1.8)
  }

  .panel(h$loss_rec,    "Reconstruction loss", "Reconstruction")
  .panel(h$loss_shrink, "Shrinkage loss",      "Shrinkage")
  .panel(h$loss_gate,   "Gate penalty",        "Gate regularisation")

  # Validation loss -- support both old (val_rmse) and new (val_loss) columns
  val_col <- if ("val_loss" %in% names(h)) "val_loss" else "val_rmse"
  .panel(h[[val_col]],  "Validation loss",     "Validation")
}


# -- Gate bar plot ------------------------------------------------------------

.plot_gates <- function(x, ...) {
  gates <- x$calibrated_gates
  if (is.null(gates)) {
    stop("No calibrated gates available. Run fit_pigauto() with validation ",
         "splits to enable gate calibration.")
  }

  trait_map <- x$trait_map
  if (is.null(trait_map)) {
    # Legacy: no trait map, plot raw gate vector
    op <- par(mar = c(6, 4.5, 3, 1), las = 2, cex.lab = 1.1)
    on.exit(par(op), add = TRUE)
    nms <- x$latent_names %||% paste0("col", seq_along(gates))
    barplot(gates, names.arg = nms, col = .pigauto_colours$gnn,
            ylab = "Calibrated gate", main = "Per-column blend gates",
            border = NA, ...)
    return(invisible(NULL))
  }

  # Aggregate gates to trait level (mean over latent columns per trait)
  n_traits <- length(trait_map)
  trait_names <- vapply(trait_map, "[[", character(1), "name")
  trait_types <- vapply(trait_map, "[[", character(1), "type")
  gate_vals <- numeric(n_traits)
  for (i in seq_len(n_traits)) {
    lc <- trait_map[[i]]$latent_cols
    gate_vals[i] <- mean(gates[lc])
  }

  colours <- vapply(trait_types, function(tp) {
    .pigauto_colours[[tp]] %||% "#999999"
  }, character(1))

  op <- par(mar = c(7, 4.5, 3, 1), las = 2, cex.lab = 1.1)
  on.exit(par(op), add = TRUE)

  bp <- barplot(gate_vals, names.arg = trait_names, col = colours,
                ylab = "Calibrated gate value", border = NA,
                main = "Blend gate per trait",
                ylim = c(0, max(gate_vals, 0.01) * 1.15), ...)

  # Reference line at gate_cap
  cap <- x$model_config$gate_cap %||% 0.8
  abline(h = cap, lty = 2, col = "#636363", lwd = 0.8)
  mtext(paste0("cap = ", cap), side = 4, at = cap, las = 1,
        cex = 0.7, col = "#636363", line = 0.3)

  # Legend
  present_types <- unique(trait_types)
  legend_cols <- vapply(present_types, function(tp) {
    .pigauto_colours[[tp]] %||% "#999999"
  }, character(1))
  legend("topright", legend = present_types, fill = legend_cols,
         border = NA, bty = "n", cex = 0.85)
}


# -- Conformal scores bar plot -----------------------------------------------

.plot_conformal_fit <- function(x, ...) {
  cs <- x$conformal_scores
  if (is.null(cs)) {
    stop("No conformal scores available. Ensure fit_pigauto() was run ",
         "with validation splits.")
  }

  cs <- cs[!is.na(cs)]
  if (length(cs) == 0L) {
    stop("All conformal scores are NA.")
  }

  # Colour by trait type
  trait_map <- x$trait_map
  if (!is.null(trait_map)) {
    type_map <- stats::setNames(
      vapply(trait_map, "[[", character(1), "type"),
      vapply(trait_map, "[[", character(1), "name")
    )
    colours <- vapply(names(cs), function(nm) {
      tp <- type_map[nm]
      if (is.na(tp)) return("#999999")
      .pigauto_colours[[tp]] %||% "#999999"
    }, character(1))
  } else {
    colours <- rep(.pigauto_colours$gnn, length(cs))
  }

  op <- par(mar = c(7, 4.5, 3, 1), las = 2, cex.lab = 1.1)
  on.exit(par(op), add = TRUE)

  bp <- barplot(cs, names.arg = names(cs), col = colours,
                ylab = "Conformal score (latent scale)", border = NA,
                main = "Conformal prediction scores (95%)",
                ylim = c(0, max(cs) * 1.15), ...)

  # Reference line at median
  med <- stats::median(cs)
  abline(h = med, lty = 2, col = "#636363", lwd = 0.8)
  mtext(paste0("median = ", round(med, 3)), side = 4, at = med,
        las = 1, cex = 0.7, col = "#636363", line = 0.3)
}


# -- 2. plot.pigauto_pred ----------------------------------------------------

#' Plot predictions from a pigauto model
#'
#' S3 \code{plot} method for \code{pigauto_pred} objects, using base R
#' graphics.  Three plot types are available: scatter plots of observed
#' vs predicted values, interval plots per species, and probability
#' charts for discrete traits.
#'
#' @param x An object of class \code{"pigauto_pred"}.
#' @param data Optional numeric matrix or data.frame of observed values
#'   in original scale.  When supplied, observed values are overlaid on
#'   scatter and interval plots.
#' @param splits Optional list (output of \code{\link{make_missing_splits}}).
#'   Used to identify which cells were held out for evaluation.
#' @param type Character.
#'   \code{"scatter"} (default): grid of observed-vs-predicted scatter
#'     plots for continuous and count traits, with 1:1 line and
#'     conformal interval bands where available.
#'   \code{"intervals"}: for each trait, species sorted by predicted
#'     value with conformal-interval ribbon; observed values shown in
#'     red where available.
#'   \code{"probabilities"}: for binary traits, boxplot of predicted
#'     probabilities grouped by true class; for categorical traits,
#'     stacked bar chart of mean predicted probabilities by true class.
#' @param trait Character vector of trait names to plot.  If \code{NULL}
#'   (default), all traits appropriate for the chosen \code{type} are
#'   plotted.
#' @param ... Additional arguments passed to base plotting functions.
#' @return Invisible \code{NULL}. Called for its side effect (plotting).
#' @examples
#' \dontrun{
#' pred <- predict(fit)
#' plot(pred)
#' plot(pred, type = "scatter", data = observed_df)
#' plot(pred, type = "intervals", trait = "Beak.Length_Culmen")
#' plot(pred, type = "probabilities", data = observed_df)
#' }
#' @importFrom graphics par plot points lines abline barplot axis text
#'   legend mtext segments arrows rect polygon boxplot hist
#' @importFrom grDevices adjustcolor
#' @export
plot.pigauto_pred <- function(x, data = NULL, splits = NULL,
                              type = "scatter", trait = NULL, ...) {
  type <- match.arg(type, c("scatter", "intervals", "probabilities"))

  if (type == "scatter") {
    .plot_scatter(x, data = data, splits = splits, traits = trait, ...)
  } else if (type == "intervals") {
    .plot_intervals(x, data = data, traits = trait, ...)
  } else {
    .plot_probabilities(x, data = data, traits = trait, ...)
  }
  invisible(NULL)
}


# -- Scatter: observed vs predicted -------------------------------------------

.plot_scatter <- function(x, data = NULL, splits = NULL, traits = NULL, ...) {
  trait_map <- x$trait_map

  # Identify plottable traits (continuous and count)
  if (!is.null(trait_map)) {
    plot_traits <- Filter(function(tm) tm$type %in% c("continuous", "count"),
                          trait_map)
    if (!is.null(traits)) {
      plot_traits <- Filter(function(tm) tm$name %in% traits, plot_traits)
    }
    if (length(plot_traits) == 0L) {
      stop("No continuous or count traits available for scatter plots.")
    }
    trait_names <- vapply(plot_traits, "[[", character(1), "name")
  } else {
    trait_names <- x$trait_names
    if (!is.null(traits)) trait_names <- intersect(trait_names, traits)
    if (length(trait_names) == 0L) {
      stop("No traits available for scatter plots.")
    }
  }

  # Imputed values (original scale)
  imp <- x$imputed

  # Layout grid
  nt <- length(trait_names)
  nc <- ceiling(sqrt(nt))
  nr <- ceiling(nt / nc)

  op <- par(mfrow = c(nr, nc), mar = c(4, 4.5, 2.5, 1), cex.lab = 1.0,
            cex.main = 1.0)
  on.exit(par(op), add = TRUE)

  col_pt   <- adjustcolor(.pigauto_colours$gnn, alpha.f = 0.4)
  col_band <- adjustcolor(.pigauto_colours$gnn, alpha.f = 0.10)

  for (nm in trait_names) {
    pred_vals <- if (is.data.frame(imp)) as.numeric(imp[[nm]]) else imp[, nm]

    # Observed values
    obs_vals <- NULL
    if (!is.null(data)) {
      obs_vals <- if (is.data.frame(data)) {
        as.numeric(data[[nm]])
      } else if (nm %in% colnames(data)) {
        as.numeric(data[, nm])
      }
    }

    if (!is.null(obs_vals) && length(obs_vals) == length(pred_vals)) {
      ok <- is.finite(obs_vals) & is.finite(pred_vals)
      if (sum(ok) == 0L) {
        plot.new()
        mtext(nm, side = 3, line = 0.5, cex = 0.9)
        next
      }

      lim <- range(c(obs_vals[ok], pred_vals[ok]))
      pad <- diff(lim) * 0.05
      lim <- lim + c(-pad, pad)

      plot(obs_vals[ok], pred_vals[ok], type = "n",
           xlim = lim, ylim = lim,
           xlab = "Observed", ylab = "Predicted", main = nm,
           las = 1, bty = "l", ...)

      # Conformal interval band around 1:1 line
      if (!is.null(x$conformal_lower) && !is.null(x$conformal_upper) &&
          nm %in% colnames(x$conformal_lower)) {
        # Approximate band half-width from conformal bounds
        mid_idx <- which.min(abs(pred_vals - mean(lim)))
        if (length(mid_idx) > 0L) {
          half_w <- (x$conformal_upper[mid_idx, nm] -
                       x$conformal_lower[mid_idx, nm]) / 2
          usr <- par("usr")
          xs <- seq(usr[1], usr[2], length.out = 200)
          polygon(c(xs, rev(xs)),
                  c(xs - half_w, rev(xs + half_w)),
                  col = col_band, border = NA)
        }
      }

      abline(0, 1, lty = 2, col = "#636363", lwd = 0.8)
      points(obs_vals[ok], pred_vals[ok], pch = 16, cex = 0.7, col = col_pt)

      # Correlation annotation
      if (sum(ok) > 2L) {
        r_val <- stats::cor(obs_vals[ok], pred_vals[ok])
        mtext(paste0("r = ", round(r_val, 3)), side = 3, adj = 1,
              cex = 0.75, col = "#636363", line = 0.2)
      }
    } else {
      # No observed data -- show predicted distribution as a ranked curve
      plot(seq_along(pred_vals), sort(pred_vals), type = "n",
           xlab = "Species (rank)", ylab = "Predicted", main = nm,
           las = 1, bty = "l", ...)
      lines(seq_along(pred_vals), sort(pred_vals),
            col = .pigauto_colours$gnn, lwd = 1.5)
    }
  }
}


# -- Intervals: species sorted by prediction with ribbon ---------------------

.plot_intervals <- function(x, data = NULL, traits = NULL, ...) {
  trait_map <- x$trait_map

  # Select traits with numeric predictions
  if (!is.null(trait_map)) {
    avail <- Filter(function(tm) {
      tm$type %in% c("continuous", "count", "ordinal")
    }, trait_map)
    avail_names <- vapply(avail, "[[", character(1), "name")
  } else {
    avail_names <- x$trait_names
  }

  if (!is.null(traits)) {
    avail_names <- intersect(avail_names, traits)
  }
  if (length(avail_names) == 0L) {
    stop("No traits with interval information available.")
  }

  imp <- x$imputed
  has_conformal <- !is.null(x$conformal_lower) && !is.null(x$conformal_upper)

  nt <- length(avail_names)
  nc <- ceiling(sqrt(nt))
  nr <- ceiling(nt / nc)

  op <- par(mfrow = c(nr, nc), mar = c(4, 4.5, 2.5, 1), cex.lab = 1.0,
            cex.main = 1.0)
  on.exit(par(op), add = TRUE)

  col_line <- .pigauto_colours$gnn
  col_band <- adjustcolor(col_line, alpha.f = 0.20)
  col_obs  <- .pigauto_colours$worse

  for (nm in avail_names) {
    pred_vals <- if (is.data.frame(imp)) as.numeric(imp[[nm]]) else imp[, nm]
    n <- length(pred_vals)
    ord <- order(pred_vals)
    pred_sorted <- pred_vals[ord]
    xi <- seq_len(n)

    ylim <- range(pred_sorted, na.rm = TRUE)

    # Expand y-limits for conformal intervals
    if (has_conformal && nm %in% colnames(x$conformal_lower)) {
      lo <- x$conformal_lower[ord, nm]
      hi <- x$conformal_upper[ord, nm]
      ylim <- range(c(ylim, lo, hi), na.rm = TRUE)
    }

    pad <- diff(ylim) * 0.05
    if (pad == 0) pad <- 0.5
    ylim <- ylim + c(-pad, pad)

    plot(xi, pred_sorted, type = "n",
         xlab = "Species (sorted by prediction)", ylab = nm,
         main = nm, ylim = ylim, las = 1, bty = "l", ...)

    # Conformal ribbon
    if (has_conformal && nm %in% colnames(x$conformal_lower)) {
      polygon(c(xi, rev(xi)), c(lo, rev(hi)),
              col = col_band, border = NA)
    }

    lines(xi, pred_sorted, col = col_line, lwd = 1.5)

    # Overlay observed values
    if (!is.null(data)) {
      obs_vals <- if (is.data.frame(data)) {
        as.numeric(data[[nm]])
      } else if (nm %in% colnames(data)) {
        as.numeric(data[, nm])
      } else {
        NULL
      }

      if (!is.null(obs_vals) && length(obs_vals) == n) {
        obs_sorted <- obs_vals[ord]
        ok <- is.finite(obs_sorted)
        if (any(ok)) {
          points(xi[ok], obs_sorted[ok], pch = 16, cex = 0.5, col = col_obs)
        }
      }
    }
  }
}


# -- Probabilities: visualisation for discrete traits -------------------------

.plot_probabilities <- function(x, data = NULL, traits = NULL, ...) {
  trait_map <- x$trait_map
  probs <- x$probabilities

  if (is.null(probs) || length(probs) == 0L) {
    stop("No probability predictions available (binary/categorical traits).")
  }

  avail_names <- names(probs)
  if (!is.null(traits)) {
    avail_names <- intersect(avail_names, traits)
  }
  if (length(avail_names) == 0L) {
    stop("No matching probability traits found.")
  }

  nt <- length(avail_names)
  nc <- ceiling(sqrt(nt))
  nr <- ceiling(nt / nc)

  op <- par(mfrow = c(nr, nc), mar = c(5, 4.5, 3, 1), cex.lab = 1.0,
            cex.main = 1.0)
  on.exit(par(op), add = TRUE)

  for (nm in avail_names) {
    prob_data <- probs[[nm]]

    # Identify trait type
    tp <- "binary"
    tm_info <- NULL
    if (!is.null(trait_map)) {
      for (tm in trait_map) {
        if (tm$name == nm) { tm_info <- tm; tp <- tm$type; break }
      }
    }

    if (tp == "binary") {
      .plot_binary_probs(nm, prob_data, data, tm_info, ...)
    } else if (tp == "categorical") {
      .plot_categorical_probs(nm, prob_data, data, tm_info, ...)
    }
  }
}


.plot_binary_probs <- function(nm, prob_vec, data, tm_info, ...) {
  # Get true classes if available
  truth <- NULL
  if (!is.null(data)) {
    truth <- if (is.data.frame(data)) data[[nm]] else {
      if (nm %in% colnames(data)) data[, nm] else NULL
    }
  }

  if (!is.null(truth) && !is.null(tm_info)) {
    lvls <- tm_info$levels
    truth_fac <- factor(truth, levels = lvls)

    col_box <- c(adjustcolor(.pigauto_colours$worse, alpha.f = 0.4),
                 adjustcolor(.pigauto_colours$gnn, alpha.f = 0.4))

    ok <- !is.na(truth_fac)
    if (any(ok)) {
      boxplot(prob_vec[ok] ~ truth_fac[ok],
              col = col_box,
              border = c(.pigauto_colours$worse, .pigauto_colours$gnn),
              ylab = paste0("P(", lvls[2], ")"),
              xlab = "True class", main = nm,
              las = 1, bty = "l", outline = FALSE, ...)
      abline(h = 0.5, lty = 2, col = "#636363", lwd = 0.8)
      # Jittered points
      for (k in seq_along(lvls)) {
        idx <- which(truth_fac[ok] == lvls[k])
        if (length(idx) > 0L) {
          jit <- jitter(rep(k, length(idx)), amount = 0.15)
          points(jit, prob_vec[ok][idx], pch = 16, cex = 0.4,
                 col = adjustcolor("#333333", alpha.f = 0.4))
        }
      }
      return(invisible(NULL))
    }
  }

  # No truth available: histogram of predicted probabilities
  hist(prob_vec, breaks = 20,
       col = adjustcolor(.pigauto_colours$gnn, 0.4),
       border = .pigauto_colours$gnn, main = nm,
       xlab = "Predicted probability", ylab = "Frequency", las = 1, ...)
  abline(v = 0.5, lty = 2, col = "#636363", lwd = 0.8)
}


.plot_categorical_probs <- function(nm, prob_mat, data, tm_info, ...) {
  if (!is.matrix(prob_mat) || ncol(prob_mat) < 2L) {
    plot.new()
    mtext(paste(nm, "-- not enough categories"), side = 3, cex = 0.9)
    return(invisible(NULL))
  }

  lvls <- colnames(prob_mat)
  if (is.null(lvls)) lvls <- paste0("C", seq_len(ncol(prob_mat)))
  K <- ncol(prob_mat)

  # Generate K distinct colours
  pal <- if (K <= 5L) {
    c("#2b8cbe", "#31a354", "#e6ab02", "#e34a33", "#7570b3")[seq_len(K)]
  } else {
    grDevices::hcl.colors(K, palette = "Dark 3")
  }

  # Group by true class if available
  truth <- NULL
  if (!is.null(data)) {
    truth <- if (is.data.frame(data)) data[[nm]] else {
      if (nm %in% colnames(data)) data[, nm] else NULL
    }
  }

  if (!is.null(truth)) {
    truth_fac <- factor(truth, levels = lvls)
    ok <- !is.na(truth_fac)

    if (any(ok)) {
      # Average predicted probabilities by true class
      avg_probs <- matrix(0, K, K, dimnames = list(lvls, lvls))
      for (k in seq_len(K)) {
        idx <- which(truth_fac[ok] == lvls[k])
        if (length(idx) > 0L) {
          avg_probs[k, ] <- colMeans(prob_mat[which(ok)[idx], , drop = FALSE])
        }
      }

      barplot(t(avg_probs), beside = FALSE, col = pal,
              names.arg = lvls, ylab = "Mean predicted probability",
              xlab = "True class", main = nm, las = 1,
              border = NA, ...)
      legend("topright", legend = rev(lvls), fill = rev(pal),
             border = NA, bty = "n", cex = 0.75, title = "Predicted")
      return(invisible(NULL))
    }
  }

  # No truth: stacked bar of all predictions sorted by dominant class
  dominant <- apply(prob_mat, 1, which.max)
  ord <- order(dominant, -apply(prob_mat, 1, max))
  # Subsample if too many species to keep bars legible
  if (length(ord) > 100L) {
    idx_sub <- round(seq(1, length(ord), length.out = 100))
    ord <- ord[idx_sub]
  }

  barplot(t(prob_mat[ord, ]), col = pal, border = NA,
          ylab = "Predicted probability", xlab = "Species",
          main = nm, las = 1, space = 0, ...)
  legend("topright", legend = rev(lvls), fill = rev(pal),
         border = NA, bty = "n", cex = 0.75)
}


# -- 3. plot_comparison -------------------------------------------------------

#' Forest-plot style comparison of benchmark results
#'
#' Takes a data.frame of benchmark results and produces a forest-plot
#' style figure comparing baseline and GNN performance per trait.
#' Each trait appears on the y-axis; points mark RMSE (continuous/count)
#' or accuracy (binary/categorical) for the two methods, connected by
#' a horizontal line colour-coded by whether the GNN improves over the
#' baseline.
#'
#' When both RMSE and accuracy metrics are present and \code{metric} is
#' \code{NULL}, a two-panel figure is produced automatically.
#'
#' @param results A data.frame with columns \code{trait}, \code{type},
#'   \code{metric}, \code{method}, and \code{value}.  Typically the
#'   output of benchmark scripts or reshaped output from
#'   \code{\link{evaluate_imputation}}.
#' @param metric Character.  Which metric to compare.
#'   Default \code{NULL} auto-selects: \code{"rmse"} for
#'   continuous/count traits and \code{"accuracy"} for
#'   binary/categorical traits.  If a single metric name is supplied,
#'   only traits with that metric are shown.
#' @param methods Character vector of length 2: the baseline method
#'   name and the GNN method name.  Default
#'   \code{c("BM_baseline", "pigauto_GNN")}.
#' @param ... Additional arguments passed to \code{plot()}.
#' @return Invisible \code{NULL}. Called for its side effect (plotting).
#' @examples
#' \dontrun{
#' results <- read.csv("benchmark_results.csv")
#' plot_comparison(results)
#' plot_comparison(results, metric = "rmse")
#' }
#' @importFrom graphics par plot points lines abline barplot axis text
#'   legend mtext segments arrows rect
#' @importFrom grDevices adjustcolor
#' @export
plot_comparison <- function(results, metric = NULL,
                            methods = c("BM_baseline", "pigauto_GNN"),
                            ...) {
  if (!is.data.frame(results)) {
    stop("'results' must be a data.frame.")
  }
  required <- c("trait", "type", "metric", "method", "value")
  missing_cols <- setdiff(required, names(results))
  if (length(missing_cols) > 0L) {
    stop("Missing columns in results: ", paste(missing_cols, collapse = ", "))
  }

  # Filter to the two methods
  d <- results[results$method %in% methods, ]
  if (nrow(d) == 0L) {
    stop("No rows match methods: ", paste(methods, collapse = ", "))
  }

  # Auto-select metric when NULL
  if (is.null(metric)) {
    d_rmse <- d[d$metric == "rmse", ]
    d_acc  <- d[d$metric == "accuracy", ]

    if (nrow(d_rmse) > 0L && nrow(d_acc) > 0L) {
      # Two-panel layout: RMSE on the left, accuracy on the right
      op <- par(mfrow = c(1, 2), mar = c(4, 9, 3, 2), cex.lab = 1.0)
      on.exit(par(op), add = TRUE)
      .forest_panel(d_rmse, "rmse", methods, lower_is_better = TRUE, ...)
      .forest_panel(d_acc, "accuracy", methods, lower_is_better = FALSE, ...)
      return(invisible(NULL))
    } else if (nrow(d_rmse) > 0L) {
      metric <- "rmse"
    } else if (nrow(d_acc) > 0L) {
      metric <- "accuracy"
    } else {
      stop("No rmse or accuracy data found in results.")
    }
  }

  d <- d[d$metric == metric, ]
  if (nrow(d) == 0L) {
    stop("No data for metric '", metric, "'.")
  }

  lower_is_better <- metric %in% c("rmse", "mae", "brier")

  op <- par(mar = c(4, 9, 3, 2), cex.lab = 1.0)
  on.exit(par(op), add = TRUE)
  .forest_panel(d, metric, methods, lower_is_better, ...)
  invisible(NULL)
}


.forest_panel <- function(d, metric, methods, lower_is_better, ...) {
  # Aggregate across replicates (mean)
  agg <- stats::aggregate(value ~ trait + type + method, data = d, FUN = mean)

  traits <- unique(agg$trait)
  n_traits <- length(traits)

  if (n_traits == 0L) {
    plot.new()
    mtext(paste("No traits for", metric), side = 3, cex = 0.9)
    return(invisible(NULL))
  }

  # Build comparison table
  comp <- data.frame(trait = traits, stringsAsFactors = FALSE)
  comp$type <- vapply(traits, function(tr) {
    tp <- unique(agg$type[agg$trait == tr])
    tp[1]
  }, character(1))
  comp$baseline <- vapply(traits, function(tr) {
    v <- agg$value[agg$trait == tr & agg$method == methods[1]]
    if (length(v) == 0L) NA_real_ else mean(v)
  }, numeric(1))
  comp$gnn <- vapply(traits, function(tr) {
    v <- agg$value[agg$trait == tr & agg$method == methods[2]]
    if (length(v) == 0L) NA_real_ else mean(v)
  }, numeric(1))

  # Sort by baseline value
  comp <- comp[order(comp$baseline), ]

  # Determine whether GNN is better
  if (lower_is_better) {
    comp$gnn_better <- comp$gnn < comp$baseline
  } else {
    comp$gnn_better <- comp$gnn > comp$baseline
  }

  col_better <- .pigauto_colours$better
  col_worse  <- .pigauto_colours$worse
  col_gnn    <- .pigauto_colours$gnn
  col_base   <- .pigauto_colours$baseline

  # Y positions

  yi <- seq_len(n_traits)
  xlim <- range(c(comp$baseline, comp$gnn), na.rm = TRUE)
  pad <- diff(xlim) * 0.08
  if (pad == 0) pad <- 0.1
  xlim <- xlim + c(-pad, pad)

  plot(NA, xlim = xlim, ylim = c(0.5, n_traits + 0.5),
       xlab = toupper(metric), ylab = "", yaxt = "n",
       main = paste0("Method comparison: ", metric),
       las = 1, bty = "l", ...)

  axis(2, at = yi, labels = comp$trait, las = 1, tick = FALSE,
       cex.axis = 0.85)

  # Connecting segments colour-coded by direction
  for (i in yi) {
    seg_col <- if (is.na(comp$gnn_better[i])) {
      "#cccccc"
    } else if (comp$gnn_better[i]) {
      col_better
    } else {
      col_worse
    }
    segments(comp$baseline[i], i, comp$gnn[i], i,
             col = seg_col, lwd = 2)
  }

  # Method points
  points(comp$baseline, yi, pch = 16, cex = 1.2, col = col_base)
  points(comp$gnn, yi, pch = 17, cex = 1.2, col = col_gnn)

  # Light horizontal grid
  abline(h = yi, col = adjustcolor("#cccccc", 0.3), lwd = 0.5)

  legend("bottomright",
         legend = c(methods[1], methods[2], "GNN better", "GNN worse"),
         pch = c(16, 17, NA, NA),
         lty = c(NA, NA, 1, 1),
         lwd = c(NA, NA, 2, 2),
         col = c(col_base, col_gnn, col_better, col_worse),
         bty = "n", cex = 0.8)
}
