#' k-fold cross-validation for pigauto trait imputation
#'
#' Performs stratified k-fold cross-validation by rotating which cells
#' serve as the test set.  Returns per-fold and aggregated metrics.
#'
#' @param data pigauto_data object (output of \code{\link{preprocess_traits}}).
#' @param tree phylo object.
#' @param k integer.  Number of folds (default 5).
#' @param seeds integer vector.  Seeds for replicate runs.
#' @param epochs integer.  Training epochs per fold (default 500).
#' @param verbose logical.  Print progress (default TRUE).
#' @param ... Additional arguments passed to \code{\link{fit_pigauto}}
#'   (e.g. \code{hidden_dim}, \code{k_eigen}, \code{use_attention}).
#' @return A list of class \code{"pigauto_cv"} with:
#'   \describe{
#'     \item{results}{Data.frame with columns: fold, rep, trait, type,
#'       metric, value.}
#'     \item{summary}{Data.frame with mean and sd across folds/reps for
#'       each trait + metric.}
#'     \item{conformal_coverage}{Data.frame of coverage per trait across
#'       folds (if available).}
#'     \item{k}{Number of folds.}
#'     \item{n_reps}{Number of replicates.}
#'   }
#' @examples
#' \dontrun{
#' data(avonet300, tree300, package = "pigauto")
#' traits <- avonet300; rownames(traits) <- traits$Species_Key
#' traits$Species_Key <- NULL
#' pd <- preprocess_traits(traits, tree300)
#' cv <- cross_validate(pd, tree300, k = 5L, seeds = 1:3, epochs = 500L)
#' print(cv)
#' summary(cv)
#' }
#' @export
cross_validate <- function(data, tree, k = 5L, seeds = 1:3,
                           epochs = 500L, verbose = TRUE, ...) {
  if (!inherits(data, "pigauto_data")) {
    stop("'data' must be a pigauto_data object (output of preprocess_traits).")
  }
  if (!inherits(tree, "phylo")) stop("'tree' must be a phylo object.")
  k <- as.integer(k)
  if (k < 2L) stop("'k' must be >= 2.")

  X_scaled  <- data$X_scaled
  trait_map <- data$trait_map
  n         <- nrow(X_scaled)
  p         <- ncol(X_scaled)

  # Pre-compute graph (depends only on tree, reuse across all folds)
  dots <- list(...)
  k_eigen <- dots$k_eigen %||% 8L
  if (verbose) message("Building phylogenetic graph (reused across folds)...")
  graph <- build_phylo_graph(tree, k_eigen = k_eigen)

  # Collect per-fold results

  all_results   <- list()
  all_conformal <- list()

  for (rep_i in seq_along(seeds)) {
    seed <- seeds[rep_i]

    # Create k non-overlapping folds from observed cells
    folds <- make_cv_folds(X_scaled, k, seed, trait_map)

    for (fold_i in seq_len(k)) {
      if (verbose) {
        message(sprintf("Rep %d/%d, Fold %d/%d (seed=%d)...",
                        rep_i, length(seeds), fold_i, k, seed))
      }

      # Build splits: fold_i = test, fold_i+1 (mod k) = val, rest = train
      val_fold <- ((fold_i %% k) + 1L)  # next fold wraps around
      test_idx <- folds[[fold_i]]
      val_idx  <- folds[[val_fold]]

      splits <- list(
        val_idx  = val_idx,
        test_idx = test_idx,
        n        = n,
        p        = p
      )

      # Fit pigauto for this fold
      fit <- fit_pigauto(
        data    = data,
        tree    = tree,
        splits  = splits,
        graph   = graph,
        epochs  = as.integer(epochs),
        seed    = seed,
        verbose = FALSE,
        ...
      )

      # Predict on test fold
      pred <- predict(fit, return_se = TRUE)

      # Evaluate using evaluate_imputation (test split)
      eval_df <- evaluate_imputation(
        pred   = pred,
        truth  = X_scaled,
        splits = splits
      )

      # Keep only test-set rows
      eval_test <- eval_df[eval_df$split == "test", , drop = FALSE]

      if (nrow(eval_test) > 0L) {
        # Pivot to long format: one row per (fold, rep, trait, type, metric)
        fold_results <- pivot_eval_to_long(eval_test, fold_i, rep_i)
        all_results[[length(all_results) + 1L]] <- fold_results
      }

      # Conformal coverage (if available)
      if (!is.null(fit$conformal_scores) && !is.null(trait_map)) {
        cov_df <- compute_conformal_coverage_fold(
          pred, X_scaled, splits$test_idx, trait_map, fold_i, rep_i
        )
        if (!is.null(cov_df)) {
          all_conformal[[length(all_conformal) + 1L]] <- cov_df
        }
      }
    }
  }

  # Combine results
  if (length(all_results) == 0L) {
    results <- data.frame(
      fold = integer(0), rep = integer(0), trait = character(0),
      type = character(0), metric = character(0), value = double(0),
      stringsAsFactors = FALSE
    )
  } else {
    results <- do.call(rbind, all_results)
    rownames(results) <- NULL
  }

  # Aggregate: mean and sd per trait + metric
  summary_df <- aggregate_cv_results(results)

  # Conformal coverage
  conformal_coverage <- NULL
  if (length(all_conformal) > 0L) {
    conformal_coverage <- do.call(rbind, all_conformal)
    rownames(conformal_coverage) <- NULL
  }

  structure(
    list(
      results            = results,
      summary            = summary_df,
      conformal_coverage = conformal_coverage,
      k                  = k,
      n_reps             = length(seeds)
    ),
    class = "pigauto_cv"
  )
}


# ---- Internal: create k non-overlapping folds from observed cells -----------
#
# For each trait, the observed cells are randomly assigned to k folds.
# This ensures every fold has cells from every trait (when possible).
# Returns a list of k integer vectors (linear column-major indices into
# the n x p matrix).

make_cv_folds <- function(X_scaled, k, seed, trait_map = NULL) {
  set.seed(seed)
  n <- nrow(X_scaled)
  p <- ncol(X_scaled)

  if (!is.null(trait_map)) {
    # Trait-level folding: assign observed (species, trait) cells to folds,
    # then expand to latent indices.
    n_traits <- length(trait_map)
    folds_latent <- vector("list", k)
    for (i in seq_len(k)) folds_latent[[i]] <- integer(0)

    for (ti in seq_along(trait_map)) {
      tm <- trait_map[[ti]]
      lc <- tm$latent_cols

      # Observed rows for this trait (check first latent column)
      obs_rows <- which(!is.na(X_scaled[, lc[1]]))
      if (length(obs_rows) == 0L) next

      # Randomly assign observed rows to folds
      fold_assign <- sample(rep(seq_len(k), length.out = length(obs_rows)))

      for (fi in seq_len(k)) {
        rows_fi <- obs_rows[fold_assign == fi]
        if (length(rows_fi) == 0L) next
        # Expand to latent column indices (column-major)
        for (lj in lc) {
          folds_latent[[fi]] <- c(folds_latent[[fi]],
                                  rows_fi + (lj - 1L) * n)
        }
      }
    }
    return(folds_latent)
  }

  # Fallback: no trait_map, fold at individual cell level
  observed <- which(!is.na(X_scaled))
  if (length(observed) == 0L) {
    stop("No observed cells in X_scaled.")
  }
  fold_assign <- sample(rep(seq_len(k), length.out = length(observed)))
  lapply(seq_len(k), function(i) observed[fold_assign == i])
}


# ---- Internal: pivot evaluate_imputation output to long format ---------------

pivot_eval_to_long <- function(eval_df, fold, rep) {
  metric_cols <- c("rmse", "pearson_r", "coverage_95", "mae",
                   "spearman_rho", "accuracy", "brier")
  rows <- list()

  for (i in seq_len(nrow(eval_df))) {
    for (mc in metric_cols) {
      val <- eval_df[[mc]][i]
      if (!is.na(val)) {
        rows[[length(rows) + 1L]] <- data.frame(
          fold   = fold,
          rep    = rep,
          trait  = eval_df$trait[i],
          type   = eval_df$type[i],
          metric = mc,
          value  = val,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}


# ---- Internal: aggregate CV results to mean +/- sd --------------------------

aggregate_cv_results <- function(results) {
  if (nrow(results) == 0L) {
    return(data.frame(
      trait = character(0), type = character(0), metric = character(0),
      mean = double(0), sd = double(0), n_folds = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  keys <- unique(results[, c("trait", "type", "metric")])
  rows <- vector("list", nrow(keys))

  for (i in seq_len(nrow(keys))) {
    mask <- results$trait  == keys$trait[i] &
            results$type   == keys$type[i] &
            results$metric == keys$metric[i]
    vals <- results$value[mask]

    rows[[i]] <- data.frame(
      trait   = keys$trait[i],
      type    = keys$type[i],
      metric  = keys$metric[i],
      mean    = mean(vals, na.rm = TRUE),
      sd      = stats::sd(vals, na.rm = TRUE),
      n_folds = length(vals),
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


# ---- Internal: conformal coverage on test fold cells -------------------------

compute_conformal_coverage_fold <- function(pred, X_scaled, test_idx,
                                            trait_map, fold, rep) {
  if (is.null(pred$conformal_lower) || is.null(pred$conformal_upper)) {
    return(NULL)
  }

  n <- nrow(X_scaled)
  rows <- list()

  for (tm in trait_map) {
    nm <- tm$name
    lc <- tm$latent_cols

    if (!(tm$type %in% c("continuous", "count", "ordinal"))) next

    # Find test rows for this trait's first latent column
    row_i <- ((test_idx - 1L) %% n) + 1L
    col_j <- ((test_idx - 1L) %/% n) + 1L
    keep  <- which(col_j == lc[1])
    ri    <- row_i[keep]
    if (length(ri) == 0L) next

    # Truth in original scale
    truth_latent <- X_scaled[ri, lc[1]]
    ok <- is.finite(truth_latent)
    if (sum(ok) == 0L) next

    # Back-transform truth to original scale for comparison with intervals
    if (tm$type == "continuous") {
      truth_orig <- truth_latent[ok] * tm$sd + tm$mean
      if (isTRUE(tm$log_transform)) truth_orig <- exp(truth_orig)
    } else if (tm$type == "count") {
      truth_orig <- expm1(truth_latent[ok] * tm$sd + tm$mean)
    } else {
      truth_orig <- truth_latent[ok] * tm$sd + tm$mean
    }

    lower <- pred$conformal_lower[ri[ok], nm]
    upper <- pred$conformal_upper[ri[ok], nm]

    covered <- (truth_orig >= lower) & (truth_orig <= upper)
    coverage <- mean(covered, na.rm = TRUE)

    rows[[length(rows) + 1L]] <- data.frame(
      fold     = fold,
      rep      = rep,
      trait    = nm,
      type     = tm$type,
      coverage = coverage,
      n        = sum(ok),
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}


# ---- Print and summary methods -----------------------------------------------

#' @export
print.pigauto_cv <- function(x, ...) {
  cat(sprintf("%d-fold cross-validation (%d replicate%s)\n",
              x$k, x$n_reps, if (x$n_reps > 1L) "s" else ""))
  rule <- paste(rep("\u2500", 50), collapse = "")
  cat(rule, "\n")


  if (nrow(x$summary) == 0L) {
    cat("No results.\n")
    return(invisible(x))
  }

  # Build display table
  traits <- unique(x$summary$trait)
  for (tr in traits) {
    sub <- x$summary[x$summary$trait == tr, , drop = FALSE]
    tp  <- sub$type[1]
    cat(sprintf("%-25s [%s]\n", tr, tp))

    for (i in seq_len(nrow(sub))) {
      cat(sprintf("  %-18s %7.4f +/- %.4f  (n=%d)\n",
                  sub$metric[i], sub$mean[i], sub$sd[i], sub$n_folds[i]))
    }
  }

  # Conformal coverage summary
  if (!is.null(x$conformal_coverage) && nrow(x$conformal_coverage) > 0L) {
    cat("\nConformal coverage (95% target):\n")
    cov_agg <- stats::aggregate(
      coverage ~ trait + type,
      data = x$conformal_coverage,
      FUN = function(v) c(mean = mean(v), sd = stats::sd(v))
    )
    cov_flat <- cbind(
      cov_agg[, c("trait", "type")],
      as.data.frame(cov_agg$coverage)
    )
    for (i in seq_len(nrow(cov_flat))) {
      cat(sprintf("  %-25s %.3f +/- %.3f\n",
                  cov_flat$trait[i], cov_flat$mean[i], cov_flat$sd[i]))
    }
  }

  invisible(x)
}

#' @export
summary.pigauto_cv <- function(object, ...) {
  cat(sprintf("%d-fold cross-validation (%d replicate%s)\n",
              object$k, object$n_reps,
              if (object$n_reps > 1L) "s" else ""))
  rule <- paste(rep("\u2500", 70), collapse = "")
  cat(rule, "\n")

  if (nrow(object$summary) == 0L) {
    cat("No results.\n")
    return(invisible(object$summary))
  }

  # Wide-format display with key metrics per trait type
  traits <- unique(object$summary$trait)

  # Header
  cat(sprintf("%-25s %-12s %15s %15s %15s\n",
              "Trait", "Type", "RMSE (mean+/-sd)",
              "r (mean+/-sd)", "Acc (mean+/-sd)"))
  cat(paste(rep("\u2500", 85), collapse = ""), "\n")

  for (tr in traits) {
    sub <- object$summary[object$summary$trait == tr, , drop = FALSE]
    tp  <- sub$type[1]

    # Extract key metrics
    fmt_metric <- function(metric_name) {
      row <- sub[sub$metric == metric_name, , drop = FALSE]
      if (nrow(row) == 0L) return("\u2500")
      sprintf("%.3f +/- %.3f", row$mean[1], row$sd[1])
    }

    rmse_str <- fmt_metric("rmse")
    r_str    <- if (tp == "continuous") {
      fmt_metric("pearson_r")
    } else if (tp == "ordinal") {
      s <- fmt_metric("spearman_rho")
      if (s != "\u2500") paste0(s, " (rho)") else s
    } else {
      "\u2500"
    }
    acc_str <- fmt_metric("accuracy")

    cat(sprintf("%-25s %-12s %15s %15s %15s\n",
                tr, tp, rmse_str, r_str, acc_str))
  }

  invisible(object$summary)
}
