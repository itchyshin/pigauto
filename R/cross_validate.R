#' k-fold cross-validation for pigauto trait imputation
#'
#' Evaluates trait-imputation performance by rotating held-out values through
#' k test folds. Hold out trait cells, whole species, or rooted clades.
#' Returns per-fold and aggregated metrics.
#'
#' @param data pigauto_data object (output of \code{\link{preprocess_traits}}).
#' @param tree phylo object.
#' @param k integer.  Number of folds (default 5).
#' @param seeds optional integer vector. Seeds for replicate runs. The default
#'   \code{NULL} performs one run using the current RNG stream.
#' @param epochs integer.  Training epochs per fold (default 500).
#' @param verbose logical.  Print progress (default TRUE).
#' @param holdout One of \code{"cell"} (default), \code{"species"},
#'   or \code{"clade"}. Species and clade modes mask every trait and
#'   every observation of each held-out species. They require at least three
#'   folds and retain the full tree, including held-out tips.
#' @param ... Additional arguments passed to \code{\link{fit_pigauto}}
#'   (e.g. \code{hidden_dim}, \code{k_eigen}, \code{use_attention}).
#' @details
#' \code{cross_validate()} precomputes one graph shared across folds and uses
#' \code{k_eigen = 8L} for that graph when \code{k_eigen} is omitted. An
#' explicit \code{k_eigen} is used for the shared graph and forwarded to
#' \code{fit_pigauto()}. Cell folds divide observed cells within each trait.
#' Species folds assign entire species at random; clade folds partition the
#' rooted tree into k disjoint rooted subtrees by repeatedly splitting the
#' largest eligible subtree. A rooted \code{ape::phylo} representation is
#' required, as determined by \code{ape::is.rooted()}; a root polytomy may
#' require a root edge. Clade folds can be unequal in size. Multifurcating
#' trees can defeat this greedy partition even when a partition exists. None
#' of the modes are
#' stratified by outcome class. The next fold is used for validation, and the
#' rest for training. Grouped folds require at least two training observations
#' and two training species for each trait.
#'
#' In species and clade modes, trait centering, scaling, observed ranges, and
#' other trait-map statistics are recalculated using the training rows only.
#' The trait types, factor levels and log-transform choices from the supplied
#' preprocessed data remain fixed across folds. Specify these in advance when
#' they would otherwise be inferred from held-out values. The phylogeny and
#' any supplied covariates remain available for held-out species. Grouped
#' modes require data freshly constructed by the current preprocess_traits().
#' Cell-mode
#' preprocessing retains its existing whole-data behaviour.
#' Interpret each metric on its documented scale in
#' \code{\link{evaluate_imputation}}; continuous and count RMSE use latent
#' scales. Conformal coverage is calculated separately after predictions are
#' decoded. Clade-structured holdouts need not be exchangeable with the
#' validation cells, so nominal conformal coverage is not guaranteed.
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
#'     \item{holdout}{The chosen holdout mode.}
#'     \item{fold_assignments}{For grouped modes, a data frame identifying the
#'       held-out species for each fold and replicate; NULL for cell mode.}
#'   }
#' @examples
#' \donttest{
#' data(avonet300, tree300, package = "pigauto")
#' tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(12L)])
#' traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
#'                      c("Mass", "Wing.Length"), drop = FALSE]
#' rownames(traits) <- tree$tip.label
#' cv <- cross_validate(preprocess_traits(traits, tree), tree, k = 3L,
#'                      seeds = 1L, epochs = 1L, verbose = FALSE,
#'                      k_eigen = 4L, hidden_dim = 8L,
#'                      use_transformer_blocks = FALSE)
#' print(cv)
#' summary(cv)
#' # Whole-species validation, with the full tree retained:
#' cv_species <- cross_validate(preprocess_traits(traits, tree), tree,
#'                              k = 3L, seeds = 1L, holdout = "species",
#'                              gnn = FALSE, verbose = FALSE)
#' cv_species$fold_assignments
#' }
#' @export
cross_validate <- function(data, tree, k = 5L, seeds = NULL,
                           epochs = 500L, verbose = TRUE,
                           holdout = c("cell", "species", "clade"), ...) {
  holdout <- match.arg(holdout)
  if (!inherits(data, "pigauto_data")) {
    stop("'data' must be a pigauto_data object (output of preprocess_traits).")
  }
  if (!inherits(tree, "phylo")) stop("'tree' must be a phylo object.")
  k <- as.integer(k)
  if (k < 2L) stop("'k' must be >= 2.")
  if (holdout != "cell" && k < 3L) stop("Grouped holdouts require k >= 3.", call. = FALSE)
  if (is.null(seeds)) {
    seeds <- list(NULL)
  } else {
    if (!is.numeric(seeds) || anyNA(seeds) || any(!is.finite(seeds))) {
      stop("`seeds` must be NULL or a finite numeric vector.", call. = FALSE)
    }
    seeds <- as.list(as.integer(seeds))
  }

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
  all_assignments <- list()

  for (rep_i in seq_along(seeds)) {
    seed <- seeds[[rep_i]]

    # Create k non-overlapping folds from observed cells
    folds <- if (holdout == "cell") make_cv_folds(X_scaled, k, seed, trait_map) else
      make_grouped_cv_folds(data, tree, k, seed, holdout)
    if (holdout != "cell") {
      species_folds <- attr(folds, "species_folds")
      for (i in seq_len(k)) {
        spp <- data$species_names[sort(species_folds[[i]])]
        all_assignments[[length(all_assignments) + 1L]] <- data.frame(
          rep = rep_i, fold = i, species = spp, stringsAsFactors = FALSE)
      }
    }

    for (fold_i in seq_len(k)) {
      if (verbose) {
        message(sprintf("Rep %d/%d, Fold %d/%d%s...",
                        rep_i, length(seeds), fold_i, k,
                        if (is.null(seed)) "" else paste0(" (seed=", seed, ")")))
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

      # Re-estimate trait scaling on training rows for grouped holdouts.
      fold_data <- if (holdout == "cell") data else prepare_grouped_cv_data(data, splits)
      fold_truth <- fold_data$X_scaled
      fold_trait_map <- fold_data$trait_map

      # Fit pigauto for this fold
      fit <- fit_pigauto(
        data    = fold_data,
        tree    = tree,
        splits  = splits,
        graph   = graph,
        epochs  = as.integer(epochs),
        seed    = seed,
        verbose = FALSE,
        ...
      )

      # Predict on test fold with validation/test cells hidden from DAE
      # context, matching the held-out input used during fit/calibration.
      pred <- predict(fit, return_se = TRUE,
                      .mask_observed_idx = c(splits$val_idx, splits$test_idx))

      # Evaluate using evaluate_imputation (test split)
      eval_df <- evaluate_imputation(
        pred   = pred,
        truth  = fold_truth,
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
          pred, fold_truth, splits$test_idx, fold_trait_map, fold_i, rep_i
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
      n_reps             = length(seeds),
      holdout            = holdout,
      fold_assignments   = if (length(all_assignments)) do.call(rbind, all_assignments) else NULL
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
  if (!is.null(seed)) set.seed(seed)
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
                   "spearman_rho", "accuracy", "brier", "zero_accuracy",
                   "aitchison", "rmse_clr", "simplex_mae")
  metric_cols <- intersect(metric_cols, names(eval_df))
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

    if (!(tm$type %in% c("continuous", "count", "ordinal", "proportion"))) next

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
    } else if (tm$type == "proportion") {
      truth_orig <- stats::plogis(truth_latent[ok] * tm$sd + tm$mean)
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

#' Print a cross-validation summary.
#'
#' Displays the number of folds and replicates, then reports each trait's
#' mean and standard deviation for every metric. Returns \code{x} invisibly.
#' Additional arguments are ignored by this print method.
#'
#' @param x A \code{pigauto_cv} object returned by \code{cross_validate()}.
#' @rdname cross_validate
#' @section Printing:
#' \code{print.pigauto_cv()} reports \code{k} folds and the replicate count,
#' followed by each trait's metric means and SDs. Displayed \code{n} is the
#' number of represented fold-by-replicate metric rows, not just the fold
#' count. This count includes rows with missing metric values; means and
#' SDs omit those values. When available, a separate conformal block shows
#' mean coverage and its SD against a nominal 95% target; the target is
#' not a guarantee. Empty summaries print \dQuote{No results.} The method
#' ignores \code{...} and returns \code{x} invisibly.
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

#' Print a compact summary of cross-validation results
#'
#' @rdname cross_validate
#' @aliases summary.pigauto_cv
#' @usage \method{summary}{pigauto_cv}(object, ...)
#' @param object A \code{pigauto_cv} object returned by
#'   \code{cross_validate()}.
#' @param ... Additional arguments forwarded by \code{cross_validate()} to \code{fit_pigauto()}; ignored by the print and summary methods.
#' @return Invisibly returns the summary data frame stored in
#'   \code{object$summary}.
#' @section Summary method:
#' \code{summary.pigauto_cv()} prints a compact per-trait table with RMSE,
#' Pearson \eqn{r} for continuous traits, Spearman \eqn{\rho} for ordinal
#' traits, and accuracy when returned (binary/categorical accuracy or
#' multi-proportion dominant-component accuracy). A dash means the metric
#' row is absent. A present but nonfinite metric may display \code{NA} or
#' \code{NaN}. An empty summary prints \dQuote{No results.} The full
#' \code{object$summary} data frame, which can contain metrics omitted from
#' this display, is returned invisibly. The method ignores \code{...}.
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

# Assign complete species or complete rooted subtrees to folds. Every latent
# column and every observation belonging to a selected species is masked.
make_grouped_cv_folds <- function(data, tree, k, seed, holdout) {
  if (!is.null(seed)) set.seed(seed)
  if (!identical(data$species_names, tree$tip.label)) {
    stop("Data species must match tree tips in the same order.", call. = FALSE)
  }
  n_tip <- length(tree$tip.label)
  if (k > n_tip) stop("k cannot exceed the number of species.", call. = FALSE)
  if (holdout == "species") {
    assignment <- sample(rep(seq_len(k), length.out = n_tip))
    species_folds <- lapply(seq_len(k), function(i) which(assignment == i))
  } else {
    if (!ape::is.rooted(tree)) stop("A rooted tree is required for clade holdouts.", call. = FALSE)
    edge <- tree$edge
    roots <- setdiff(unique(edge[, 1L]), edge[, 2L])
    if (length(roots) != 1L) stop("A rooted tree is required for clade holdouts.", call. = FALSE)
    children <- split(edge[, 2L], edge[, 1L])
    descendants <- function(node) {
      stack <- integer(nrow(edge) + 1L)
      stack[1L] <- node
      top <- 1L
      tips <- integer(n_tip)
      n_found <- 0L
      while (top > 0L) {
        current <- stack[top]
        top <- top - 1L
        if (current <= n_tip) {
          n_found <- n_found + 1L
          tips[n_found] <- current
        } else {
          ch <- children[[as.character(current)]]
          stack[top + seq_along(ch)] <- ch
          top <- top + length(ch)
        }
      }
      tips[seq_len(n_found)]
    }
    nodes <- roots
    while (length(nodes) < k) {
      possible <- which(vapply(nodes, function(node) {
        ch <- children[[as.character(node)]]
        length(ch) >= 2L && length(nodes) + length(ch) - 1L <= k
      }, logical(1)))
      if (!length(possible)) stop("The greedy clade partition did not find k rooted clades; choose another k or resolve polytomies.", call. = FALSE)
      sizes <- vapply(nodes[possible], function(node) length(descendants(node)), integer(1))
      choice <- possible[which.max(sizes)]
      nodes <- c(nodes[-choice], children[[as.character(nodes[choice])]])
    }
    species_folds <- lapply(nodes, descendants)
  }
  n <- nrow(data$X_scaled)
  p <- ncol(data$X_scaled)
  row_species <- if (isTRUE(data$multi_obs)) data$obs_to_species else seq_len(n)
  folds <- lapply(species_folds, function(species) {
    rows <- which(row_species %in% species)
    idx <- as.integer(as.vector(outer(rows, seq_len(p) - 1L, function(r, j) r + j * n)))
    idx <- idx[!is.na(data$X_scaled[idx])]
    if (!length(idx)) {
      stop("A grouped fold has no observed trait values; use fewer folds or a different seed.", call. = FALSE)
    }
    idx
  })
  attr(folds, "species_folds") <- species_folds
  folds
}

# Holdout rows are hidden when estimating trait-specific transformations.
# The original trait types, factor levels, and log-transform choices are
# treated as an input schema; the numerical scaling is learned per fold.
prepare_grouped_cv_data <- function(data, splits) {
  if (is.null(data$X_original) || is.null(data$preprocess_center) || is.null(data$preprocess_scale)) {
    stop("Grouped cross-validation requires a freshly preprocessed data object; rerun preprocess_traits().", call. = FALSE)
  }
  n <- nrow(data$X_scaled)
  held_rows <- unique(((c(splits$val_idx, splits$test_idx) - 1L) %% n) + 1L)
  training <- data$X_original
  training[held_rows, ] <- NA
  tmap <- data$trait_map
  for (tm in tmap) {
    cols <- if (tm$type == "multi_proportion") tm$input_cols else tm$name
    observed <- stats::complete.cases(training[, cols, drop = FALSE])
    if (sum(observed) < 2L) {
      stop("Trait ", sQuote(tm$name), " has fewer than two training observations in a grouped fold.", call. = FALSE)
    }
    row_species <- if (isTRUE(data$multi_obs)) data$obs_to_species else seq_len(n)
    if (length(unique(row_species[observed])) < 2L) {
      stop("Trait ", sQuote(tm$name), " has fewer than two training species in a grouped fold.", call. = FALSE)
    }
  }
  standalone <- tmap[vapply(tmap, function(x) x$type != "multi_proportion", logical(1))]
  types <- vapply(standalone, `[[`, character(1), "type")
  groups <- lapply(tmap[vapply(tmap, function(x) x$type == "multi_proportion", logical(1))], `[[`, "input_cols")
  logs <- vapply(standalone, function(x) isTRUE(x$log_transform), logical(1))
  new_map <- build_trait_map(training, types, names(logs)[logs],
                             data$preprocess_center, data$preprocess_scale,
                             multi_proportion_groups = groups)
  if (!identical(vapply(tmap, `[[`, integer(1), "n_latent"),
                 vapply(new_map, `[[`, integer(1), "n_latent"))) {
    stop("Grouped fold changed the trait schema.", call. = FALSE)
  }
  data$trait_map <- new_map
  data$X_scaled <- encode_to_latent(data$X_original, new_map)$X
  data$means <- vapply(new_map, function(x) if (length(x$mean) == 1L) x$mean else NA_real_, numeric(1))
  data$sds <- vapply(new_map, function(x) if (length(x$sd) == 1L) x$sd else NA_real_, numeric(1))
  data
}
