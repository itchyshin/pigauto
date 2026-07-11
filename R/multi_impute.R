#' Generate experimental stochastic completion datasets
#'
#' Run pigauto's full imputation pipeline and return `M` stochastic
#' completions of the trait matrix instead of a single point estimate.
#' The conformal-width and Brownian/MC-dropout draws returned here are
#' experimental prediction-diagnostic draws. A 3,000-fit known-DGP campaign
#' found that neither method passed any of the 12 downstream fixed-effect
#' gate cells. Do not use these datasets for downstream inference or Rubin
#' pooling. A separate analysis-aware backend, [multi_impute_analysis()], has
#' passed its package-level fixed-effect gate for a narrow set of analyses.
#'
#' @section When to use this:
#'
#' This function is useful for comparing stochastic prediction behavior from
#' one tree. It is not the analysis-aware inferential backend.
#'
#' [multi_impute_trees()] provides an experimental posterior-tree sensitivity
#' path, but tree uncertainty is not supported by [multi_impute_analysis()].
#'
#' @param traits data.frame with species as rownames and trait columns.
#'   Same input format as [impute()]. Supported column types are
#'   numeric, integer, factor, ordered factor, and logical.
#' @param tree object of class `phylo` aligned with `traits`.
#' @param m integer. Number of stochastic completion datasets to generate
#'   (default `100`). Observed cells are identical across all `M`
#'   datasets; only originally-missing cells vary.
#' @param draws_method character. How stochastic draws are generated for
#'   missing cells. One of:
#'   \describe{
#'     \item{`"conformal"`}{(default) Run the model once, then sample each
#'       originally-missing cell from a Normal distribution centred on the
#'       point estimate with SD = conformal_score / 1.96. Converting a
#'       split-conformal residual quantile to a Normal scale is a heuristic;
#'       the conformal coverage guarantee does not establish that these draws
#'       are proper multiple imputations. Falls back to BM-SE-based Normal sampling
#'       when conformal scores are unavailable, and to Bernoulli / Categorical
#'       draws for discrete traits.}
#'     \item{`"mc_dropout"`}{Run `M` stochastic GNN forward passes in training
#'       mode (dropout active) on top of stochastic Brownian-motion baseline
#'       draws. Brownian draws still contribute between-draw variation
#'       when a calibrated GNN gate is zero.}
#'   }
#' @param species_col character or `NULL`. If set, marks the column
#'   in `traits` containing species identifiers and enables multiple
#'   observations per species. See [impute()] for details.
#' @param trait_types named character vector overriding auto-detected
#'   trait types for specific columns. Required for `"proportion"` and
#'   `"zi_count"`. See [impute()] and [preprocess_traits()]. Default
#'   `NULL` (auto-detect).
#' @param multi_proportion_groups named list declaring compositional
#'   trait groups (rows summing to 1), e.g.
#'   `list(diet = c("plant", "invert", "vert"))`. Forwarded to
#'   [impute()] / [preprocess_traits()]. Default `NULL`.
#' @param log_transform logical. Auto-log positive continuous columns
#'   (default `TRUE`).
#' @param missing_frac numeric. Fraction of observed cells held out for
#'   validation/test during training (default `0.25`). Passed through
#'   to [impute()].
#' @param covariates data.frame or matrix of environmental covariates
#'   (fully observed, numeric). Passed through to [impute()].
#'   Default `NULL` (no covariates).
#' @param epochs integer. Maximum GNN training epochs (default `2000`).
#' @param verbose logical. Print progress (default `TRUE`).
#' @param seed integer. Random seed (default `1`).
#' @param ... additional arguments forwarded to [fit_pigauto()] via
#'   [impute()]. See [fit_pigauto()] for the full list; the "Safety
#'   floor" section below describes the relevant new v0.9.1.9002
#'   argument.
#'
#' @return An object of class `"pigauto_mi"` with components:
#'   \describe{
#'     \item{`datasets`}{A list of length `m`. Each element is a
#'       data.frame with the same shape and column types as the input
#'       `traits`; observed cells are preserved and missing cells are
#'       filled with the corresponding stochastic draw. These datasets are
#'       for prediction diagnostics, not downstream inference.}
#'     \item{`m`}{Number of stochastic completion datasets.}
#'     \item{`pooled_point`}{A single data.frame whose missing cells
#'       are replaced by the MC-averaged point estimate. Convenient for
#'       reporting but does *not* provide a valid downstream MI analysis.}
#'     \item{`se`}{Matrix of per-cell uncertainty summaries combining the
#'       baseline SE and the between-draw standard deviation.}
#'     \item{`imputed_mask`}{Logical matrix; `TRUE` where a cell was
#'       originally missing.}
#'     \item{`fit`}{The underlying [`pigauto_fit`][fit_pigauto()]
#'       object, retained for diagnostics and for calls to [predict()]
#'       on new data.}
#'     \item{`data`}{The [`pigauto_data`][preprocess_traits()] object.}
#'     \item{`tree`}{The input phylogeny.}
#'     \item{`species_col`}{Passed-through species-column name or
#'       `NULL`.}
#'   }
#'
#' @details
#' These draws do not condition on a declared substantive analysis model.
#' Consequently, stochastic variation alone does not make them proper or
#' congenial multiple imputations. The analysis-aware backend requires the
#' analysis model before generating draws and dispatches only across its
#' documented supported model classes.
#'
#' **`draws_method = "conformal"` (default)**: Run the model once; missing
#' cells are sampled from
#' \eqn{x_{ij}^{(k)} \sim \mathrm{N}(\hat\mu_{ij},\; q_{j}/1.96)}
#' where \eqn{q_j} is the trait-level split-conformal residual quantile.
#' Dividing this quantile by 1.96 is a pragmatic Normal-scale construction,
#' not a consequence of the conformal coverage guarantee. For discrete traits (binary,
#' categorical) it uses Bernoulli / categorical draws from the estimated
#' probability vector. For \code{multi_proportion} groups it draws the
#' K CLR latent columns with their BM latent SEs, projects back to
#' sum-zero CLR space, and decodes to the simplex.
#'
#' **`draws_method = "mc_dropout"`**: Run `M` GNN forward passes in
#' training mode (dropout active) on top of stochastic BM baseline draws.
#' When `r_cal = 0`, the GNN-dropout term disappears but the BM draw still
#' contributes between-draw variance.
#'
#' Nakagawa & Freckleton (2008, 2011) review the consequences of
#' ignoring missing data in ecological and comparative analyses and
#' argue for multiple imputation as the default.
#'
#' @references
#' Rubin DB (1987). *Multiple Imputation for Nonresponse in Surveys.*
#' Wiley.
#'
#' Nakagawa S, Freckleton RP (2008). "Missing inaction: the dangers of
#' ignoring missing data." *Trends in Ecology & Evolution* 23(11):
#' 592-596.
#'
#' Nakagawa S, Freckleton RP (2011). "Model averaging, missing data and
#' multiple imputation: a case study for behavioural ecology."
#' *Behavioral Ecology and Sociobiology* 65(1): 103-116.
#'
#' @section Safety floor (v0.9.1.9002+):
#'   When \code{fit_pigauto()} was called with \code{safety_floor = TRUE}
#'   (the default since v0.9.1.9002), the 3-way blend
#'   \code{r_BM * BM + r_GNN * GNN + r_MEAN * MEAN} propagates through
#'   every imputation draw automatically via the updated
#'   \code{predict.pigauto_fit()}.  For \code{draws_method = "mc_dropout"}
#'   the mean term contributes no between-draw variance (it is a
#'   deterministic scalar per column); between-draw variance comes from the
#'   BM-draw and GNN-dropout terms
#'   only.  For \code{draws_method = "conformal"} the blend centre is the
#'   3-way prediction and conformal scores remain calibrated on the
#'   blended residuals.
#'
#' @seealso [impute()] for point imputation and [multi_impute_analysis()]
#'   for the narrow analysis-aware inferential backend.
#'
#' @examples
#' \dontrun{
#' library(pigauto)
#' data(avonet300, tree300)
#' df <- avonet300; rownames(df) <- df$Species_Key; df$Species_Key <- NULL
#'
#' # Generate 100 complete datasets
#' mi <- multi_impute(df, tree300, m = 100)
#' print(mi)
#'
#' # Inspect stochastic prediction sensitivity only.
#' lapply(mi$datasets, head)
#' }
#'
#' @export
multi_impute <- function(traits, tree, m = 100L,
                         draws_method = c("conformal", "mc_dropout"),
                         species_col = NULL,
                         trait_types = NULL,
                         multi_proportion_groups = NULL,
                         log_transform = TRUE,
                         missing_frac = 0.25,
                         covariates = NULL,
                         epochs = 2000L, verbose = TRUE, seed = 1L, ...) {

  draws_method <- match.arg(draws_method)
  m <- as.integer(m)
  if (!is.finite(m) || m < 2L) {
    stop("`m` must be an integer >= 2 (stochastic diagnostics need at least ",
         "two draws). Got m = ", m, ".", call. = FALSE)
  }

  if (draws_method == "mc_dropout") {
    # ---- MC dropout: M stochastic GNN forward passes (training mode) --------
    # Run the pipeline with n_imputations = m so predict() runs the model M
    # times with dropout active. Each pass yields a different latent matrix;
    # the M decoded data.frames are returned in pred$imputed_datasets.
    res <- impute(
      traits        = traits,
      tree          = tree,
      species_col   = species_col,
      trait_types   = trait_types,
      multi_proportion_groups = multi_proportion_groups,
      log_transform = log_transform,
      missing_frac  = missing_frac,
      n_imputations = m,
      covariates    = covariates,
      epochs        = as.integer(epochs),
      verbose       = verbose,
      seed          = as.integer(seed),
      ...
    )

    pred <- res$prediction
    if (is.null(pred$imputed_datasets) || length(pred$imputed_datasets) != m) {
      stop("predict.pigauto_fit() did not return ", m,
           " imputed datasets. This is an internal error -- please report.",
           call. = FALSE)
    }
    # Pass `input_row_order` so build_completed re-aligns predictions back
    # to the user's input row order (multi-obs reordering fix, 2026-04-26).
    input_row_order <- res$data$input_row_order
    datasets <- lapply(pred$imputed_datasets, function(imp_df) {
      build_completed(traits, imp_df, species_col,
                       input_row_order = input_row_order)$completed
    })

  } else {
    # ---- Conformal: single pass + conformal-width Normal sampling -----------
    # Run once to get point estimates, conformal scores, and probabilities.
    res <- impute(
      traits        = traits,
      tree          = tree,
      species_col   = species_col,
      trait_types   = trait_types,
      multi_proportion_groups = multi_proportion_groups,
      log_transform = log_transform,
      missing_frac  = missing_frac,
      n_imputations = 1L,
      covariates    = covariates,
      epochs        = as.integer(epochs),
      verbose       = verbose,
      seed          = as.integer(seed),
      ...
    )

    pred      <- res$prediction
    trait_map <- res$fit$trait_map
    imask     <- res$imputed_mask

    input_row_order <- res$data$input_row_order
    datasets <- lapply(seq_len(m), function(i) {
      # `imask` is in user-input row order (built by build_completed from
      # the original traits data.frame).  But pred$imputed / pred$se /
      # pred$probabilities are in INTERNAL (tree-tip) order from the
      # model.  .sample_conformal_draw must therefore map input-order
      # mask indices to internal-order before perturbing imp; otherwise
      # it perturbs the wrong cells and build_completed re-aligns the
      # unperturbed (internal) values back into the user's missing cells,
      # producing zero between-imputation variance.  See
      # adversarial_review_opus.md (HIGH severity finding C.1).
      imp_df <- .sample_conformal_draw(pred, imask, trait_map,
                                       seed_i = as.integer(seed) + i,
                                       input_row_order = input_row_order)
      build_completed(traits, imp_df, species_col,
                       input_row_order = input_row_order)$completed
    })
  }

  structure(
    list(
      datasets        = datasets,
      m               = m,
      draws_method    = draws_method,
      pooled_point    = res$completed,
      se              = pred$se,
      # Exposed for downstream coverage / interval-width / Brier analyses
      # without having to round-trip through `mi$fit` + a fresh predict().
      # `conformal_lower` / `conformal_upper` are the trait-level 95%
      # prediction intervals (back-transformed to user scale when the
      # trait was log-transformed); `probabilities` is the per-discrete-
      # trait probability list (binary: numeric vector, categorical: K-
      # column matrix). See predict.pigauto_fit() for shape and scale
      # conventions. Either may be NULL when the underlying conformal
      # scores were unavailable (e.g. zero validation cells for a trait).
      conformal_lower = pred$conformal_lower,
      conformal_upper = pred$conformal_upper,
      probabilities   = pred$probabilities,
      imputed_mask    = res$imputed_mask,
      fit             = res$fit,
      data            = res$data,
      tree            = tree,
      species_col     = species_col,
      evaluation      = res$evaluation
    ),
    class = c("pigauto_mi", "list")
  )
}


# ---- Internal: one conformal-width draw -------------------------------------
# Samples missing cells from N(mu, conformal_score / 1.96) for continuous
# types (on the appropriate transformed scale), and from Bernoulli /
# Categorical for discrete types. Falls back to BM SE when conformal score
# is not available for a trait.
.sample_conformal_draw <- function(pred, imputed_mask, trait_map, seed_i,
                                    input_row_order = NULL) {
  set.seed(seed_i)
  imp    <- pred$imputed
  probs  <- pred$probabilities
  cscores <- pred$conformal_scores  # named vector, NA for discrete traits

  # Build input-row -> internal-row index map.  `imputed_mask` rows are in
  # user-input order; `pred$imputed` / `pred$se` / `pred$probabilities` are
  # in internal (tree-tip) order.  `input_row_order[k] = i` means internal
  # row k holds original input row i, so the inverse map is
  # `match(seq_along(input_row_order), input_row_order)`.  When NULL or
  # identity, internal order == input order and the conversion is a no-op.
  input_to_internal <- if (!is.null(input_row_order)) {
    match(seq_along(input_row_order), input_row_order)
  } else NULL

  to_internal <- function(rows_input) {
    if (is.null(input_to_internal)) return(rows_input)
    out <- input_to_internal[rows_input]
    out[!is.na(out)]    # drop input rows with no internal counterpart
  }

  for (tm in trait_map) {
    nm   <- tm$name
    if (tm$type == "multi_proportion") {
      comp_cols <- if (!is.null(tm$input_cols)) tm$input_cols else tm$levels
      if (!all(comp_cols %in% names(imp))) next
      if (!all(comp_cols %in% colnames(imputed_mask))) next

      missing_rows <- rowSums(imputed_mask[, comp_cols, drop = FALSE]) > 0
      rows <- to_internal(which(missing_rows))
      if (length(rows) == 0L) next

      lc <- tm$latent_cols
      latent_mu <- pred$imputed_latent[rows, lc, drop = FALSE]
      s_latent <- pred$se_latent[rows, lc, drop = FALSE]
      s_latent[!is.finite(s_latent)] <- 0

      latent_draw <- latent_mu +
        matrix(stats::rnorm(length(rows) * length(lc)), nrow = length(rows)) *
          s_latent

      clr_mat <- latent_draw
      for (k in seq_along(lc)) {
        clr_mat[, k] <- clr_mat[, k] * tm$sd[k] + tm$mean[k]
      }
      clr_mat <- clr_mat - rowMeans(clr_mat)
      prop_mat <- softmax_rows(clr_mat)

      for (k in seq_along(comp_cols)) {
        imp[[comp_cols[k]]][rows] <- prop_mat[, k]
      }
      next
    }

    if (!(nm %in% names(imp))) next
    if (!(nm %in% colnames(imputed_mask))) next
    # Map input-order mask indices to internal-order indices for indexing
    # into pred$imputed / pred$se / pred$probabilities.
    rows <- to_internal(which(imputed_mask[, nm]))
    if (length(rows) == 0L) next
    N    <- length(rows)

    # Conformal half-width → approximate 1-sigma SD
    cs <- if (!is.null(cscores) && nm %in% names(cscores) &&
               is.finite(cscores[nm]))
            cscores[nm] / 1.96
          else
            NULL

    if (tm$type == "continuous") {
      mu       <- imp[[nm]][rows]           # original-scale point estimate
      # s_latent = SE in z-score (latent) space
      s_latent <- if (!is.null(cs)) cs else pred$se[rows, nm] / tm$sd
      if (isTRUE(tm$log_transform)) {
        # Draw on the z-score scale then back-transform to avoid the
        # delta-method approximation error (s_orig/mu is tiny when mu >> 1).
        z_mu  <- (log(pmax(mu, .Machine$double.eps)) - tm$mean) / tm$sd
        z_drw <- rnorm(N, z_mu, s_latent)
        imp[[nm]][rows] <- exp(z_drw * tm$sd + tm$mean)
      } else {
        # Non-log: draw directly in original scale (s_orig = s_latent * sd)
        imp[[nm]][rows] <- rnorm(N, mu, s_latent * tm$sd)
      }

    } else if (tm$type == "count") {
      # log1p transform: draw on z-score scale, back via expm1
      mu       <- as.numeric(imp[[nm]][rows])
      s_latent <- if (!is.null(cs)) cs else pred$se[rows, nm] / tm$sd
      z_mu  <- (log1p(pmax(mu, 0)) - tm$mean) / tm$sd
      z_drw <- rnorm(N, z_mu, s_latent)
      draw  <- pmax(round(expm1(z_drw * tm$sd + tm$mean)), 0L)
      imp[[nm]][rows] <- as.integer(draw)

    } else if (tm$type == "ordinal") {
      K        <- length(tm$levels)
      int_mu   <- as.integer(imp[[nm]][rows]) - 1L
      s_latent <- if (!is.null(cs)) cs else pred$se[rows, nm] / tm$sd
      s_orig   <- s_latent * tm$sd   # ordinal is integer-scale, no transform
      draw_i   <- pmin(pmax(round(rnorm(N, int_mu, s_orig)), 0L), K - 1L)
      imp[[nm]][rows] <- factor(tm$levels[as.integer(draw_i) + 1L],
                                levels = tm$levels, ordered = TRUE)

    } else if (tm$type == "proportion") {
      # logit transform: draw on z-score scale, back via plogis
      mu       <- imp[[nm]][rows]
      s_latent <- if (!is.null(cs)) cs else pred$se[rows, nm] / tm$sd
      p_mu     <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
      logit_mu <- stats::qlogis(p_mu)
      z_mu     <- (logit_mu - tm$mean) / tm$sd
      z_drw    <- rnorm(N, z_mu, s_latent)
      imp[[nm]][rows] <- stats::plogis(z_drw * tm$sd + tm$mean)

    } else if (tm$type == "binary") {
      p   <- probs[[nm]][rows]
      idx <- rbinom(N, 1L, pmin(pmax(p, 0), 1)) + 1L
      imp[[nm]][rows] <- factor(tm$levels[idx], levels = tm$levels)

    } else if (tm$type == "categorical") {
      pm  <- probs[[nm]][rows, , drop = FALSE]
      idx <- apply(pm, 1L, function(p) {
        p <- pmax(p, 0); p <- p / sum(p)
        sample.int(length(p), 1L, prob = p)
      })
      imp[[nm]][rows] <- factor(tm$levels[idx], levels = tm$levels)

    } else if (tm$type == "zi_count") {
      p_nz    <- pmin(pmax(probs[[nm]][rows], 0), 1)
      gate    <- rbinom(N, 1L, p_nz)
      mu      <- as.numeric(imp[[nm]][rows])
      s       <- pred$se[rows, nm]
      cond_mu <- ifelse(p_nz > 0.01, mu / p_nz, mu)
      draw_c  <- as.integer(pmax(round(rnorm(N, cond_mu, s)), 0L))
      draw_c[gate == 0L] <- 0L
      imp[[nm]][rows] <- draw_c
    }
  }
  imp
}


#' @export
print.pigauto_mi <- function(x, ...) {
  n_sp  <- length(x$data$species_names)
  traits <- vapply(x$data$trait_map, "[[", character(1), "name")
  p <- length(traits)

  total_cells <- length(x$imputed_mask)
  n_imp_cells <- sum(x$imputed_mask)
  pct <- if (total_cells > 0) 100 * n_imp_cells / total_cells else 0

  cat("pigauto experimental stochastic completion diagnostics\n")
  method_label <- switch(x$draws_method %||% "mc_dropout",
    mc_dropout = "MC dropout",
    conformal  = "conformal-width sampling",
    x$draws_method
  )
  cat(sprintf("  M        : %d completion draws (%s)\n", x$m, method_label))
  cat(sprintf("  Species  : %d\n", n_sp))
  cat(sprintf("  Traits   : %d -- %s\n", p,
              paste(traits, collapse = ", ")))
  cat(sprintf("  Cells    : %d imputed / %d total (%.1f%%)\n",
              n_imp_cells, total_cells, pct))

  cat("\n  Access diagnostic draws:  mi$datasets[[i]]\n")
  cat("  Downstream inference:     unsupported for these draws\n")
  cat("  Analysis-aware MI:        multi_impute_analysis(...)\n")
  invisible(x)
}
