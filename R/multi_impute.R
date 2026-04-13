#' Generate M complete datasets for multiple imputation
#'
#' Run pigauto's full imputation pipeline and return `M` stochastic
#' completions of the trait matrix instead of a single point estimate.
#' The `M` datasets are the input needed for the classical multiple
#' imputation workflow: fit a downstream model on each dataset, then
#' pool the results with Rubin's rules via [pool_mi()]. This is the
#' standard way to propagate imputation uncertainty into phylogenetic
#' comparative analyses (PGLS, PGLMM, etc.) rather than treating
#' imputed cells as if they were observed.
#'
#' @param traits data.frame with species as rownames and trait columns.
#'   Same input format as [impute()]. Supported column types are
#'   numeric, integer, factor, ordered factor, and logical.
#' @param tree object of class `phylo` aligned with `traits`.
#' @param m integer. Number of imputation datasets to generate
#'   (default `100`). Observed cells are identical across all `M`
#'   datasets; only originally-missing cells vary.
#' @param draws_method character. How stochastic draws are generated for
#'   missing cells. One of:
#'   \describe{
#'     \item{`"conformal"`}{(default) Run the model once, then sample each
#'       originally-missing cell from a Normal distribution centred on the
#'       point estimate with SD = conformal_score / 1.96. The conformal score
#'       is the empirical 97.5th percentile of held-out absolute residuals,
#'       so the draw width is calibrated against actual prediction error —
#'       not a model assumption. Falls back to BM-SE-based Normal sampling
#'       when conformal scores are unavailable, and to Bernoulli / Categorical
#'       draws for discrete traits. **Preferred default for pigauto** because
#'       MC dropout gives zero variance whenever the calibrated gate is zero
#'       (i.e. whenever the BM baseline already fits well), which is common
#'       for continuous traits with strong phylogenetic signal.}
#'     \item{`"mc_dropout"`}{Run `M` stochastic GNN forward passes in training
#'       mode (dropout active). Useful when the calibrated gate is open
#'       (r_cal > 0, i.e. GNN meaningfully corrects the BM baseline). When
#'       all gates are zero — as is typical for continuous traits on datasets
#'       with strong phylogenetic signal — MC dropout is deterministic and
#'       falls back silently to the BM-only point estimate for every draw.
#'       Check `mi$fit$calibrated_gates` before using this method.}
#'   }
#' @param species_col character or `NULL`. If set, marks the column
#'   in `traits` containing species identifiers and enables multiple
#'   observations per species. See [impute()] for details.
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
#'   [impute()].
#'
#' @return An object of class `"pigauto_mi"` with components:
#'   \describe{
#'     \item{`datasets`}{A list of length `m`. Each element is a
#'       data.frame with the same shape and column types as the input
#'       `traits`; observed cells are preserved and missing cells are
#'       filled with the corresponding imputation draw. Pass this list
#'       to [with_imputations()] to fit downstream models.}
#'     \item{`m`}{Number of imputations.}
#'     \item{`pooled_point`}{A single data.frame whose missing cells
#'       are replaced by the MC-averaged point estimate. Convenient for
#'       reporting but does *not* propagate imputation uncertainty --
#'       use `datasets` + [pool_mi()] for inference.}
#'     \item{`se`}{Matrix of per-cell standard errors combining the
#'       baseline SE and the between-imputation standard deviation.}
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
#' Multiple imputation is a method for doing *downstream analysis*
#' under missing data, not an end in itself. Plugging a single
#' point-estimate imputation into a regression underestimates standard
#' errors because it treats imputed cells as if they were observed.
#' The standard remedy, due to Rubin (1987), is to generate `M`
#' stochastic completions, fit the downstream model on each, and pool
#' the results. `multi_impute()` + [with_imputations()] + [pool_mi()]
#' implement this workflow end to end.
#'
#' **`draws_method = "conformal"` (default)**: Run the model once; missing
#' cells are sampled from
#' \eqn{x_{ij}^{(k)} \sim \mathrm{N}(\hat\mu_{ij},\; q_{j}/1.96)}
#' where \eqn{q_j} is the trait-level conformal score (the empirical
#' 97.5th percentile of held-out absolute residuals, in latent z-score
#' units back-transformed to the original scale). The draw width is
#' therefore calibrated against actual prediction error regardless of
#' whether the BM or GNN term dominates. For discrete traits (binary,
#' categorical) it uses Bernoulli / categorical draws from the estimated
#' probability vector. This is the preferred default for pigauto.
#'
#' **`draws_method = "mc_dropout"`**: Run `M` GNN forward passes in
#' training mode (dropout active). **Caution**: when the per-trait
#' calibrated gate `r_cal = 0` (which happens whenever the BM baseline
#' already fits well, typically for continuous traits with strong
#' phylogenetic signal), every MC pass is identical to the BM point
#' estimate and draws have zero between-imputation variance. Check
#' `mi$fit$calibrated_gates` after fitting — if all gates for the traits
#' of interest are zero, use `draws_method = "conformal"` instead.
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
#' @seealso [impute()] for single-point imputation, [with_imputations()]
#'   for applying a model-fitting function across the `M` datasets,
#'   [pool_mi()] for Rubin's rules pooling of the resulting fits.
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
#' # Downstream analysis: phylogenetic GLS via nlme, pooled with Rubin's rules
#' fits <- with_imputations(mi, function(d) {
#'   d$species <- rownames(d)
#'   nlme::gls(
#'     log(Mass) ~ log(Wing.Length),
#'     correlation = ape::corBrownian(phy = tree300, form = ~species),
#'     data = d, method = "ML"
#'   )
#' })
#' pool_mi(fits)
#' }
#'
#' @export
multi_impute <- function(traits, tree, m = 100L,
                         draws_method = c("conformal", "mc_dropout"),
                         species_col = NULL,
                         log_transform = TRUE,
                         missing_frac = 0.25,
                         covariates = NULL,
                         epochs = 2000L, verbose = TRUE, seed = 1L, ...) {

  draws_method <- match.arg(draws_method)
  m <- as.integer(m)
  if (!is.finite(m) || m < 2L) {
    stop("`m` must be an integer >= 2 (multiple imputation needs at least ",
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
    datasets <- lapply(pred$imputed_datasets, function(imp_df) {
      build_completed(traits, imp_df, species_col)$completed
    })

  } else {
    # ---- Conformal: single pass + conformal-width Normal sampling -----------
    # Run once to get point estimates, conformal scores, and probabilities.
    res <- impute(
      traits        = traits,
      tree          = tree,
      species_col   = species_col,
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

    datasets <- lapply(seq_len(m), function(i) {
      imp_df <- .sample_conformal_draw(pred, imask, trait_map,
                                       seed_i = as.integer(seed) + i)
      build_completed(traits, imp_df, species_col)$completed
    })
  }

  structure(
    list(
      datasets     = datasets,
      m            = m,
      draws_method = draws_method,
      pooled_point = res$completed,
      se           = pred$se,
      imputed_mask = res$imputed_mask,
      fit          = res$fit,
      data         = res$data,
      tree         = tree,
      species_col  = species_col,
      evaluation   = res$evaluation
    ),
    class = c("pigauto_mi", "list")
  )
}


# ---- Internal: one conformal-width draw -------------------------------------
# Samples missing cells from N(mu, conformal_score / 1.96) for continuous
# types (on the appropriate transformed scale), and from Bernoulli /
# Categorical for discrete types. Falls back to BM SE when conformal score
# is not available for a trait.
.sample_conformal_draw <- function(pred, imputed_mask, trait_map, seed_i) {
  set.seed(seed_i)
  imp    <- pred$imputed
  probs  <- pred$probabilities
  cscores <- pred$conformal_scores  # named vector, NA for discrete traits

  for (tm in trait_map) {
    nm   <- tm$name
    if (!(nm %in% names(imp))) next
    if (!(nm %in% colnames(imputed_mask))) next
    rows <- which(imputed_mask[, nm])
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

  cat("pigauto multiple imputation\n")
  method_label <- switch(x$draws_method %||% "mc_dropout",
    mc_dropout = "MC dropout",
    conformal  = "conformal-width sampling",
    x$draws_method
  )
  cat(sprintf("  M        : %d imputations (%s)\n", x$m, method_label))
  cat(sprintf("  Species  : %d\n", n_sp))
  cat(sprintf("  Traits   : %d -- %s\n", p,
              paste(traits, collapse = ", ")))
  cat(sprintf("  Cells    : %d imputed / %d total (%.1f%%)\n",
              n_imp_cells, total_cells, pct))

  cat("\n  Access imputation draws:  mi$datasets[[i]]\n")
  cat("  Fit downstream models:    with_imputations(mi, fit_fun)\n")
  cat("  Pool with Rubin's rules:  pool_mi(fits)\n")
  invisible(x)
}
