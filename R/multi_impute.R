#' Generate experimental stochastic completion datasets
#'
#' Run pigauto's full imputation pipeline and return `M` stochastic
#' completions of the trait matrix instead of a single point estimate.
#' The conformal-width and Brownian/MC-dropout draws returned here are
#' experimental prediction-diagnostic draws. Do not use these datasets for
#' downstream inference or Rubin pooling. For continuous traits,
#' `draws_method = "posterior"` instead returns proper Bayesian posterior
#' imputations that can be passed to [with_imputations()] and [pool_mi()]
#' for the analyses described in "Posterior draws" below. The separate
#' analysis-aware backend, [multi_impute_analysis()], covers its own
#' documented narrow regime.
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
#'       a nominal held-out conformal diagnostic does not establish that these
#'       draws are proper multiple imputations. Falls back to BM-SE-based Normal sampling
#'       when conformal scores are unavailable, and to Bernoulli / Categorical
#'       draws for discrete traits.}
#'     \item{`"mc_dropout"`}{Run `M` stochastic GNN forward passes in training
#'       mode (dropout active) on top of stochastic Brownian-motion baseline
#'       draws. Brownian draws still contribute between-draw variation
#'       when a calibrated GNN gate is zero.}
#'     \item{`"posterior"`}{Continuous traits only, one row per species, no
#'       covariates. Fit a Bayesian multivariate phylogenetic mixed model by
#'       MCMC and return `m` completions drawn from the posterior predictive
#'       distribution of the missing cells. No GNN is fitted. See "Posterior
#'       draws" below.}
#'   }
#'   Conformal and MC-dropout draws perturb missing cells around a point
#'   prediction rather than drawing them jointly from their conditional
#'   distribution given the observed data. In a 16-regime simulation of a
#'   downstream phylogenetic regression (PGLS slope, two traits, 120
#'   replicates per regime), conformal draws biased the pooled slope by
#'   -0.20 to -0.46, and the pooled 95% intervals covered the truth in 0 to
#'   17% of replicates, in all 16 regimes; MC-dropout draws biased it by
#'   -0.03 to -0.38. For downstream inference on continuous traits use
#'   `"posterior"`.
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
#' @param seed optional integer. When supplied, makes fitting and imputation
#'   draws reproducible; the default `NULL` uses the current RNG stream.
#' @param gnn logical. Passed through to [impute()] / [fit_pigauto()]. When
#'   `TRUE` (default), the usual GNN correction is trained. When `FALSE`,
#'   no GNN is used (baseline-only fit; see [fit_pigauto()]). With
#'   `draws_method = "mc_dropout"` this degrades gracefully to BM-posterior
#'   draws (there is no dropout to run) and a one-time message is printed;
#'   `draws_method = "conformal"` is unaffected. User `covariates` are
#'   ignored under `gnn = FALSE` (with a warning), as in [impute()].
#' @param lambda_mode character. Pagel-lambda mode for the BM baseline,
#'   forwarded to [impute()] / [fit_pigauto()]. `"estimate"` (default) fits
#'   a per-trait Pagel's lambda on each continuous-family (BM-eligible)
#'   latent column; discrete traits stay at lambda = 1. `"fixed_1"`
#'   preserves the pre-lambda Brownian correlation matrix everywhere;
#'   `"cv"` and `"bayes"` are alternative per-column estimators. See
#'   [fit_pigauto()] for the full contract, including the
#'   `predict_method = "exact"` / `joint_refine_iter > 0` interaction with
#'   `lambda_block`.
#' @param posterior_control list of settings for `draws_method =
#'   "posterior"` (ignored otherwise). Elements, with defaults:
#'   \describe{
#'     \item{`n_chains`}{`4L`. Number of MCMC chains. Chain 1 starts at the
#'       per-trait REML Pagel's lambda (up to 2,000 tips; above that, or if
#'       the REML fit fails, at lambda = 0.9), and the other chains start
#'       from values dispersed around it. The values used are returned in
#'       `mi$posterior$start`.}
#'     \item{`n_iter`}{`5000L`. Sweeps per chain after burn-in.}
#'     \item{`burnin`}{`1000L`. Burn-in sweeps per chain (also used to tune
#'       the Metropolis step sizes).}
#'     \item{`thin`}{Default: `n_iter` integer-divided by
#'       `ceiling(keep_draws / n_chains)`, so
#'       that the chains together keep `keep_draws` sweeps.}
#'     \item{`keep_draws`}{`1000L`. Target number of kept posterior
#'       predictive draws across chains, used for the per-cell intervals.
#'       Must be at least `m`. The default `thin` keeps at least this many
#'       whenever `n_iter` is at least `ceiling(keep_draws / n_chains)`. If
#'       a user-set `n_iter` or `thin` keeps fewer, a warning is raised and
#'       the intervals use the sweeps that were kept.}
#'     \item{`param_uncertainty`}{`"full"` (default), `"none"` or `"both"`.
#'       `"none"` is an improper plug-in mode for validation only: the
#'       covariance matrices are fixed at their posterior means from a full
#'       run, and only the missing cells (and the trait means) are drawn.
#'       Its result carries `mi_workflow =
#'       "pigauto_posterior_plugin_diagnostic"`, and [with_imputations()]
#'       and [pool_mi()] refuse it.
#'       `"both"` is also for validation only: it runs the sampler once and
#'       returns the proper results exactly as `"full"` does, plus the
#'       plug-in draws from the same run in `mi$posterior_improper`.}
#'     \item{`seed`}{Integer or `NULL`. Defaults to `seed`.}
#'   }
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
#'       baseline SE and the between-draw standard deviation. A
#'       \strong{descriptive spread for reporting and ranking only} — it is
#'       \strong{not} a Rubin's-rules pooled standard error (it contains no
#'       within-imputation variance component and no small-sample df
#'       correction) and must not be used for downstream pooling.}
#'     \item{`mi_workflow`}{`"pigauto_diagnostic_mi"`, recording that these
#'       are prediction-diagnostic completions and cannot be passed to
#'       [with_imputations()] or [pool_mi()].}
#'     \item{`imputed_mask`}{Logical matrix; `TRUE` where a cell was
#'       originally missing.}
#'     \item{`fit`}{The underlying [`pigauto_fit`][fit_pigauto()]
#'       object, retained for diagnostics and for calls to [predict()]
#'       on new data.}
#'     \item{`data`}{The [`pigauto_data`][preprocess_traits()] object.}
#'     \item{`tree`}{The input phylogeny.}
#'     \item{`species_col`}{Passed-through species-column name or
#'       `NULL`.}
#'     \item{`posterior`}{`draws_method = "posterior"` only. A list with
#'       `cell_interval` (data.frame: `row` (row of `traits`), `trait`,
#'       `lower`, `upper`, `median`; 95% posterior predictive interval on
#'       the original scale), `diagnostics` (data.frame: `parameter`, `rhat`,
#'       `ess_bulk`, with attribute `"converged"`), `params` (list:
#'       `Sigma_P` and `Sigma_E` as `K x K x draws` arrays, `lambda` and `mu`
#'       as `draws x K` matrices, on the latent scale), `converged`,
#'       `control`, `hyper` (priors), `start` (REML lambda and chain-1 start
#'       lambda), `draw_index` (kept sweeps used for the `m` datasets),
#'       `wall_s` and `sweeps`. For this method `fit` is `NULL`, `se` holds
#'       the posterior predictive SD of each imputed cell, the class is
#'       `c("pigauto_posterior_mi", "pigauto_mi", "list")` and
#'       `mi_workflow` is `"pigauto_posterior_mi_v1"`
#'       (`"pigauto_posterior_plugin_diagnostic"` with
#'       `posterior_control$param_uncertainty = "none"`, which
#'       [with_imputations()] and [pool_mi()] refuse).}
#'     \item{`posterior_improper`}{Only with
#'       `posterior_control$param_uncertainty = "both"` (validation only). A
#'       list with `datasets` (`m` completed data.frames drawn with the
#'       covariance matrices fixed at their posterior means, from the same
#'       chain run as `datasets`), `cell_interval` (same columns as
#'       `posterior$cell_interval`), and the fixed `Sigma_P` and `Sigma_E`.
#'       Not for downstream inference: it ignores parameter uncertainty.}
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
#' not an inference consequence of a nominal held-out conformal diagnostic. For discrete traits (binary,
#' categorical) it uses Bernoulli / categorical draws from the estimated
#' probability vector. For \code{zi_count} the gate is Bernoulli from
#' P(nonzero) and the conditional magnitude is drawn on the log1p-z
#' latent using the conformal score (or the magnitude latent SE when
#' scores are unavailable) — not the reported SE of expected count. For
#' \code{multi_proportion} groups it draws the
#' K CLR latent columns with their BM latent SEs, projects back to
#' sum-zero CLR space, and decodes to the simplex.
#'
#' **`draws_method = "mc_dropout"`**: Run `M` GNN forward passes in
#' training mode (dropout active) on top of stochastic BM baseline draws.
#' When `r_cal = 0`, the GNN-dropout term disappears but the BM draw still
#' contributes between-draw variance.
#'
#' **`draws_method = "posterior"`**: see the "Posterior draws" section.
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
#' Gelman A (2006). "Prior distributions for variance parameters in
#' hierarchical models." *Bayesian Analysis* 1(3): 515-534.
#'
#' Hadfield JD (2010). "MCMC methods for multi-response generalized linear
#' mixed models: the MCMCglmm R package." *Journal of Statistical Software*
#' 33(2): 1-22.
#'
#' Hadfield JD, Nakagawa S (2010). "General quantitative genetic methods for
#' comparative biology: phylogenies, taxonomies and multi-trait models for
#' continuous and categorical characters." *Journal of Evolutionary Biology*
#' 23(3): 494-508.
#'
#' Vehtari A, Gelman A, Simpson D, Carpenter B, Buerkner P-C (2021).
#' "Rank-normalization, folding, and localization: an improved R-hat for
#' assessing convergence of MCMC." *Bayesian Analysis* 16(2): 667-718.
#'
#' Nakagawa S, Freckleton RP (2011). "Model averaging, missing data and
#' multiple imputation: a case study for behavioural ecology."
#' *Behavioral Ecology and Sociobiology* 65(1): 103-116.
#'
#' @section Posterior draws:
#' With `draws_method = "posterior"`, the latent (z-scored, optionally
#' log-transformed) traits follow the multivariate phylogenetic mixed model
#' \deqn{\mathrm{vec}(Y) \sim \mathrm{N}(1\mu^\top,\;
#'   \Sigma_P \otimes R + \Sigma_E \otimes I_n),}
#' with \eqn{R} the phylogenetic correlation matrix of the tree, full
#' \eqn{K \times K} phylogenetic (\eqn{\Sigma_P}) and residual
#' (\eqn{\Sigma_E}) covariance matrices, and a flat prior on \eqn{\mu}.
#' The implied Pagel's lambda of trait \eqn{k} is
#' \eqn{\Sigma_P[k,k] / (\Sigma_P[k,k] + \Sigma_E[k,k])}.
#'
#' Priors: \eqn{\Sigma_P = \mathrm{diag}(\alpha)\,\Sigma_W\,
#' \mathrm{diag}(\alpha)} with \eqn{\Sigma_W \sim \mathrm{IW}(K+1, I_K)} and
#' \eqn{\alpha \sim \mathrm{N}(0, 1000\, I_K)} (parameter expansion as in
#' MCMCglmm; Hadfield 2010; Gelman 2006);
#' \eqn{\Sigma_E \sim \mathrm{IW}(K+1,\, 0.01\,\mathrm{diag}(s^2))}, where
#' \eqn{s^2} are the observed variances of the latent traits; and a flat
#' prior on \eqn{\mu}. Here \eqn{\mathrm{IW}(\nu, S)} has density
#' proportional to
#' \eqn{|\Sigma|^{-(\nu+K+1)/2}\exp\{-\mathrm{tr}(S\Sigma^{-1})/2\}}.
#' MCMCglmm parameterises the inverse-Wishart by \eqn{(V, \nu)} with scale
#' matrix \eqn{\nu V}, so its equivalent of the \eqn{\Sigma_W} prior is
#' \eqn{V = I_K/(K+1)}, \eqn{\nu = K+1}, not MCMCglmm's default prior. The
#' prior settings used are returned in `mi$posterior$hyper`.
#'
#' The sampler draws the phylogenetic effects, the trait means and the missing
#' cells as one block with a sparse Cholesky factor of the Hadfield and
#' Nakagawa (2010) node precision, and updates the covariance matrices by
#' parameter-expanded Gibbs steps plus Metropolis moves with the
#' phylogenetic effects integrated out. Convergence is checked with
#' rank-normalised split R-hat and bulk effective sample size (Vehtari et
#' al. 2021); `mi$posterior$diagnostics` reports them, and a warning is
#' raised when any R-hat is at least 1.05 or any ESS is at most 400.
#'
#' The `m` datasets are posterior predictive draws spaced evenly across the
#' kept sweeps of all chains. Per-cell 95% intervals in
#' `mi$posterior$cell_interval` come from all kept sweeps (`keep_draws` of
#' them by default), never from the `m` datasets.
#'
#' These draws come from a linear-Gaussian model of the imputed traits on
#' the scale where they were imputed (the log scale for traits
#' log-transformed by `log_transform`). They are proper imputations for
#' analyses that are linear in the imputed traits on that scale and whose
#' variables are all among the imputed traits, for example a phylogenetic
#' regression of one imputed trait on others (on the log scale for
#' log-transformed traits). Not covered: covariates from outside the imputed
#' traits, because the imputation model does not contain them; nonlinear
#' terms (such as squares) or interactions among imputed traits; and
#' analysing a log-transformed trait on its raw scale. The GNN is not used,
#' and `gnn`, `epochs`, `missing_frac`, `lambda_mode` and other fitting
#' arguments are ignored (a message lists any that were supplied).
#' Non-continuous traits, multiple observations per species (or any
#' `species_col`), `covariates`, and input with no missing cells are errors.
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
#' \donttest{
#' library(pigauto)
#' data(avonet300, tree300)
#' tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
#' df <- avonet300[match(tree$tip.label, avonet300$Species_Key),
#'                 c("Mass", "Wing.Length"), drop = FALSE]
#' rownames(df) <- tree$tip.label
#' df$Mass[seq_len(3L)] <- NA_real_
#' mi <- multi_impute(df, tree, m = 2L, epochs = 5L, verbose = FALSE)
#' print(mi)
#' lapply(mi$datasets, head)
#' }
#'
#' @export
multi_impute <- function(traits, tree, m = 100L,
                         draws_method = c("conformal", "mc_dropout",
                                          "posterior"),
                         species_col = NULL,
                         trait_types = NULL,
                         multi_proportion_groups = NULL,
                         log_transform = TRUE,
                         missing_frac = 0.25,
                         covariates = NULL,
                         epochs = 2000L, verbose = TRUE, seed = NULL,
                         gnn = TRUE,
                         lambda_mode = c("estimate", "fixed_1", "cv", "bayes"),
                         posterior_control = list(),
                         ...) {

  # Recorded before match.arg() reassigns lambda_mode (used only by
  # draws_method = "posterior" to report ignored fitting arguments).
  supplied_fit_args <- c(if (!missing(gnn)) "gnn",
                         if (!missing(epochs)) "epochs",
                         if (!missing(missing_frac)) "missing_frac",
                         if (!missing(lambda_mode)) "lambda_mode")
  draws_method <- match.arg(draws_method)
  lambda_mode <- match.arg(lambda_mode)
  m <- as.integer(m)
  if (!is.finite(m) || m < 2L) {
    stop("`m` must be an integer >= 2 (stochastic diagnostics need at least ",
         "two draws). Got m = ", m, ".", call. = FALSE)
  }

  if (draws_method == "posterior") {
    # ---- Posterior: Bayesian phylogenetic mixed model (R/mi_posterior.R) ----
    # No GNN is fitted; arguments that only configure the GNN fit are
    # reported as ignored.
    ignored <- c(supplied_fit_args, names(list(...)))
    return(.multi_impute_posterior(
      traits = traits, tree = tree, m = m, species_col = species_col,
      trait_types = trait_types,
      multi_proportion_groups = multi_proportion_groups,
      log_transform = log_transform, covariates = covariates,
      verbose = verbose, seed = seed, control = posterior_control,
      ignored_args = ignored))
  }

  if (draws_method == "mc_dropout") {
    # gnn = FALSE: there is no GNN to run dropout on, so predict() draws
    # from the BM posterior instead (see fit_pigauto()'s gnn = FALSE
    # contract). Not an error -- just a different (still valid) source of
    # between-imputation variance.
    if (isFALSE(gnn)) {
      message("gnn = FALSE: mc_dropout draws are BM-posterior draws (no dropout).")
    }
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
      seed          = if (is.null(seed)) NULL else as.integer(seed),
      gnn           = gnn,
      lambda_mode   = lambda_mode,
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
      seed          = if (is.null(seed)) NULL else as.integer(seed),
      gnn           = gnn,
      lambda_mode   = lambda_mode,
      ...
    )

    pred      <- res$prediction
    trait_map <- res$fit$trait_map
    imask     <- res$imputed_mask

    input_row_order <- res$data$input_row_order
    # `imask` is in user-input row order (built by build_completed from
    # the original traits data.frame).  But pred$imputed / pred$se /
    # pred$probabilities are in INTERNAL (tree-tip) order from the
    # model.  .sample_conformal_draw (called inside .conformal_draws)
    # therefore maps input-order mask indices to internal-order before
    # perturbing imp; otherwise it perturbs the wrong cells and
    # build_completed re-aligns the unperturbed (internal) values back
    # into the user's missing cells, producing zero between-imputation
    # variance.  See adversarial_review_opus.md (HIGH severity finding
    # C.1).  P1-11: this sampling block is shared with
    # multi_impute_trees(draws_method = "conformal") via .conformal_draws().
    datasets <- .conformal_draws(traits, pred, imask, trait_map, m,
                                 species_col = species_col, seed = seed,
                                 input_row_order = input_row_order)
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
      mi_workflow     = "pigauto_diagnostic_mi",
      evaluation      = res$evaluation
    ),
    class = c("pigauto_diagnostic_mi", "pigauto_mi", "list")
  )
}


# ---- Internal: M conformal-width completed datasets from one prediction ----
# Shared sampling block (P1-11): given a single-pass prediction (`pred`),
# draws `m` conformal-width stochastic completions and re-aligns each to
# the user's input row order via build_completed(). Both multi_impute()
# (draws_method = "conformal") and multi_impute_trees()
# (draws_method = "conformal") call this so the sampling mechanism is
# identical in both places -- only how `pred` / `imputed_mask` / `trait_map`
# are obtained differs (single-tree impute() vs per-tree impute() /
# predict()).
#
# @param traits           original traits data.frame (user input).
# @param pred              a `pigauto_pred` object from a single-pass
#                          (`n_imputations = 1`) predict() / impute() call.
# @param imputed_mask      logical matrix in user-input row order, TRUE for
#                          originally-missing cells (e.g. `res$imputed_mask`).
# @param trait_map         trait map from the fit that produced `pred`.
# @param m                 integer, number of draws to generate.
# @param species_col       character or NULL, forwarded to build_completed().
# @param seed              optional integer base seed; draw `i` uses
#                          `seed + i`. NULL uses the current RNG stream.
# @param input_row_order    forwarded to .sample_conformal_draw() /
#                          build_completed() for row re-alignment.
# @return list of `m` completed data.frames.
.conformal_draws <- function(traits, pred, imputed_mask, trait_map, m,
                             species_col = NULL, seed = NULL,
                             input_row_order = NULL) {
  lapply(seq_len(m), function(i) {
    imp_df <- .sample_conformal_draw(
      pred, imputed_mask, trait_map,
      seed_i = if (is.null(seed)) NULL else as.integer(seed) + i,
      input_row_order = input_row_order
    )
    build_completed(traits, imp_df, species_col,
                    input_row_order = input_row_order)$completed
  })
}


# ---- Internal: one conformal-width draw -------------------------------------
# Samples missing cells from N(mu, conformal_score / 1.96) for continuous
# types (on the appropriate transformed scale), and from Bernoulli /
# Categorical for discrete types. zi_count: Bernoulli gate + Normal draw
# on the magnitude log1p-z latent (conformal score or se_latent mag),
# never pred$se of E[X]. Falls back to BM / latent SE when the conformal
# score is not available for a trait.
.sample_conformal_draw <- function(pred, imputed_mask, trait_map, seed_i,
                                    input_row_order = NULL) {
  if (!is.null(seed_i)) set.seed(seed_i)
  imp    <- pred$imputed
  probs  <- pred$probabilities
  cscores <- pred$conformal_scores  # named; NA for binary/categorical

  # Build input-row -> internal-row index map.  `imputed_mask` rows are in
  # user-input order; `pred$imputed` / `pred$se` / `pred$probabilities` are
  # in internal (tree-tip) order.  `input_row_order[k] = i` means internal
  # row k holds original input row i, so the inverse map is over the full
  # input mask: `match(seq_len(nrow(imputed_mask)), input_row_order)`.  This
  # preserves NA entries for data-only rows after filtering. When NULL or
  # identity, internal order == input order and the conversion is a no-op.
  input_to_internal <- if (!is.null(input_row_order)) {
    match(seq_len(nrow(imputed_mask)), input_row_order)
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
      # Gate: Bernoulli from P(nonzero). Magnitude: log1p-z latent draw
      # using the conformal score (or se_latent of the mag column). Do
      # not use pred$se — that slot is the delta-method SE of E[X].
      p_nz <- pmin(pmax(probs[[nm]][rows], 0), 1)
      gate <- stats::rbinom(N, 1L, p_nz)
      lc   <- tm$latent_cols
      z_mu <- as.numeric(pred$imputed_latent[rows, lc[2L]])
      if (!is.null(cs)) {
        s_latent <- cs
      } else if (!is.null(pred$se_latent) &&
                 ncol(pred$se_latent) >= lc[2L]) {
        s_latent <- pred$se_latent[rows, lc[2L]]
        s_latent[!is.finite(s_latent)] <- 0
      } else {
        s_latent <- 0
      }
      z_drw  <- stats::rnorm(N, z_mu, s_latent)
      draw_c <- as.integer(pmax(round(expm1(z_drw * tm$sd + tm$mean)), 0L))
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

  if (identical(x$draws_method, "posterior")) {
    return(.print_posterior_mi(x, n_sp, traits, p, n_imp_cells, total_cells,
                               pct))
  }
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


# ---- Internal: print method body for draws_method = "posterior" ------------
.print_posterior_mi <- function(x, n_sp, traits, p, n_imp_cells, total_cells,
                                pct) {
  post <- x$posterior
  ctl <- post$control
  dg <- post$diagnostics
  cat("pigauto posterior multiple imputation (phylogenetic mixed model)\n")
  cat(sprintf("  M        : %d completed datasets (posterior predictive)\n",
              x$m))
  cat(sprintf("  Species  : %d\n", n_sp))
  cat(sprintf("  Traits   : %d -- %s\n", p, paste(traits, collapse = ", ")))
  cat(sprintf("  Cells    : %d imputed / %d total (%.1f%%)\n",
              n_imp_cells, total_cells, pct))
  cat(sprintf("  MCMC     : %d chains x (%d burn-in + %d sweeps, thin %d)%s\n",
              ctl$n_chains, ctl$burnin, ctl$n_iter, ctl$thin,
              switch(ctl$param_uncertainty,
                     none = "; covariances fixed (plug-in mode)",
                     both = "; plug-in draws also in mi$posterior_improper",
                     "")))
  lam <- colMeans(post$params$lambda)
  cat(sprintf("  lambda   : %s (posterior means)\n",
              paste(sprintf("%s = %.2f", names(lam), lam), collapse = ", ")))
  cat(sprintf(paste0("  Converged: %s (max R-hat %.3f, needs < 1.05; ",
                     "min bulk ESS %.0f, needs > 400)\n"),
              if (isTRUE(post$converged)) "yes" else "NO",
              max(dg$rhat, na.rm = TRUE), min(dg$ess_bulk, na.rm = TRUE)))
  if (!isTRUE(post$converged)) {
    cat("             Increase posterior_control$n_iter before using these draws.\n")
  }
  cat("\n  Access datasets:          mi$datasets[[i]]\n")
  cat("  Per-cell 95% intervals:   mi$posterior$cell_interval\n")
  if (identical(x$mi_workflow, "pigauto_posterior_plugin_diagnostic")) {
    cat("  Downstream inference:     unsupported for these draws (plug-in\n")
    cat("                            validation mode; rerun with\n")
    cat("                            param_uncertainty = \"full\")\n")
  } else {
    cat("  Downstream inference:     with_imputations(mi, f) then pool_mi()\n")
  }
  invisible(x)
}
