#' Fit the phylogenetic baseline
#'
#' Dispatches to pigauto's phylogenetic baseline machinery and returns imputed
#' latent-scale means and standard errors for every species.
#'
#' @details
#' When \code{splits} is supplied the val and test cells are masked to
#' \code{NA} before fitting, so the baseline is evaluated under the same
#' conditions as \code{\link{fit_pigauto}}.
#'
#' Continuous-family columns use Brownian-motion conditional MVN baselines on
#' the phylogenetic correlation matrix, either independently or through the
#' joint MVN path when the data and optional dependencies support it. Binary,
#' ordinal, categorical, and zero-inflated gate columns use the appropriate
#' label-propagation or threshold/liability baseline candidates, with
#' per-column fallbacks when a joint path is not available.
#'
#' \strong{Covariates and the joint baseline (P1-8)}: \code{data$covariates}
#' is only used by the per-column BM path (\code{bm_impute_col_with_cov()}).
#' The joint MVN and threshold-joint (Rphylopars) baselines do not accept a
#' covariate design matrix, so when a joint path is selected (BM-eligible
#' columns >= 2, or binary/ordinal cols present, with Rphylopars available)
#' any supplied covariates are ignored for the BASELINE and a warning is
#' emitted; covariates still reach the GNN correction via
#' \code{\link{fit_pigauto}} regardless of which baseline path fires.
#'
#' \strong{Per-type lambda dispatch (arc/lambda-per-type; corrected in S4/S6,
#' feat/joint-lambda-default)}: \code{lambda_mode} only ever governs the
#' baseline for CONTINUOUS-FAMILY columns (continuous, count, proportion,
#' zi_count magnitude) -- NOT ordinal. Binary, ordinal, categorical, and
#' zero-inflated gate columns always stay at lambda = 1 in EVERY path
#' (threshold-joint, OVR-categorical, label propagation, and the per-trait
#' ordinal path-selection fallback below), regardless of \code{lambda_mode}
#' -- there is no discrete-trait analogue of Pagel's lambda, and previously
#' forcing these columns onto label propagation any time \code{lambda_mode
#' != "fixed_1"} cost 19pp of Trophic.Level accuracy on AVONET (0.789 ->
#' 0.600; see \code{docs/dev-log/2026-08-16-external-comparison-results.md}).
#' When \code{lambda_mode != "fixed_1"} and the threshold-joint baseline
#' fires for a dataset with binary/ordinal AND continuous-family columns,
#' the joint liability fit's own continuous-column OUTPUT is now USED
#' (each continuous-family column gets its own \code{lambda_k} via the
#' joint solver's \code{lambda_cols}), rather than being discarded for a
#' separate lambda-aware per-column BM re-fit as in the pre-S4 design.
#' Ordinal liability columns inside that same joint fit, and the
#' BM-via-MVN / K-class-OVR-LP alternatives the "Per-trait ordinal path
#' selection" block below compares against it, all stay at lambda = 1.
#'
#' @param data object of class \code{"pigauto_data"}.
#' @param tree object of class \code{"phylo"}.
#' @param splits list (output of \code{\link{make_missing_splits}}) or
#'   \code{NULL}.
#' @param model character. Evolutionary model: \code{"BM"} (default) or
#'   \code{"OU"}.
#' @param graph optional list returned by \code{\link{build_phylo_graph}}.
#'   When supplied, \code{graph$D} (cophenetic distances) is reused for
#'   label propagation and \code{graph$R_phy} (phylogenetic correlation
#'   matrix) is reused for BM imputation, avoiding duplicate \eqn{O(n^2)}
#'   allocations. When \code{NULL} (default), both matrices are computed
#'   here.
#' @param multi_obs_aggregation character. How to aggregate multiple
#'   observations per species before the Level-C joint baseline:
#'   \code{"hard"} (default) thresholds binary proportions at 0.5 and uses
#'   argmax for categorical, matching Phase 10 behaviour.  \code{"soft"}
#'   preserves species-level proportions and dispatches the truncated-Gaussian
#'   soft E-step (\code{estep_liability_binary_soft}) so that intermediate
#'   class frequencies contribute fractional liability evidence.  Only
#'   relevant for multi-obs data with binary or categorical traits when the
#'   Level-C joint baseline is active.
#' @param lambda_mode character. Pagel-lambda mode for the CONTINUOUS-FAMILY
#'   baseline (continuous, count, proportion, zi_count magnitude columns;
#'   NOT ordinal, which the threshold-joint path treats as a liability
#'   column via \code{estep_liability_ordinal()} -- see \dQuote{Per-type
#'   lambda dispatch} in Details). \code{"estimate"} (default, S4) fits
#'   each continuous-family column's own Pagel's lambda via profile REML;
#'   when the joint MVN or threshold-joint baseline fires, this now runs
#'   INSIDE that joint fit (\code{R/joint_mvn_solver.R}'s per-column
#'   \code{lambda_cols}) rather than being discarded in favour of a
#'   separate per-column re-fit. \code{"fixed_1"} preserves the classic
#'   Brownian correlation matrix (lambda = 1 everywhere). \code{"cv"} and
#'   \code{"bayes"} have no joint analogue and force the per-column BM path
#'   for continuous-family columns (as before); binary/ordinal/categorical/
#'   zi_gate columns are unaffected by \code{lambda_mode} in every case and
#'   keep the threshold-joint / OVR-categorical baseline at lambda = 1.
#'   \strong{Covariate caveat}: when \code{data$covariates} is supplied,
#'   the per-column path switches to \code{bm_impute_col_with_cov()},
#'   which accepts a numeric lambda or \code{"estimate"} but not
#'   \code{"cv"} / \code{"bayes"}; those two modes are silently ignored
#'   (fit at lambda = 1) for BM-eligible columns and a warning is emitted.
#' @param lambda_fixed optional named numeric vector (names = latent column
#'   names, i.e. \code{colnames(data$X_scaled)}) giving a FIXED lambda per
#'   continuous-family column, overriding \code{lambda_mode} entirely for
#'   those columns (spec 4.5: reproduce a previous fit's exact per-trait
#'   lambda at predict time without re-estimating -- typically supply that
#'   fit's own \code{$lambda_per_trait}). Columns not present in
#'   \code{lambda_fixed} keep their lambda = 1 default. \code{NULL}
#'   (default) means every continuous-family column follows
#'   \code{lambda_mode} normally.
#' @param em_iterations integer. Number of Phase 6 EM iterations for the
#'   threshold-joint baseline (binary + ordinal + OVR categorical). Default
#'   \code{0L} disables the EM loop and preserves v0.9.1 output byte-for-byte.
#'   When \code{>= 1}, the BM rate \eqn{\Sigma} learned by the in-house
#'   joint solver (\code{R/joint_mvn_solver.R}) at iteration \eqn{k} is fed back as the
#'   per-trait prior SD at iteration \eqn{k+1}, up to \code{em_iterations}
#'   times or until \code{em_tol} convergence.  \code{em_iterations = 1L} is
#'   a degenerate single-pass run and produces the same baseline output as
#'   \code{0L}; \code{>= 2L} is needed for actual iteration. Only affects
#'   the threshold-joint path (continuous-only traits pass through the
#'   existing joint MVN path unchanged).
#' @param em_tol numeric. Relative-Frobenius convergence tolerance for the
#'   Phase 6 / 7 EM loop. Early-stops when
#'   \eqn{||\Sigma_k - \Sigma_{k-1}||_F / ||\Sigma_{k-1}||_F < }
#'   \code{em_tol}.  Default \code{1e-3}.
#' @param em_offdiag logical. Phase 7 opt-in: when \code{TRUE} AND
#'   \code{em_iterations >= 2L}, each liability cell's prior at iteration
#'   \eqn{k+1} is the conditional-MVN \eqn{(\mu, sd)} given the posterior
#'   liability of other traits at iteration \eqn{k}, using the full off-
#'   diagonal entries of \eqn{\Sigma}. Binary + ordinal only (OVR categorical
#'   stays on Phase 6 diagonal). Default \code{FALSE} preserves Phase 6
#'   behaviour.
#' @param joint_solver character. Which solver estimates the joint
#'   Sigma / posterior for the joint MVN, threshold-joint, and OVR
#'   categorical baselines. \code{"inhouse"} (default) uses the
#'   single-pass in-house solver (\code{R/joint_mvn_solver.R}); under
#'   \code{lambda_mode = "fixed_1"} it is byte-identical to prior
#'   releases, but \code{lambda_mode = "estimate"} (the current default)
#'   is new behaviour, not a byte-compatibility guarantee.
#'   \code{"rphylopars"} delegates to \code{Rphylopars::phylopars()};
#'   under \code{"estimate"} this calls it with \code{model = "lambda"}
#'   (substantially slower than \code{model = "BM"}), with automatic
#'   fallback to \code{"inhouse"} on failure or implausible output (a
#'   plausibility guard that also fires under \code{lambda_mode =
#'   "fixed_1"} -- see NEWS). Only affects the joint MVN / threshold-joint
#'   / OVR categorical paths above; ignored when those paths don't fire.
#'   \code{lambda_mode} does NOT disable the continuous-only joint MVN
#'   path: both \code{fit_joint_mvn_baseline()} and
#'   \code{fit_joint_threshold_baseline()} accept a \code{lambda_mode} /
#'   \code{lambda_fixed} argument and estimate lambda inside the joint fit
#'   via \code{lambda_cols} -- see \dQuote{Per-type lambda dispatch} in
#'   Details. Only \code{"cv"} / \code{"bayes"} force continuous-family
#'   columns off the joint path entirely (no joint analogue for those two
#'   modes).
#' @param predict_method character. Prediction route for the in-house joint
#'   solver. \code{"auto"} (default, S5b/S5c) fits the baseline once with
#'   the \code{"exact"} route and once with the \code{"per_column"} route,
#'   then for each TRAIT (each \code{trait_map} entry; a categorical
#'   trait's K latent columns and a zi_count trait's gate + magnitude
#'   columns are chosen together) picks whichever route has the lower loss
#'   on that trait's cells: mean squared error on the z-scored latent
#'   scale for continuous/count/proportion/ordinal/zi_count magnitude,
#'   mean log-loss of \code{plogis(mu)} against the observed 0/1 truth for
#'   binary/zi_count gate, and mean multinomial log-loss for categorical.
#'   \strong{Which cells the choice sees}: a trait's held-out validation
#'   rows (species in multi-observation data) are split into two halves
#'   (seeded via \code{seed}, or the ambient RNG state when
#'   \code{seed = NULL}) when the trait has at least 38 of them and the two
#'   candidate fits predict it differently. The route choice sees only half
#'   A; half B is reserved so \code{fit_pigauto()}/\code{impute()} calibrate
#'   the GNN gate and compute conformal scores on rows that never informed
#'   the choice. Below 38 rows each half would keep fewer than the 19 a 95\%
#'   split-conformal interval needs, so the same rows both choose the route
#'   and calibrate; when both fits agree, no choice is made and all rows
#'   calibrate. With fewer than 5 validation cells, or on a genuine tie, a
#'   trait keeps \code{"exact"}, as it does when \code{splits} is
#'   \code{NULL} (no validation cells at all makes every trait
#'   \code{"exact"}). The chosen route per trait is returned as
#'   \code{$predict_method_by_trait}. \code{"exact"} uses the full
#'   cross-trait conditional
#'   mean and variance of \code{vec(L) ~ MVN(0, Sigma \%x\% R(lambda_block))}
#'   in the sparse precision form (Hadfield & Nakagawa, 2010), with each
#'   column centred at its own GLS phylogenetic mean at \code{lambda_block}
#'   before the solve (mean-model consistency; see
#'   \code{docs/dev-log/exact-default/}). Discrete liability columns
#'   (binary, zi gate, ordinal-via-OVR synthetic columns) share
#'   \code{lambda_block} under \code{"exact"} rather than staying fixed at
#'   lambda = 1, so the whole joint fit uses one internally-consistent
#'   \code{R(lambda_block)}. Falls back to \code{"per_column"} above
#'   roughly 20000 unknown cells (roughly 4000 species at 5 traits), on a
#'   singular/unusable Sigma, or when fewer than 2 joint columns or no
#'   Henderson sparse precision are available; the fallback prints a
#'   one-time \code{message()} (not a warning) per R session when
#'   \code{predict_method} was left at its default, or a \code{warning()}
#'   every time when \code{"exact"} was requested explicitly.
#'   \code{"per_column"} retains the original per-column conditional
#'   prediction route (each column's own posterior, no cross-trait
#'   borrowing in the prediction step; \code{se} is that column's own
#'   conditional SE). Neither option changes covariance estimation or the
#'   \code{"rphylopars"} solver.
#' @param joint_refine_iter integer, default \code{0L}. Enables
#'   cross-trait refinement of the joint baseline's cell imputations
#'   using the estimated Sigma (the in-house solver's \code{max_iter}
#'   EM cell-refinement; \code{R/joint_mvn_solver.R}). \code{0L} preserves
#'   current behaviour byte-for-byte. The refinement is guarded: the
#'   Sigma step must shrink each iteration, or the loop rolls back to the
#'   last good iterate and sets \code{$diverged}. \strong{Has no effect on
#'   any trait predicted by the \code{"exact"} route}: the exact
#'   conditional solve returns its analytic posterior directly and never
#'   enters this EM loop. Under the default \code{predict_method =
#'   "auto"}, \code{joint_refine_iter} therefore only ever applies to
#'   traits that \code{"auto"} routes to \code{"per_column"} (see
#'   \code{$predict_method_by_trait} below).
#' @param predict_route optional named character vector (names = trait
#'   names, i.e. \code{vapply(data$trait_map, function(tm) tm$name,
#'   character(1))}; values \code{"exact"} or \code{"per_column"}). Forces
#'   the prediction route for each named trait, bypassing the
#'   \code{"auto"} validation-loss comparison entirely (and overriding
#'   \code{predict_method}, if a value other than \code{"auto"} was also
#'   supplied). Traits not named default to \code{"exact"}. A name that
#'   matches no trait in \code{data$trait_map} is ignored, with a warning
#'   (once per call). Used internally to replay a validation-split
#'   \code{"auto"} choice onto a production refit that has no validation
#'   cells (e.g. \code{fit_pigauto()}'s \code{baseline_full}, fit with
#'   \code{splits = NULL}) so that refit reuses the SAME per-trait route
#'   rather than re-deciding with zero validation evidence. \code{NULL}
#'   (default) disables forcing; \code{predict_method} governs normally.
#' @param seed optional integer. Seeds the deterministic half-A (route
#'   choice) / half-B (reserved for \code{fit_pigauto()}/\code{impute()}'s
#'   gate calibration and conformal scoring) split of each trait's
#'   validation cells under \code{predict_method = "auto"} with real
#'   validation cells (see \dQuote{"auto"} above). \code{NULL} (default)
#'   leaves the split to whatever the ambient RNG state is at call time,
#'   matching this package's other optional \code{seed} arguments.
#'   Ignored when \code{predict_method} is \code{"exact"} /
#'   \code{"per_column"}, when \code{predict_route} is supplied, or when
#'   there are no validation cells.
#' @return A list with:
#'   \describe{
#'     \item{mu}{Numeric matrix (n_species x p_latent), baseline means in
#'       latent scale.}
#'     \item{se}{Numeric matrix (n_species x p_latent), standard errors.}
#'     \item{path}{Named character vector, one entry per \code{trait_map}
#'       entry (name = trait/group name), recording which dispatch branch
#'       produced that trait's baseline: \code{"joint_mvn"} (Phase 2
#'       continuous-only joint), \code{"threshold_joint"} (Phase 3/5 binary +
#'       ordinal + continuous-family liability joint), \code{"ovr_categorical"}
#'       (Phase 6 one-vs-rest K-fits), \code{"per_column_bm"} (per-column
#'       Brownian-motion conditional MVN, including the lambda-aware path and
#'       the per-trait ordinal fallback that beats threshold-joint on held-out
#'       val MSE), \code{"multi_proportion_bm"} (per-component
#'       \code{per_column_bm} on the CLR columns of a multi_proportion group
#'       -- same kernel as \code{"per_column_bm"}, labelled separately because
#'       multi_proportion is architecturally distinct), \code{"label_propagation"}
#'       (phylogenetic label propagation for binary/categorical/ordinal, and
#'       the per-trait ordinal fallback that beats the other two candidates),
#'       \code{"zi_gate_lp"} (label propagation for a zero-inflated count's
#'       gate column), or \code{"zi_mag_constant"} (global mean/sd fallback
#'       for a zero-inflated count's magnitude column when fewer than 5
#'       non-zero observations exist). For \code{zi_count} traits (2 latent
#'       columns: gate then magnitude), \code{path} reports the GATE column's
#'       dispatch -- the magnitude column can independently land on
#'       \code{"joint_mvn"}, \code{"threshold_joint"}, \code{"per_column_bm"},
#'       or \code{"zi_mag_constant"} and is not separately surfaced here.}
#'     \item{lambda_per_trait}{Named numeric vector (length \code{p_latent},
#'       names = \code{colnames(data$X_scaled)}), the lambda actually used
#'       for each latent column: 1 for every column not eligible for lambda
#'       estimation (all discrete columns; continuous-family columns under
#'       \code{lambda_mode \%in\% c("cv", "bayes")}), else the estimated or
#'       fixed value.}
#'     \item{lambda_block}{Numeric scalar, the shared lambda used internally
#'       by whichever joint fit ran for its Sigma M-step / opt-in exact
#'       conditional / opt-in EM refine (see \code{R/joint_mvn_solver.R}'s
#'       \code{$lambda_block}); \code{NA} when no joint fit ran.}
#'     \item{lambda_mode}{Character, echoes the resolved \code{lambda_mode}
#'       argument.}
#'     \item{predict_method_used}{Character scalar, \code{"exact"},
#'       \code{"per_column"}, or \code{"auto"}. Under a concrete
#'       \code{predict_method} (\code{"exact"} / \code{"per_column"}),
#'       aggregated across every joint fit that ran: \code{"exact"} only if
#'       every joint fit achieved exact; \code{"per_column"} if any fell
#'       back, if none ran (a pure per-column baseline), or if every joint
#'       fit used \code{joint_solver = "rphylopars"} (no exact/per_column
#'       dichotomy there). Under \code{predict_method = "auto"} (or a
#'       forced \code{predict_route}), always \code{"auto"}; see
#'       \code{$predict_method_by_trait} for the resolved per-trait routes.}
#'     \item{predict_method_by_trait}{Named character vector (names =
#'       trait names), one entry per \code{trait_map} entry, giving the
#'       route (\code{"exact"} or \code{"per_column"}) actually used for
#'       that trait's cells. Under a concrete \code{predict_method}, every
#'       entry equals \code{predict_method_used}. Under \code{"auto"} (or a
#'       forced \code{predict_route}), the per-trait validation choice (see
#'       \code{predict_method}'s \code{"auto"} entry above).}
#'   }
#' @examples
#' \donttest{
#' data(avonet300, tree300, package = "pigauto")
#' tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
#' traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
#'                      c("Mass", "Wing.Length"), drop = FALSE]
#' rownames(traits) <- tree$tip.label
#' pd     <- preprocess_traits(traits, tree)
#' splits <- make_missing_splits(pd$X_scaled, trait_map = pd$trait_map)
#' bl     <- fit_baseline(pd, tree, splits)
#' }
#' @importFrom stats complete.cases rnorm rbinom
#' @export
fit_baseline <- function(data, tree, splits = NULL, model = "BM",
                         graph = NULL,
                         multi_obs_aggregation = c("hard", "soft"),
                         lambda_mode = c("estimate", "fixed_1", "cv", "bayes"),
                         lambda_fixed = NULL,
                         em_iterations = 0L,
                         em_tol = 1e-3,
                         em_offdiag = FALSE,
                         joint_solver = c("inhouse", "rphylopars"),
                         predict_method = c("auto", "exact", "per_column"),
                         joint_refine_iter = 0L,
                         predict_route = NULL,
                         seed = NULL) {
  # S5c (rose-review.md, required change 4): `predict_method_explicit` used
  # to be a documented-but-"internal use only" argument of this exported
  # signature. It is resolved here, from THIS call's own missing() check,
  # and threaded to the unexported .fit_baseline_dispatch() instead --
  # fit_pigauto()/impute() (the only callers that need to propagate a
  # DIFFERENT caller's missing() verdict) call .fit_baseline_dispatch()
  # directly rather than this exported wrapper.
  predict_method_explicit <- !missing(predict_method)
  .fit_baseline_dispatch(
    data = data, tree = tree, splits = splits, model = model, graph = graph,
    multi_obs_aggregation = multi_obs_aggregation, lambda_mode = lambda_mode,
    lambda_fixed = lambda_fixed, em_iterations = em_iterations, em_tol = em_tol,
    em_offdiag = em_offdiag, joint_solver = joint_solver,
    predict_method = predict_method, joint_refine_iter = joint_refine_iter,
    predict_method_explicit = predict_method_explicit,
    predict_route = predict_route, seed = seed)
}

# S5b/S5c: "auto" and `predict_route` are dispatcher-level concerns, handled
# entirely by this unexported dispatcher -- .fit_baseline_core() (the
# pre-S5b fit_baseline() body) only ever sees a concrete "exact" /
# "per_column". Takes `predict_method_explicit` directly (unlike the
# exported fit_baseline(), which infers it from its OWN missing() check)
# so fit_pigauto()/impute() can propagate whether THEIR OWN caller
# explicitly requested `predict_method`, without exposing that plumbing on
# the public signature.
#
#' @noRd
.fit_baseline_dispatch <- function(data, tree, splits = NULL, model = "BM",
                         graph = NULL,
                         multi_obs_aggregation = c("hard", "soft"),
                         lambda_mode = c("estimate", "fixed_1", "cv", "bayes"),
                         lambda_fixed = NULL,
                         em_iterations = 0L,
                         em_tol = 1e-3,
                         em_offdiag = FALSE,
                         joint_solver = c("inhouse", "rphylopars"),
                         predict_method = c("auto", "exact", "per_column"),
                         joint_refine_iter = 0L,
                         predict_method_explicit = NULL,
                         predict_route = NULL,
                         seed = NULL) {
  if (is.null(predict_method_explicit)) {
    predict_method_explicit <- !missing(predict_method)
  }
  predict_method <- match.arg(predict_method)

  if (!is.null(predict_route)) {
    if (!is.character(predict_route) || is.null(names(predict_route)) ||
        any(!nzchar(names(predict_route))) ||
        !all(predict_route %in% c("exact", "per_column"))) {
      stop("'predict_route' must be a named character vector (names = ",
           "trait names) with values 'exact' or 'per_column'.",
           call. = FALSE)
    }
    return(.fit_baseline_route(
      data = data, tree = tree, splits = splits, model = model, graph = graph,
      multi_obs_aggregation = multi_obs_aggregation, lambda_mode = lambda_mode,
      lambda_fixed = lambda_fixed, em_iterations = em_iterations,
      em_tol = em_tol, em_offdiag = em_offdiag, joint_solver = joint_solver,
      joint_refine_iter = joint_refine_iter, predict_route = predict_route))
  }

  if (identical(predict_method, "auto")) {
    return(.fit_baseline_auto(
      data = data, tree = tree, splits = splits, model = model, graph = graph,
      multi_obs_aggregation = multi_obs_aggregation, lambda_mode = lambda_mode,
      lambda_fixed = lambda_fixed, em_iterations = em_iterations,
      em_tol = em_tol, em_offdiag = em_offdiag, joint_solver = joint_solver,
      joint_refine_iter = joint_refine_iter, seed = seed))
  }

  .fit_baseline_core(
    data = data, tree = tree, splits = splits, model = model, graph = graph,
    multi_obs_aggregation = multi_obs_aggregation, lambda_mode = lambda_mode,
    lambda_fixed = lambda_fixed, em_iterations = em_iterations, em_tol = em_tol,
    em_offdiag = em_offdiag, joint_solver = joint_solver,
    predict_method = predict_method, joint_refine_iter = joint_refine_iter,
    predict_method_explicit = predict_method_explicit)
}

#' @noRd
.fit_baseline_core <- function(data, tree, splits = NULL, model = "BM",
                         graph = NULL,
                         multi_obs_aggregation = c("hard", "soft"),
                         lambda_mode = c("estimate", "fixed_1", "cv", "bayes"),
                         lambda_fixed = NULL,
                         em_iterations = 0L,
                         em_tol = 1e-3,
                         em_offdiag = FALSE,
                         joint_solver = c("inhouse", "rphylopars"),
                         predict_method = c("exact", "per_column"),
                         joint_refine_iter = 0L,
                         predict_method_explicit = NULL) {
  # S3 default flip (docs/dev-log/exact-default/S3-default-report.md):
  # resolve BEFORE match.arg() reassigns `predict_method` (reassignment
  # does not retroactively change what missing() reports, but this still
  # has to run first). A concrete TRUE/FALSE passed in by a caller further
  # up the stack (fit_pigauto(), impute()) overrides the local missing()
  # check, because those callers always forward a resolved concrete value
  # here, which would otherwise make missing() wrongly report FALSE
  # (explicit) regardless of what the ORIGINAL top-level caller did.
  if (is.null(predict_method_explicit)) {
    predict_method_explicit <- !missing(predict_method)
  }
  multi_obs_aggregation <- match.arg(multi_obs_aggregation)
  soft_aggregate <- identical(multi_obs_aggregation, "soft")
  lambda_mode <- match.arg(lambda_mode)
  joint_solver <- match.arg(joint_solver)
  predict_method <- match.arg(predict_method)
  if (!is.numeric(joint_refine_iter) || length(joint_refine_iter) != 1L ||
      !is.finite(joint_refine_iter) ||
      joint_refine_iter != as.integer(joint_refine_iter) ||
      joint_refine_iter < 0L) {
    stop("'joint_refine_iter' must be a non-negative integer scalar.",
         call. = FALSE)
  }
  joint_refine_iter <- as.integer(joint_refine_iter)
  # Translate dispatcher mode to the kernel-layer lambda argument.
  # "fixed_1"  -> lambda = 1.0   (back-compat; bit-identical to v0.9.x)
  # "estimate" -> lambda = "estimate" (per-column ML)
  # "cv"       -> per-column CV-selected lambda, computed below
  bm_lambda <- switch(lambda_mode,
    "fixed_1"  = 1.0,
    "estimate" = "estimate",
    "cv"       = "cv",
    "bayes"    = "bayes",
    1.0
  )
  if (!is.null(lambda_fixed)) {
    if (!is.numeric(lambda_fixed) || is.null(names(lambda_fixed)) ||
        any(!nzchar(names(lambda_fixed)))) {
      stop("'lambda_fixed' must be a named numeric vector (names = latent ",
           "column names).", call. = FALSE)
    }
    if (!all(is.finite(lambda_fixed)) || any(lambda_fixed < 0) ||
        any(lambda_fixed > 1)) {
      stop("'lambda_fixed' entries must be finite and in [0, 1].",
           call. = FALSE)
    }
  }
  em_iterations <- as.integer(em_iterations)
  em_offdiag    <- isTRUE(em_offdiag)
  if (!is.finite(em_iterations) || em_iterations < 0L) {
    stop("'em_iterations' must be a non-negative integer.", call. = FALSE)
  }
  if (em_offdiag && em_iterations < 2L) {
    # Silent: em_offdiag has no effect at em=0 (no EM at all) or em=1
    # (plug-in path; no previous Σ to condition on).
    em_offdiag <- FALSE
  }
  if (!inherits(data, "pigauto_data")) {
    stop("'data' must be a pigauto_data object (output of preprocess_traits).")
  }
  if (!inherits(tree, "phylo")) stop("'tree' must be a phylo object.")

  X   <- data$X_scaled
  p   <- ncol(X)

  multi_obs <- isTRUE(data$multi_obs)
  if (multi_obs) {
    # Multi-obs: X is n_obs x p, species_names is n_species
    n_obs     <- data$n_obs
    n_species <- data$n_species
    spp       <- data$species_names          # unique species (n_species)
    obs_spp   <- data$obs_species            # species per obs (n_obs)
    obs_to_sp <- data$obs_to_species         # integer mapping (n_obs)
  } else {
    n_obs     <- nrow(X)
    n_species <- n_obs
    spp       <- data$species_names
    obs_spp   <- spp
    obs_to_sp <- NULL
  }

  # ---- Phylogenetic similarity for discrete-trait label propagation ------
  # Reuse build_phylo_graph()'s cached D if supplied; otherwise compute.
  # At n = 10,000 each cophenetic() call is ~15 seconds and ~800 MB of
  # allocation, so caching through `graph` is a meaningful speedup even
  # though this stage is not the dominant scaling bottleneck.
  if (!is.null(graph) && !is.null(graph$D)) {
    D_phylo <- graph$D
  } else {
    D_phylo <- ape::cophenetic.phylo(tree)
  }
  # Reorder to match species order
  D_phylo  <- D_phylo[spp, spp]
  sigma_lp <- stats::median(D_phylo) * 0.5
  sim_phylo <- exp(-(D_phylo^2) / (2 * sigma_lp^2))
  diag(sim_phylo) <- 0  # exclude self for label propagation

  # Mask val + test cells before fitting
  if (!is.null(splits)) {
    X[splits$val_idx]  <- NA
    X[splits$test_idx] <- NA
  }

  trait_map <- data$trait_map

  # Output at species level (n_species x p)
  mu <- matrix(0, nrow = n_species, ncol = p)
  se <- matrix(0, nrow = n_species, ncol = p)
  dimnames(mu) <- list(spp, data$latent_names)
  dimnames(se) <- list(spp, data$latent_names)

  # Dispatch tracker (S1): one entry per latent column, filled in at the
  # same site that writes mu[, col] / se[, col] for that column, so it can
  # never drift from the branch that actually ran. Reduced to a per-trait
  # `path` vector (first latent column of each trait_map entry) at the end.
  col_path <- rep(NA_character_, p)

  # S5c (rose-review.md, required change 3): per-column "exact"/"per_column"
  # route actually used, tracked at the same sites that set `col_path` to a
  # joint-dispatch branch. Unlike `predict_method_used` (a single aggregate
  # over every joint fit in this call), this lets `predict_method_by_trait`
  # report each trait's OWN outcome even when only SOME joint fits in a
  # concrete "exact"/"per_column" request fell back (e.g. one OVR class
  # fit). Left NA_character_ until the corresponding dispatch site sets it;
  # any column that never goes through a joint call (label propagation,
  # per-column BM, multi_proportion, zi magnitude constant) is filled with
  # "per_column" right before `path` / `predict_method_by_trait` are built
  # below.
  trait_route <- rep(NA_character_, p)

  # S4 lambda diagnostics: one entry per latent column, 1 for every column
  # never touched by lambda estimation (every discrete column, and every
  # continuous-family column under lambda_mode = "fixed_1" / "cv" / "bayes").
  # Filled in at the same sites that populate mu[, col] / se[, col] for a
  # lambda-aware column, from either a joint fit's own $lambda_per_trait_fit
  # / $lambda_per_trait (mapped back by column NAME) or a per-column
  # bm_impute_col() / bm_impute_col_with_cov() call's $lambda_hat (or the
  # known fixed numeric value when no $lambda_hat is returned).
  # `lambda_block_out` is the block value from whichever joint fit ran
  # (NA if none did); see R/joint_mvn_solver.R's $lambda_block.
  lambda_per_trait <- stats::setNames(rep(1, p), colnames(X))
  lambda_block_out <- NA_real_

  # ---- Identify BM-eligible columns (continuous in latent space) -----------
  bm_cols <- integer(0)
  has_multi_proportion <- FALSE  # track for joint-dispatch guard
  mp_cols <- integer(0)          # multi_proportion latent cols (path labelling only)
  zi_mag_fallback <- integer(0)  # ZI magnitude cols with too few non-zero obs
  for (tm in trait_map) {
    if (tm$type %in% c("continuous", "count", "ordinal", "proportion",
                       "multi_proportion")) {
      # multi_proportion: K independent BM fits, one per CLR column
      if (tm$type == "multi_proportion") {
        has_multi_proportion <- TRUE
        mp_cols <- c(mp_cols, tm$latent_cols)
      }
      bm_cols <- c(bm_cols, tm$latent_cols)
    } else if (tm$type == "zi_count") {
      # Magnitude column (col 2) is BM-eligible if enough non-zero obs
      mag_col <- tm$latent_cols[2]
      n_finite <- sum(is.finite(X[, mag_col]))
      if (n_finite >= 5L) {
        bm_cols <- c(bm_cols, mag_col)
      } else {
        # Fallback: constant imputation (global mean of non-zero values)
        zi_mag_fallback <- c(zi_mag_fallback, mag_col)
        finite_vals <- X[is.finite(X[, mag_col]), mag_col]
        mu[, mag_col] <- if (length(finite_vals) > 0) mean(finite_vals) else 0
        se[, mag_col] <- if (length(finite_vals) > 1) stats::sd(finite_vals) else 0
        col_path[mag_col] <- "zi_mag_constant"
      }
    }
  }

  # ---- Level-C Phase 2, 3, 4 & 5: joint baseline dispatch -------------------
  # Binary and ZI-count gate cols both take the threshold-joint path (both
  # are binary-like observations with a truncated-Gaussian E-step).
  binary_cols  <- integer(0)
  zi_gate_cols <- integer(0)
  cat_cols     <- integer(0)
  ordinal_cols <- integer(0)
  for (tm in trait_map) {
    if (tm$type == "binary") {
      binary_cols <- c(binary_cols, tm$latent_cols)
    } else if (tm$type == "zi_count") {
      zi_gate_cols <- c(zi_gate_cols, tm$latent_cols[1])
    } else if (tm$type == "categorical") {
      cat_cols <- c(cat_cols, tm$latent_cols)
    } else if (tm$type == "ordinal") {
      ordinal_cols <- c(ordinal_cols, tm$latent_cols)
    }
  }
  # ZI gates join binary for dispatch purposes.
  binary_cols <- c(binary_cols, zi_gate_cols)

  # Phase 4 scope: threshold-joint fires when (binary) + BM are present.
  # Categorical-only (no binary) datasets fall back to Phase 2 MVN + LP:
  # Rphylopars has numerical instability with multi-categorical liability
  # matrices (the rank-(K-1) drop + multiple cat groups combine badly).
  # Phase 6 EM will refine this once Sigma is estimated stably.
  #
  # arc/lambda-per-type (2026-08) + S4 (feat/joint-lambda-default, 2026-09):
  # `lambda_mode` governs ONLY where CONTINUOUS-FAMILY columns get their
  # baseline mu/se, not whether the discrete-trait joint machinery runs at
  # all. Originally `lambda_mode != "fixed_1"` set `force_per_column <-
  # TRUE`, which also disabled `use_threshold_joint` and the OVR-categorical
  # loop below -- i.e. binary/ordinal/categorical traits were pushed onto
  # plain label propagation any time a user asked for lambda estimation on
  # their continuous traits. Measured cost: on AVONET data this dropped
  # Trophic.Level (categorical) accuracy from 0.789 to 0.600 (19pp) while
  # lambda_mode="bayes" only ever improves continuous traits -- discrete
  # traits have no lambda concept and were always fit at lambda = 1
  # regardless (see docs/dev-log/2026-08-16-external-comparison-results.md
  # and NEWS.md). `force_per_column` below therefore gates ONLY
  # `use_continuous_joint` (the continuous-only joint MVN path); the
  # threshold-joint / OVR-categorical dispatch further down is
  # unconditionally eligible whenever its other preconditions hold,
  # independent of `lambda_mode`.
  #
  # S4 update: `force_per_column` is now TRUE only for "cv" / "bayes" --
  # they have no joint-solver analogue (R/joint_mvn_solver.R's `lambda`
  # argument only understands "fixed_1", "estimate", or a numeric value;
  # see docs/dev-log/2026-09-22-joint-lambda-alignment.md section 5).
  # "estimate" now stays on the joint path: the in-house solver estimates a
  # per-trait Pagel's lambda for the continuous-family columns INSIDE the
  # joint fit (`lambda_cols`, S2/S4), so there is no need to discard that
  # fit's continuous-column output and re-fit it separately any more (see
  # the `cont_idx` block a few lines down, which used to gate on
  # `identical(lambda_mode, "fixed_1")` and now gates on
  # `lambda_mode %in% c("fixed_1", "estimate")`). Binary/ordinal/
  # categorical/zi_gate columns keep the joint/OVR baseline, fit at
  # lambda = 1 as they always have -- lambda_mode never touches them
  # (section 7 B iii cut of the alignment note: this is a spec non-goal,
  # and `tests/testthat/test-lambda-per-type.R` locks it in as a
  # regression gate).
  force_per_column <- lambda_mode %in% c("cv", "bayes")
  # S4: what to hand the joint helpers below. "cv"/"bayes" have no joint
  # analogue (per the alignment note section 5) so they always run the
  # joint fit's OWN internal lambda at "fixed_1" -- irrelevant to their
  # output anyway, since `force_per_column` disables use_continuous_joint
  # and the threshold-joint `cont_idx` populate block for those two modes,
  # so the joint fit's continuous-column output is always discarded for
  # "cv"/"bayes" regardless of what lambda it ran at internally.
  lambda_mode_joint <- if (lambda_mode %in% c("fixed_1", "estimate")) {
    lambda_mode
  } else {
    "fixed_1"
  }
  use_threshold_joint <- (length(binary_cols) + length(ordinal_cols)) >= 1L &&
    length(bm_cols) >= 1L &&
    !has_multi_proportion &&
    joint_mvn_available()

  use_continuous_joint <- !use_threshold_joint &&
    length(bm_cols) >= 2L &&
    !has_multi_proportion &&
    !force_per_column &&
    joint_mvn_available()

  # P1-8: the joint MVN / threshold-joint (Rphylopars) paths have no covariate
  # design matrix, so user covariates never reach the BASELINE when a joint
  # path fires. They still reach the GNN correction in fit_pigauto(). Warn
  # once rather than silently ignoring them.
  if (!is.null(data$covariates) && (use_threshold_joint || use_continuous_joint)) {
    warning(
      "'covariates' are ignored by the joint ",
      if (use_threshold_joint) "threshold-joint" else "MVN",
      " baseline (Rphylopars has no covariate design). The baseline is ",
      "covariate-free; covariates still reach the GNN correction. To get a ",
      "covariate-aware baseline, force the per-column path (e.g. supply a ",
      "single BM-eligible trait, or uninstall/mask Rphylopars).",
      call. = FALSE
    )
  }

  # Aggregated across whichever joint dispatcher(s) run below, into the
  # $predict_method_used field on this function's return value (S3 default
  # flip). "exact" only if every joint fit that ran achieved exact;
  # "per_column" if any fell back, failed, or if no joint fit ran at all
  # (i.e. the whole baseline is per-column); NA_character_ only if every
  # joint fit used joint_solver = "rphylopars" (where the dichotomy does
  # not apply) and none fell back to in-house.
  predict_method_used_all <- character(0)

  if (use_threshold_joint) {
    jt <- if (em_iterations >= 1L) {
      fit_joint_threshold_baseline_em(data, tree, splits = splits,
                                       graph = graph,
                                       soft_aggregate = soft_aggregate,
                                       em_iterations = em_iterations,
                                       em_tol = em_tol,
                                       em_offdiag = em_offdiag,
                                       joint_solver = joint_solver, predict_method = predict_method,
                                       joint_refine_iter = joint_refine_iter,
                                       lambda_mode = lambda_mode_joint,
                                       lambda_fixed = lambda_fixed,
                                       predict_method_explicit = predict_method_explicit)
    } else {
      fit_joint_threshold_baseline(data, tree, splits = splits,
                                    graph = graph,
                                    soft_aggregate = soft_aggregate,
                                    joint_solver = joint_solver, predict_method = predict_method,
                                    joint_refine_iter = joint_refine_iter,
                                    lambda_mode = lambda_mode_joint,
                                    lambda_fixed = lambda_fixed,
                                    predict_method_explicit = predict_method_explicit)
    }
    if (!is.null(jt$predict_method_used) && !is.na(jt$predict_method_used)) {
      predict_method_used_all <- c(predict_method_used_all, jt$predict_method_used)
    }

    populated_cols <- integer(0)

    # Continuous-family passthrough (mu_liab on z-score scale).
    # Excludes binary (needs logit decode) and ordinal (needs threshold decode).
    #
    # S4 (was the "arc/lambda-per-type hybrid discard"): the joint liability
    # fit above (`jt`) always USES these continuous-family columns internally
    # to estimate the joint Sigma that the binary/ordinal posteriors condition
    # on -- that's unavoidable and desirable (it's the whole point of the
    # joint model). Previously its continuous-column mu/se were always fit at
    # lambda = 1 (neither Rphylopars nor the pre-S2 in-house solver had a
    # lambda argument), so under `lambda_mode != "fixed_1"` that lambda = 1
    # OUTPUT was discarded here and these columns fell through to `bm_cols`
    # for a separate lambda-aware per-column re-fit a few hundred lines down.
    # As of S2/S4 the in-house solver estimates a per-trait Pagel's lambda for
    # exactly these columns INSIDE the joint fit (`lambda_mode_joint` /
    # `lambda_fixed` threaded into `fit_joint_threshold_baseline()` above, via
    # its own `lambda_cols` mapping -- see that function's roxygen). The joint
    # output is therefore already lambda-aware under "estimate" too, so there
    # is nothing left to discard: keep it for both "fixed_1" and "estimate".
    # "cv"/"bayes" have no joint analogue (section 5 of the alignment note)
    # and still discard, falling through to the per-column path exactly as
    # before -- `lambda_mode_joint` forced the joint fit itself to "fixed_1"
    # for those two modes, but the OUTPUT is unused either way since this
    # `if` excludes them.
    cont_idx <- which(!(jt$liab_types %in% c("binary", "categorical", "ordinal")))
    if (lambda_mode %in% c("fixed_1", "estimate")) {
      for (idx in cont_idx) {
        col <- jt$liab_cols[idx]
        if (any(!is.na(jt$mu_liab[, idx]))) {
          mu[, col] <- jt$mu_liab[, idx]
          se[, col] <- jt$se_liab[, idx]
          populated_cols <- c(populated_cols, col)
          col_path[col] <- "threshold_joint"
        }
      }
      if (!is.null(jt$lambda_per_trait_fit)) {
        lambda_per_trait[names(jt$lambda_per_trait_fit)] <- jt$lambda_per_trait_fit
      }
      if (!is.null(jt$lambda_block) && is.finite(jt$lambda_block)) {
        lambda_block_out <- jt$lambda_block
      }
    }

    # Binary -> logit(P)
    bin_idx <- which(jt$liab_types == "binary")
    for (idx in bin_idx) {
      col <- jt$liab_cols[idx]
      if (all(is.na(jt$mu_liab[, idx]))) next
      dec <- decode_binary_liability(mu_liab = jt$mu_liab[, idx],
                                      se_liab = jt$se_liab[, idx])
      mu[, col] <- dec$mu_logit
      se[, col] <- 0
      populated_cols <- c(populated_cols, col)
      col_path[col] <- "threshold_joint"
    }

    # Ordinal -> z-scored integer class via threshold decode
    ord_idx <- which(jt$liab_types == "ordinal")
    ordinal_threshold_populated <- integer(0)   # subset for path-selection below
    for (idx in ord_idx) {
      col <- jt$liab_cols[idx]
      if (all(is.na(jt$mu_liab[, idx]))) next
      # Find the trait_map entry for this ordinal col
      tm_ord <- NULL
      for (tm in trait_map) {
        if (tm$type == "ordinal" && col %in% tm$latent_cols) {
          tm_ord <- tm; break
        }
      }
      if (is.null(tm_ord)) next
      dec <- decode_ordinal_liability(mu_liab = jt$mu_liab[, idx],
                                        se_liab = jt$se_liab[, idx],
                                        tm = tm_ord)
      mu[, col] <- dec$mu_z
      se[, col] <- 0
      populated_cols <- c(populated_cols, col)
      ordinal_threshold_populated <- c(ordinal_threshold_populated, col)
      col_path[col] <- "threshold_joint"
    }

    # S5c: attribute this joint fit's own exact/per_column verdict to every
    # column it actually populated (may be downgraded below for an ordinal
    # column that the path-selection block re-routes to bm_mvn/lp).
    if (!is.null(jt$predict_method_used) && !is.na(jt$predict_method_used)) {
      trait_route[populated_cols] <- jt$predict_method_used
    }

    bm_cols      <- setdiff(bm_cols,      populated_cols)
    binary_cols  <- setdiff(binary_cols,  populated_cols)
    ordinal_cols <- setdiff(ordinal_cols, populated_cols)

    # ---- Per-trait ordinal path selection (Opus #6, 2026-04-30) -----------
    # The threshold-joint path is theoretically more flexible than per-
    # column BM-via-MVN on z-scored integer class for ordinal traits, but
    # at small K (especially K=3, e.g. AVONET Migration) the K-1 thresholds
    # are pinned to a narrow band by phylopars EM and produce systematically
    # worse predictions than a Gaussian conditional MVN on z-scored
    # integers.  See `useful/MEMO_2026-04-29_phase6_migration_bisect.md`
    # for the bisect localising the regression to commit a541dbd.
    #
    # Rather than ship a `K <= 3 -> LP` heuristic, we compute BOTH paths
    # for each populated ordinal trait and pick the lower-val-MSE one
    # against the held-out val cells in `data$X_scaled`.  Single-obs only
    # for now (multi-obs would require species aggregation of the
    # alternative path; out of scope for this fix).
    ordinal_path_chosen <- character(0)
    if (length(ordinal_threshold_populated) > 0L &&
        !is.null(splits) && !multi_obs) {
      R_phy_local <- if (!is.null(graph) && !is.null(graph$R_phy)) {
        graph$R_phy[spp, spp]
      } else {
        phylo_cor_matrix(tree)[spp, spp]
      }
      # Linear-index decode helpers for splits$val_idx (integer indices
      # into the original n_obs x p_latent matrix).
      n_rows_sp <- nrow(data$X_scaled)
      val_idx   <- splits$val_idx
      val_col   <- ((val_idx - 1L) %/% n_rows_sp) + 1L
      val_row   <- ((val_idx - 1L) %% n_rows_sp) + 1L
      truth_full <- data$X_scaled

      for (col in ordinal_threshold_populated) {
        val_rows_j <- val_row[val_col == col]
        if (length(val_rows_j) == 0L) {
          ordinal_path_chosen[as.character(col)] <- "threshold_joint"
          next
        }
        # Threshold-joint prediction (currently in mu, possibly
        # species-level; for single-obs n_species == n_obs).
        tj_pred <- mu[, col]
        # BM-via-MVN alternative on the masked z-scored ordinal column.
        # Ordinal columns are documented to stay at lambda = 1 regardless
        # of `lambda_mode` (there is no discrete-trait analogue of Pagel's
        # lambda) -- see "Per-type lambda dispatch" in the roxygen Details
        # above. Do NOT substitute `bm_lambda` here: it tracks the
        # continuous-family setting and would leak an estimated lambda into
        # this candidate even though `lambda_per_trait` keeps reporting 1
        # for ordinal columns (Rose review, 2026-09-23).
        bm_res <- bm_impute_col(X[, col], R_phy_local, lambda = 1.0)
        # Val MSE for both paths.
        truth_j  <- truth_full[val_rows_j, col]
        finite_t <- is.finite(truth_j)
        if (!any(finite_t)) {
          ordinal_path_chosen[as.character(col)] <- "threshold_joint"
          next
        }
        tj_diff  <- tj_pred[val_rows_j[finite_t]] - truth_j[finite_t]
        bm_diff  <- bm_res$mu[val_rows_j[finite_t]] - truth_j[finite_t]
        tj_mse   <- if (any(is.finite(tj_diff))) {
                      mean(tj_diff[is.finite(tj_diff)]^2)
                    } else NA_real_
        bm_mse   <- if (any(is.finite(bm_diff))) {
                      mean(bm_diff[is.finite(bm_diff)]^2)
                    } else NA_real_

        # A third candidate, label propagation via K-class OVR (Phase F,
        # a3b89e6), was removed on 2026-09-26: it decoded the 0..K-1 ordinal
        # coding as 1..K and so dropped every lowest-class observation, and
        # with that fixed it lowered ordinal accuracy on the 18 core cells
        # (docs/dev-log/ordinal-route/).

        # Pick the lower finite val MSE of the two options.
        mses <- c(threshold_joint = tj_mse,
                  bm_mvn          = bm_mse)
        mses <- mses[is.finite(mses)]
        if (length(mses) == 0L) {
          ordinal_path_chosen[as.character(col)] <- "threshold_joint"
        } else {
          chosen <- names(mses)[which.min(mses)]
          if (chosen == "bm_mvn") {
            mu[, col] <- bm_res$mu
            se[, col] <- bm_res$se
            col_path[col] <- "per_column_bm"
            # S5c: this column's own outcome is now a per-column route,
            # regardless of what the threshold-joint fit as a whole did.
            trait_route[col] <- "per_column"
          } # else: threshold_joint already in mu/se, col_path, trait_route
          ordinal_path_chosen[as.character(col)] <- chosen
        }
      }
    }

  } else if (use_continuous_joint) {
    joint <- fit_joint_mvn_baseline(data, tree, splits = splits, graph = graph,
                                     soft_aggregate = soft_aggregate,
                                     joint_solver = joint_solver, predict_method = predict_method,
                                     joint_refine_iter = joint_refine_iter,
                                     lambda_mode = lambda_mode_joint,
                                     lambda_fixed = lambda_fixed,
                                     predict_method_explicit = predict_method_explicit)
    mu[, bm_cols] <- joint$mu[, bm_cols]
    se[, bm_cols] <- joint$se[, bm_cols]
    col_path[bm_cols] <- "joint_mvn"
    if (!is.null(joint$predict_method_used) && !is.na(joint$predict_method_used)) {
      trait_route[bm_cols] <- joint$predict_method_used
    }
    if (!is.null(joint$lambda_per_trait)) {
      lambda_per_trait[names(joint$lambda_per_trait)] <- joint$lambda_per_trait
    }
    if (!is.null(joint$lambda_block) && is.finite(joint$lambda_block)) {
      lambda_block_out <- joint$lambda_block
    }
    if (!is.null(joint$predict_method_used) && !is.na(joint$predict_method_used)) {
      predict_method_used_all <- c(predict_method_used_all, joint$predict_method_used)
    }
    bm_cols <- integer(0)
  }

  # ---- Categorical -> K independent OVR fits (Phase 6) -------------------
  # Each categorical trait gets K separate threshold-joint fits, one per
  # class. This sidesteps the rank-(K-1) phylopars instability of the
  # single-fit approach and fires regardless of whether the threshold-joint
  # / continuous-joint dispatchers above ran. If phylopars is unavailable
  # OR a fit fails for any reason, the per-trait result falls through to LP
  # below.
  #
  # arc/lambda-per-type: NOT gated on `force_per_column` -- categorical
  # traits have no lambda concept (OVR fits are threshold-joint calls,
  # always at lambda = 1), so `lambda_mode` never forces them onto label
  # propagation. See the dispatch comment above `use_threshold_joint`.
  if (length(cat_cols) > 0L && joint_mvn_available()) {
    for (tm in trait_map) {
      if (tm$type != "categorical") next
      k_cols <- tm$latent_cols
      if (!all(k_cols %in% cat_cols)) next  # already handled
      # Extract the trait name from the "<name>=<level>" column names
      col_name_1 <- colnames(data$X_scaled)[k_cols[1]]
      trait_name <- sub("=.*$", "", col_name_1)
      probs <- tryCatch(
        if (em_iterations >= 1L) {
          fit_ovr_categorical_fits_em(data, tree, trait_name = trait_name,
                                       splits = splits, graph = graph,
                                       soft_aggregate = soft_aggregate,
                                       em_iterations = em_iterations,
                                       em_tol = em_tol,
                                       joint_solver = joint_solver, predict_method = predict_method,
                                       joint_refine_iter = joint_refine_iter,
                                       predict_method_explicit = predict_method_explicit)
        } else {
          fit_ovr_categorical_fits(data, tree, trait_name = trait_name,
                                    splits = splits, graph = graph,
                                    soft_aggregate = soft_aggregate,
                                    joint_solver = joint_solver, predict_method = predict_method,
                                    joint_refine_iter = joint_refine_iter,
                                    predict_method_explicit = predict_method_explicit)
        },
        error = function(e) NULL
      )
      if (is.null(probs)) next
      # If OVR came back all-NA (every class's fit failed), leave for LP.
      if (all(is.na(probs))) next
      pmu <- attr(probs, "predict_method_used")
      if (!is.null(pmu) && !is.na(pmu)) {
        predict_method_used_all <- c(predict_method_used_all, pmu)
        # S5c: this categorical trait's OWN outcome, not the aggregate
        # across every other joint fit in this call.
        trait_route[k_cols] <- pmu
      }
      log_probs <- decode_ovr_categorical(probs)
      mu[, k_cols] <- log_probs
      se[, k_cols] <- 0
      col_path[k_cols] <- "ovr_categorical"
      cat_cols <- setdiff(cat_cols, k_cols)
    }
  }

  # ---- Internal BM imputation on BM-eligible columns -----------------------
  if (model == "OU") {
    message("OU not yet supported by the internal BM baseline; using BM. ",
            "Install Rphylopars for OU support.")
  }

  if (length(bm_cols) > 0) {
    # Retrieve or compute the phylogenetic correlation matrix
    if (!is.null(graph) && !is.null(graph$R_phy)) {
      R_phy <- graph$R_phy
    } else {
      R_phy <- phylo_cor_matrix(tree)
    }
    R_phy <- R_phy[spp, spp]

    # Aggregate multi-obs to species-level means for BM imputation
    if (multi_obs) {
      X_sp <- matrix(NA_real_, n_species, length(bm_cols))
      colnames(X_sp) <- colnames(X)[bm_cols]
      for (j in seq_along(bm_cols)) {
        col_vals <- X[, bm_cols[j]]
        sp_means <- tapply(col_vals, obs_spp, function(v) {
          v <- v[!is.na(v)]
          if (length(v) == 0L) NA_real_ else mean(v)
        })
        X_sp[match(names(sp_means), spp), j] <- as.numeric(sp_means)
      }
    } else {
      X_sp <- X[, bm_cols, drop = FALSE]
    }

    # ---- Covariate-aware design matrix (Fix G, 2026-04-25) ------------------
    # When `data$covariates` is non-NULL, the BM baseline becomes a GLS
    # regression on covariates: y = X*beta + u, u ~ MVN(0, sigma^2 * R).
    # This puts linear covariate effects into the BASELINE so the GNN's
    # delta only has to learn the nonlinear / interactive residuals.
    #
    # Without this, the GNN had to re-derive linear cov effects from
    # scratch through several non-linear layers + a regularised gate.
    # Empirically that converged to a much worse solution than direct
    # GLS regression (see useful/GNN_ARCHITECTURE_EXPLAINED.md).
    cov_design <- NULL
    if (!is.null(data$covariates)) {
      cov_mat <- as.matrix(data$covariates)
      if (multi_obs) {
        # Aggregate covariates to species level (mean across obs per species)
        cov_sp <- matrix(NA_real_, n_species, ncol(cov_mat))
        colnames(cov_sp) <- colnames(cov_mat)
        for (j in seq_len(ncol(cov_mat))) {
          sp_means <- tapply(cov_mat[, j], obs_spp, mean, na.rm = TRUE)
          cov_sp[match(names(sp_means), spp), j] <- as.numeric(sp_means)
        }
        cov_design <- cbind(intercept = 1, cov_sp)
      } else {
        cov_design <- cbind(intercept = 1, cov_mat)
      }
      # Replace any residual NAs (defensive): mean impute by column
      for (j in seq_len(ncol(cov_design))) {
        bad <- !is.finite(cov_design[, j])
        if (any(bad)) cov_design[bad, j] <- mean(cov_design[!bad, j])
      }
      # S4: `bm_impute_col_with_cov()` gained a `lambda` argument (S3,
      # feat/joint-lambda-default) that accepts a numeric scalar or
      # `"estimate"`, so `lambda_mode %in% c("fixed_1", "estimate")` (or a
      # numeric `lambda_fixed` value) now reaches it. It does NOT accept
      # `"cv"` / `"bayes"` -- those two modes still fall back to lambda = 1
      # with a warning here (not per column).
      if (lambda_mode %in% c("cv", "bayes") && is.null(lambda_fixed)) {
        warning(
          "lambda_mode = '", lambda_mode, "' is not supported by the ",
          "covariate-aware BM baseline; bm_impute_col_with_cov() only ",
          "accepts a numeric lambda or \"estimate\". Pagel's lambda is ",
          "ignored (fit at lambda = 1) for BM-eligible columns in this fit.",
          call. = FALSE
        )
      }
    }

    # Impute each BM-eligible column (covariate-aware when cov_design supplied)
    for (j in seq_along(bm_cols)) {
      col_name_j <- colnames(X_sp)[j]
      # lambda_fixed (spec 4.5 predict-time rebuild) overrides lambda_mode
      # entirely for continuous-family columns; falls back to 1 for a name
      # lambda_fixed doesn't cover.
      lam_j <- if (!is.null(lambda_fixed)) {
        val <- unname(lambda_fixed[col_name_j])
        if (is.na(val)) 1.0 else val
      } else if (!is.null(cov_design) && lambda_mode %in% c("cv", "bayes")) {
        1.0   # bm_impute_col_with_cov() has no cv/bayes concept
      } else {
        bm_lambda
      }
      if (is.null(cov_design)) {
        res_j <- bm_impute_col(X_sp[, j], R_phy, lambda = lam_j)
      } else {
        res_j <- bm_impute_col_with_cov(X_sp[, j], cov_design, R_phy, lambda = lam_j)
      }
      mu[, bm_cols[j]] <- res_j$mu
      se[, bm_cols[j]] <- res_j$se
      lambda_per_trait[bm_cols[j]] <- if (!is.null(res_j$lambda_hat)) {
        res_j$lambda_hat
      } else if (is.numeric(lam_j)) {
        lam_j
      } else if (identical(lam_j, "estimate") && sum(!is.na(X_sp[, j])) >= 10L) {
        # bm_impute_col() returns early for a fully observed column (nothing to
        # impute) without estimating lambda. Report the estimate anyway so the
        # stored lambda_per_trait describes the trait, not the missingness.
        ml_lambda_for_col(X_sp[, j], R_phy)
      } else {
        1
      }
      col_path[bm_cols[j]] <- if (bm_cols[j] %in% mp_cols) {
        "multi_proportion_bm"
      } else {
        "per_column_bm"
      }
    }
  }

  # ---- Binary baseline: phylogenetic label propagation -------------------
  for (tm in trait_map) {
    if (tm$type != "binary") next
    lc   <- tm$latent_cols
    # If the Phase 3 threshold-joint populated this col, skip LP.
    # binary_cols is the set of UNpopulated binary latent cols after the
    # threshold dispatch; if our `lc` is not in binary_cols, the joint
    # fit already handled it.
    if (!all(lc %in% binary_cols)) next
    col_path[lc] <- "label_propagation"

    # Get species-level observations
    if (multi_obs) {
      sp_vals <- tapply(X[, lc], obs_spp, function(v) {
        v <- v[!is.na(v)]
        if (length(v) == 0) NA_real_ else mean(v)
      })
      vals_species <- rep(NA_real_, n_species)
      names(vals_species) <- spp
      vals_species[names(sp_vals)] <- as.numeric(sp_vals)
    } else {
      vals_species <- X[, lc]
      names(vals_species) <- spp
    }

    observed <- !is.na(vals_species)
    if (sum(observed) == 0) {
      mu[, lc] <- logit(0.5)
      se[, lc] <- 0
      next
    }

    # Phylo-weighted probability for each species
    sim_obs <- sim_phylo[, observed, drop = FALSE]
    row_weights <- rowSums(sim_obs)
    row_weights[row_weights < 1e-10] <- 1e-10
    probs <- as.numeric(sim_obs %*% vals_species[observed]) / row_weights
    probs <- pmin(pmax(probs, 0.01), 0.99)  # clip for stability

    mu[, lc] <- logit(probs)
    se[, lc] <- 0
  }

  # ---- Categorical baseline: phylogenetic label propagation ---------------
  for (tm in trait_map) {
    if (tm$type != "categorical") next
    K    <- tm$n_latent
    lc   <- tm$latent_cols
    if (!all(lc %in% cat_cols)) next  # handled by threshold-joint path
    col_path[lc] <- "label_propagation"
    oh   <- X[, lc, drop = FALSE]  # n_obs x K one-hot (with NAs)

    # Get species-level one-hot observations
    if (multi_obs) {
      # Average one-hot within species (handles multiple obs)
      oh_species <- matrix(NA_real_, n_species, K)
      for (s in seq_len(n_species)) {
        rows <- which(obs_to_sp == s)
        obs_rows <- rows[complete.cases(oh[rows, , drop = FALSE])]
        if (length(obs_rows) > 0) {
          oh_species[s, ] <- colMeans(oh[obs_rows, , drop = FALSE])
        }
      }
    } else {
      oh_species <- oh
    }

    # Which species have observed values
    observed <- complete.cases(oh_species)
    if (sum(observed) == 0) {
      freqs <- rep(1 / K, K)
      log_freqs <- log(freqs)
      for (k in seq_len(K)) mu[, lc[k]] <- log_freqs[k]
      se[, lc] <- 0
      next
    }

    # Phylo-weighted category probabilities per species
    sim_obs <- sim_phylo[, observed, drop = FALSE]
    row_weights <- rowSums(sim_obs)
    row_weights[row_weights < 1e-10] <- 1e-10

    # Weighted category probs: (n_species x K)
    weighted_probs <- (sim_obs %*% oh_species[observed, , drop = FALSE]) /
      row_weights
    # Add small floor and renormalise
    weighted_probs <- pmax(weighted_probs, 1e-6)
    weighted_probs <- weighted_probs / rowSums(weighted_probs)

    for (k in seq_len(K)) {
      mu[, lc[k]] <- log(weighted_probs[, k])
    }
    se[, lc] <- 0
  }

  # ---- ZI count gate baseline: phylogenetic label propagation ---------------
  for (tm in trait_map) {
    if (tm$type != "zi_count") next
    lc_gate <- tm$latent_cols[1]
    # Phase 5: if threshold-joint handled this gate, skip LP
    if (!(lc_gate %in% binary_cols)) next
    col_path[lc_gate] <- "zi_gate_lp"

    # Get species-level gate values (0 = zero, 1 = non-zero)
    if (multi_obs) {
      sp_vals <- tapply(X[, lc_gate], obs_spp, function(v) {
        v <- v[!is.na(v)]
        if (length(v) == 0) NA_real_ else mean(v)
      })
      vals_species <- rep(NA_real_, n_species)
      names(vals_species) <- spp
      vals_species[names(sp_vals)] <- as.numeric(sp_vals)
    } else {
      vals_species <- X[, lc_gate]
      names(vals_species) <- spp
    }

    observed <- !is.na(vals_species)
    if (sum(observed) == 0) {
      mu[, lc_gate] <- logit(0.5)
      se[, lc_gate] <- 0
      next
    }

    # Phylo-weighted non-zero probability
    sim_obs <- sim_phylo[, observed, drop = FALSE]
    row_weights <- rowSums(sim_obs)
    row_weights[row_weights < 1e-10] <- 1e-10
    probs <- as.numeric(sim_obs %*% vals_species[observed]) / row_weights
    probs <- pmin(pmax(probs, 0.01), 0.99)

    mu[, lc_gate] <- logit(probs)
    se[, lc_gate] <- 0
  }

  # ---- Dispatch record (S1) ------------------------------------------------
  # One entry per trait_map entry, named after that trait/group. Multi-column
  # traits (categorical, zi_count, multi_proportion) get all their latent
  # columns from the same dispatch branch EXCEPT zi_count, whose gate and
  # magnitude columns can go through different branches (e.g. gate via
  # "zi_gate_lp", magnitude via "per_column_bm"); `path` reports the GATE
  # column's dispatch there (latent_cols[1]), matching the "zi_gate_lp" name.
  trait_names_path <- vapply(trait_map, function(tm) tm$name, character(1))
  path <- vapply(trait_map, function(tm) col_path[tm$latent_cols[1]],
                 character(1))
  names(path) <- trait_names_path

  # predict_method_used (S3 default flip): "exact" only if every joint fit
  # that ran achieved exact; "per_column" if any joint fit fell back or
  # failed, if no joint fit ran at all (a pure per-column baseline is
  # genuinely "per_column", whether or not the caller asked for "exact"),
  # or if every joint fit used the rphylopars solver (where the
  # exact/per_column dichotomy does not apply and fit_joint_solver()
  # reports NA_character_, filtered out of predict_method_used_all above).
  predict_method_used <- if (length(predict_method_used_all) == 0L) {
    "per_column"
  } else if (all(predict_method_used_all == "exact")) {
    "exact"
  } else {
    "per_column"
  }

  # S4 fix (root cause A, docs/dev-log/exact-default/S4-fixes-report.md):
  # attach the block value used internally as an attribute on the returned
  # $lambda_per_trait vector, so a caller that replays it as a later
  # fit_baseline(lambda_fixed = ...) call reproduces this fit's mu/se
  # exactly rather than approximating lambda_block via mean(lambda_vec)
  # (see .mvn_resolve_lambda()'s numeric-vector branch in
  # R/joint_mvn_solver.R). NA when no joint fit ran (nothing to carry).
  if (is.finite(lambda_block_out)) {
    attr(lambda_per_trait, "lambda_block") <- lambda_block_out
  }
  # S5c (rose-review.md, required change 3): report each trait's OWN
  # route, read off `trait_route` (tracked at every dispatch site above),
  # not the single aggregate `predict_method_used`. Previously this was
  # `rep(predict_method_used, ...)`, which meant one fallen-back joint fit
  # (e.g. a single OVR class) mislabelled every OTHER trait as
  # "per_column" too, even traits whose own joint fit achieved exact.  Any
  # column that never went through a tracked joint call (label
  # propagation, per-column BM, multi_proportion, zi magnitude constant)
  # is genuinely "per_column"; the final NA-fill below covers those.
  trait_route[is.na(trait_route)] <- "per_column"
  predict_method_by_trait <- vapply(
    trait_map, function(tm) trait_route[tm$latent_cols[1]], character(1))
  names(predict_method_by_trait) <- trait_names_path

  out <- list(mu = mu, se = se, path = path,
              lambda_per_trait = lambda_per_trait,
              lambda_block = lambda_block_out,
              lambda_mode = lambda_mode,
              predict_method_used = predict_method_used,
              predict_method_by_trait = predict_method_by_trait)
  if (exists("ordinal_path_chosen", inherits = FALSE) &&
      length(ordinal_path_chosen) > 0L) {
    out$ordinal_path_chosen <- ordinal_path_chosen
  }
  out
}

# ---- S5b: "auto" per-trait route selection --------------------------------
# See docs/dev-log/exact-default/S5b-route-choice-report.md and the
# `predict_method` / `predict_route` roxygen on fit_baseline() above.

#' @noRd
.pigauto_route_val_loss <- function(tm, val_row, species_row, X_truth, mu) {
  type <- tm$type
  if (type %in% c("continuous", "count", "ordinal", "proportion")) {
    col   <- tm$latent_cols[1]
    truth <- X_truth[val_row, col]
    pred  <- mu[species_row, col]
    ok    <- is.finite(truth) & is.finite(pred)
    n     <- sum(ok)
    loss  <- if (n > 0L) mean((pred[ok] - truth[ok])^2) else NA_real_
    return(list(n = n, loss = loss))
  }
  if (type == "binary") {
    col   <- tm$latent_cols[1]
    truth <- X_truth[val_row, col]
    prob  <- stats::plogis(mu[species_row, col])
    ok    <- is.finite(truth) & is.finite(prob)
    n     <- sum(ok)
    if (n == 0L) return(list(n = 0L, loss = NA_real_))
    p    <- pmin(pmax(prob[ok], 1e-6), 1 - 1e-6)
    y    <- truth[ok]
    loss <- mean(-(y * log(p) + (1 - y) * log1p(-p)))
    return(list(n = n, loss = loss))
  }
  if (type == "categorical") {
    k_cols <- tm$latent_cols
    oh     <- X_truth[val_row, k_cols, drop = FALSE]
    logp   <- mu[species_row, k_cols, drop = FALSE]
    ok_row <- stats::complete.cases(oh)
    n      <- sum(ok_row)
    if (n == 0L) return(list(n = 0L, loss = NA_real_))
    oh_ok   <- oh[ok_row, , drop = FALSE]
    logp_ok <- logp[ok_row, , drop = FALSE]
    true_k  <- max.col(oh_ok, ties.method = "first")
    ll <- logp_ok[cbind(seq_len(n), true_k)]
    ok2 <- is.finite(ll)
    n2 <- sum(ok2)
    loss <- if (n2 > 0L) mean(-ll[ok2]) else NA_real_
    return(list(n = n2, loss = loss))
  }
  if (type == "zi_count") {
    gate_col <- tm$latent_cols[1]
    mag_col  <- tm$latent_cols[2]

    truth_g <- X_truth[val_row, gate_col]
    prob_g  <- stats::plogis(mu[species_row, gate_col])
    ok_g    <- is.finite(truth_g) & is.finite(prob_g)
    n_g     <- sum(ok_g)
    loss_g  <- if (n_g > 0L) {
      p <- pmin(pmax(prob_g[ok_g], 1e-6), 1 - 1e-6)
      y <- truth_g[ok_g]
      mean(-(y * log(p) + (1 - y) * log1p(-p)))
    } else NA_real_

    truth_m <- X_truth[val_row, mag_col]
    pred_m  <- mu[species_row, mag_col]
    ok_m    <- is.finite(truth_m) & is.finite(pred_m)
    n_m     <- sum(ok_m)
    loss_m  <- if (n_m > 0L) mean((pred_m[ok_m] - truth_m[ok_m])^2) else NA_real_

    parts <- c(loss_g, loss_m)
    parts <- parts[is.finite(parts)]
    loss  <- if (length(parts) > 0L) sum(parts) else NA_real_
    return(list(n = n_g + n_m, loss = loss))
  }
  # multi_proportion (never joint-dispatched; both routes are identical) or
  # any other type with no exact/per_column distinction: nothing to compare.
  list(n = 0L, loss = NA_real_)
}

#' @noRd
.pigauto_choose_predict_route <- function(data, splits, fit_exact, fit_pc) {
  trait_map   <- data$trait_map
  trait_names <- vapply(trait_map, function(tm) tm$name, character(1))
  # S5c (rose-review.md, "Edge case B"): default each trait to whatever
  # fit_exact ITSELF actually ran for that trait (from its own, now
  # per-trait-accurate $predict_method_by_trait -- see required change 3),
  # not a hardcoded "exact" label. When fit_exact's own exact route was not
  # usable for a trait (too large, singular Sigma, etc.) it already fell
  # back to per_column internally; reporting "exact" here on a tie would
  # misdescribe what actually ran. Falls back to a plain "exact" default
  # only if fit_exact predates this field.
  route <- fit_exact$predict_method_by_trait
  if (is.null(route)) {
    route <- stats::setNames(rep("exact", length(trait_names)), trait_names)
  }

  X_truth <- data$X_scaled
  n_obs   <- nrow(X_truth)
  val_idx <- splits$val_idx
  val_col <- ((val_idx - 1L) %/% n_obs) + 1L
  val_row <- ((val_idx - 1L) %% n_obs) + 1L
  multi_obs   <- isTRUE(data$multi_obs)
  species_row <- if (multi_obs) data$obs_to_species[val_row] else val_row

  for (tm in trait_map) {
    keep <- val_col %in% tm$latent_cols
    if (!any(keep)) next
    vr <- val_row[keep]
    sr <- species_row[keep]
    res_exact <- .pigauto_route_val_loss(tm, vr, sr, X_truth, fit_exact$mu)
    res_pc    <- .pigauto_route_val_loss(tm, vr, sr, X_truth, fit_pc$mu)
    n_total <- max(res_exact$n, res_pc$n)
    if (n_total < 5L || !is.finite(res_exact$loss) || !is.finite(res_pc$loss) ||
        isTRUE(all.equal(res_exact$loss, res_pc$loss))) {
      next  # tie / insufficient evidence -> whatever fit_exact actually ran
    }
    route[[tm$name]] <- if (res_pc$loss < res_exact$loss) "per_column" else "exact"
  }
  route
}

# ---- per-trait route (half A) / calibration+conformal (half B) -----------
# validation split for "auto" (rose-review.md required change 1; unit and
# no-choice rules from rose-review-2.md B1/B2). Choosing the route that
# minimises validation error makes those same residuals optimistic, so a
# trait's validation UNITS are split into a ROUTE half (used only by
# .pigauto_choose_predict_route()) and a SCORE half (the cells
# fit_pigauto()/impute() restrict gate calibration + conformal scoring to,
# under "auto" only).
#
# The unit is the held-out trait row, not the latent cell: make_missing_splits()
# holds out all K cells of a categorical row (2 for zi_count), and
# calibrate_gates() keeps a row only by its first latent column, so a
# cell-level split leaked rows across halves and shrank calibration. In
# multi-obs data the unit is the species, because the route loss is scored
# on species-level baseline predictions.
#
# A trait is split only when (a) both candidate fits give different
# predictions on its validation cells (otherwise no route is chosen and
# halving would only widen its conformal interval: single-trait fits,
# multi_proportion, traits the joint path cannot reach), and (b) it has at
# least 38 units, so each half keeps the 19 a 95% split-conformal interval
# needs. Benchmark round 3 split every trait with >= 10 cells and lost up to
# 0.032 coverage at n <= 300. Otherwise all its cells both choose the route
# and calibrate.
#
#' @noRd
.pigauto_split_route_score <- function(val_idx, n_obs, trait_map, seed,
                                        mu_exact = NULL, mu_pc = NULL,
                                        unit_of_row = NULL) {
  trait_names <- vapply(trait_map, function(tm) tm$name, character(1))
  val_col <- ((val_idx - 1L) %/% n_obs) + 1L
  val_row <- ((val_idx - 1L) %% n_obs) + 1L
  val_unit <- if (is.null(unit_of_row)) val_row else unit_of_row[val_row]
  base_row <- if (is.null(unit_of_row)) val_row else unit_of_row[val_row]

  route_idx <- integer(0)
  score_idx <- integer(0)
  n_route <- stats::setNames(integer(length(trait_map)), trait_names)
  n_score <- stats::setNames(integer(length(trait_map)), trait_names)

  for (i in seq_along(trait_map)) {
    tm    <- trait_map[[i]]
    keep  <- val_col %in% tm$latent_cols
    idx_j <- val_idx[keep]
    units <- unique(val_unit[keep])
    n_j   <- length(units)

    differs <- TRUE
    if (!is.null(mu_exact) && !is.null(mu_pc)) {
      cells <- cbind(base_row[keep], val_col[keep])
      d <- abs(mu_exact[cells] - mu_pc[cells])
      differs <- any(is.finite(d) & d > 1e-10)
    }

    if (!differs || n_j < 38L) {
      route_idx <- c(route_idx, idx_j)
      score_idx <- c(score_idx, idx_j)
      n_route[[tm$name]] <- n_j
      n_score[[tm$name]] <- n_j
      next
    }
    if (!is.null(seed)) set.seed(seed + 29L + i)
    route_u <- sample(units, ceiling(n_j / 2))
    in_route <- val_unit[keep] %in% route_u
    route_idx <- c(route_idx, idx_j[in_route])
    score_idx <- c(score_idx, idx_j[!in_route])
    n_route[[tm$name]] <- length(route_u)
    n_score[[tm$name]] <- n_j - length(route_u)
  }

  list(route_idx = route_idx, score_idx = score_idx,
       n_route = n_route, n_score = n_score)
}

# S5c (required change 2): resolve which `lambda_fixed` vector a given
# route should use, consulting a `lambda_by_route` attribute when present.
# `.fit_baseline_auto()` / `.fit_baseline_route()` attach this attribute
# (a list with $exact / $per_column, each a full per-column lambda vector
# from the fit that actually produced that route) to the $lambda_per_trait
# they return, precisely so that a later `fit_baseline(..., lambda_fixed =
# bl$lambda_per_trait)` rebuild of a MIXED-route "auto"/`predict_route` fit
# replays each route's OWN lambda values -- not a single merged vector
# that would, for example, pin a per_column-routed discrete trait's lambda
# at 1 inside what would otherwise be a fresh "exact" candidate fit,
# perturbing that fit's Sigma estimate (and therefore every other exact-
# routed trait's mu, through shared cross-trait covariance) relative to
# the original. Falls back to using `lambda_fixed` as-is when the
# attribute is absent (a bare vector, or a single-route fit).
#
#' @noRd
.pigauto_lambda_fixed_for_route <- function(lambda_fixed, route) {
  if (is.null(lambda_fixed)) return(lambda_fixed)
  by_route <- attr(lambda_fixed, "lambda_by_route")
  if (!is.null(by_route) && !is.null(by_route[[route]])) {
    return(by_route[[route]])
  }
  lambda_fixed
}

#' @noRd
.fit_baseline_auto <- function(data, tree, splits, model, graph,
                                multi_obs_aggregation, lambda_mode, lambda_fixed,
                                em_iterations, em_tol, em_offdiag, joint_solver,
                                joint_refine_iter, seed = NULL) {
  trait_map   <- data$trait_map
  trait_names <- vapply(trait_map, function(tm) tm$name, character(1))

  no_val <- is.null(splits) || length(splits$val_idx) == 0L
  if (no_val) {
    out <- .fit_baseline_core(
      data = data, tree = tree, splits = splits, model = model, graph = graph,
      multi_obs_aggregation = multi_obs_aggregation, lambda_mode = lambda_mode,
      lambda_fixed = .pigauto_lambda_fixed_for_route(lambda_fixed, "exact"),
      em_iterations = em_iterations,
      em_tol = em_tol, em_offdiag = em_offdiag, joint_solver = joint_solver,
      predict_method = "exact", joint_refine_iter = joint_refine_iter,
      predict_method_explicit = FALSE)
    # No validation cells at all: "auto" cannot compare routes, so every
    # trait is REQUESTED at "exact" (spec: "No splits ... use exact for
    # every trait"). `predict_method_by_trait` keeps the core fit's own
    # actual per-trait outcome (it can still legitimately show
    # "per_column" for a trait whose exact route was not usable, e.g. too
    # large or a singular Sigma -- the exact-fallback machinery already
    # handles that below fit_baseline(), unrelated to "auto").
    out$predict_method_used <- "auto"
    out$score_val_idx <- if (is.null(splits)) integer(0) else splits$val_idx
    out$route_val_n <- stats::setNames(integer(length(trait_names)), trait_names)
    out$score_val_n <- stats::setNames(integer(length(trait_names)), trait_names)
    return(out)
  }

  fit_exact <- .fit_baseline_core(
    data = data, tree = tree, splits = splits, model = model, graph = graph,
    multi_obs_aggregation = multi_obs_aggregation, lambda_mode = lambda_mode,
    lambda_fixed = .pigauto_lambda_fixed_for_route(lambda_fixed, "exact"),
    em_iterations = em_iterations,
    em_tol = em_tol, em_offdiag = em_offdiag, joint_solver = joint_solver,
    predict_method = "exact", joint_refine_iter = joint_refine_iter,
    predict_method_explicit = FALSE)
  fit_pc <- .fit_baseline_core(
    data = data, tree = tree, splits = splits, model = model, graph = graph,
    multi_obs_aggregation = multi_obs_aggregation, lambda_mode = lambda_mode,
    lambda_fixed = .pigauto_lambda_fixed_for_route(lambda_fixed, "per_column"),
    em_iterations = em_iterations,
    em_tol = em_tol, em_offdiag = em_offdiag, joint_solver = joint_solver,
    predict_method = "per_column", joint_refine_iter = joint_refine_iter,
    predict_method_explicit = FALSE)

  # S5c required change 1: choose the route on a RANDOM HALF of each
  # trait's validation cells (`route_idx`); the other half (`score_idx`) is
  # reserved for fit_pigauto()/impute()'s gate calibration + conformal
  # scoring, so those residuals are not post-selected on the same cells
  # that picked the route.
  rsplit <- .pigauto_split_route_score(
    splits$val_idx, nrow(data$X_scaled), trait_map, seed,
    mu_exact = fit_exact$mu, mu_pc = fit_pc$mu,
    unit_of_row = if (isTRUE(data$multi_obs)) data$obs_to_species else NULL)
  splits_route <- splits
  splits_route$val_idx <- rsplit$route_idx
  route <- .pigauto_choose_predict_route(data, splits_route, fit_exact, fit_pc)

  mu   <- fit_exact$mu
  se   <- fit_exact$se
  path <- fit_exact$path
  lambda_per_trait <- fit_exact$lambda_per_trait
  for (tm in trait_map) {
    if (identical(route[[tm$name]], "per_column")) {
      mu[, tm$latent_cols]   <- fit_pc$mu[, tm$latent_cols, drop = FALSE]
      se[, tm$latent_cols]   <- fit_pc$se[, tm$latent_cols, drop = FALSE]
      path[tm$name]          <- fit_pc$path[tm$name]
      lambda_per_trait[tm$latent_cols] <- fit_pc$lambda_per_trait[tm$latent_cols]
    }
  }
  lb_attr <- attr(fit_exact$lambda_per_trait, "lambda_block")
  if (!is.null(lb_attr)) attr(lambda_per_trait, "lambda_block") <- lb_attr
  attr(lambda_per_trait, "lambda_by_route") <- list(
    exact = fit_exact$lambda_per_trait,
    per_column = fit_pc$lambda_per_trait
  )

  out <- list(mu = mu, se = se, path = path,
              lambda_per_trait = lambda_per_trait,
              lambda_block = fit_exact$lambda_block,
              lambda_mode = fit_exact$lambda_mode,
              predict_method_used = "auto",
              predict_method_by_trait = route,
              score_val_idx = rsplit$score_idx,
              route_val_n = rsplit$n_route,
              score_val_n = rsplit$n_score)
  if (!is.null(fit_exact$ordinal_path_chosen)) {
    out$ordinal_path_chosen <- fit_exact$ordinal_path_chosen
  }
  out
}

#' @noRd
.fit_baseline_route <- function(data, tree, splits, model, graph,
                                 multi_obs_aggregation, lambda_mode, lambda_fixed,
                                 em_iterations, em_tol, em_offdiag, joint_solver,
                                 joint_refine_iter, predict_route) {
  trait_map   <- data$trait_map
  trait_names <- vapply(trait_map, function(tm) tm$name, character(1))
  resolved <- stats::setNames(rep("exact", length(trait_names)), trait_names)
  known <- intersect(names(predict_route), trait_names)
  resolved[known] <- predict_route[known]
  # S5c (required change 4): warn once per call on a name that matches no
  # trait in this dataset -- previously silently ignored, so a misspelled
  # trait name gave no feedback at all.
  unknown_names <- setdiff(names(predict_route), trait_names)
  if (length(unknown_names) > 0L) {
    warning(
      "'predict_route' has name(s) that match no trait in this dataset ",
      "(ignored): ", paste(unknown_names, collapse = ", "), call. = FALSE
    )
  }

  needed <- unique(resolved)
  fits <- list()
  for (r in needed) {
    fits[[r]] <- .fit_baseline_core(
      data = data, tree = tree, splits = splits, model = model, graph = graph,
      multi_obs_aggregation = multi_obs_aggregation, lambda_mode = lambda_mode,
      lambda_fixed = .pigauto_lambda_fixed_for_route(lambda_fixed, r),
      em_iterations = em_iterations,
      em_tol = em_tol, em_offdiag = em_offdiag, joint_solver = joint_solver,
      predict_method = r, joint_refine_iter = joint_refine_iter,
      predict_method_explicit = FALSE)
  }

  if (length(needed) == 1L) {
    out <- fits[[needed]]
    out$predict_method_by_trait <- resolved
    out$predict_method_used <- "auto"
    return(out)
  }

  base <- fits[[needed[1]]]
  mu   <- base$mu
  se   <- base$se
  path <- base$path
  lambda_per_trait <- base$lambda_per_trait
  for (tm in trait_map) {
    r <- resolved[[tm$name]]
    f <- fits[[r]]
    mu[, tm$latent_cols] <- f$mu[, tm$latent_cols, drop = FALSE]
    se[, tm$latent_cols] <- f$se[, tm$latent_cols, drop = FALSE]
    path[tm$name] <- f$path[tm$name]
    lambda_per_trait[tm$latent_cols] <- f$lambda_per_trait[tm$latent_cols]
  }
  lb_attr <- attr(base$lambda_per_trait, "lambda_block")
  if (!is.null(lb_attr)) attr(lambda_per_trait, "lambda_block") <- lb_attr
  attr(lambda_per_trait, "lambda_by_route") <- lapply(fits, function(f) f$lambda_per_trait)

  out <- list(mu = mu, se = se, path = path,
              lambda_per_trait = lambda_per_trait,
              lambda_block = base$lambda_block,
              lambda_mode = base$lambda_mode,
              predict_method_used = "auto",
              predict_method_by_trait = resolved)
  if (!is.null(base$ordinal_path_chosen)) {
    out$ordinal_path_chosen <- base$ordinal_path_chosen
  }
  out
}
