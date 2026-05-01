#' Fit a pigauto model for trait imputation
#'
#' Trains a pigauto model: a gated ensemble of a phylogenetic baseline
#' and an attention-based graph neural network correction, implemented
#' as an internal torch module (\code{ResidualPhyloDAE}; "Residual" here
#' refers to the ResNet-style skip connections in the GNN layers, not
#' to a statistical residual). For continuous, count, and ordinal traits
#' the baseline is Brownian motion (phylogenetic correlation matrix); for binary
#' and categorical traits it is phylogenetic label propagation. Supports
#' all five trait types via a unified latent space.
#'
#' @details
#' **Blend formulation:**
#' The prediction is \eqn{\hat{x} = (1-r)\mu + r\delta}, where \eqn{\mu} is
#' the BM baseline, \eqn{\delta} is the model's direct prediction, and
#' \eqn{r = \sigma(\rho) \times \mathrm{cap}} is a per-column learnable gate
#' bounded in \eqn{(0, \mathrm{gate\_cap})}.  When \eqn{r = 0}, the prediction
#' collapses to the baseline.
#' The gate is regularised toward zero via the shrinkage penalty on
#' \eqn{\delta - \mu}, so the model defaults to the baseline unless the
#' GNN's correction demonstrably helps on the validation set.
#'
#' **Training objective** (per epoch):
#' \enumerate{
#'   \item A random subset of observed cells is corrupted with a learnable
#'     mask token.
#'   \item The model predicts \eqn{\delta} from graph context.
#'   \item Loss = type-specific reconstruction on corrupted cells +
#'     \code{lambda_shrink} * MSE(\eqn{\delta - \mu}) +
#'     \code{lambda_gate} * MSE(\eqn{r}).
#' }
#'
#' The gate penalty on \eqn{r} is necessary because when \eqn{\delta =
#' \mu} (the BM-optimal solution for observed cells), the reconstruction
#' and shrinkage losses both equal zero regardless of \eqn{r}, leaving no
#' gradient to close the gate.  The explicit penalty ensures gates default
#' toward zero when the GNN correction provides no benefit.
#'
#' **Type-specific losses:**
#' \describe{
#'   \item{continuous/count/ordinal}{MSE}
#'   \item{binary}{BCE with logits}
#'   \item{categorical}{cross-entropy over K latent columns}
#' }
#'
#' @param data object of class \code{"pigauto_data"}.
#' @param tree object of class \code{"phylo"}.
#' @param splits list (output of \code{\link{make_missing_splits}}) or
#'   \code{NULL}.
#' @param graph list (output of \code{\link{build_phylo_graph}}) or
#'   \code{NULL}.
#' @param baseline list (output of \code{\link{fit_baseline}}) or \code{NULL}.
#' @param hidden_dim integer. Hidden layer width (default \code{64}).
#' @param k_eigen integer. Number of spectral node features (default
#'   \code{8}).
#' @param dropout numeric. Dropout rate (default \code{0.10}).
#' @param lr numeric. AdamW learning rate (default \code{0.003}).
#' @param weight_decay numeric. AdamW weight decay (default \code{1e-4}).
#' @param epochs integer. Maximum training epochs (default \code{3000}).
#' @param n_gnn_layers integer. Number of graph message-passing layers
#'   (default \code{2}).  Each layer has its own learnable alpha gate,
#'   layer normalisation, and ResNet-style skip connection.
#' @param gate_cap numeric. Upper bound for the per-column blend gate
#'   (default \code{0.8}).  Safety comes from regularisation, not the cap.
#' @param use_attention logical. Use attention in the GNN layers (default
#'   \code{TRUE}).
#' @param use_transformer_blocks logical. Replace the legacy attention
#'   stack with pre-norm transformer-encoder blocks (multi-head attention
#'   + FFN + two residual skips). Default \code{TRUE}. Set \code{FALSE}
#'   to reconstruct pre-v0.9.0 fits (single-head attention with a
#'   learnable alpha gate per layer).
#' @param n_heads integer. Number of attention heads when
#'   \code{use_transformer_blocks = TRUE} (default \code{4}). Each head
#'   learns its own phylogenetic bandwidth (B2 rate-aware attention).
#' @param ffn_mult integer. Feed-forward width multiplier inside each
#'   transformer block (default \code{4}, giving \code{hidden_dim * 4}).
#' @param corruption_rate numeric. Final corruption fraction if
#'   \code{corruption_ramp > 0}; otherwise the fixed corruption rate per
#'   epoch (default \code{0.55}).
#' @param corruption_start numeric. Initial corruption fraction for the
#'   curriculum schedule (default \code{0.20}).  Ignored if
#'   \code{corruption_ramp = 0}.
#' @param corruption_ramp integer. Epochs over which corruption linearly
#'   ramps from \code{corruption_start} to \code{corruption_rate}
#'   (default \code{500}).  Set to \code{0} for fixed corruption.
#' @param refine_steps integer. Iterative refinement steps at inference
#'   (default \code{8}).
#' @param lambda_shrink numeric. Weight on the shrinkage penalty
#'   \code{||delta - baseline||^2} that keeps the GNN correction close
#'   to the phylogenetic baseline (default \code{0.03}).
#' @param lambda_gate numeric. Weight on the gate regularisation penalty
#'   that pushes learnable gates toward zero.  Prevents gates from staying
#'   open when the GNN provides no useful correction (default \code{0.01}).
#' @param warmup_epochs integer. Linear learning-rate warmup over the first
#'   N epochs (default \code{200}).  After warmup, a cosine schedule decays
#'   the LR to \code{1e-5}.
#' @param edge_dropout numeric. Fraction of adjacency edges randomly zeroed
#'   each training epoch for graph regularisation (default \code{0.1}).
#'   Set to \code{0} to disable.
#' @param eval_every integer. Evaluate on val every N epochs (default
#'   \code{100}).
#' @param patience integer. Early-stopping patience in eval cycles (default
#'   \code{10}).
#' @param clip_norm numeric. Gradient clip norm (default \code{1.0}).
#' @param conformal_method character. How the conformal prediction score
#'   is estimated from validation residuals.  \code{"split"} (default,
#'   backward-compatible) takes a single sample quantile; \code{"bootstrap"}
#'   takes \code{conformal_bootstrap_B} bootstrap resamples of the val
#'   residuals and averages the per-resample quantiles.  Empirically the
#'   bootstrap variant reduces the conformal-score variance across seeds
#'   by ~30% at small \code{n_val}, but on the simulation design used to
#'   evaluate it (n=150 species, 35\% trait_MAR) it did not translate
#'   into better 95\% coverage stability (coverage SD 0.094 vs 0.102
#'   across 10 seeds).  Ship as an opt-in experimental knob; defaults to
#'   \code{"split"}.
#' @param conformal_bootstrap_B integer. Bootstrap resamples used when
#'   \code{conformal_method = "bootstrap"}; default \code{500}.  Ignored
#'   otherwise.
#' @param conformal_split_val logical. Default \code{FALSE} (pre-2026-04-28
#'   single-set behaviour, retained because forcing the split regresses
#'   the AVONET300 / OVR-categorical / BIEN safety-floor smoke benches by
#'   2-26%). When \code{TRUE}, the validation set is split per latent
#'   column into a calibration half (used to pick the calibrated blend
#'   gate) and a conformal half (used to compute the conformal residual
#'   quantile). This restores split-conformal exchangeability — without
#'   the split, the gate is selected to minimise residual MSE on the very
#'   cells whose residuals drive the conformal quantile, producing
#'   systematic undercoverage (most visible at small \code{n_val}; see
#'   the \code{coverage_investigation} memo). Use this when accurate 95%
#'   coverage matters more than the bench-grade RMSE; the split only
#'   activates per-column when that column has at least
#'   \code{2 * min_val_cells} val cells (smaller columns silently fall
#'   back to the single-set path to keep \code{calibrate_gates()}'s
#'   half-A / half-B cross-check stable).
#' @param gate_method character. How the per-trait calibrated gate is
#'   chosen.  \code{"single_split"} (default, backward-compatible) runs
#'   the grid search on a single random half-A / half-B split of the val
#'   rows; \code{"median_splits"} repeats the whole procedure for
#'   \code{gate_splits_B} random splits and takes the median
#'   \code{best_g}.  \code{"cv_folds"} (2026-04-30) partitions val cells
#'   into \code{gate_cv_folds} (default 5) deterministic non-overlapping
#'   folds and runs the grid + half-B-verify procedure once per fold
#'   (training set = K-1 folds, held-out = remaining fold), taking the
#'   componentwise median of K winning weight vectors.  \code{"cv_folds"}
#'   uses larger training sets per split (K-1/K vs 1/2 in
#'   \code{median_splits}) and has a standard cross-validation
#'   interpretation, motivated by the open val→test drift observed on
#'   4/32 binary cells in the discrete-bench memo.
#'   \code{"median_splits"} slightly reduces gate
#'   bimodality at small \code{n_val} (SD 0.406 → 0.360 across 10 seeds
#'   on the evaluation sim) with a small coverage-SD improvement
#'   (0.094 → 0.086).  Negligible runtime cost (B × cheap grid searches).
#' @param gate_splits_B integer. Random splits used when
#'   \code{gate_method = "median_splits"}; default \code{31} (odd so the
#'   median is well-defined).
#' @param gate_cv_folds integer. Number of CV folds when
#'   \code{gate_method = "cv_folds"}; default \code{5}, must be
#'   \code{>= 2}.  Capped at \code{n_val} per trait so each fold has at
#'   least 1 cell.  When effective K \code{< 2} (e.g. \code{n_val = 1}),
#'   the code falls back to a single split.
#' @param safety_floor logical. When \code{TRUE} (default), post-training
#'   calibration searches a 3-way simplex \code{r_BM * BM + r_GNN * GNN
#'   + r_MEAN * MEAN} so the grand mean is always in the candidate set,
#'   guaranteeing \code{pigauto_val_RMSE <= mean_val_RMSE} by
#'   construction. When \code{FALSE}, the v0.9.1 1-D calibration is used
#'   exactly (\code{r_MEAN = 0}).
#' @param phylo_signal_gate logical. When \code{TRUE} (default since
#'   v0.9.1.9003), compute per-trait Pagel's \eqn{\lambda} on
#'   training-observed cells before fitting; for traits with
#'   \code{lambda < phylo_signal_threshold}, force
#'   \code{(r_cal_bm = 0, r_cal_gnn = 0, r_cal_mean = 1)} directly
#'   and skip BM + GNN training on those traits.  Requires the
#'   \code{phytools} package.  Falls back to safety-floor-only
#'   behaviour (\code{phylo_signal_gate = FALSE} effective) when
#'   \code{phytools} is absent.
#' @param phylo_signal_threshold numeric, default \code{0.2}.  Traits
#'   with Pagel's \eqn{\lambda} below this value are routed to the
#'   grand-mean corner of the safety-floor simplex.
#' @param phylo_signal_method character, currently only \code{"lambda"}
#'   is fully implemented.  Reserved \code{"blomberg_k"} path returns
#'   Blomberg's K via \code{phytools::phylosig()} but uses the same
#'   threshold --- which is NOT dimensionally comparable; users
#'   selecting K must supply a K-appropriate threshold.
#' @param min_val_cells integer. Warn at fit time if any trait has fewer
#'   than \code{min_val_cells} validation cells available for gate
#'   calibration and conformal-score estimation.  Default \code{10}: the
#'   floor of pathological territory, where the conformal quantile
#'   collapses to \code{max(val_residuals)} and gate calibration becomes
#'   essentially a coin flip between \code{0} and \code{gate_cap}.
#'   Recommended operational target is \code{n_val >= 20-30} per trait;
#'   achieve this by increasing \code{missing_frac} or collecting more
#'   species.  See \strong{Calibration at small n} below.
#' @param verbose logical. Print training progress (default \code{TRUE}).
#' @param seed integer. Random seed (default \code{1}).
#' @section Calibration at small n:
#'
#' pigauto's 95\% intervals are \emph{conformal}, not parametric: the
#' interval half-width for each trait is the empirical \code{(1 - alpha)}
#' quantile of \eqn{|y - \hat y|} on held-out validation cells.  When the
#' number of validation cells per trait (\code{n_val}) is large
#' (\eqn{\gtrsim 30}), split-conformal gives near-exact marginal coverage
#' under mild exchangeability assumptions.
#'
#' At small \code{n_val} (\eqn{<} 20, and especially \eqn{<} 10) two
#' things degrade at once:
#'
#' \itemize{
#'   \item The conformal quantile clamps to \code{max(residuals)} because
#'     the required quantile level exceeds 1.  The score is the single
#'     largest val residual and has substantial sampling variance across
#'     fits.
#'   \item The gate calibration's half-A / half-B split ends up with only
#'     a few cells per half, so the grid-search winner is essentially
#'     random between \code{0} and \code{gate_cap}.
#' }
#'
#' Empirically (pigauto simulation harness, n=150 species, 35\%
#' trait_MAR, 10 random seeds), default behaviour produces 92\% mean
#' coverage with per-fit coverage ranging \code{[0.73, 1.00]}.  The mean
#' is close to the 95\% target; the per-fit variance is the issue.
#' \code{gate_method = "median_splits"} and \code{conformal_method =
#' "bootstrap"} reduce their respective estimator variances but
#' empirically do not meaningfully narrow the coverage distribution.
#' \strong{Treat the 95\% intervals as approximate in the small-\code{n}
#' regime}; increase \code{missing_frac} (more held-out data for
#' calibration) or collect more species if tight per-fit coverage is
#' required.
#'
#' The \code{min_val_cells} warning fires at fit time whenever any trait
#' falls below this threshold, so you know when to interpret the
#' intervals cautiously.
#'
#' @return An object of class \code{"pigauto_fit"}.
#' @importFrom torch torch_tensor torch_float torch_bool torch_long
#' @importFrom torch optim_adam with_no_grad nn_utils_clip_grad_norm_
#' @importFrom torch nnf_mse_loss torch_rand_like torch_where torch_zeros_like
#' @export
fit_pigauto <- function(
    data,
    tree,
    splits            = NULL,
    graph             = NULL,
    baseline          = NULL,
    hidden_dim        = 64L,
    k_eigen           = "auto",
    n_gnn_layers      = 2L,
    gate_cap          = 0.8,
    use_attention     = TRUE,
    use_transformer_blocks = TRUE,
    n_heads           = 4L,
    ffn_mult          = 4L,
    dropout           = 0.10,
    lr                = 0.003,
    weight_decay      = 1e-4,
    epochs            = 3000L,
    corruption_rate   = 0.55,
    corruption_start  = 0.20,
    corruption_ramp   = 500L,
    refine_steps      = 8L,
    lambda_shrink     = 0.03,
    lambda_gate       = 0.01,
    warmup_epochs     = 200L,
    edge_dropout      = 0.1,
    eval_every        = 100L,
    patience          = 10L,
    clip_norm         = 1.0,
    conformal_method  = c("split", "bootstrap"),
    conformal_bootstrap_B = 500L,
    conformal_split_val = FALSE,
    gate_method       = c("single_split", "median_splits", "cv_folds"),
    gate_splits_B     = 31L,
    gate_cv_folds     = 5L,
    safety_floor      = TRUE,
    phylo_signal_gate = TRUE,
    phylo_signal_threshold = 0.2,
    phylo_signal_method = c("lambda", "blomberg_k"),
    min_val_cells     = 20L,
    verbose           = TRUE,
    seed              = 1L
) {
  conformal_method    <- match.arg(conformal_method)
  gate_method         <- match.arg(gate_method)
  phylo_signal_method <- match.arg(phylo_signal_method)
  if (!inherits(data, "pigauto_data")) {
    stop("'data' must be a pigauto_data object.")
  }

  set.seed(seed)
  torch::torch_manual_seed(seed)

  device <- get_device()
  if (verbose) message("Using device: ", as.character(device))

  # ---- Phylogenetic-signal gate -----------------------------------------------
  # Compute per-trait Pagel's lambda on training-observed cells. Traits with
  # lambda < phylo_signal_threshold are routed to the grand-mean corner (0,0,1)
  # of the safety-floor simplex after calibration.
  if (phylo_signal_gate) {
    phylo_signal_per_trait <- compute_phylo_signal_per_trait(
      data = data, tree = tree,
      method = phylo_signal_method, min_tips = 20L)
  } else {
    phylo_signal_per_trait <- rep(NA_real_, length(data$trait_map))
    names(phylo_signal_per_trait) <- vapply(data$trait_map,
                                              function(tm) tm$name,
                                              character(1L))
  }
  phylo_gate_triggered <- phylo_signal_gate &
    !is.na(phylo_signal_per_trait) &
    (phylo_signal_per_trait < phylo_signal_threshold)
  names(phylo_gate_triggered) <- names(phylo_signal_per_trait)
  # Build a latent-column-level index vector for downstream overrides.
  gated_latent_cols <- integer(0L)
  if (any(phylo_gate_triggered)) {
    gated_trait_names <- names(phylo_gate_triggered)[phylo_gate_triggered]
    for (tm in data$trait_map) {
      if (tm$name %in% gated_trait_names) {
        gated_latent_cols <- c(gated_latent_cols, tm$latent_cols)
      }
    }
  }

  # ---- Graph ----------------------------------------------------------------
  if (is.null(graph)) {
    if (verbose) message("Computing phylogenetic graph...")
    graph <- build_phylo_graph(tree, k_eigen = k_eigen)
  }
  # Use the actual k_eigen from the graph (may differ from input if "auto")
  k_eigen <- ncol(graph$coords)

  # ---- Baseline -------------------------------------------------------------
  if (is.null(baseline)) {
    if (verbose) message("Fitting baseline...")
    # Pass graph through so fit_baseline can reuse graph$D instead of
    # calling ape::cophenetic.phylo() a second time on the same tree.
    baseline <- fit_baseline(data, tree, splits = splits, graph = graph)
  }

  # ---- Trait map ------------------------------------------------------------
  trait_map <- data$trait_map
  has_trait_map <- !is.null(trait_map)

  # ---- Multi-obs handling --------------------------------------------------
  multi_obs  <- isTRUE(data$multi_obs)
  n_obs      <- nrow(data$X_scaled)
  n_species  <- data$n_species %||% n_obs
  obs_to_sp  <- data$obs_to_species  # NULL when single-obs

  # ---- Data preparation ----------------------------------------------------
  n <- n_obs                           # n = n_obs (rows in X)
  p <- ncol(data$X_scaled)

  X_truth <- data$X_scaled             # n_obs x p_latent

  # Baseline mu is at species level (n_species x p).
  # For multi-obs, expand to observation level.
  MU_species <- baseline$mu            # n_species x p_latent
  if (multi_obs) {
    MU <- MU_species[obs_to_sp, , drop = FALSE]
    rownames(MU) <- NULL
  } else {
    MU <- MU_species
  }

  # Fill missing cells with baseline predictions
  X_fill <- X_truth
  X_fill[is.na(X_fill)] <- MU[is.na(X_fill)]

  DELTA_tgt <- X_fill - MU            # observed - baseline (legacy
                                       # tensor t_DELTA below; kept for
                                       # backward compatibility, not
                                       # used as a training target)

  # Observed mask (TRUE = actually observed, not just filled)
  M_obs_mat <- !is.na(X_truth)
  if (!is.null(splits)) M_obs_mat[c(splits$val_idx, splits$test_idx)] <- FALSE

  # Held-out-masked input: like X_fill but with val/test cells replaced
  # with the baseline prediction.  This is what the model SHOULD see at
  # val/test evaluation time and at gate calibration time, so that it
  # does not have direct access to the held-out truth values.  Without
  # this masking, val cells' true values leak into the model input and
  # the calibration signal becomes meaningless.
  X_fill_heldout <- X_fill
  if (!is.null(splits)) {
    hold_idx <- c(splits$val_idx, splits$test_idx)
    X_fill_heldout[hold_idx] <- MU[hold_idx]
  }

  # ---- Trait-level corruption mask expansion --------------------------------
  # For categorical traits: corrupt all K latent columns together
  if (has_trait_map) {
    # Build mapping from latent column -> trait index (for grouped corruption)
    latent_to_trait <- integer(p)
    for (i in seq_along(trait_map)) {
      latent_to_trait[trait_map[[i]]$latent_cols] <- i
    }
    n_orig_traits <- length(trait_map)
  }

  # ---- Tensors ---------------------------------------------------------------
  t_X      <- torch::torch_tensor(X_fill,    dtype = torch::torch_float(),
                                  device = device)
  t_X_eval <- torch::torch_tensor(X_fill_heldout, dtype = torch::torch_float(),
                                  device = device)
  t_MU     <- torch::torch_tensor(MU,        dtype = torch::torch_float(),
                                  device = device)
  t_DELTA  <- torch::torch_tensor(DELTA_tgt, dtype = torch::torch_float(),
                                  device = device)
  t_M_obs  <- torch::torch_tensor(M_obs_mat, dtype = torch::torch_bool(),
                                  device = device)
  t_adj    <- torch::torch_tensor(graph$adj,    dtype = torch::torch_float(),
                                  device = device)
  t_coords <- torch::torch_tensor(graph$coords, dtype = torch::torch_float(),
                                  device = device)

  # Squared cophenetic distances for B2 rate-aware attention.
  # graph$D_sq is set by build_phylo_graph() and is NOT freed in impute.R's
  # cleanup block (only D and R_phy are freed there), so it is available here.
  # Falls back to NULL if the graph was built before B2.1 (old cache).
  t_D_sq <- if (!is.null(graph$D_sq)) {
    torch::torch_tensor(graph$D_sq, dtype = torch::torch_float(),
                        device = device)
  } else {
    NULL
  }

  # Drop the cophenetic distance matrix from fit_pigauto's local view of
  # the graph now that we have the torch tensors we need. This keeps the
  # fit object slim on disk (predict() only reads $adj and $coords), but
  # does NOT free memory at the caller's scope -- R's copy-on-modify
  # means the caller still holds the original list. Callers that want
  # to free memory during training should set graph$D <- NULL before
  # calling fit_pigauto(); impute() already does this, and the
  # bundled benchmark scripts do the same.
  graph$D <- NULL

  # Observation-to-species mapping for multi-obs (1-indexed for torch-R)
  if (multi_obs) {
    t_obs_to_sp <- torch::torch_tensor(
      as.integer(obs_to_sp),
      dtype = torch::torch_long(), device = device
    )
  } else {
    t_obs_to_sp <- NULL
  }

  # Val / test masks
  if (!is.null(splits)) {
    val_mat <- matrix(FALSE, n, p); val_mat[splits$val_idx]  <- TRUE
    t_val   <- torch::torch_tensor(val_mat, dtype = torch::torch_bool(),
                                   device = device)
    t_truth <- torch::torch_tensor(X_truth, dtype = torch::torch_float(),
                                   device = device)
    # Replace NaN in truth with 0 (won't affect masked evaluation)
    t_truth_safe <- t_truth$clone()
    t_truth_safe[t_truth$isnan()] <- 0
    has_val <- TRUE
  } else {
    has_val <- FALSE
  }

  # ---- Covariates (environmental conditioners) ------------------------------
  has_covariates <- !is.null(data$covariates)
  n_cov_cols     <- if (has_covariates) ncol(data$covariates) else 0L
  t_covariates   <- NULL
  if (has_covariates) {
    t_covariates <- torch::torch_tensor(
      data$covariates, dtype = torch::torch_float(), device = device
    )
  }

  # ---- Model ----------------------------------------------------------------
  cov_dim <- p + 1L + n_cov_cols
  # n_user_cov tells the model how many user covariates exist (the last
  # n_cov_cols columns of the covs tensor). The model's obs_refine MLP
  # re-injects these covariates after species-level GNN message passing
  # via a residual connection so the GNN's delta is conditioned on the
  # raw covariates rather than only on the species-level hidden state
  # (which has been propagated through phylogeny-only graph layers).
  #
  # Historical note (pre-2026-04-25): n_user_cov was forced to 0 in
  # single-obs mode, which meant user covariates entered the model only
  # via the input encoder and were diluted through the GNN layers. The
  # GNN had no path to extract nonlinear covariate-effect signal beyond
  # what a single linear projection could carry. Fixed by feeding the
  # raw user covariates back into the obs_refine MLP in single-obs mode
  # too, so single-obs imputation can finally exploit nonlinear
  # covariate structure.
  n_user_cov <- n_cov_cols
  model <- ResidualPhyloDAE(
    input_dim              = p,
    hidden_dim             = as.integer(hidden_dim),
    coord_dim              = as.integer(k_eigen),
    cov_dim                = as.integer(cov_dim),
    per_column_rs          = TRUE,
    n_gnn_layers           = as.integer(n_gnn_layers),
    gate_cap               = gate_cap,
    use_attention          = use_attention,
    n_user_cov             = as.integer(n_user_cov),
    dropout                = dropout,
    use_transformer_blocks = use_transformer_blocks,
    n_heads                = as.integer(n_heads),
    ffn_mult               = as.integer(ffn_mult)
  )

  # Type-aware gate init: set res_raw so effective gate starts at a
  # non-trivial value for all trait types.  Discrete traits previously
  # started near 0 (disc_target = 0.001), which killed gradient flow
  # through the sigmoid and prevented the GNN from learning discrete
  # corrections.  Starting at 0.10 gives a healthy gradient in both
  # directions.  logit(target / gate_cap) gives the correct res_raw.
  cont_target <- 0.135
  disc_target <- 0.10
  cont_raw <- log(cont_target / gate_cap / (1 - cont_target / gate_cap))
  disc_raw <- log(disc_target / gate_cap / (1 - disc_target / gate_cap))

  if (has_trait_map) {
    init_vals <- rep(cont_raw, p)
    for (tm in trait_map) {
      if (tm$type %in% c("binary", "categorical", "ordinal")) {
        init_vals[tm$latent_cols] <- disc_raw
      } else if (tm$type == "zi_count") {
        # Gate column (col 1): like binary (start closed)
        init_vals[tm$latent_cols[1]] <- disc_raw
        # Magnitude column (col 2): like count (start slightly open)
        init_vals[tm$latent_cols[2]] <- cont_raw
      }
      # proportion: uses cont_raw (default), no special case needed
    }
  } else {
    init_vals <- rep(cont_raw, p)
  }
  torch::with_no_grad({
    model$res_raw$copy_(
      torch::torch_tensor(init_vals, dtype = torch::torch_float())
    )
  })

  model$to(device = device)
  opt <- torch::optim_adamw(model$parameters, lr = lr,
                             weight_decay = weight_decay)

  gpu_mem_checkpoint("after model$to(device) + optimizer")

  # ---- Training loop --------------------------------------------------------
  best_val      <- Inf
  best_state    <- NULL
  patience_left <- patience
  history       <- data.frame(
    epoch = integer(0), loss_rec = double(0),
    loss_shrink = double(0), loss_gate = double(0),
    val_loss = double(0), lr = double(0)
  )

  for (epoch in seq_len(epochs)) {
    model$train()
    opt$zero_grad()

    # ---- Learning rate schedule (warmup + cosine decay) --------------------
    scheduled_lr <- cosine_lr(epoch, warmup_epochs, epochs, lr)
    opt$param_groups[[1]]$lr <- scheduled_lr

    # ---- Adaptive corruption rate ------------------------------------------
    if (corruption_ramp > 0L) {
      eff_corruption <- corruption_start +
        (corruption_rate - corruption_start) *
        min(epoch / corruption_ramp, 1.0)
    } else {
      eff_corruption <- corruption_rate
    }

    # ---- Edge dropout (graph regularisation) --------------------------------
    if (edge_dropout > 0 && model$training) {
      edge_mask <- torch::torch_rand_like(t_adj) > edge_dropout
      t_adj_train <- t_adj * edge_mask$to(dtype = torch::torch_float())
      # Renormalise rows to maintain expected message scale
      row_sums <- t_adj_train$sum(dim = 2L, keepdim = TRUE)$clamp(min = 1e-8)
      orig_sums <- t_adj$sum(dim = 2L, keepdim = TRUE)$clamp(min = 1e-8)
      t_adj_train <- t_adj_train * (orig_sums / row_sums)
    } else {
      t_adj_train <- t_adj
    }

    # ---- Corruption masking -------------------------------------------------
    if (has_trait_map) {
      # Corrupt at trait level, then expand to latent columns
      u_trait <- matrix(stats::runif(n * n_orig_traits), n, n_orig_traits)
      corrupt_trait <- u_trait < eff_corruption
      # Only corrupt where observed
      obs_trait <- matrix(FALSE, n, n_orig_traits)
      for (i in seq_along(trait_map)) {
        lc <- trait_map[[i]]$latent_cols
        obs_trait[, i] <- as.logical(M_obs_mat[, lc[1]])
      }
      corrupt_trait <- corrupt_trait & obs_trait
      # Expand to latent columns
      corrupt_latent <- matrix(FALSE, n, p)
      for (i in seq_along(trait_map)) {
        for (lc in trait_map[[i]]$latent_cols) {
          corrupt_latent[, lc] <- corrupt_trait[, i]
        }
      }
      masked_bool <- torch::torch_tensor(corrupt_latent, dtype = torch::torch_bool(),
                                          device = device)
    } else {
      u           <- torch::torch_rand_like(t_X)
      masked_bool <- (u < eff_corruption) & t_M_obs
    }

    if (as.numeric(masked_bool$sum()$cpu()$item()) == 0L) next

    Mask_t  <- masked_bool$to(dtype = torch::torch_float())
    mask_ind <- Mask_t$any(dim = 2L, keepdim = TRUE)$to(
      dtype = torch::torch_float()
    )

    covs_t <- make_covs_tensor(t_MU, mask_ind, t_covariates)

    tok   <- model$mask_token$expand(c(n, p))
    X_in  <- torch::torch_where(masked_bool, tok, t_X)

    out   <- model(X_in, t_coords, covs_t, t_adj_train, t_obs_to_sp,
                   D_sq = t_D_sq)
    delta <- out$delta
    rs    <- out$rs

    # Direct prediction loss: train delta to predict the truth directly.
    # This gives delta a full-strength gradient signal on every corrupted
    # cell, unlike the earlier blend formulation where the effective
    # gradient on delta was scaled by rs (approx 0.1 to 0.2), causing
    # delta to stay close to the baseline and calibration to fall back
    # to gate = 0.
    if (has_trait_map) {
      loss_rec <- compute_mixed_loss(delta, t_X, masked_bool, trait_map)
    } else {
      loss_rec <- torch::nnf_mse_loss(
        delta[masked_bool], t_X[masked_bool]
      )
    }

    # Auxiliary blend loss: (1-rs)*mu + rs*delta against the truth.
    # This gives rs a gradient signal so that the model can learn which
    # columns benefit from the GNN, but with a small weight so it doesn't
    # dominate the direct supervision on delta.  Post-hoc calibration on
    # the validation set still has the final say over the gate.
    pred_blend <- (1 - rs) * t_MU + rs * delta
    if (has_trait_map) {
      loss_blend <- compute_mixed_loss(pred_blend, t_X, masked_bool, trait_map)
    } else {
      loss_blend <- torch::nnf_mse_loss(
        pred_blend[masked_bool], t_X[masked_bool]
      )
    }

    # Shrinkage: keep delta from drifting to extreme values on cells with
    # very weak supervision.  Lighter than before because direct
    # supervision handles most of the regularisation.
    deviation <- (delta - t_MU)[masked_bool]
    loss_shrink <- torch::torch_mean(deviation$pow(2))

    # Gate regularisation pushes rs toward zero by default, so the model
    # only opens a gate when the direct delta signal is clearly helpful.
    loss_gate <- torch::torch_mean(rs$pow(2))

    loss <- loss_rec + 0.1 * loss_blend +
            lambda_shrink * loss_shrink + lambda_gate * loss_gate
    loss$backward()

    if (!all_grads_finite(model$parameters)) {
      opt$zero_grad(); next
    }
    torch::nn_utils_clip_grad_norm_(model$parameters, max_norm = clip_norm)
    opt$step()

    # ---- Evaluation ---------------------------------------------------------
    if (epoch %% eval_every == 0L) {
      val_loss_val <- NA_real_
      if (has_val) {
        model$eval()
        torch::with_no_grad({
          mask_ind0 <- torch::torch_zeros(c(n, 1L), device = device)
          covs0     <- make_covs_tensor(t_MU, mask_ind0, t_covariates)
          # Use t_X_eval (val/test cells replaced with baseline) so
          # that held-out truth does not leak into the model input.
          out0      <- model(t_X_eval, t_coords, covs0, t_adj, t_obs_to_sp,
                             D_sq = t_D_sq)
          pred0     <- (1 - out0$rs) * t_MU + out0$rs * out0$delta
        })
        if (has_trait_map) {
          val_loss_val <- composite_val_loss(pred0, t_truth_safe, t_val,
                                             trait_map)
        } else {
          val_loss_val <- as.numeric(
            rmse_torch(pred0, t_truth_safe, t_val)$cpu()$item()
          )
        }
      }

      history <- rbind(history, data.frame(
        epoch       = epoch,
        loss_rec    = as.numeric(loss_rec$item()),
        loss_shrink = as.numeric(loss_shrink$item()),
        loss_gate   = as.numeric(loss_gate$item()),
        val_loss    = if (is.na(val_loss_val)) NA_real_ else val_loss_val,
        lr          = scheduled_lr
      ))

      if (verbose) {
        msg <- sprintf(
          "epoch %4d | rec %.4f | shrink %.4f | gate %.4f | lr %.1e",
          epoch, as.numeric(loss_rec$item()),
          as.numeric(loss_shrink$item()),
          as.numeric(loss_gate$item()),
          scheduled_lr
        )
        if (!is.na(val_loss_val)) {
          msg <- paste0(msg, sprintf(" | val %.4f | best %.4f | pat %d",
                                     val_loss_val, best_val, patience_left))
        }
        message(msg)
      }

      if (!is.na(val_loss_val)) {
        if (val_loss_val + 1e-6 < best_val) {
          best_val      <- val_loss_val
          best_state    <- model$state_dict()
          patience_left <- patience
        } else {
          patience_left <- patience_left - 1L
        }
        if (patience_left <= 0L) {
          if (verbose) message("Early stopping at epoch ", epoch)
          break
        }
      }
    }
  }

  gpu_mem_checkpoint("end of training loop (before restoring best state)")
  if (!is.null(best_state)) model$load_state_dict(best_state)
  gpu_mem_checkpoint("after model$load_state_dict(best_state)")

  # ---- Post-training gate calibration (validation set) --------------------
  calibrated_gates <- NULL
  calibrated_gates_list <- NULL
  mean_baseline_per_col <- NULL
  if (has_val && has_trait_map) {
    if (verbose) message("Calibrating gates on validation set...")
    gpu_mem_checkpoint("entering gate calibration")
    model$eval()

    # Forward pass on val set (use t_X_eval so held-out truths do not leak).
    torch::with_no_grad({
      mask_ind0 <- torch::torch_zeros(c(n, 1L), device = device)
      covs0     <- make_covs_tensor(t_MU, mask_ind0, t_covariates)
      out_cal   <- model(t_X_eval, t_coords, covs0, t_adj, t_obs_to_sp,
                         D_sq = t_D_sq)
      delta_cal <- as.matrix(out_cal$delta$cpu())
      mu_cal    <- as.matrix(t_MU$cpu())
    })

    X_truth_r    <- X_truth
    val_mask_mat <- matrix(FALSE, n, p)
    val_mask_mat[splits$val_idx] <- TRUE

    # Fix C.3 (Opus 2026-04-28): split val into a CALIBRATION half (used by
    # `calibrate_gates()` to pick per-trait blend weights) and a CONFORMAL
    # half (used by `compute_conformal_scores()` to estimate residual
    # quantiles). Reusing the same val cells for both breaks split-conformal
    # exchangeability — the gate is selected to minimise val residual MSE,
    # which post-selects the very residuals that drive the conformal
    # quantile, producing systematic undercoverage. The split is per
    # latent column to avoid wasting cells on traits with already-tiny
    # val sets.
    #
    # Min-val safeguard: only split a column when its val set has at
    # least `2 * min_val_cells` cells.  Below that threshold the half-A /
    # half-B cross-check inside `calibrate_gates()` (which itself halves
    # the input) becomes too noisy and the calibrated gate degrades in
    # ways that show up as a +2-25% RMSE regression on benches like the
    # AVONET300 safety-floor smoke test.  Columns under the threshold
    # silently fall back to the legacy single-set behaviour for that
    # column (calibration uses the full val cells; conformal also uses
    # them but accepts the modest undercoverage risk).  Set
    # `conformal_split_val = FALSE` to disable splitting everywhere.
    split_threshold <- 2L * as.integer(min_val_cells)
    if (isTRUE(conformal_split_val)) {
      set.seed(seed + 23L)
      val_mask_cal  <- matrix(FALSE, n, p)
      val_mask_conf <- matrix(FALSE, n, p)
      for (jcol in seq_len(p)) {
        idx_j <- which(val_mask_mat[, jcol])
        n_j   <- length(idx_j)
        if (n_j == 0L) next
        if (n_j < split_threshold) {
          # Below threshold: do not split — both halves get the full val
          # cells.  Documented downside: conformal scores for this column
          # are post-selected on the gate calibration cells, undercovering
          # by an amount bounded by the size of the gate-grid search
          # (small in practice).
          val_mask_cal[idx_j, jcol]  <- TRUE
          val_mask_conf[idx_j, jcol] <- TRUE
          next
        }
        cal_n   <- ceiling(n_j / 2)
        cal_idx <- sample(idx_j, cal_n)
        val_mask_cal[cal_idx, jcol]                        <- TRUE
        val_mask_conf[setdiff(idx_j, cal_idx), jcol]       <- TRUE
      }
    } else {
      val_mask_cal  <- val_mask_mat
      val_mask_conf <- val_mask_mat
    }

    # Safety-floor: compute per-latent-column grand mean on training-observed
    # cells only. Excludes val + test hold-out to prevent leakage. Works in
    # both single-obs and multi-obs mode (X_scaled is obs-level either way).
    mean_baseline_per_col <- if (safety_floor) {
      mb        <- numeric(p)
      latent_nm <- colnames(data$X_scaled)
      if (!is.null(latent_nm)) names(mb) <- latent_nm

      # Build per-column trait-type lookup by walking trait_map.
      # multi_proportion CLR columns are treated as "continuous" on the latent
      # scale (each component is a z-scored CLR value).
      col_type <- character(p)
      for (tm_entry in data$trait_map) {
        for (idx in seq_along(tm_entry$latent_cols)) {
          lc <- tm_entry$latent_cols[idx]
          col_type[lc] <- if (tm_entry$type == "zi_count") {
            if (idx == 1L) "binary" else "zi_mag"
          } else if (tm_entry$type == "multi_proportion") {
            "continuous"
          } else {
            tm_entry$type
          }
        }
      }

      # Build test mask to exclude test cells (val already in val_mask_mat).
      test_mask_mat_sf <- matrix(FALSE, n, p)
      if (!is.null(splits) && length(splits$test_idx) > 0L) {
        test_mask_mat_sf[splits$test_idx] <- TRUE
      }
      # Training-observed = not val, not test, not NA.
      train_mask_mat <- !val_mask_mat & !test_mask_mat_sf & !is.na(data$X_scaled)

      for (j in seq_len(p)) {
        mb[j] <- mean_baseline_scalar(
          x_col      = data$X_scaled[, j],
          train_mask = train_mask_mat[, j],
          trait_type = col_type[j]
        )
      }
      # Guard: replace any NA that slipped through with 0 (e.g. columns with
      # no training observations at all — degenerate but possible).
      mb[is.na(mb)] <- 0
      mb
    } else {
      NULL
    }

    calibrated_gates_list <- calibrate_gates(
      trait_map             = trait_map,
      mu_cal                = mu_cal,
      delta_cal             = delta_cal,
      X_truth_r             = X_truth_r,
      val_mask_mat          = val_mask_cal,
      gate_grid             = seq(0, gate_cap, length.out = 9L),
      gate_cap              = gate_cap,
      gate_method           = gate_method,
      gate_splits_B         = gate_splits_B,
      gate_cv_folds         = gate_cv_folds,
      safety_floor          = safety_floor,
      mean_baseline_per_col = mean_baseline_per_col,
      simplex_step          = 0.05,
      min_val_cells         = min_val_cells,
      seed                  = seed,
      latent_names          = data$latent_names,
      verbose               = verbose
    )
    calibrated_gates <- calibrated_gates_list$r_cal_gnn   # legacy scalar slot

    # Apply phylo-signal gate override: gated latent columns get (0, 0, 1)
    # regardless of what the 3-way simplex calibrator picked.
    if (length(gated_latent_cols) > 0L && !is.null(calibrated_gates_list)) {
      calibrated_gates_list$r_cal_bm[gated_latent_cols]   <- 0
      calibrated_gates_list$r_cal_gnn[gated_latent_cols]  <- 0
      calibrated_gates_list$r_cal_mean[gated_latent_cols] <- 1
      calibrated_gates <- calibrated_gates_list$r_cal_gnn
    }
  }

  # ---- Conformal prediction scores (validation set) -----------------------
  conformal_scores <- NULL
  if (has_val && has_trait_map && !is.null(calibrated_gates)) {
    gpu_mem_checkpoint("entering compute_conformal_scores")
    conformal_scores <- compute_conformal_scores(
      trait_map        = trait_map,
      calibrated_gates = calibrated_gates,
      mu_cal           = mu_cal,
      delta_cal        = delta_cal,
      X_truth_r        = X_truth_r,
      val_mask_mat     = val_mask_conf,
      method           = conformal_method,
      bootstrap_B      = conformal_bootstrap_B,
      verbose          = verbose
    )
  }

  # ---- Final test loss -------------------------------------------------------
  test_loss <- NA_real_
  if (!is.null(splits) && length(splits$test_idx) > 0L) {
    test_mat <- matrix(FALSE, n, p); test_mat[splits$test_idx] <- TRUE
    t_test   <- torch::torch_tensor(test_mat, dtype = torch::torch_bool(),
                                    device = device)
    model$eval()
    torch::with_no_grad({
      mask_ind0 <- torch::torch_zeros(c(n, 1L), device = device)
      covs0     <- make_covs_tensor(t_MU, mask_ind0, t_covariates)
      # Use t_X_eval so held-out test cells are not leaked to the model.
      out0      <- model(t_X_eval, t_coords, covs0, t_adj, t_obs_to_sp,
                         D_sq = t_D_sq)
      pred0     <- (1 - out0$rs) * t_MU + out0$rs * out0$delta
      t_truth_f <- torch::torch_tensor(X_truth, dtype = torch::torch_float(),
                                       device = device)
      t_truth_f[t_truth_f$isnan()] <- 0
    })
    if (has_trait_map) {
      test_loss <- composite_val_loss(pred0, t_truth_f, t_test, trait_map)
    } else {
      test_loss <- as.numeric(
        rmse_torch(pred0, t_truth_f, t_test)$cpu()$item()
      )
    }
  }

  # ---- Build result object ---------------------------------------------------
  model_config <- list(
    hidden_dim             = hidden_dim,
    k_eigen                = k_eigen,
    n_gnn_layers           = n_gnn_layers,
    gate_cap               = gate_cap,
    use_attention          = use_attention,
    use_transformer_blocks = use_transformer_blocks,
    n_heads                = as.integer(n_heads),
    ffn_mult               = as.integer(ffn_mult),
    dropout                = dropout,
    refine_steps           = refine_steps,
    cov_dim                = cov_dim,
    input_dim              = p,
    per_column_rs          = TRUE,
    n_user_cov             = n_user_cov,
    seed                   = as.integer(seed)
  )

  # Move model state to CPU before returning. Otherwise the returned
  # fit object holds GPU tensor references and the training-time
  # activations / optimizer state cannot be released by torch's CUDA
  # caching allocator. At n >= 5000 this keeps ~40 GB live on GPU, and
  # the downstream predict() path then OOMs on any card < 80 GB. See
  # GPU bundle jobs 4741993/4741994/4742888/4742889 on 2026-04-21 for
  # the Vulcan L40S reproducer.
  gpu_mem_checkpoint("before state_dict CPU move")
  model_state_cpu <- lapply(model$state_dict(),
                             function(t) t$detach()$cpu())
  gpu_mem_checkpoint("after state_dict CPU move (before rm + empty_cache)")

  # Drop all local GPU tensor refs, then force torch's caching
  # allocator to reclaim by calling cuda_empty_cache().  R's `rm`
  # removes the bindings; the underlying tensors become unreachable
  # and torch can release them at the next empty_cache() call.
  local_gpu_tensors <- grep("^t_", ls(), value = TRUE)
  if (length(local_gpu_tensors)) rm(list = local_gpu_tensors)
  rm(model)
  invisible(gc(full = TRUE, verbose = FALSE))
  if (torch::cuda_is_available()) {
    try(torch::cuda_empty_cache(), silent = TRUE)
  }
  gpu_mem_checkpoint("end of fit_pigauto (after rm + gc + empty_cache)")

  # Backward-compat: store val_rmse and test_rmse names
  # (graph$D was stripped earlier, right after tensor creation, so the
  # graph stored here is already the slim version without the cophenetic
  # distance matrix.)
  structure(
    list(
      model_state    = model_state_cpu,
      model_config   = model_config,
      graph          = graph,
      baseline       = baseline,
      norm           = list(
        means         = data$means,
        sds           = data$sds,
        log_transform = data$log_transform
      ),
      species_names  = data$species_names,
      obs_species    = data$obs_species,
      obs_to_species = data$obs_to_species,
      # Phase G' (2026-05-01): retain X_scaled so predict.pigauto_fit
      # can recover original-units observed values for PMM
      # (match_observed = "pmm").  The matrix is n_obs x p_latent with
      # NA at originally-missing cells; PMM uses non-NA cells as the
      # donor pool.  Adds ~8 * n_obs * p_latent bytes to the fit object;
      # at AVONET full (n=10k, p~7) that's ~600 KB -- negligible.
      X_scaled       = data$X_scaled,
      n_species      = n_species,
      n_obs          = n_obs,
      multi_obs      = multi_obs,
      trait_names    = data$trait_names,
      latent_names   = data$latent_names,
      trait_map      = trait_map,
      splits         = splits,
      history        = history,
      val_rmse         = best_val,
      test_rmse        = test_loss,
      calibrated_gates = calibrated_gates,
      r_cal      = if (!is.null(calibrated_gates_list))
                     calibrated_gates_list$r_cal_gnn else NULL,
      r_cal_bm   = if (!is.null(calibrated_gates_list))
                     calibrated_gates_list$r_cal_bm  else NULL,
      r_cal_gnn  = if (!is.null(calibrated_gates_list))
                     calibrated_gates_list$r_cal_gnn else NULL,
      r_cal_mean = if (!is.null(calibrated_gates_list))
                     calibrated_gates_list$r_cal_mean else NULL,
      mean_baseline_per_col  = mean_baseline_per_col,
      safety_floor           = safety_floor,
      phylo_signal_per_trait = phylo_signal_per_trait,
      phylo_gate_triggered   = phylo_gate_triggered,
      phylo_signal_method    = phylo_signal_method,
      phylo_signal_threshold = phylo_signal_threshold,
      conformal_scores = conformal_scores,
      covariates       = data$covariates,
      cov_means        = data$cov_means,
      cov_sds          = data$cov_sds,
      cov_names        = data$cov_names
    ),
    class = "pigauto_fit"
  )
}
