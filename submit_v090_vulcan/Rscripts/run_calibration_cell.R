#!/usr/bin/env Rscript
#
# Rscripts/run_calibration_cell.R
#
# Pigauto calibration-grid sim, one array cell per SLURM task.
#
# Purpose
#   Probe pigauto's 95% coverage and per-trait NRMSE/accuracy across a
#   (phylogenetic-signal x missingness-mechanism x trait-type) grid, using
#   Dan Noble's sim_bace() harness from the BACE package. Replicates the
#   Phase-1 coverage finding on pigauto at scale.
#
# Grid (60 cells)
#   4 signal scenarios: all_high, all_moderate, all_low, mixed
#   3 mechanisms:       phylo_MAR, trait_MAR, trait_MNAR
#   5 trait types:      gaussian, binary, multinomial, poisson, ordinal
#   -> 4 x 3 x 5 = 60 cells (matches submit_calibration_grid.sh --array=1-60)
#
# Per-task runtime
#   N_REPS x per-fit wall (~3 min at n=150 single-obs, 500 epochs)
#   Default N_REPS = 50  ->  ~2.5 hr per task
#
# Usage
#   Rscript run_calibration_cell.R <SLURM_ARRAY_TASK_ID>   # 1..60
#
# Output (written to results/ relative to CWD at submit time)
#   results/cell_<TASK_ID>.rds   tidy data.frame; one row per (rep, trait)
#                                 metric; carries scenario, mechanism, and
#                                 focal trait_type tags.
#
# Dependencies
#   pigauto  (load-path is $R_LIBS_USER from submit_calibration_grid.sh)
#   BACE     (install on login node: pak::pak("local::~/BACE_src"))
#   ape, MASS
#
# Resume
#   The script checkpoints the partial RDS after every rep. Re-running with
#   the same TASK_ID resumes from the last rep_id already recorded.

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  library(MASS)
  library(BACE)
  library(pigauto)
})

# ---------------------------------------------------------------------------
# Arg parsing
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript run_calibration_cell.R <SLURM_ARRAY_TASK_ID (1..60)>",
       call. = FALSE)
}
TASK_ID <- suppressWarnings(as.integer(args[1]))
if (!is.finite(TASK_ID) || TASK_ID < 1L) {
  stop("TASK_ID must be a positive integer. Got: ", args[1], call. = FALSE)
}

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

CONFIG <- list(
  n_reps     = 50L,    # reps per cell; matches ~2.5 hr SLURM wall budget
  n_species  = 150L,   # single-obs tree size (matches Dan's BACE harness)
  miss_rate  = 0.35,   # target missing fraction per variable
  epochs     = 500L,
  m_impute   = 20L,    # MI draws per dataset for 95% coverage estimation
  alpha      = 0.05,   # coverage target = 1 - alpha
  dep_strength = 1.5,  # MAR/MNAR coefficient (Dan's default)
  base_seed  = 2026L
)

SIGNAL   <- list(high = 0.90, moderate = 0.60, low = 0.20)

CELL_GRID <- expand.grid(
  scenario   = c("all_high", "all_moderate", "all_low", "mixed"),
  mechanism  = c("phylo_MAR", "trait_MAR", "trait_MNAR"),
  trait_type = c("gaussian", "binary", "multinomial", "poisson", "ordinal"),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
CELL_GRID$cell_id <- seq_len(nrow(CELL_GRID))
stopifnot(nrow(CELL_GRID) == 60L)

if (TASK_ID > nrow(CELL_GRID)) {
  stop("TASK_ID ", TASK_ID, " exceeds grid size ", nrow(CELL_GRID),
       call. = FALSE)
}

cell       <- CELL_GRID[TASK_ID, ]
scenario   <- cell$scenario
mechanism  <- cell$mechanism
trait_type <- cell$trait_type

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_rds <- file.path(out_dir, sprintf("cell_%02d.rds", TASK_ID))

# ---------------------------------------------------------------------------
# Helpers (vendored inline so the script is self-contained; these used to
# live in dev/coverage_sim/sim_bace_source.R but are too small to warrant
# a separate file in this bundle.)
# ---------------------------------------------------------------------------

# Translate scenario label -> length-5 phylo_signal vector for sim_bace().
# Order matches (response, x1, x2, x3, x4) = (focal, four predictors).
make_phylo_signal <- function(scenario) {
  switch(scenario,
    all_high     = rep(SIGNAL$high,     5L),
    all_moderate = rep(SIGNAL$moderate, 5L),
    all_low      = rep(SIGNAL$low,      5L),
    mixed        = c(SIGNAL$high, SIGNAL$moderate, SIGNAL$low,
                     SIGNAL$moderate, SIGNAL$high),
    stop("Unknown scenario: ", scenario, call. = FALSE)
  )
}

# Coerce a column (factor / ordered / numeric) to a numeric vector suitable
# for use as a MAR / MNAR linear predictor.
trait_to_numeric <- function(x) {
  if (is.factor(x)) return(as.numeric(x) - 1)
  if (is.logical(x)) return(as.numeric(x))
  as.numeric(x)
}

# Root-find an intercept so that plogis(c_hat + lp) has mean = target rate.
# Fall back to the empirical quantile if the bracketing search fails.
calibrate_intercept <- function(lp, target_rate) {
  obj <- function(c_hat) mean(plogis(c_hat + lp)) - target_rate
  res <- tryCatch(
    stats::uniroot(obj, interval = c(-20, 20), extendInt = "yes")$root,
    error = function(e) NA_real_
  )
  if (is.na(res)) res <- stats::quantile(-lp, probs = target_rate)
  res
}

# Single-observation missingness injector. Operates per-variable on a
# complete data.frame with species rownames. Returns both the masked
# data.frame and the corresponding TRUE = missing mask.
inject_singleobs <- function(cd, tree, mech, rate,
                             vars = c("y", "x1", "x2", "x3", "x4")) {
  n  <- nrow(cd)
  mm <- as.data.frame(matrix(FALSE, n, length(vars),
                             dimnames = list(NULL, vars)))
  md <- cd

  # Phylogenetic covariance for phylo_MAR (avoid ape::vcv per-variable cost)
  Sp <- if (mech == "phylo_MAR") ape::vcv(tree, corr = TRUE) else NULL

  # Target-map for trait_MAR: variable v's missingness depends on tmc[[v]]
  tmc <- c(y = "x3", x1 = "y", x2 = "y", x3 = "x1", x4 = "y")

  for (v in vars) {
    lp <- switch(mech,
      phylo_MAR = {
        z <- MASS::mvrnorm(1L, mu = rep(0, nrow(Sp)), Sigma = Sp)
        names(z) <- rownames(Sp)
        -CONFIG$dep_strength * z[rownames(cd)]
      },
      trait_MAR  = -CONFIG$dep_strength *
        as.numeric(scale(trait_to_numeric(cd[[tmc[[v]]]]))),
      trait_MNAR = -CONFIG$dep_strength *
        as.numeric(scale(trait_to_numeric(cd[[v]]))),
      stop("Unknown mechanism: ", mech, call. = FALSE)
    )
    c_hat <- calibrate_intercept(lp, rate)
    mm[[v]] <- stats::rbinom(n, 1L, plogis(c_hat + lp)) == 1L
    md[[v]][mm[[v]]] <- NA
  }

  list(miss_data = md, miss_mask = mm)
}

# 95% coverage and point-metric extraction from a pigauto_mi object.
#
# For each held-out cell on the focal variable, compute:
#   - coverage95: truth ∈ [q_{alpha/2}, q_{1-alpha/2}] of MI draws.
#   - point_rmse (continuous), point_acc (discrete): pooled point estimate.
#
# For discrete traits the MI credible region is the support of draws (modal
# prediction's confidence). coverage95 there = "truth is in the top-95% mass
# of the MI draw distribution", which collapses to:
#   - binary:    truth ∈ {modal class} whenever P(truth) >= 0.025.
#   - categorical / ordinal: ranking-based, see compute_coverage_discrete().
compute_cell_metrics <- function(truth_df, mi, miss_mask, trait_type,
                                 focal_var, alpha = 0.05) {
  masked_idx <- which(miss_mask[[focal_var]])
  if (length(masked_idx) == 0L) {
    return(data.frame(metric = character(0), value = numeric(0),
                      n_cells = integer(0)))
  }

  draws <- vapply(mi$datasets, function(d) d[masked_idx, focal_var],
                  FUN.VALUE = truth_df[masked_idx, focal_var])
  if (!is.matrix(draws)) draws <- matrix(draws, ncol = length(mi$datasets))

  truth <- truth_df[masked_idx, focal_var]

  if (trait_type %in% c("gaussian", "poisson")) {
    # Continuous / count: quantile-based conformal CI from MI draws
    lo <- apply(draws, 1L, stats::quantile, probs = alpha / 2,     na.rm = TRUE)
    hi <- apply(draws, 1L, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)
    truth_num <- trait_to_numeric(truth)
    cov_hits  <- truth_num >= lo & truth_num <= hi
    point_est <- rowMeans(matrix(as.numeric(draws),
                                 nrow = nrow(draws)), na.rm = TRUE)
    rmse <- sqrt(mean((truth_num - point_est)^2, na.rm = TRUE))
    data.frame(
      metric  = c("coverage95", "rmse", "nrmse"),
      value   = c(mean(cov_hits, na.rm = TRUE), rmse,
                  rmse / (stats::sd(truth_num, na.rm = TRUE) + 1e-12)),
      n_cells = length(masked_idx)
    )
  } else {
    # Discrete (binary / multinomial / ordinal): bootstrap-style coverage
    # on the MI draw distribution. For each masked cell the MI gives a
    # distribution over classes; coverage95 = P(truth in smallest set whose
    # mass >= 1 - alpha), i.e. the 95% credible set from the draws.
    truth_chr <- as.character(truth)
    draw_chr  <- matrix(as.character(draws), nrow = nrow(draws))
    cov_hits  <- vapply(seq_len(nrow(draw_chr)), function(i) {
      tab <- prop.table(table(draw_chr[i, ], useNA = "no"))
      tab <- sort(tab, decreasing = TRUE)
      cumm <- cumsum(tab)
      k95  <- which(cumm >= 1 - alpha)[1]
      if (is.na(k95)) k95 <- length(tab)
      truth_chr[i] %in% names(tab)[seq_len(k95)]
    }, logical(1))
    # Point accuracy: modal MI prediction == truth
    modes <- vapply(seq_len(nrow(draw_chr)), function(i) {
      tab <- table(draw_chr[i, ])
      names(tab)[which.max(tab)]
    }, character(1))
    acc <- mean(modes == truth_chr, na.rm = TRUE)
    data.frame(
      metric  = c("coverage95", "accuracy"),
      value   = c(mean(cov_hits, na.rm = TRUE), acc),
      n_cells = length(masked_idx)
    )
  }
}

# Recode the focal trait type in the simulated data.frame so pigauto
# infers the right class automatically.
recast_focal <- function(cd, focal_var, trait_type) {
  if (trait_type == "binary") {
    # sim_bace returns 0/1 numeric -> factor with 2 levels
    cd[[focal_var]] <- factor(as.character(cd[[focal_var]]), levels = c("0", "1"))
  } else if (trait_type == "ordinal") {
    # Treat the focal as ordered factor for pigauto's ordinal path
    lv <- sort(unique(stats::na.omit(cd[[focal_var]])))
    cd[[focal_var]] <- factor(as.character(cd[[focal_var]]),
                              levels = as.character(lv), ordered = TRUE)
  } else if (trait_type == "multinomial") {
    cd[[focal_var]] <- factor(as.character(cd[[focal_var]]))
  } else if (trait_type == "poisson") {
    cd[[focal_var]] <- as.integer(cd[[focal_var]])
  }
  cd
}

# ---------------------------------------------------------------------------
# Main: one-cell loop over reps
# ---------------------------------------------------------------------------

log_line <- function(...) {
  cat(sprintf("[cell=%02d %s/%s/%s] ",
              TASK_ID, scenario, mechanism, trait_type),
      ..., "\n", sep = "")
  flush.console()
}

log_line(sprintf("starting: n_reps=%d, n_species=%d, miss_rate=%.2f, epochs=%d",
                 CONFIG$n_reps, CONFIG$n_species, CONFIG$miss_rate,
                 CONFIG$epochs))

# Map trait_type -> sim_bace response/predictor encoding. We run sim_bace
# with response_type = focal; the other 4 variables are a fixed mix of
# binary + multinomial + poisson + ordinal (threshold) so that the
# covariance structure and MAR target-map tmc stay interpretable across
# cells. This keeps comparability across trait_type cells while still
# isolating the focal type's coverage.
sim_response_type <- switch(trait_type,
  gaussian    = "gaussian",
  binary      = "binary",
  multinomial = "multinomial3",
  poisson     = "poisson",
  ordinal     = "threshold3",
  stop("Unknown trait_type: ", trait_type, call. = FALSE)
)
# Fixed predictor ensemble (Dan's setup). Focal is the response; the four
# predictors stay constant so (scenario, mechanism) are what changes.
sim_predictor_types <- c("binary", "multinomial3", "poisson", "threshold3")
phylo_vec           <- make_phylo_signal(scenario)
focal_var           <- "y"  # sim_bace names the response "y" by default

# Resume: if an RDS exists for this TASK_ID, pick up from the last rep_id.
existing <- if (file.exists(out_rds)) {
  tryCatch(readRDS(out_rds), error = function(e) NULL)
} else NULL
done_reps <- if (!is.null(existing) && nrow(existing)) {
  sort(unique(existing$rep))
} else integer(0)
if (length(done_reps)) {
  log_line(sprintf("resuming; %d reps already recorded", length(done_reps)))
}

rows <- if (is.null(existing)) list() else list(existing)

for (rep_id in setdiff(seq_len(CONFIG$n_reps), done_reps)) {
  t0      <- Sys.time()
  rep_seed <- as.integer(CONFIG$base_seed + 101L * rep_id +
                         17L * TASK_ID + as.integer(rep_id))
  set.seed(rep_seed)

  # ---- Simulate -----------------------------------------------------------
  sim <- tryCatch(
    sim_bace(
      response_type   = sim_response_type,
      predictor_types = sim_predictor_types,
      var_names       = c("y", "x1", "x2", "x3", "x4"),
      phylo_signal    = phylo_vec,
      n_cases         = CONFIG$n_species,
      n_species       = CONFIG$n_species,
      beta_sparsity   = 0.3,
      missingness     = rep(0, 5L)
    ),
    error = function(e) e
  )
  if (inherits(sim, "error")) {
    log_line(sprintf("rep %d: sim_bace FAILED: %s", rep_id,
                     conditionMessage(sim)))
    next
  }

  cd   <- sim$data
  tree <- sim$tree
  rownames(cd) <- cd$species
  cd$species  <- NULL
  cd          <- recast_focal(cd, focal_var, trait_type)

  # Keep a pristine copy of the COMPLETE data for coverage evaluation
  truth_df <- cd

  # ---- Inject missingness -------------------------------------------------
  miss <- tryCatch(
    inject_singleobs(cd, tree, mechanism, CONFIG$miss_rate),
    error = function(e) e
  )
  if (inherits(miss, "error")) {
    log_line(sprintf("rep %d: inject_singleobs FAILED: %s", rep_id,
                     conditionMessage(miss)))
    next
  }

  # ---- multi_impute -------------------------------------------------------
  mi <- tryCatch(
    multi_impute(
      traits        = miss$miss_data,
      tree          = tree,
      m             = CONFIG$m_impute,
      species_col   = NULL,
      draws_method  = "conformal",
      log_transform = FALSE,
      missing_frac  = 0.25,
      epochs        = CONFIG$epochs,
      verbose       = FALSE,
      seed          = rep_seed
    ),
    error = function(e) e
  )
  if (inherits(mi, "error")) {
    log_line(sprintf("rep %d: multi_impute FAILED: %s", rep_id,
                     conditionMessage(mi)))
    next
  }

  wall_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # ---- Metrics ------------------------------------------------------------
  metrics <- tryCatch(
    compute_cell_metrics(
      truth_df    = truth_df,
      mi          = mi,
      miss_mask   = miss$miss_mask,
      trait_type  = trait_type,
      focal_var   = focal_var,
      alpha       = CONFIG$alpha
    ),
    error = function(e) e
  )
  if (inherits(metrics, "error")) {
    log_line(sprintf("rep %d: metrics FAILED: %s", rep_id,
                     conditionMessage(metrics)))
    next
  }

  metrics$scenario   <- scenario
  metrics$mechanism  <- mechanism
  metrics$trait_type <- trait_type
  metrics$cell_id    <- TASK_ID
  metrics$rep        <- rep_id
  metrics$seed       <- rep_seed
  metrics$wall_s     <- wall_s

  rows[[length(rows) + 1L]] <- metrics

  # Checkpoint after every rep so SLURM kills leave usable partials
  combined <- do.call(rbind, rows)
  saveRDS(combined, out_rds)

  cov_row <- metrics[metrics$metric == "coverage95", ]
  log_line(sprintf("rep %d done: cov95=%.2f  n_cells=%d  wall=%.0fs",
                   rep_id,
                   if (nrow(cov_row)) cov_row$value[1] else NA_real_,
                   if (nrow(cov_row)) cov_row$n_cells[1] else NA_integer_,
                   wall_s))

  rm(sim, cd, tree, truth_df, miss, mi, metrics); gc(verbose = FALSE)
}

# ---------------------------------------------------------------------------
# Finalise
# ---------------------------------------------------------------------------

if (length(rows) == 0L) {
  log_line("no successful reps; writing empty RDS")
  saveRDS(data.frame(), out_rds)
} else {
  combined <- do.call(rbind, rows)
  saveRDS(combined, out_rds)
  log_line(sprintf("done: %d reps recorded in %s",
                   length(unique(combined$rep)), out_rds))
}
