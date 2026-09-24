# script/mi_gls/regimes.R
#
# Regime grid for the MI-GLS-attenuation campaign (arc/mi-gls-attenuation)
# AND its posterior-MI extension (arc/mi-posterior). Regimes 1-16 are
# UNCHANGED from the original campaign (same ids, same seeds -- see
# script/mi_gls/dgp_v2.R -- so v1 (01_cell.R) and v2 (01_cell_v2.R) results
# join by regime_id + rep). Regimes 17-24 extend the grid to validate
# multi_impute(draws_method = "posterior") (docs/dev-log/mi-posterior/
# design.md) when BOTH traits are missing and have DIFFERENT phylogenetic
# signal, which the shared-lambda tree-transform DGP of regimes 1-16 cannot
# express.
#
# rho (trait correlation between x and y) is fixed at 0.7 throughout
# regimes 1-16 -- it is not part of that grid.
#
# `lambda` is a Pagel lambda applied to the tree FOR SIMULATION ONLY (via
# the internal transform_tree_pagel(), R/pagel_lambda.R): traits are
# generated on the lambda-transformed tree, but every downstream step
# (imputation and both lm/gls fits) uses the ORIGINAL (lambda = 1) tree.
# At lambda = 0.5 this deliberately misspecifies the phylogenetic model
# used for imputation/inference relative to the true DGP.
#
# `missing`:
#   - "x_only": 30% of x cells missing (drawn per the `mechanism` below);
#     y always fully observed.
#   - "both": 30% of x cells AND (independently) 30% of y cells missing,
#     each drawn separately per `mechanism`.
#
# `mechanism`:
#   - "MCAR": each candidate cell missing independently with p = 0.30.
#   - "MAR_phylo": missingness concentrated in 2 randomly chosen clades
#     (15-35% of tips each), as in script/mondrian_confirmation/06_mi_se_sim.R
#     (itself following ~/pigauto_regime_map/mech_cell.R) -- clade tips get
#     7x the baseline miss probability, rescaled so the expected count is
#     30% of n. "both" draws a FRESH pair of clades independently for y.

regimes <- expand.grid(
  lambda    = c(1, 0.5),
  n         = c(300L, 1000L),
  mechanism = c("MCAR", "MAR_phylo"),
  missing   = c("x_only", "both"),
  KEEP.OUT.ATTRS  = FALSE,
  stringsAsFactors = FALSE
)
regimes$rho <- 0.7
regimes$regime_id <- seq_len(nrow(regimes))
regimes <- regimes[, c("regime_id", "lambda", "n", "mechanism", "missing", "rho")]

# Oracle proper-MI is only run when lambda == 1 (correctly-specified model)
# AND missing == "x_only" (the oracle formula in 01_cell.R conditions on
# ALL of y being observed; see 13_mi_gls_attenuation_diag.R's derivation).
regimes$run_oracle <- regimes$lambda == 1 & regimes$missing == "x_only"

## ---- regimes 17-24: two-trait Sigma_P / Sigma_E Kronecker DGP -----------
##
## docs/dev-log/mi-posterior/design.md section 1:
##   vec(Y) ~ N(mu, Sigma_P %x% R + Sigma_E %x% I_n),  R = cov2cor(vcv(tree))
##   lambda_k = Sigma_P[k,k] / (Sigma_P[k,k] + Sigma_E[k,k])
## Sigma_P, Sigma_E chosen with unit total variance per trait
## (Sigma_P[k,k] + Sigma_E[k,k] = 1), so:
##   Sigma_P[k,k] = lambda_k,  Sigma_E[k,k] = 1 - lambda_k
##   Sigma_P[1,2] = corr_phylo * sqrt(lambda_x * lambda_y)
##   Sigma_E[1,2] = corr_resid * sqrt((1-lambda_x) * (1-lambda_y))
##
## 17-20: per-trait lambda differ (x = 0.3, y = 0.9); phylogenetic and
##   residual correlation both 0.7. MCAR (17, 18), MAR_phylo (19, 20);
##   n alternates 300/1000.
## 21-24: lambda = 0.7 for both traits; phylogenetic correlation 0.7,
##   residual correlation 0. MCAR (21, 22), MAR_phylo (23, 24); n
##   alternates 300/1000.
## Both traits missing (missing = "both") in all eight; oracle proper-MI
## (v1's true-model formula) does not generalise to this two-lambda DGP, so
## run_oracle = FALSE throughout.

# lambda_x, lambda_y, phylogenetic correlation, residual correlation ->
# Sigma_P, Sigma_E (2 x 2), unit total variance per trait. Shared with
# script/mi_gls/dgp_v2.R (which sources this file) and the parameter
# self-check below.
build_sigma_pair <- function(lambda_x, lambda_y, corr_phylo, corr_resid) {
  cp_p <- corr_phylo * sqrt(lambda_x * lambda_y)
  cp_e <- corr_resid * sqrt((1 - lambda_x) * (1 - lambda_y))
  Sigma_P <- matrix(c(lambda_x, cp_p, cp_p, lambda_y), 2, 2)
  Sigma_E <- matrix(c(1 - lambda_x, cp_e, cp_e, 1 - lambda_y), 2, 2)
  list(Sigma_P = Sigma_P, Sigma_E = Sigma_E)
}

# Marginal (single-tip) regression coefficient of y on x: R_ii = 1 for
# every tip, so this does not depend on the tree. It is the OLS estimand
# ONLY. In regimes 17-24 Cov(x) = Sp11 R + Se11 I and Cov(y, x) =
# Sp12 R + Se12 I are not proportional, so E[y | x] is not a scalar times
# x and neither analysis model targets this value: corBrownian GLS and
# phylolm(lambda) weight tips by the tree and converge elsewhere (review
# finding dgp#0; complete-data coverage of this value is near 0 there).
# It is therefore DESCRIPTIVE only. The downstream-coverage truth for
# regimes 17-24 is a pseudo-truth, the mean complete-data slope over the
# expected reps of that (regime, analysis model), computed in
# 03_summarise_v2.R (design.md section 5c, decision D1). Regimes 1-16
# keep rho = 0.7, which is exact for any analysis model there because the
# DGP covariance is proportional (E[y | x] = 0.7 x).
true_beta_pop_kron <- function(Sigma_P, Sigma_E) unname(Sigma_P[1, 2] + Sigma_E[1, 2])

extra <- expand.grid(
  lambda_pair = c("differ", "same"),
  mechanism   = c("MCAR", "MAR_phylo"),
  n           = c(300L, 1000L),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
# Order to match the task spec exactly: 17,18 = MCAR/differ (n=300,1000);
# 19,20 = MAR_phylo/differ; 21,22 = MCAR/same; 23,24 = MAR_phylo/same.
extra <- extra[order(match(extra$lambda_pair, c("differ", "same")),
                     match(extra$mechanism, c("MCAR", "MAR_phylo")),
                     match(extra$n, c(300L, 1000L))), ]
extra$lambda_x    <- ifelse(extra$lambda_pair == "differ", 0.3, 0.7)
extra$lambda_y    <- ifelse(extra$lambda_pair == "differ", 0.9, 0.7)
extra$corr_phylo  <- 0.7
extra$corr_resid  <- ifelse(extra$lambda_pair == "differ", 0.7, 0.0)
extra$missing     <- "both"
extra$lambda      <- NA_real_   # not used by the Kronecker DGP (see lambda_x/lambda_y)
extra$rho         <- NA_real_   # not used by the Kronecker DGP (see corr_phylo/corr_resid)
extra$run_oracle  <- FALSE
extra$regime_id   <- 16L + seq_len(nrow(extra))
extra$true_beta_pop <- vapply(seq_len(nrow(extra)), function(i) {
  sg <- build_sigma_pair(extra$lambda_x[i], extra$lambda_y[i],
                         extra$corr_phylo[i], extra$corr_resid[i])
  true_beta_pop_kron(sg$Sigma_P, sg$Sigma_E)
}, numeric(1))

regimes$lambda_x <- NA_real_
regimes$lambda_y <- NA_real_
regimes$corr_phylo <- NA_real_
regimes$corr_resid <- NA_real_
regimes$true_beta_pop <- regimes$rho   # well-defined and identical to rho for 1-16 (see 01_cell.R header)

extra <- extra[, names(regimes)]
regimes <- rbind(regimes, extra)
rownames(regimes) <- NULL

## ---- expected campaign grid (arc/mi-posterior, 2026-09-24) ---------------
##
## One source for which rows the v2 harness must produce, read by
## 01_cell_v2.R (which methods to run), 03_summarise_v2.R (which summary
## rows to write) and the G6/G7 gates (which rows must exist). The grid is
## always built from `regimes`, never from the files on disk (design.md
## section 5c, decisions D3 and D6).

# Replicates per regime in the approved plan.
mi_gls_v2_planned_reps <- 200L

# Posterior MI methods run in a regime. posterior_none (parameters fixed)
# exists only for G6 rule 4 (proper vs improper SE ratio in the
# both-missing regimes), so it runs only where both traits are missing.
mi_gls_v2_methods <- function(missing) {
  if (identical(missing, "both")) c("posterior_full", "posterior_none") else "posterior_full"
}

# Traits with masked cells in a regime (per-cell coverage keys).
mi_gls_v2_traits <- function(missing) {
  if (identical(missing, "both")) c("x", "y") else "x"
}

# "1-3,7" -> c(1L, 2L, 3L, 7L)
mi_gls_v2_parse_ids <- function(spec) {
  parts <- trimws(strsplit(spec, ",", fixed = TRUE)[[1L]])
  parts <- parts[nzchar(parts)]
  ids <- unlist(lapply(parts, function(p) {
    if (grepl("^[0-9]+-[0-9]+$", p)) {
      ab <- as.integer(strsplit(p, "-", fixed = TRUE)[[1L]])
      seq(ab[1L], ab[2L])
    } else if (grepl("^[0-9]+$", p)) {
      as.integer(p)
    } else {
      stop("cannot parse regime spec '", spec, "' (use e.g. 1-24 or 1,21)", call. = FALSE)
    }
  }))
  sort(unique(ids))
}

# Expected grid for the summariser and the gates. Env MI_REGIMES (default:
# every regime in `regimes`) and MI_N_REPS (default
# mi_gls_v2_planned_reps). `full` is TRUE only for the planned grid (all
# regimes, at least the planned rep count); the gates print their pass
# token only when `full` is TRUE.
mi_gls_v2_expected <- function() {
  spec <- Sys.getenv("MI_REGIMES", "")
  ids <- if (nzchar(spec)) mi_gls_v2_parse_ids(spec) else regimes$regime_id
  bad <- setdiff(ids, regimes$regime_id)
  if (length(bad)) {
    stop("MI_REGIMES names unknown regime id(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }
  n_reps <- suppressWarnings(as.integer(Sys.getenv("MI_N_REPS", as.character(mi_gls_v2_planned_reps))))
  if (!is.finite(n_reps) || n_reps < 1L) {
    stop("MI_N_REPS must be a positive integer; got '", Sys.getenv("MI_N_REPS"), "'", call. = FALSE)
  }
  list(regime_ids = ids, n_reps = n_reps,
       full = setequal(ids, regimes$regime_id) && n_reps >= mi_gls_v2_planned_reps,
       label = sprintf("regimes=%s (%d) reps=%d", if (nzchar(spec)) spec else "all",
                       length(ids), n_reps))
}

if (sys.nframe() == 0L) {
  print(regimes)

  cat("\n---- regime 17-24 parameter check (large-n draw, R = identity) ----\n")
  cat("Single-tip marginal Var/Cov does not depend on the phylogenetic\n")
  cat("correlation matrix R (R_ii = 1 for every tip regardless of tree\n")
  cat("shape), so drawing i.i.d. (a, e) pairs (R = I) from Sigma_P and\n")
  cat("Sigma_E is a valid, tree-free check of the Kronecker construction.\n\n")

  check_kron_params <- function(lambda_x, lambda_y, corr_phylo, corr_resid,
                                n_draw = 200000L) {
    sg <- build_sigma_pair(lambda_x, lambda_y, corr_phylo, corr_resid)
    La <- chol(sg$Sigma_P); Le <- chol(sg$Sigma_E)
    A <- matrix(stats::rnorm(n_draw * 2), n_draw, 2) %*% La
    E <- matrix(stats::rnorm(n_draw * 2), n_draw, 2) %*% Le
    Y <- A + E
    data.frame(
      lambda_x_target = lambda_x, lambda_x_emp = stats::var(A[, 1]) / stats::var(Y[, 1]),
      lambda_y_target = lambda_y, lambda_y_emp = stats::var(A[, 2]) / stats::var(Y[, 2]),
      corr_phylo_target = corr_phylo, corr_phylo_emp = stats::cor(A[, 1], A[, 2]),
      corr_resid_target = corr_resid, corr_resid_emp = stats::cor(E[, 1], E[, 2]),
      var_x_emp = stats::var(Y[, 1]), var_y_emp = stats::var(Y[, 2])
    )
  }

  set.seed(1L)
  for (rid in 17:24) {
    r <- regimes[regimes$regime_id == rid, ]
    chk <- check_kron_params(r$lambda_x, r$lambda_y, r$corr_phylo, r$corr_resid)
    cat(sprintf("regime %2d: lambda_x %.2f->%.4f  lambda_y %.2f->%.4f  corr_P %.2f->%.4f  corr_E %.2f->%.4f  var_x=%.4f var_y=%.4f  true_beta_pop=%.4f\n",
               rid, chk$lambda_x_target, chk$lambda_x_emp,
               chk$lambda_y_target, chk$lambda_y_emp,
               chk$corr_phylo_target, chk$corr_phylo_emp,
               chk$corr_resid_target, chk$corr_resid_emp,
               chk$var_x_emp, chk$var_y_emp, r$true_beta_pop))
  }
}
