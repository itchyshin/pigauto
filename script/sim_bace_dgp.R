#!/usr/bin/env Rscript
# script/sim_bace_dgp.R
#
# Thin DGP wrapper around BACE::sim_bace() for the comprehensive
# pigauto-vs-baselines simulation study.
#
# This is sourced (not run directly) by:
#   - script/bench_sim_bace_pigauto.R  (sweep driver)
#
# Goals:
#  * Standardise the output format so every bench / analysis script
#    consumes a consistent list shape.
#  * Apply MCAR mask to the response trait only (predictors must be
#    fully observed since pigauto's `covariates` argument doesn't
#    accept NAs).
#  * Generate cleanly with full reproducibility (seed control).

suppressPackageStartupMessages({
  if (!requireNamespace("BACE", quietly = TRUE)) {
    stop("BACE package required. Install via: ",
         "devtools::install('BACE') from the pigauto root.")
  }
  library(BACE)
  library(ape)
})

#' Generate a simulated comparative-data scenario via BACE::sim_bace().
#'
#' @param n_species integer. Number of species (tree tips).
#' @param multi_obs_ratio integer. n_cases = multi_obs_ratio * n_species.
#'   1 = single-obs; >1 = multi-obs (replicates per species).
#' @param phylo_signal numeric in [0, 1). Phylogenetic signal of the
#'   response trait. Predictors get fixed phylo_signal = 0.3 (they're
#'   ecological covariates with moderate phylo structure).
#' @param beta_resp_strength numeric. Strength of the response-on-predictor
#'   effect. 0 = no covariate signal; 0.5 = moderate; 1.0 = strong.
#' @param response_type character. One of "gaussian", "binary", "thresholdK"
#'   (e.g. "threshold4").
#' @param n_predictors integer. Number of predictor traits (default 3).
#' @param miss_frac numeric in (0, 1). MCAR fraction of response held out.
#' @param seed integer. Random seed (controls tree, traits, mask).
#' @param birth, death numeric. Tree sim rates.
#' @return A list with:
#'   - tree:           phylo, with tips matching the data's species column
#'   - df_complete:    data.frame, full data (no missingness)
#'   - df_observed:    data.frame, response NA-masked at held-out cells
#'   - mask:           logical vector (length n_cases), TRUE = held-out
#'   - response_name:  character (default "y")
#'   - predictor_names: character vector ("x1", ..., "xK")
#'   - response_type:  character (echoes input for convenience)
#'   - n_species:      integer
#'   - n_cases:        integer
#'   - meta:           list of parameter values used
sim_bace_dgp <- function(n_species,
                          multi_obs_ratio = 1L,
                          phylo_signal = 0.4,
                          beta_resp_strength = 0.5,
                          response_type = "gaussian",
                          n_predictors = 3L,
                          miss_frac = 0.30,
                          seed = 2026L,
                          birth = 0.8,
                          death = 0.4) {

  stopifnot(
    is.numeric(n_species), n_species >= 30L,
    multi_obs_ratio >= 1L,
    phylo_signal >= 0, phylo_signal < 1,
    beta_resp_strength >= 0,
    miss_frac > 0, miss_frac < 1,
    n_predictors >= 1L
  )

  n_cases <- as.integer(n_species * multi_obs_ratio)

  # Predictor types: mix of gaussian + binary so both pigauto baselines
  # are exercised (cov-aware GLS on continuous, label-prop on binary).
  predictor_types <- if (n_predictors == 3L) {
    c("gaussian", "gaussian", "binary")
  } else if (n_predictors == 5L) {
    c("gaussian", "gaussian", "gaussian", "binary", "binary")
  } else {
    rep("gaussian", n_predictors)
  }

  # Phylo signal vector: response + predictors. Predictors get moderate
  # signal (0.3) to mimic real ecological covariates (climate, etc.).
  phylo_vec <- c(phylo_signal, rep(0.3, n_predictors))

  # beta_resp: scale by beta_resp_strength.  When strength = 0, response
  # has zero coefficient on every predictor (covariates carry no signal).
  beta_resp <- if (beta_resp_strength == 0) {
    rep(0, n_predictors)
  } else {
    # Spread the strength across predictors with mild heterogeneity
    base_betas <- c(0.7, 0.4, 0.5, 0.3, 0.6)[seq_len(n_predictors)]
    beta_resp_strength * base_betas
  }

  # Missingness: only response gets MCAR mask.  Predictors fully observed
  # because pigauto's `covariates` argument requires no NAs.
  miss_vec <- c(miss_frac, rep(0, n_predictors))

  # Run BACE simulator
  set.seed(seed)
  sim <- BACE::sim_bace(
    response_type   = response_type,
    predictor_types = predictor_types,
    beta_resp       = beta_resp,
    phylo_signal    = phylo_vec,
    n_cases         = n_cases,
    n_species       = as.integer(n_species),
    missingness     = miss_vec,
    birth           = birth,
    death           = death
  )

  df_complete <- sim$complete_data
  df_observed <- sim$data
  tree        <- sim$tree

  # Identify response and predictor column names
  response_name <- "y"
  predictor_names <- grep("^x[0-9]+$", names(df_complete), value = TRUE)
  if (length(predictor_names) != n_predictors) {
    # BACE may name them differently; fall back to all-non-y, non-species
    predictor_names <- setdiff(names(df_complete),
                                 c("species", "animal", response_name))
  }

  # Mask: cells in df_observed where response is NA
  mask <- is.na(df_observed[[response_name]])

  # Sanity: tree tips should match all species labels in df
  sp_labels <- as.character(unique(df_complete$species))
  if (!all(sp_labels %in% tree$tip.label)) {
    warning("Some species in df not in tree$tip.label")
  }

  list(
    tree           = tree,
    df_complete    = df_complete,
    df_observed    = df_observed,
    mask           = mask,
    response_name  = response_name,
    predictor_names = predictor_names,
    response_type  = response_type,
    n_species      = as.integer(n_species),
    n_cases        = n_cases,
    meta = list(
      n_species          = as.integer(n_species),
      multi_obs_ratio    = multi_obs_ratio,
      phylo_signal       = phylo_signal,
      beta_resp_strength = beta_resp_strength,
      response_type      = response_type,
      n_predictors       = n_predictors,
      miss_frac          = miss_frac,
      seed               = seed,
      predictor_types    = predictor_types
    )
  )
}
