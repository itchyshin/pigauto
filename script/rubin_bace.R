# script/rubin_bace.R
#
# BACE arms for the Rubin's-rules MI-arm comparison lane (arc/rubin-freq-bace,
# docs/dev-log/arc/2026-09-24-rubin-freq-bace-plan.md). S3 owns this file and its test only -- no
# edits to R/ or BACE/ (BACE is Dan Noble's code, installed as a package, read-only). Assumes
# script/campaign_gnn_off_lib.R has already been sourced (for make_cell(), run_bace()'s data-prep
# conventions it reuses, and bace_diagnostics()); rubin_cell.R sources it first.
#
# Two arms on top of one BACE::bace() fit:
#   mi_bace_shipped() -- BACE's own n_final imputed datasets, unmodified ("BACE as shipped")
#   mi_bace_resid()   -- BACE's shipped datasets plus an INDEPENDENT per-cell residual draw on
#                        gaussian-modelled missing cells ("BACE + residual draw (post hoc)")
#
# IMPORTANT premise correction (2026-09-24): the brief for this arm was written against
# BACE/R/model_functions.R as it reads on disk (gitignored, local-only, and STALE relative to what
# is actually installed -- confirmed by deparsing the installed namespace: the on-disk
# bace_final_imp() computes `levels <- dat_prep$levels`, a call into a two-element UNNAMED list
# that is always NULL and would throw on every categorical trait; the INSTALLED bace_final_imp()
# has clearly been patched past that and past several other things). On disk, .predict_bace()'s
# gaussian branch takes the posterior MEAN (`pred_prob <- .pred_cont(model) * sd_val + mean_val`),
# so the stale source supports "BACE's continuous imputations carry no residual draw." The
# INSTALLED bace_final_imp() -- the one BACE::bace() actually runs -- calls
# `.predict_bace(model, dat_prep, response_var, type, sample = TRUE, formula = ..., data_full =
# ..., cluster_col = ...)` for every trait, hardcoded, not exposed as a bace()-level argument. For
# type == "gaussian" the installed .predict_bace()'s sample = TRUE branch draws ONE posterior
# iteration's fitted value and ADDS `rnorm(n_obs, 0, sqrt(sigma2_units))` before back-transforming
# by `* sd_val + mean_val` -- i.e. BACE's shipped continuous imputations, as this package is
# actually installed and run here, ALREADY are posterior-predictive draws with a residual
# component, not posterior means. Confirmed empirically too (script/tests-rubin dev log): spread
# across n_final shipped datasets at one missing cell was the same order of magnitude as
# sqrt(mean VCV[,"units"]) on the z scale.
#
# Consequence: mi_bace_resid(), built exactly to spec below, adds a SECOND independent residual
# draw on top of the one BACE's own installed final-imputation step already adds -- it does not
# "add back a missing predictive draw," it roughly doubles the injected residual variance for
# gaussian traits. Both arms are still well-defined and testable as written; whether "BACE +
# residual draw" remains the right label, or whether the plan doc's Q1 needs revisiting, is a
# design call outside this file (flagged to main).
#
# Review update (Meng, 2026-09-24, B2 and N8): with the installed build this arm is a NEGATIVE CONTROL
# (a doubled residual), not a fix. A further property matters more: every final run starts from the
# same converged dataset (bace_final_imp: .one_final_run(run, last_data)), whose fills came from
# sample = FALSE runs, so predictors later in the formula order are identical across the M datasets.
# A BACE arm that chains final runs is the candidate replacement; that design call is Shinichi's.

#' Fit BACE with n_final = M, reusing run_bace()'s data-prep conventions (script/campaign_gnn_off_lib.R
#' L601-653: Species column from rownames, one formula per trait regressed on all the others,
#' n_cores = 1L, skip_conv = TRUE, ovr_categorical = TRUE).
#'
#' @param cell a make_cell() result (uses cell$tree, cell$df_miss)
#' @param M number of BACE final imputation datasets (BACE's n_final)
#' @param nitt,burnin,thin,runs passed straight to BACE::bace()
#' @return list(outb = the "bace_complete" object, wall_s = fit wall time in seconds,
#'   converged = BACE's own convergence verdict, ess_med = median effective sample size over the
#'   final fits' fixed effects, diag = the full bace_diagnostics() list -- same extraction G9b/v1 use)
fit_bace_mi <- function(cell, M = 20L, nitt, burnin, thin, runs) {
  tree_b <- cell$tree
  if (any(tree_b$edge.length == 0)) tree_b$edge.length[tree_b$edge.length == 0] <- 1e-8
  df_b <- cell$df_miss; df_b$Species <- rownames(cell$df_miss)
  all_traits <- setdiff(names(df_b), "Species")
  fixformula <- lapply(all_traits, function(v)
    paste0(v, " ~ ", paste(setdiff(all_traits, v), collapse = " + ")))

  t0 <- Sys.time()
  outb <- BACE::bace(fixformula = fixformula, ran_phylo_form = "~1|Species", phylo = tree_b,
                     data = df_b, nitt = nitt, burnin = burnin, thin = thin, runs = runs,
                     n_final = M, n_cores = 1L, verbose = FALSE, skip_conv = TRUE,
                     ovr_categorical = TRUE)
  wall_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  d <- bace_diagnostics(outb)   # same convergence/ESS extraction v1/G9b use (campaign_gnn_off_lib.R)
  list(outb = outb, wall_s = wall_s, converged = d$converged, ess_med = d$ess_med, diag = d,
       fixformula = fixformula, tree_b = tree_b, mcmc = c(nitt = nitt, burnin = burnin, thin = thin))
}

#' Chained BACE (the arm Meng's review B2 proposed; Shinichi approved it 2026-09-24). BACE's own final
#' step starts all n_final runs from the same converged dataset, so predictors later in the formula order
#' are identical across the M imputations. Here final run m starts from final dataset m - 1 instead
#' (run 1 from BACE's converged dataset), using BACE's own bace_final_imp() with n_final = 1 and the same
#' MCMC settings; BACE's code is called, never edited. It reuses the shipped fit's convergence phase
#' (fb$outb$initial_results), so the shipped and chained arms are paired on one fit. Cost: M extra
#' sweeps, run sequentially.
#'
#' @param fb a fit_bace_mi() result
#' @param df_miss the cell's missing-data data.frame
#' @param M number of chained datasets
#' @return list(datasets = M data.frames aligned like mi_bace_shipped(), models = per-run model lists,
#'   wall_s = chain wall time)
mi_bace_chain <- function(fb, df_miss, M = 20L) {
  obj <- fb$outb$initial_results
  k <- length(obj$data)
  raw <- vector("list", M); models <- vector("list", M)
  t0 <- Sys.time()
  for (m in seq_len(M)) {
    r <- BACE:::bace_final_imp(obj, fixformula = fb$fixformula, ran_phylo_form = "~1|Species",
                               phylo = fb$tree_b, nitt = fb$mcmc[["nitt"]], thin = fb$mcmc[["thin"]],
                               burnin = fb$mcmc[["burnin"]], n_final = 1L, species = FALSE,
                               verbose = FALSE, n_cores = 1L, ovr_categorical = TRUE)
    raw[[m]] <- r$all_datasets[[1]]; models[[m]] <- r$all_models[[1]]
    obj$data[[k]] <- raw[[m]]   # the next run starts from this draw
  }
  wall_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(datasets = mi_bace_shipped(list(imputed_datasets = raw), df_miss), models = models, wall_s = wall_s)
}

#' BACE's own n_final imputed datasets, unchanged, aligned to df_miss's row order/column
#' names/classes (BACE's internal Species column dropped; storage class only, never values, is
#' coerced back to df_miss's -- e.g. BACE returns a rounded double for count traits, df_miss has
#' them as integer).
#'
#' @param outb a "bace_complete" object (fit_bace_mi()$outb)
#' @param df_miss the cell's missing-data data.frame (make_cell()$df_miss)
#' @return list of length M, each a data.frame shaped like df_miss
mi_bace_shipped <- function(outb, df_miss) {
  lapply(outb$imputed_datasets, function(d) {
    d <- d[rownames(df_miss), names(df_miss), drop = FALSE]
    for (v in names(df_miss)) {
      if (is.integer(df_miss[[v]]) && !is.integer(d[[v]])) {
        d[[v]] <- as.integer(round(d[[v]]))
      } else if (is.factor(df_miss[[v]]) && !identical(levels(d[[v]]), levels(df_miss[[v]]))) {
        d[[v]] <- factor(as.character(d[[v]]), levels = levels(df_miss[[v]]),
                         ordered = is.ordered(df_miss[[v]]))
      }
    }
    d
  })
}

#' sd_val BACE used to z-score gaussian trait v: .extract_gaussian_attrs() runs on the response with its
#' missing rows reset to NA, so it is the sd of the OBSERVED values, the same in every final run (Meng
#' review N8, docs/dev-log/arc/2026-09-24-rubin-review.md; the earlier run-chaining reconstruction was
#' wrong because every final run starts from the same converged dataset).
#'
#' @param df_miss the cell's missing-data data.frame
#' @param v a gaussian-modelled trait name
#' @return numeric scalar, sd_val on the data scale
.bace_gaussian_sd <- function(df_miss, v) stats::sd(df_miss[[v]], na.rm = TRUE)

#' BACE's shipped datasets (mi_bace_shipped()) plus an independent residual draw added to every
#' MISSING cell of every trait BACE modelled as "gaussian" (identified from
#' outb$final_results$types, not hard-coded -- in this DGP that is c1, c2, prp: all numeric and
#' non-integer, so BACE's is.numeric & !is.integer rule classes them gaussian; cnt is integer-typed
#' and goes through BACE's AIC(normal) vs AIC(poisson) test, which favours poisson for these count
#' rates). Observed cells are never touched. Label: "BACE + residual draw (post hoc)" -- and see
#' the file header for why, against the INSTALLED BACE package, this adds a SECOND residual draw
#' rather than adding one back. The draw does not propagate through BACE's chained imputation (a
#' trait imputed downstream of v in the same final run still only sees v's un-perturbed BACE
#' prediction, because the perturbation is applied here, after BACE has already returned), and it
#' ignores residual correlation between the per-trait model fits (each trait's noise is drawn
#' independently, even though BACE fits per-trait models on a shared, correlated DGP).
#'
#' @param outb a "bace_complete" object (fit_bace_mi()$outb)
#' @param df_miss the cell's missing-data data.frame
#' @param seed RNG seed (one call is reproducible; different seeds give independent draws)
#' @return list of length M, each a data.frame shaped like df_miss
mi_bace_resid <- function(outb, df_miss, seed) {
  out <- mi_bace_shipped(outb, df_miss)
  types <- outb$final_results$types
  gaussian_traits <- names(types)[vapply(types, identical, logical(1), y = "gaussian")]
  gaussian_traits <- intersect(gaussian_traits, names(df_miss))
  M <- length(out)

  set.seed(seed + 7919L)   # distinct from the DGP seed (Meng N8)
  for (i in seq_len(M)) {
    for (v in gaussian_traits) {
      miss_idx <- which(is.na(df_miss[[v]]))
      if (!length(miss_idx)) next
      fit_v <- outb$final_results$all_models[[i]][[v]]
      if (is.null(fit_v)) next   # defensive: gaussian traits with missing data are always fit
      sd_val <- .bace_gaussian_sd(df_miss, v)
      units <- as.matrix(fit_v$VCV)[, "units"]
      sigma2_units_i <- units[sample.int(length(units), 1L)]   # NOT sample(units, 1): sample() on
      # a length-1 numeric would sample from 1:units instead of returning units itself
      noise <- stats::rnorm(length(miss_idx), mean = 0, sd = sqrt(sigma2_units_i) * sd_val)
      out[[i]][miss_idx, v] <- out[[i]][miss_idx, v] + noise
    }
  }
  out
}
