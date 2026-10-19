# Frozen scenario design for the 2026-10-18 validation. No model fitting here.
# campaign_sim_design.R and mi_gls/regimes.R remain the source of their grids.
validation_manifest <- function(root = ".") {
  root <- normalizePath(root, mustWork = TRUE)
  design <- file.path(root, "script/campaign_sim_design.R")
  if (!file.exists(design)) stop("missing campaign_sim_design.R")
  campaign <- function(stage) {
    out <- suppressWarnings(system2("Rscript", c(shQuote(design), "--stage", stage),
                                    stdout = TRUE, stderr = FALSE))
    if (!is.null(attr(out, "status")) && attr(out, "status") != 0L)
      stop("campaign design failed: ", stage)
    read.csv(text = paste(out, collapse = "\n"), stringsAsFactors = FALSE)
  }
  rows <- list()
  add <- function(slice, adapter, scenario, n = NA_integer_, regime_id = NA_integer_,
                  model = NA_character_, mechanism = NA_character_, lambda = NA_real_,
                  rho = NA_real_, evo = NA_character_, frac = NA_real_, reps,
                  claim_tier = "descriptive") {
    d <- data.frame(slice, adapter, n, scenario, regime_id, model, mechanism,
                    lambda, rho, evo, frac, reps, claim_tier,
                    stringsAsFactors = FALSE)
    d$scenario_id <- sprintf("%s_%03d", slice, seq_len(nrow(d)))
    d <- d[, c("scenario_id", "slice", "adapter", "n", "scenario",
               "regime_id", "model", "mechanism", "lambda", "rho", "evo",
               "frac", "reps", "claim_tier")]
    rows[[length(rows) + 1L]] <<- d
  }
  for (stage in c("core", "factorial")) {
    x <- campaign(stage)
    add(paste0("prediction_", stage), "prediction", x$cell, n = x$n,
        mechanism = x$miss, lambda = x$lambda, rho = x$rho, evo = x$evo,
        frac = x$frac, reps = 200L, claim_tier = "recovery")
  }
  low <- expand.grid(n = c(100L, 300L), lambda = c(0, 0.1),
                     rho = c(0, 0.5), KEEP.OUT.ATTRS = FALSE)
  add("prediction_low_signal", "prediction",
      sprintf("types_mixed_BM_l%s_r%s_mcar0.3_n%d", low$lambda, low$rho, low$n),
      n = low$n, mechanism = "mcar", lambda = low$lambda, rho = low$rho,
      evo = "BM", frac = 0.3, reps = 200L, claim_tier = "recovery")
  add("zi_count", "zi_count",
      c("zf_0.2", "zf_0.4", "zf_0.6", "zf_0.8",
        "mean_nz_5", "mean_nz_20", "mean_nz_100"),
      n = 300L, mechanism = "mcar", frac = 0.25, reps = 200L,
      claim_tier = "recovery")
  add("multi_proportion", "multi_proportion",
      c("signal_0.2", "signal_0.4", "signal_0.6", "signal_0.8",
        "signal_1.0", "K_3", "K_8", "K_12"),
      n = 300L, mechanism = "mcar", frac = 0.25, reps = 200L,
      claim_tier = "recovery")
  mo <- expand.grid(lambda = c(0.5, 0.9), beta = c(0, 0.5, 1),
                    frac = c(0.5, 0.8), KEEP.OUT.ATTRS = FALSE)
  add("multiobs", "multiobs", sprintf("ctmax_l%s_b%s_f%s", mo$lambda, mo$beta, mo$frac),
      n = 300L, mechanism = "species_missing", lambda = mo$lambda,
      frac = mo$frac, reps = 100L, claim_tier = "recovery")
  add("multiobs_mixed", "multiobs_mixed",
      c("high_phylo_correlated", "mod_phylo_correlated", "low_phylo_correlated"),
      n = 150L, mechanism = "species_missing", lambda = c(1, 0.5, 0.2),
      rho = c(0.7, 0.6, 0.6), frac = 0.3, reps = 100L,
      claim_tier = "recovery")
  cov <- expand.grid(scenario = c("high_phylo_no_env", "high_phylo_mod_env",
                                  "mod_phylo_strong_env", "low_phylo_strong_env"),
                     frac = c(0.2, 0.4, 0.6), KEEP.OUT.ATTRS = FALSE,
                     stringsAsFactors = FALSE)
  lambda_cov <- c(0.9, 0.7, 0.3, 0.1)[match(cov$scenario,
                    c("high_phylo_no_env", "high_phylo_mod_env",
                      "mod_phylo_strong_env", "low_phylo_strong_env"))]
  add("covariate", "covariate", paste0(cov$scenario, "_f", cov$frac),
      n = 300L, mechanism = "mcar", lambda = lambda_cov, frac = cov$frac,
      reps = 100L, claim_tier = "recovery")
  add("tree", "tree", paste0("avonet_f", c(0.2, 0.4, 0.6)), n = 300L,
      mechanism = "mcar", frac = c(0.2, 0.4, 0.6), reps = 20L,
      claim_tier = "descriptive")
  e <- new.env(parent = baseenv())
  sys.source(file.path(root, "script/mi_gls/regimes.R"), envir = e)
  ids <- c(25L, 26L, 27L, 28L, 33L, 34L, 35L, 38L)
  pr <- e$regimes[match(ids, e$regimes$regime_id), ]
  if (anyNA(pr$regime_id) || !all(pr$gated)) stop("posterior regime source mismatch")
  add("posterior", "posterior", paste0("regime_", pr$regime_id),
      n = pr$n, regime_id = pr$regime_id, model = "gls",
      mechanism = pr$mechanism, lambda = pr$lambda, rho = pr$rho,
      frac = 0.3, reps = 500L, claim_tier = "inferential")
  an <- expand.grid(model = c("lm", "glm", "lmer"), n = c(100L, 300L),
                    mechanism = c("MCAR", "MAR"), KEEP.OUT.ATTRS = FALSE,
                    stringsAsFactors = FALSE)
  add("analysis", "analysis", sprintf("%s_n%d_%s", an$model, an$n, an$mechanism),
      n = an$n, model = an$model, mechanism = an$mechanism, frac = 0.3,
      reps = 500L, claim_tier = "inferential")
  add("real", "real", c("pantheria", "avonet", "bien", "fishbase"),
      mechanism = "observed_mask", frac = 0.3, reps = 20L,
      claim_tier = "descriptive")
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  validate_manifest(out)
  out
}

validate_manifest <- function(manifest) {
  needed <- c("scenario_id", "slice", "adapter", "n", "scenario", "regime_id",
              "model", "mechanism", "lambda", "rho", "evo", "frac", "reps",
              "claim_tier")
  if (!is.data.frame(manifest) || !identical(names(manifest), needed) ||
      nrow(manifest) < 1L) stop("manifest schema invalid")
  if (anyNA(manifest$scenario_id) || any(!nzchar(manifest$scenario_id)) ||
      anyDuplicated(manifest$scenario_id) ||
      anyNA(manifest$scenario) || any(!nzchar(manifest$scenario)) ||
      anyNA(manifest$adapter) || any(!nzchar(manifest$adapter)))
    stop("manifest identity invalid")
  if (!is.numeric(manifest$reps) || anyNA(manifest$reps) ||
      any(manifest$reps < 1L) || any(manifest$reps != as.integer(manifest$reps)))
    stop("manifest replication invalid")
  if (any(!manifest$claim_tier %in% c("inferential", "recovery", "descriptive")))
    stop("manifest claim tier invalid")
  if (any(!is.na(manifest$frac) & (manifest$frac <= 0 | manifest$frac >= 1)))
    stop("manifest missingness fraction invalid")
  TRUE
}

validation_jobs <- function(manifest) {
  validate_manifest(manifest)
  idx <- rep(seq_len(nrow(manifest)), manifest$reps)
  out <- manifest[idx, , drop = FALSE]
  out$rep <- sequence(manifest$reps)
  out$job_id <- sprintf("%s__%04d", out$scenario_id, out$rep)
  out$seed <- 20261018L + as.integer(idx) * 10000L + out$rep
  out <- out[, c("job_id", names(manifest), "rep", "seed")]
  rownames(out) <- NULL
  if (anyDuplicated(out$job_id) || anyDuplicated(out$seed))
    stop("job identity collision")
  out
}
