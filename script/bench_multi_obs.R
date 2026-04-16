#!/usr/bin/env Rscript
#
# script/bench_multi_obs.R
#
# Benchmark: multi-observation-per-species imputation with covariates
#
# Purpose
#   Exercise pigauto's multi-obs code path (species_col) and covariate
#   conditioning on a simulated CTmax-like dataset.  Tests whether the
#   GNN + covariate pipeline recovers species-level and observation-level
#   trait values better than a naive species-mean baseline, and whether
#   covariates improve accuracy when a within-species predictor exists.
#
# Data-generating process
#   CTmax_ij = mu + phylo_i + beta * acclim_temp_ij + epsilon_ij
#
#   phylo_i        ~ BM(tree, sigma = sqrt(lambda))
#   acclim_temp_ij ~ N(20, 5)       (observation-level covariate)
#   epsilon_ij     ~ N(0, 0.5)      (residual)
#   mu = 30                          (grand mean)
#
# Sweep
#   lambda (phylo signal):       0.5, 0.9
#   beta (acclimation slope):    0.0, 0.5, 1.0
#   species_missing_frac:        0.5, 0.8
#   n_species:                   200
#   obs_per_species:             Poisson(5), clamped [1, 20]
#   n_reps:                      2
#   epochs:                      200
#
# Methods
#   1. species_mean   -- observed species mean, grand mean for unobserved
#   2. pigauto_no_cov -- impute() with species_col, no covariates
#   3. pigauto_cov    -- impute() with species_col + acclimation temperature
#
# Metrics (on held-out observations)
#   - obs_rmse:     observation-level RMSE
#   - sp_rmse:      species-level RMSE (mean pred vs mean truth, unobserved spp)
#   - obs_pearson_r: observation-level Pearson r
#
# Output
#   script/bench_multi_obs.rds    tidy results + metadata
#   script/bench_multi_obs.md     human-readable summary
#
# Run with
#   cd pigauto && /usr/local/bin/Rscript script/bench_multi_obs.R

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/phase-10-multiobs",
    quiet = TRUE
  )
})

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "bench_multi_obs.rds")
out_md  <- file.path(here, "script", "bench_multi_obs.md")

script_start <- proc.time()[["elapsed"]]

log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start), ..., "\n",
      sep = "")
  flush.console()
}

# -------------------------------------------------------------------------
# Constants
# -------------------------------------------------------------------------

n_species         <- 200L
mu_ctmax          <- 30.0
sigma_residual    <- 0.5
acclim_mean       <- 20.0
acclim_sd         <- 5.0
obs_lambda_pois   <- 5L        # Poisson mean for obs-per-species
obs_min           <- 1L
obs_max           <- 20L
n_reps            <- 2L
epochs            <- 200L
within_obs_miss   <- 0.20      # fraction of remaining (observed) spp obs masked

# Sweep parameters
lambdas            <- c(0.5, 0.9)
betas              <- c(0.0, 0.5, 1.0)
sp_missing_fracs   <- c(0.5, 0.8)

# -------------------------------------------------------------------------
# Helpers: simulation
# -------------------------------------------------------------------------

simulate_ctmax_data <- function(tree, lambda, beta, seed) {
  set.seed(seed)
  sp <- tree$tip.label
  n_sp <- length(sp)

  # Species-level phylogenetic effect (BM with sigma = sqrt(lambda))
  phylo_vals <- ape::rTraitCont(tree, model = "BM",
                                sigma = sqrt(lambda),
                                root.value = 0)
  phylo_vals <- as.numeric(phylo_vals[sp])
  names(phylo_vals) <- sp

  # Variable number of observations per species
  n_obs_per <- rpois(n_sp, obs_lambda_pois)
  n_obs_per <- pmin(pmax(n_obs_per, obs_min), obs_max)
  names(n_obs_per) <- sp

  # Build observation-level data
  species_vec   <- rep(sp, n_obs_per)
  n_total       <- length(species_vec)
  acclim_temp   <- rnorm(n_total, mean = acclim_mean, sd = acclim_sd)
  epsilon       <- rnorm(n_total, mean = 0, sd = sigma_residual)
  phylo_per_obs <- phylo_vals[species_vec]

  ctmax <- mu_ctmax + phylo_per_obs + beta * acclim_temp + epsilon

  # True species means (averaging over within-species variation)
  # This is the expected value at the species-level acclimation mean
  sp_true_mean <- mu_ctmax + phylo_vals + beta * acclim_mean

  df <- data.frame(
    species    = species_vec,
    CTmax      = ctmax,
    acclim_temp = acclim_temp,
    stringsAsFactors = FALSE
  )

  list(
    df           = df,
    sp           = sp,
    n_obs_per    = n_obs_per,
    phylo_vals   = phylo_vals,
    sp_true_mean = sp_true_mean
  )
}


# -------------------------------------------------------------------------
# Helpers: masking
# -------------------------------------------------------------------------

mask_species <- function(df, sim, sp_missing_frac, within_miss_frac, seed) {
  # Mask entire species (all their observations become NA)
  # Plus mask a fraction of observations within the observed species
  set.seed(seed + 999L)
  sp <- sim$sp
  n_sp <- length(sp)

  n_sp_miss <- floor(sp_missing_frac * n_sp)
  missing_sp <- sample(sp, n_sp_miss)
  observed_sp <- setdiff(sp, missing_sp)

  df_masked <- df
  is_missing_sp <- df$species %in% missing_sp
  df_masked$CTmax[is_missing_sp] <- NA

  # Within observed species, mask some observations
  obs_idx <- which(!is_missing_sp)
  n_within_miss <- floor(within_miss_frac * length(obs_idx))
  if (n_within_miss > 0L) {
    within_miss_idx <- sample(obs_idx, n_within_miss)
    df_masked$CTmax[within_miss_idx] <- NA
  }

  list(
    df_masked    = df_masked,
    missing_sp   = missing_sp,
    observed_sp  = observed_sp,
    all_na_rows  = which(is.na(df_masked$CTmax))
  )
}


# -------------------------------------------------------------------------
# Helpers: species-mean baseline
# -------------------------------------------------------------------------

species_mean_impute <- function(df_masked) {
  # For each species: use observed mean if available, else grand mean
  sp_means <- tapply(df_masked$CTmax, df_masked$species,
                     function(x) mean(x, na.rm = TRUE))
  # Grand mean from all observed values
  grand_mean <- mean(df_masked$CTmax, na.rm = TRUE)

  preds <- numeric(nrow(df_masked))
  for (i in seq_len(nrow(df_masked))) {
    sp_i <- df_masked$species[i]
    sm <- sp_means[[sp_i]]
    preds[i] <- if (is.finite(sm)) sm else grand_mean
  }
  preds
}


# -------------------------------------------------------------------------
# Helpers: metrics
# -------------------------------------------------------------------------

rmse <- function(truth, pred) {
  ok <- is.finite(truth) & is.finite(pred)
  if (sum(ok) == 0L) return(NA_real_)
  sqrt(mean((truth[ok] - pred[ok])^2))
}

pearson_r <- function(truth, pred) {
  ok <- is.finite(truth) & is.finite(pred)
  if (sum(ok) < 3L) return(NA_real_)
  stats::cor(truth[ok], pred[ok])
}

compute_metrics <- function(truth_obs, pred_obs, truth_sp_mean, pred_sp_mean,
                            missing_sp) {
  # Observation-level metrics on ALL held-out observations
  obs_rmse_val <- rmse(truth_obs, pred_obs)
  obs_r_val    <- pearson_r(truth_obs, pred_obs)

  # Species-level RMSE on unobserved species only
  sp_idx <- names(truth_sp_mean) %in% missing_sp
  sp_rmse_val <- rmse(truth_sp_mean[sp_idx], pred_sp_mean[sp_idx])

  data.frame(
    obs_rmse    = obs_rmse_val,
    sp_rmse     = sp_rmse_val,
    obs_pearson_r = obs_r_val,
    stringsAsFactors = FALSE
  )
}


# -------------------------------------------------------------------------
# Helpers: extract pigauto predictions
# -------------------------------------------------------------------------

extract_pigauto_preds <- function(result, df_full, mask_info) {
  # result is a pigauto_result from impute()
  # df_full is the complete (unmasked) data
  # Returns observation-level predictions and species-level means

  completed <- result$completed

  # Observation-level predictions for held-out cells
  held_out_idx <- mask_info$all_na_rows
  pred_obs <- completed$CTmax[held_out_idx]
  truth_obs <- df_full$CTmax[held_out_idx]

  # Species-level mean predictions for missing species
  sp_preds <- tapply(completed$CTmax, completed$species, mean, na.rm = TRUE)

  list(
    pred_obs  = pred_obs,
    truth_obs = truth_obs,
    sp_preds  = sp_preds
  )
}


# -------------------------------------------------------------------------
# Main runner for a single cell
# -------------------------------------------------------------------------

run_one_cell <- function(lambda, beta, sp_missing_frac, rep_id) {
  rep_seed <- rep_id * 1000L + round(lambda * 100) + round(beta * 10) +
    round(sp_missing_frac * 100)

  cell_label <- sprintf("lambda=%.1f, beta=%.1f, sp_miss=%.1f, rep=%d",
                        lambda, beta, sp_missing_frac, rep_id)
  log_line(sprintf("=== %s (seed=%d) ===", cell_label, rep_seed))

  # 1. Simulate tree and data
  set.seed(rep_seed)
  tree <- ape::rtree(n_species)
  sim  <- simulate_ctmax_data(tree, lambda, beta, rep_seed)
  df_full <- sim$df

  log_line(sprintf("  n_obs = %d, obs/species: mean=%.1f, range=[%d, %d]",
                   nrow(df_full),
                   mean(sim$n_obs_per),
                   min(sim$n_obs_per),
                   max(sim$n_obs_per)))

  # 2. Mask species and within-species observations
  mask_info <- mask_species(df_full, sim, sp_missing_frac,
                            within_obs_miss, rep_seed)
  df_masked <- mask_info$df_masked
  n_missing <- sum(is.na(df_masked$CTmax))
  log_line(sprintf("  Missing obs: %d / %d (%.1f%%), missing species: %d / %d",
                   n_missing, nrow(df_masked),
                   100 * n_missing / nrow(df_masked),
                   length(mask_info$missing_sp), n_species))

  # Truth vectors for held-out observations
  held_out_idx <- mask_info$all_na_rows
  truth_obs <- df_full$CTmax[held_out_idx]

  # True species means for missing species
  truth_sp_mean <- sim$sp_true_mean

  results <- list()

  # ---- Method 1: species_mean -----------------------------------------------
  log_line("  [1/3] species_mean...")
  t0 <- proc.time()[["elapsed"]]
  m1 <- tryCatch({
    pred_all <- species_mean_impute(df_masked)
    pred_obs <- pred_all[held_out_idx]
    sp_pred_mean <- tapply(pred_all, df_masked$species, mean, na.rm = TRUE)
    # For missing species, the prediction is the grand mean for all obs
    metrics <- compute_metrics(truth_obs, pred_obs,
                               truth_sp_mean, sp_pred_mean[sim$sp],
                               mask_info$missing_sp)
    metrics$method <- "species_mean"
    metrics
  }, error = function(e) {
    log_line(sprintf("    ERROR: %s", e$message))
    data.frame(obs_rmse = NA_real_, sp_rmse = NA_real_,
               obs_pearson_r = NA_real_, method = "species_mean",
               stringsAsFactors = FALSE)
  })
  t1 <- proc.time()[["elapsed"]]
  log_line(sprintf("    Done in %.1fs", t1 - t0))
  results[[1L]] <- m1

  # ---- Method 2: pigauto_no_cov --------------------------------------------
  log_line("  [2/3] pigauto_no_cov...")
  t0 <- proc.time()[["elapsed"]]
  m2 <- tryCatch({
    # Prepare traits data.frame (species + CTmax only)
    traits_nocov <- df_masked[, c("species", "CTmax"), drop = FALSE]

    result_nocov <- impute(
      traits     = traits_nocov,
      tree       = tree,
      species_col = "species",
      log_transform = FALSE,
      missing_frac = 0.0,    # don't add extra splits; we mask manually
      epochs     = epochs,
      verbose    = FALSE,
      seed       = rep_seed
    )

    pred_obs <- result_nocov$completed$CTmax[held_out_idx]
    sp_pred_mean <- tapply(result_nocov$completed$CTmax,
                           result_nocov$completed$species,
                           mean, na.rm = TRUE)
    metrics <- compute_metrics(truth_obs, pred_obs,
                               truth_sp_mean, sp_pred_mean[sim$sp],
                               mask_info$missing_sp)
    metrics$method <- "pigauto_no_cov"
    metrics
  }, error = function(e) {
    log_line(sprintf("    ERROR: %s", e$message))
    data.frame(obs_rmse = NA_real_, sp_rmse = NA_real_,
               obs_pearson_r = NA_real_, method = "pigauto_no_cov",
               stringsAsFactors = FALSE)
  })
  t1 <- proc.time()[["elapsed"]]
  log_line(sprintf("    Done in %.1fs", t1 - t0))
  results[[2L]] <- m2

  # ---- Method 3: pigauto_cov -----------------------------------------------
  log_line("  [3/3] pigauto_cov...")
  t0 <- proc.time()[["elapsed"]]
  m3 <- tryCatch({
    # Prepare traits data.frame (species + CTmax)
    traits_cov <- df_masked[, c("species", "CTmax"), drop = FALSE]
    # Covariates: acclimation temperature (fully observed)
    covs <- data.frame(acclim_temp = df_masked$acclim_temp)

    result_cov <- impute(
      traits     = traits_cov,
      tree       = tree,
      species_col = "species",
      covariates  = covs,
      log_transform = FALSE,
      missing_frac = 0.0,
      epochs     = epochs,
      verbose    = FALSE,
      seed       = rep_seed
    )

    pred_obs <- result_cov$completed$CTmax[held_out_idx]
    sp_pred_mean <- tapply(result_cov$completed$CTmax,
                           result_cov$completed$species,
                           mean, na.rm = TRUE)
    metrics <- compute_metrics(truth_obs, pred_obs,
                               truth_sp_mean, sp_pred_mean[sim$sp],
                               mask_info$missing_sp)
    metrics$method <- "pigauto_cov"
    metrics
  }, error = function(e) {
    log_line(sprintf("    ERROR: %s", e$message))
    data.frame(obs_rmse = NA_real_, sp_rmse = NA_real_,
               obs_pearson_r = NA_real_, method = "pigauto_cov",
               stringsAsFactors = FALSE)
  })
  t1 <- proc.time()[["elapsed"]]
  log_line(sprintf("    Done in %.1fs", t1 - t0))
  results[[3L]] <- m3

  # Combine
  out <- do.call(rbind, results)
  out$lambda            <- lambda
  out$beta              <- beta
  out$sp_missing_frac   <- sp_missing_frac
  out$rep               <- rep_id
  out$n_obs             <- nrow(df_full)
  out$n_held_out        <- length(held_out_idx)
  out$n_missing_sp      <- length(mask_info$missing_sp)
  out
}


# -------------------------------------------------------------------------
# Build the full grid of cells
# -------------------------------------------------------------------------

cells <- expand.grid(
  lambda          = lambdas,
  beta            = betas,
  sp_missing_frac = sp_missing_fracs,
  rep_id          = seq_len(n_reps),
  stringsAsFactors = FALSE
)
cells$cell_key <- sprintf("%.1f/%.1f/%.1f/%d",
                          cells$lambda, cells$beta,
                          cells$sp_missing_frac, cells$rep_id)

n_cells <- nrow(cells)
log_line(sprintf("Total cells to run: %d", n_cells))
log_line(sprintf("Grid: %d lambda x %d beta x %d sp_miss x %d reps",
                 length(lambdas), length(betas),
                 length(sp_missing_fracs), n_reps))

# -------------------------------------------------------------------------
# Sequential execution (multi-obs + torch is not fork-safe)
# -------------------------------------------------------------------------

all_results <- vector("list", n_cells)

for (ci in seq_len(n_cells)) {
  row <- cells[ci, ]
  log_line(sprintf("\n--- Cell %d / %d ---", ci, n_cells))

  cell_result <- tryCatch(
    run_one_cell(row$lambda, row$beta, row$sp_missing_frac, row$rep_id),
    error = function(e) {
      log_line(sprintf("  CELL ERROR: %s", e$message))
      data.frame(
        obs_rmse = NA_real_, sp_rmse = NA_real_, obs_pearson_r = NA_real_,
        method = "ERROR", lambda = row$lambda, beta = row$beta,
        sp_missing_frac = row$sp_missing_frac, rep = row$rep_id,
        n_obs = NA_integer_, n_held_out = NA_integer_,
        n_missing_sp = NA_integer_,
        stringsAsFactors = FALSE
      )
    }
  )

  all_results[[ci]] <- cell_result

  # Periodic save (every 4 cells)
  if (ci %% 4L == 0L || ci == n_cells) {
    partial <- do.call(rbind, all_results[seq_len(ci)])
    rownames(partial) <- NULL
    saveRDS(list(results = partial, partial = TRUE, cells_done = ci),
            out_rds)
    log_line(sprintf("  Saved partial results (%d / %d cells)", ci, n_cells))
  }

  # Free memory between cells
  invisible(gc(verbose = FALSE))
}


# -------------------------------------------------------------------------
# Finalise
# -------------------------------------------------------------------------

results <- do.call(rbind, all_results)
rownames(results) <- NULL

total_wall <- proc.time()[["elapsed"]] - script_start

log_line(sprintf("Total wall: %.1fs (%.1f min)", total_wall, total_wall / 60))
log_line(sprintf("Total rows: %d", nrow(results)))

# Filter out error rows
good <- results[results$method != "ERROR", ]
errs <- results[results$method == "ERROR", ]
if (nrow(errs) > 0L) {
  log_line(sprintf("WARNING: %d error rows", nrow(errs)))
}

# Save final RDS
saveRDS(list(
  results          = good,
  errors           = if (nrow(errs) > 0L) errs else NULL,
  total_wall       = total_wall,
  n_species        = n_species,
  n_reps           = n_reps,
  epochs           = epochs,
  lambdas          = lambdas,
  betas            = betas,
  sp_missing_fracs = sp_missing_fracs,
  mu_ctmax         = mu_ctmax,
  sigma_residual   = sigma_residual,
  acclim_mean      = acclim_mean,
  acclim_sd        = acclim_sd,
  obs_lambda_pois  = obs_lambda_pois,
  within_obs_miss  = within_obs_miss,
  timestamp        = Sys.time(),
  commit = tryCatch(system("git rev-parse HEAD", intern = TRUE),
                    error = function(e) "unknown")
), out_rds)
log_line(sprintf("Wrote %s", out_rds))

# -------------------------------------------------------------------------
# Markdown summary
# -------------------------------------------------------------------------

machine <- tryCatch(
  sprintf("%s %s (%s), R %s",
          Sys.info()[["sysname"]], Sys.info()[["release"]],
          Sys.info()[["machine"]],
          paste(R.version$major, R.version$minor, sep = ".")),
  error = function(e) "machine info unavailable"
)

md <- character()
md <- c(md, "# Multi-observation-per-species benchmark")
md <- c(md, "")
md <- c(md, sprintf("Run on: %s", format(Sys.time())))
md <- c(md, sprintf("Machine: %s", machine))
md <- c(md, sprintf("Species: %d, reps: %d, epochs: %d",
                     n_species, n_reps, epochs))
md <- c(md, sprintf("Total wall: %.1f min", total_wall / 60))
md <- c(md, "")
md <- c(md, "## Data-generating process")
md <- c(md, "")
md <- c(md, "```")
md <- c(md, "CTmax_ij = mu + phylo_i + beta * acclim_temp_ij + epsilon_ij")
md <- c(md, sprintf("mu = %.1f, sigma_residual = %.2f", mu_ctmax, sigma_residual))
md <- c(md, sprintf("acclim_temp ~ N(%.1f, %.1f)", acclim_mean, acclim_sd))
md <- c(md, sprintf("obs_per_species ~ Pois(%d), clamped [%d, %d]",
                     obs_lambda_pois, obs_min, obs_max))
md <- c(md, "```")
md <- c(md, "")
md <- c(md, "## Methods")
md <- c(md, "")
md <- c(md, "- **species_mean**: mean of observed CTmax per species; grand mean for unobserved species.")
md <- c(md, "- **pigauto_no_cov**: `impute(traits, tree, species_col = \"species\")`")
md <- c(md, "- **pigauto_cov**: `impute(traits, tree, species_col = \"species\", covariates = acclim_temp)`")
md <- c(md, "")

if (nrow(good) > 0L) {
  # ---- Observation-level RMSE summary ---
  md <- c(md, "## Observation-level RMSE (lower is better)")
  md <- c(md, "")
  rmse_agg <- aggregate(obs_rmse ~ method + lambda + beta + sp_missing_frac,
                        data = good, FUN = mean, na.rm = TRUE)
  rmse_agg <- rmse_agg[order(rmse_agg$lambda, rmse_agg$beta,
                              rmse_agg$sp_missing_frac,
                              match(rmse_agg$method,
                                    c("species_mean", "pigauto_no_cov",
                                      "pigauto_cov"))), ]

  md <- c(md, "| method | lambda | beta | sp_miss | obs_RMSE |")
  md <- c(md, "|--------|--------|------|---------|----------|")
  for (i in seq_len(nrow(rmse_agg))) {
    md <- c(md, sprintf("| %s | %.1f | %.1f | %.1f | %.4f |",
                        rmse_agg$method[i], rmse_agg$lambda[i],
                        rmse_agg$beta[i], rmse_agg$sp_missing_frac[i],
                        rmse_agg$obs_rmse[i]))
  }

  # ---- Species-level RMSE summary ---
  md <- c(md, "")
  md <- c(md, "## Species-level RMSE (unobserved species, lower is better)")
  md <- c(md, "")
  sp_rmse_agg <- aggregate(sp_rmse ~ method + lambda + beta + sp_missing_frac,
                           data = good, FUN = mean, na.rm = TRUE)
  sp_rmse_agg <- sp_rmse_agg[order(sp_rmse_agg$lambda, sp_rmse_agg$beta,
                                    sp_rmse_agg$sp_missing_frac,
                                    match(sp_rmse_agg$method,
                                          c("species_mean", "pigauto_no_cov",
                                            "pigauto_cov"))), ]

  md <- c(md, "| method | lambda | beta | sp_miss | sp_RMSE |")
  md <- c(md, "|--------|--------|------|---------|---------|")
  for (i in seq_len(nrow(sp_rmse_agg))) {
    md <- c(md, sprintf("| %s | %.1f | %.1f | %.1f | %.4f |",
                        sp_rmse_agg$method[i], sp_rmse_agg$lambda[i],
                        sp_rmse_agg$beta[i], sp_rmse_agg$sp_missing_frac[i],
                        sp_rmse_agg$sp_rmse[i]))
  }

  # ---- Observation-level Pearson r ---
  md <- c(md, "")
  md <- c(md, "## Observation-level Pearson r (higher is better)")
  md <- c(md, "")
  r_agg <- aggregate(obs_pearson_r ~ method + lambda + beta + sp_missing_frac,
                     data = good, FUN = mean, na.rm = TRUE)
  r_agg <- r_agg[order(r_agg$lambda, r_agg$beta,
                        r_agg$sp_missing_frac,
                        match(r_agg$method,
                              c("species_mean", "pigauto_no_cov",
                                "pigauto_cov"))), ]

  md <- c(md, "| method | lambda | beta | sp_miss | pearson_r |")
  md <- c(md, "|--------|--------|------|---------|-----------|")
  for (i in seq_len(nrow(r_agg))) {
    md <- c(md, sprintf("| %s | %.1f | %.1f | %.1f | %.4f |",
                        r_agg$method[i], r_agg$lambda[i],
                        r_agg$beta[i], r_agg$sp_missing_frac[i],
                        r_agg$obs_pearson_r[i]))
  }

  # ---- Covariate lift summary ---
  md <- c(md, "")
  md <- c(md, "## Covariate lift (pigauto_cov / pigauto_no_cov RMSE ratio)")
  md <- c(md, "")
  md <- c(md, "Ratio < 1 means covariates help; ratio > 1 means they hurt.")
  md <- c(md, "")

  rmse_nocov <- rmse_agg[rmse_agg$method == "pigauto_no_cov", ]
  rmse_cov   <- rmse_agg[rmse_agg$method == "pigauto_cov", ]
  if (nrow(rmse_nocov) > 0L && nrow(rmse_cov) > 0L) {
    merged <- merge(rmse_nocov, rmse_cov,
                    by = c("lambda", "beta", "sp_missing_frac"),
                    suffixes = c("_nocov", "_cov"))
    merged$ratio <- merged$obs_rmse_cov / merged$obs_rmse_nocov

    md <- c(md, "| lambda | beta | sp_miss | RMSE_nocov | RMSE_cov | ratio |")
    md <- c(md, "|--------|------|---------|------------|----------|-------|")
    for (i in seq_len(nrow(merged))) {
      md <- c(md, sprintf("| %.1f | %.1f | %.1f | %.4f | %.4f | %.3f |",
                          merged$lambda[i], merged$beta[i],
                          merged$sp_missing_frac[i],
                          merged$obs_rmse_nocov[i],
                          merged$obs_rmse_cov[i],
                          merged$ratio[i]))
    }
  }

} else {
  md <- c(md, "## Results")
  md <- c(md, "")
  md <- c(md, "No successful results to report.")
}

md <- c(md, "")
md <- c(md, "---")
md <- c(md, sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M")))

writeLines(md, out_md)
log_line(sprintf("Wrote %s", out_md))
log_line("done")
