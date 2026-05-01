#!/usr/bin/env Rscript
# script/diag_bace_predictor_phylo.R
#
# Diagnostic: does BACE-sim multi-obs regression depend on predictor
# phylogenetic signal?
#
# Hypothesis: BACE generates predictors with `phylo_signal[i] = 0.3` baked in
# (default in sim_bace_dgp.R). The predictor's between-species structure is
# correlated with the response's phylo random effect, double-counting through
# pigauto's BM baseline. If we DROP the predictor's phylo signal to 0, pigauto
# should improve.
#
# Sweep:
#   predictor_phylo_signal in {0.0, 0.3, 0.6}
#   beta_resp_strength      in {0.0, 0.5}
#   multi_obs_ratio         = 4
#   n_species               = 100
#   response_type           = "gaussian"
#   response_phylo_signal   = 0.4 (held constant)
#   reps                    = 2
#
# Methods compared (cheap subset):
#   - column_mean
#   - pigauto_no_cov
#   - pigauto_cov_sfT
#
# Total cells: 3 * 2 * 2 = 12 cells, each ~30-60s --> 6-12 min wall.

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  pkg_path <- "/Users/z3437171/Dropbox/Github Local/pigauto"
  devtools::load_all(pkg_path, quiet = TRUE)
  source(file.path(pkg_path, "script", "sim_bace_dgp.R"))
})

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "diag_bace_predictor_phylo.rds")
out_md  <- file.path(here, "useful", "diag_bace_predictor_phylo.md")

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - t0),
                              ..., "\n", sep = "")

# Modified sim_bace_dgp wrapper to allow an explicit predictor_phylo_signal
sim_dgp_pred_phylo <- function(n_species, multi_obs_ratio, response_phylo_signal,
                                predictor_phylo_signal, beta_resp_strength,
                                seed, miss_frac = 0.30) {
  n_predictors <- 3L
  predictor_types <- c("gaussian", "gaussian", "binary")
  phylo_vec <- c(response_phylo_signal,
                 rep(predictor_phylo_signal, n_predictors))
  beta_resp <- if (beta_resp_strength == 0) {
    rep(0, n_predictors)
  } else {
    base <- c(0.7, 0.4, 0.5)
    beta_resp_strength * base
  }
  miss_vec <- c(miss_frac, rep(0, n_predictors))

  set.seed(seed)
  sim <- BACE::sim_bace(
    response_type   = "gaussian",
    predictor_types = predictor_types,
    beta_resp       = beta_resp,
    phylo_signal    = phylo_vec,
    n_cases         = as.integer(n_species * multi_obs_ratio),
    n_species       = as.integer(n_species),
    missingness     = miss_vec,
    birth           = 0.8,
    death           = 0.4
  )
  list(
    tree           = sim$tree,
    df_complete    = sim$complete_data,
    df_observed    = sim$data,
    mask           = is.na(sim$data$y),
    response_name  = "y",
    predictor_names = grep("^x[0-9]+$", names(sim$complete_data), value = TRUE)
  )
}

run_methods <- function(d, seed) {
  resp_n <- d$response_name
  pred_n <- d$predictor_names
  df_obs <- d$df_observed
  df_full <- d$df_complete
  tree <- d$tree
  mask <- d$mask
  truth <- df_full[[resp_n]][mask]

  out <- list()

  # column_mean
  cm <- mean(df_obs[[resp_n]], na.rm = TRUE)
  pred <- rep(cm, sum(mask))
  out$column_mean <- list(rmse = sqrt(mean((truth - pred)^2)), n = sum(mask))

  # pigauto_no_cov
  pred_pn <- tryCatch({
    res <- impute(traits = df_obs[, c("species", resp_n), drop = FALSE],
                   tree = tree, species_col = "species",
                   epochs = 200L, verbose = FALSE, seed = seed,
                   missing_frac = 0.0)
    res$completed[[resp_n]][mask]
  }, error = function(e) { log_line("    no_cov ERROR: ", e$message); rep(NA_real_, sum(mask)) })
  out$pigauto_no_cov <- list(rmse = sqrt(mean((truth - pred_pn)^2, na.rm = TRUE)),
                              n = sum(mask))

  # pigauto_cov_sfT (covariates = the three predictors, expanded)
  covs_df <- df_obs[, pred_n, drop = FALSE]
  # Expand factors via model.matrix
  cov_mm <- as.data.frame(model.matrix(~ . - 1, data = covs_df))
  pred_pc <- tryCatch({
    res <- impute(traits = df_obs[, c("species", resp_n), drop = FALSE],
                   tree = tree, species_col = "species",
                   covariates = cov_mm,
                   epochs = 200L, verbose = FALSE, seed = seed,
                   missing_frac = 0.0)
    res$completed[[resp_n]][mask]
  }, error = function(e) { log_line("    cov ERROR: ", e$message); rep(NA_real_, sum(mask)) })
  out$pigauto_cov_sfT <- list(rmse = sqrt(mean((truth - pred_pc)^2, na.rm = TRUE)),
                               n = sum(mask))

  out
}

cells <- expand.grid(
  predictor_phylo_signal = c(0.0, 0.3, 0.6),
  beta_resp_strength     = c(0.0, 0.5),
  rep_id                  = 1:2,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

log_line("Total cells: ", nrow(cells))

results <- list()
for (i in seq_len(nrow(cells))) {
  row <- cells[i, ]
  seed <- 7000L + as.integer(round(row$predictor_phylo_signal * 10)) * 100L +
            as.integer(round(row$beta_resp_strength * 10)) * 10L + row$rep_id
  log_line(sprintf("Cell %d/%d: pred_phylo=%.1f, beta=%.1f, rep=%d (seed=%d)",
                    i, nrow(cells), row$predictor_phylo_signal,
                    row$beta_resp_strength, row$rep_id, seed))

  d <- sim_dgp_pred_phylo(
    n_species = 100L,
    multi_obs_ratio = 4L,
    response_phylo_signal = 0.4,
    predictor_phylo_signal = row$predictor_phylo_signal,
    beta_resp_strength = row$beta_resp_strength,
    seed = seed
  )
  log_line(sprintf("  n_cases = %d, mask = %d", nrow(d$df_observed), sum(d$mask)))

  res <- run_methods(d, seed)
  for (m in names(res)) {
    results[[length(results) + 1L]] <- data.frame(
      predictor_phylo_signal = row$predictor_phylo_signal,
      beta_resp_strength      = row$beta_resp_strength,
      rep                     = row$rep_id,
      method                  = m,
      rmse                    = res[[m]]$rmse,
      n                       = res[[m]]$n,
      stringsAsFactors        = FALSE
    )
  }

  # Save partial
  partial <- do.call(rbind, results)
  saveRDS(partial, out_rds)
  invisible(gc(verbose = FALSE))
}

df <- do.call(rbind, results)
log_line(sprintf("Done. Total wall: %.1fs", proc.time()[["elapsed"]] - t0))

# Markdown summary
md <- character()
md <- c(md, "# BACE multi-obs predictor-phylo-signal diagnostic")
md <- c(md, "")
md <- c(md, sprintf("Generated: %s", format(Sys.time())))
md <- c(md, "")
md <- c(md,
  "Sweep: predictor_phylo_signal in {0.0, 0.3, 0.6} x beta in {0, 0.5} x 2 reps.",
  "Held constant: n_species=100, multi_obs_ratio=4, response_phylo_signal=0.4,",
  "response_type=gaussian, n_predictors=3 (gaus, gaus, bin), miss_frac=0.30.",
  "")

# Per-cell aggregates
agg <- aggregate(rmse ~ method + predictor_phylo_signal + beta_resp_strength,
                  data = df, FUN = mean, na.rm = TRUE)
agg <- agg[order(agg$beta_resp_strength, agg$predictor_phylo_signal,
                  match(agg$method, c("column_mean", "pigauto_no_cov",
                                      "pigauto_cov_sfT"))), ]

md <- c(md, "## Mean RMSE per cell (lower is better)")
md <- c(md, "")
md <- c(md, "| beta | pred_phylo | column_mean | pigauto_no_cov | pigauto_cov |")
md <- c(md, "|------|-----------|-------------|----------------|-------------|")

bs <- unique(agg$beta_resp_strength)
ps <- unique(agg$predictor_phylo_signal)
for (b in bs) for (p in ps) {
  sub <- agg[agg$beta_resp_strength == b & agg$predictor_phylo_signal == p, ]
  cm  <- sub$rmse[sub$method == "column_mean"]
  pn  <- sub$rmse[sub$method == "pigauto_no_cov"]
  pc  <- sub$rmse[sub$method == "pigauto_cov_sfT"]
  md <- c(md, sprintf("| %.1f | %.1f | %.4f | %.4f | %.4f |",
                       b, p, cm, pn, pc))
}

md <- c(md, "")
md <- c(md, "## Diagnostic ratios")
md <- c(md, "")
md <- c(md, "Ratios are pigauto / column_mean. < 1 means pigauto helps.")
md <- c(md, "")
md <- c(md, "| beta | pred_phylo | no_cov / col_mean | cov / col_mean |")
md <- c(md, "|------|-----------|-------------------|---------------|")
for (b in bs) for (p in ps) {
  sub <- agg[agg$beta_resp_strength == b & agg$predictor_phylo_signal == p, ]
  cm  <- sub$rmse[sub$method == "column_mean"]
  pn  <- sub$rmse[sub$method == "pigauto_no_cov"]
  pc  <- sub$rmse[sub$method == "pigauto_cov_sfT"]
  md <- c(md, sprintf("| %.1f | %.1f | %.3f | %.3f |", b, p, pn / cm, pc / cm))
}

writeLines(md, out_md)
log_line(sprintf("Wrote %s", out_md))
log_line(sprintf("Wrote %s", out_rds))
