#!/usr/bin/env Rscript
#
# script/bench_covariate_sim.R
#
# Benchmark: when do environmental covariates improve imputation?
#
# Purpose
#   Demonstrate that covariates improve imputation accuracy when
#   trait variation has a substantial environmental component beyond
#   phylogenetic signal. Complements the Delhey benchmark (where
#   high phylogenetic signal makes covariates redundant).
#
# Design
#   - Tree:       tree300 (300 species, BirdTree Hackett MCC)
#   - Traits:     3 continuous, simulated per scenario
#   - Covariates:  3 uncorrelated normal (fully observed)
#   - Scenarios:   4 (lambda, beta) combinations
#     (a) high_phylo_no_env:    lambda=0.9, beta=0.0
#     (b) high_phylo_mod_env:   lambda=0.7, beta=0.5
#     (c) mod_phylo_strong_env: lambda=0.3, beta=1.0
#     (d) low_phylo_strong_env: lambda=0.1, beta=1.5
#   - Missingness: 20%, 40%, 60% MCAR
#   - Methods:     pigauto (no covs) vs pigauto_covs (with covs)
#   - Replicates:  3 per cell
#   - Metrics:     RMSE, Pearson r on held-out cells
#
# Trait generation
#   y = beta * B %*% env + sqrt(lambda) * bm + sqrt(1-lambda) * eps
#   where B is a 3x3 coefficient matrix, bm ~ BM(tree), eps ~ N(0,1)
#
# Output
#   script/bench_covariate_sim.rds
#   script/bench_covariate_sim.md

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto",
    quiet = TRUE
  )
})

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "bench_covariate_sim.rds")
out_md  <- file.path(here, "script", "bench_covariate_sim.md")

script_start <- proc.time()[["elapsed"]]

log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start), ..., "\n",
      sep = "")
  flush.console()
}

# -------------------------------------------------------------------------
# Data
# -------------------------------------------------------------------------

data(tree300, package = "pigauto")
n <- length(tree300$tip.label)
spp <- tree300$tip.label

# -------------------------------------------------------------------------
# Simulate BM traits on the tree
# -------------------------------------------------------------------------

simulate_bm <- function(tree, n_traits, seed) {
  set.seed(seed)
  C <- ape::vcv(tree)
  R <- stats::cov2cor(C)
  L <- chol(R)
  # Each column is an independent BM realisation (unit variance)
  z <- matrix(rnorm(n * n_traits), nrow = n, ncol = n_traits)
  bm <- t(L) %*% z
  rownames(bm) <- tree$tip.label
  bm
}

# -------------------------------------------------------------------------
# Parameters
# -------------------------------------------------------------------------

scenarios <- list(
  list(name = "high_phylo_no_env",    lambda = 0.9, beta = 0.0),
  list(name = "high_phylo_mod_env",   lambda = 0.7, beta = 0.5),
  list(name = "mod_phylo_strong_env", lambda = 0.3, beta = 1.0),
  list(name = "low_phylo_strong_env", lambda = 0.1, beta = 1.5)
)

miss_fracs <- c(0.20, 0.40, 0.60)
n_reps     <- 3L
n_traits   <- 3L
n_covs     <- 3L
epochs     <- 300L

# -------------------------------------------------------------------------
# Helper: punch holes
# -------------------------------------------------------------------------

punch_holes <- function(df, frac, seed) {
  set.seed(seed)
  for (col in names(df)) {
    n_miss <- round(nrow(df) * frac)
    idx <- sample.int(nrow(df), n_miss)
    df[[col]][idx] <- NA
  }
  df
}

# -------------------------------------------------------------------------
# Helper: compute metrics on held-out cells
# -------------------------------------------------------------------------

compute_metrics <- function(truth, imputed, mask) {
  # mask: logical matrix, TRUE = was missing (imputed)
  t_vec <- truth[mask]
  p_vec <- imputed[mask]
  ok <- is.finite(t_vec) & is.finite(p_vec)
  t_vec <- t_vec[ok]
  p_vec <- p_vec[ok]
  if (length(t_vec) < 3) return(list(rmse = NA, pearson_r = NA, n = 0L))
  list(
    rmse      = sqrt(mean((t_vec - p_vec)^2)),
    pearson_r = cor(t_vec, p_vec),
    n         = length(t_vec)
  )
}

# -------------------------------------------------------------------------
# Main loop
# -------------------------------------------------------------------------

all_results <- list()
result_idx  <- 0L

for (sc in scenarios) {
  for (frac in miss_fracs) {
    for (rep_id in seq_len(n_reps)) {
      rep_seed <- rep_id * 1000L + round(frac * 100) +
        match(sc$name, vapply(scenarios, `[[`, character(1), "name")) * 10000L

      log_line(sprintf("=== %s, frac=%.2f, rep=%d ===",
                        sc$name, frac, rep_id))

      # ---- Generate data ----
      set.seed(rep_seed)

      # BM component (unit variance)
      bm <- simulate_bm(tree300, n_traits, seed = rep_seed)

      # Environmental covariates (fully observed, uncorrelated normal)
      env <- matrix(rnorm(n * n_covs), nrow = n, ncol = n_covs)
      rownames(env) <- spp
      colnames(env) <- paste0("env", seq_len(n_covs))

      # Coefficient matrix: each trait loads on covariates
      # B is n_traits x n_covs, fixed per scenario (same for all reps)
      set.seed(sc$beta * 100 + 7)
      B <- matrix(rnorm(n_traits * n_covs), nrow = n_traits, ncol = n_covs)

      # Independent noise
      eps <- matrix(rnorm(n * n_traits), nrow = n, ncol = n_traits)

      # Trait generation: y = beta * env %*% t(B) + sqrt(lambda) * bm + sqrt(1-lambda) * eps
      y <- sc$beta * (env %*% t(B)) +
           sqrt(sc$lambda) * bm +
           sqrt(1 - sc$lambda) * eps

      colnames(y) <- paste0("trait", seq_len(n_traits))
      rownames(y) <- spp

      # Full data
      df_full <- as.data.frame(y)
      rownames(df_full) <- spp

      # Covariates data frame
      cov_df <- as.data.frame(env)
      rownames(cov_df) <- spp

      # Introduce MCAR missingness
      df_miss <- punch_holes(df_full, frac, seed = rep_seed + 999)
      mask <- is.na(as.matrix(df_miss))
      n_na <- sum(mask)
      log_line(sprintf("  Missing: %d / %d (%.1f%%)",
                       n_na, length(mask), 100 * n_na / length(mask)))

      # ---- pigauto WITHOUT covariates ----
      log_line("  pigauto (no covariates)...")
      t0 <- proc.time()[["elapsed"]]
      res_no_cov <- tryCatch({
        impute(df_miss, tree300, epochs = epochs, verbose = FALSE,
               eval_every = 25L, patience = 15L, seed = rep_seed)
      }, error = function(e) {
        log_line("    ERROR: ", conditionMessage(e))
        NULL
      })
      t1 <- proc.time()[["elapsed"]]
      log_line(sprintf("    Done in %.1fs", t1 - t0))

      # ---- pigauto WITH covariates ----
      log_line("  pigauto (with covariates)...")
      t0 <- proc.time()[["elapsed"]]
      res_cov <- tryCatch({
        impute(df_miss, tree300, covariates = cov_df,
               epochs = epochs, verbose = FALSE,
               eval_every = 25L, patience = 15L, seed = rep_seed)
      }, error = function(e) {
        log_line("    ERROR: ", conditionMessage(e))
        NULL
      })
      t1 <- proc.time()[["elapsed"]]
      log_line(sprintf("    Done in %.1fs", t1 - t0))

      # ---- Compute metrics ----
      truth_mat <- as.matrix(df_full)

      for (method_name in c("pigauto", "pigauto_covs")) {
        res_obj <- if (method_name == "pigauto") res_no_cov else res_cov
        if (is.null(res_obj)) next

        # Extract completed values for trait columns
        completed <- res_obj$completed
        imputed_mat <- as.matrix(completed[, colnames(df_full), drop = FALSE])

        for (j in seq_len(n_traits)) {
          trait_name <- paste0("trait", j)
          t_j <- truth_mat[, trait_name]
          p_j <- imputed_mat[, trait_name]
          m_j <- mask[, trait_name]

          met <- compute_metrics(
            matrix(t_j, ncol = 1),
            matrix(p_j, ncol = 1),
            matrix(m_j, ncol = 1)
          )

          result_idx <- result_idx + 1L
          all_results[[result_idx]] <- data.frame(
            scenario     = sc$name,
            lambda       = sc$lambda,
            beta         = sc$beta,
            missing_frac = frac,
            rep          = rep_id,
            method       = method_name,
            trait        = trait_name,
            rmse         = met$rmse,
            pearson_r    = met$pearson_r,
            n_imputed    = met$n,
            stringsAsFactors = FALSE
          )
        }
      }

      log_line("  Cell done.\n")
    }
  }
}

# -------------------------------------------------------------------------
# Assemble results
# -------------------------------------------------------------------------

results <- do.call(rbind, all_results)
rownames(results) <- NULL

log_line("Total rows: ", nrow(results))

# -------------------------------------------------------------------------
# Save RDS
# -------------------------------------------------------------------------

meta <- list(
  n_species  = n,
  n_traits   = n_traits,
  n_covs     = n_covs,
  epochs     = epochs,
  scenarios  = scenarios,
  miss_fracs = miss_fracs,
  n_reps     = n_reps,
  timestamp  = Sys.time(),
  wall_time  = proc.time()[["elapsed"]] - script_start
)

saveRDS(list(results = results, meta = meta), out_rds)
log_line("Wrote ", out_rds)

# -------------------------------------------------------------------------
# Markdown summary
# -------------------------------------------------------------------------

md <- character()
md <- c(md, "# Covariate effectiveness simulation benchmark\n")
md <- c(md, sprintf("- Species: %d (tree300)", n))
md <- c(md, sprintf("- Traits: %d continuous (simulated)", n_traits))
md <- c(md, sprintf("- Covariates: %d (uncorrelated normal)", n_covs))
md <- c(md, sprintf("- Scenarios: %d (lambda x beta combinations)", length(scenarios)))
md <- c(md, sprintf("- Missingness: %s", paste(miss_fracs, collapse = ", ")))
md <- c(md, sprintf("- Replicates: %d", n_reps))
md <- c(md, sprintf("- Total wall time: %.1f min\n", meta$wall_time / 60))

# RMSE table
if (nrow(results) > 0) {
  md <- c(md, "## Mean RMSE by scenario and method\n")
  md <- c(md, "| scenario | lambda | beta | miss_frac | method | mean_RMSE | mean_r |")
  md <- c(md, "|----------|--------|------|-----------|--------|-----------|--------|")

  agg <- aggregate(
    cbind(rmse, pearson_r) ~ scenario + lambda + beta + missing_frac + method,
    data = results, FUN = mean, na.rm = TRUE
  )
  agg <- agg[order(agg$lambda, -agg$beta, agg$missing_frac, agg$method), ]

  for (i in seq_len(nrow(agg))) {
    md <- c(md, sprintf("| %s | %.1f | %.1f | %.2f | %s | %.4f | %.4f |",
                        agg$scenario[i], agg$lambda[i], agg$beta[i],
                        agg$missing_frac[i], agg$method[i],
                        agg$rmse[i], agg$pearson_r[i]))
  }

  # Covariate lift ratio
  md <- c(md, "\n## Covariate lift (pigauto_covs / pigauto RMSE ratio)\n")
  md <- c(md, "| scenario | lambda | beta | miss_frac | RMSE_no_cov | RMSE_covs | ratio |")
  md <- c(md, "|----------|--------|------|-----------|-------------|-----------|-------|")

  rmse_no  <- aggregate(rmse ~ scenario + lambda + beta + missing_frac,
                        data = results[results$method == "pigauto", ],
                        FUN = mean, na.rm = TRUE)
  rmse_cov <- aggregate(rmse ~ scenario + lambda + beta + missing_frac,
                        data = results[results$method == "pigauto_covs", ],
                        FUN = mean, na.rm = TRUE)
  merged <- merge(rmse_no, rmse_cov,
                  by = c("scenario", "lambda", "beta", "missing_frac"),
                  suffixes = c("_no", "_cov"))
  merged$ratio <- merged$rmse_cov / merged$rmse_no
  merged <- merged[order(merged$lambda, -merged$beta, merged$missing_frac), ]

  for (i in seq_len(nrow(merged))) {
    md <- c(md, sprintf("| %s | %.1f | %.1f | %.2f | %.4f | %.4f | %.3f |",
                        merged$scenario[i], merged$lambda[i], merged$beta[i],
                        merged$missing_frac[i], merged$rmse_no[i],
                        merged$rmse_cov[i], merged$ratio[i]))
  }
}

md <- c(md, "\n---\nGenerated:", format(Sys.time(), "%Y-%m-%d %H:%M"))
writeLines(md, out_md)
log_line("Wrote ", out_md)
log_line("done")
