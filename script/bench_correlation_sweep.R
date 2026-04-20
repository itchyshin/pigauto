#!/usr/bin/env Rscript
#
# script/bench_correlation_sweep.R
#
# Phase 8.1: cross-trait-correlation (ρ) sweep on continuous traits.
# Does pigauto's joint MVN baseline lift over per-column BM as ρ grows?
#
# Design
#   - Tree:    ape::rcoal(300).
#   - Traits:  4 continuous traits from a joint MVN with cross-trait
#              covariance ρ (and diag = 1).
#   - Pagel's λ fixed at 1.0 (pure BM on tree) — we isolate the ρ effect.
#   - ρ grid: {0, 0.2, 0.4, 0.6, 0.8}.
#   - Reps:   3.
#   - Methods: mean_baseline, pigauto_default, pigauto_em5.
#   - miss_frac: 0.30, MCAR per trait.
#
# Output
#   script/bench_correlation_sweep.{rds,md}
#
# Runtime: ~15 min local CPU.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  library(MASS)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto",
    quiet = TRUE
  )
})

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "bench_correlation_sweep.rds")
out_md  <- file.path(here, "script", "bench_correlation_sweep.md")

CONFIG <- list(
  rhos      = c(0.0, 0.2, 0.4, 0.6, 0.8),
  n_species = 300L,
  n_traits  = 4L,
  miss_frac = 0.30,
  n_reps    = 3L,
  methods   = c("mean_baseline", "pigauto_default", "pigauto_em5"),
  epochs    = 200L
)

# -------------------------------------------------------------------------
# Data generator: K continuous traits with cross-trait ρ under BM
# -------------------------------------------------------------------------

# Kronecker: Σ_total = Σ_trait ⊗ V_tree, sample via column-wise independent
# draws then linearly combine.  Simpler: for each tip, draw a K-vector from
# MVN(0, Σ_trait) with variance 1; multiply across tips by the tree
# Cholesky factor to induce phylogenetic correlation.
sim_correlated_mvn_traits <- function(tree, rho, K, seed) {
  set.seed(seed)
  n <- length(tree$tip.label)
  V_tree  <- ape::vcv(tree)
  L_tree  <- t(chol(V_tree + diag(1e-8, n)))
  Sigma_tr <- matrix(rho, K, K) + diag(1 - rho, K)
  L_trait  <- t(chol(Sigma_tr + diag(1e-8, K)))
  # Z is n × K iid standard normal. Trait-side correlation: Z %*% t(L_trait).
  # Tree-side correlation: L_tree %*% (above).
  Z    <- matrix(stats::rnorm(n * K), n, K)
  Z_tr <- Z %*% t(L_trait)
  out  <- L_tree %*% Z_tr
  rownames(out) <- tree$tip.label
  colnames(out) <- paste0("c", seq_len(K))
  as.data.frame(out)
}

inject_mcar <- function(df, miss_frac, seed) {
  set.seed(seed)
  n <- nrow(df)
  mask <- lapply(names(df), function(v) sample.int(n, ceiling(n * miss_frac)))
  names(mask) <- names(df)
  for (v in names(df)) df[mask[[v]], v] <- NA
  list(data = df, mask = mask)
}

method_mean <- function(df_miss, tree, mask) {
  out <- df_miss
  for (v in names(out)) out[[v]][mask[[v]]] <- mean(out[[v]], na.rm = TRUE)
  out
}

method_pigauto <- function(df_miss, tree, em_iter, seed) {
  res <- pigauto::impute(df_miss, tree, log_transform = FALSE,
                           missing_frac = 0.25, n_imputations = 1L,
                           epochs = CONFIG$epochs, verbose = FALSE,
                           seed = seed, em_iterations = em_iter)
  res$completed
}

eval_cell <- function(truth, completed, mask) {
  rows <- list()
  for (v in names(truth)) {
    idx <- mask[[v]]
    if (length(idx) == 0L) next
    rmse <- sqrt(mean((truth[[v]][idx] - completed[[v]][idx])^2, na.rm = TRUE))
    pear <- suppressWarnings(stats::cor(truth[[v]][idx], completed[[v]][idx],
                                          use = "complete.obs"))
    rows[[length(rows) + 1L]] <- data.frame(trait = v,
                                              metric = c("rmse", "pearson_r"),
                                              value  = c(rmse, pear))
  }
  do.call(rbind, rows)
}

cat(sprintf("Phase 8.1 ρ sweep: %d ρ × %d reps × %d methods\n",
            length(CONFIG$rhos), CONFIG$n_reps, length(CONFIG$methods)))
results <- list()
t_start <- proc.time()[["elapsed"]]

for (rep_i in seq_len(CONFIG$n_reps)) {
  set.seed(rep_i)
  tree <- ape::rcoal(CONFIG$n_species,
                      tip.label = paste0("sp", seq_len(CONFIG$n_species)))
  for (rho in CONFIG$rhos) {
    sim_seed <- 1000L * rep_i + as.integer(rho * 10)
    truth <- sim_correlated_mvn_traits(tree, rho, CONFIG$n_traits, sim_seed)
    miss  <- inject_mcar(truth, CONFIG$miss_frac, sim_seed + 1)

    for (meth in CONFIG$methods) {
      t0 <- proc.time()[["elapsed"]]
      completed <- tryCatch(switch(meth,
        mean_baseline   = method_mean(miss$data, tree, miss$mask),
        pigauto_default = method_pigauto(miss$data, tree, 0L, sim_seed + 2),
        pigauto_em5     = method_pigauto(miss$data, tree, 5L, sim_seed + 2)
      ), error = function(e) { message("  ", meth, " @ ρ=", rho, ": ",
                                          conditionMessage(e)); NULL })
      wall <- proc.time()[["elapsed"]] - t0
      if (is.null(completed)) next

      ev <- eval_cell(truth, completed, miss$mask)
      ev$method <- meth; ev$rho <- rho; ev$rep <- rep_i; ev$wall_s <- wall
      results[[length(results) + 1L]] <- ev
      cat(sprintf("  rep=%d ρ=%.1f %-18s %.1fs\n", rep_i, rho, meth, wall))
      saveRDS(list(results = do.call(rbind, results),
                    config = CONFIG,
                    script_wall = proc.time()[["elapsed"]] - t_start), out_rds)
    }
  }
}

all_res <- do.call(rbind, results)
saveRDS(list(results = all_res, config = CONFIG,
              script_wall = proc.time()[["elapsed"]] - t_start), out_rds)

summary_tbl <- aggregate(value ~ method + rho + metric,
                          data = all_res, FUN = mean)
md <- c(
  "# Phase 8.1: cross-trait-correlation (ρ) sweep",
  "",
  sprintf("n=%d, K=%d continuous traits, miss_frac=%.2f, n_reps=%d",
          CONFIG$n_species, CONFIG$n_traits, CONFIG$miss_frac, CONFIG$n_reps),
  sprintf("Total wall: %.1f min",
          (proc.time()[["elapsed"]] - t_start) / 60),
  "",
  "## Means per (method, ρ, metric), averaged over traits and reps",
  "",
  "```",
  capture.output(print(summary_tbl, row.names = FALSE)),
  "```"
)
writeLines(md, out_md)
cat("\n=== DONE === ", out_rds, " ", out_md, "\n")
