#!/usr/bin/env Rscript
#
# script/bench_evo_model_sweep.R
#
# Phase 8.2: evolutionary-model sweep. Does pigauto gracefully degrade
# as the data-generating model diverges from BM?
#
# Design
#   - Tree:    ape::rcoal(300).
#   - Traits:  4 continuous (via pigauto's internal simulate_non_bm()).
#   - Models:  BM, OU (α=2), regime_shift (Δ=2 SD), nonlinear.
#   - Reps:    3.
#   - Methods: mean_baseline, pigauto_default, pigauto_em5.
#   - miss_frac: 0.30 MCAR.
#
# Output
#   script/bench_evo_model_sweep.{rds,md}
#
# Runtime: ~15 min.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto",
    quiet = TRUE
  )
})

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "bench_evo_model_sweep.rds")
out_md  <- file.path(here, "script", "bench_evo_model_sweep.md")

CONFIG <- list(
  models    = c("BM", "OU", "regime_shift", "nonlinear"),
  n_species = 300L,
  n_traits  = 4L,
  miss_frac = 0.30,
  n_reps    = 3L,
  methods   = c("mean_baseline", "pigauto_default", "pigauto_em5"),
  epochs    = 200L
)

# -------------------------------------------------------------------------
# Data generator — BM via Cholesky, else delegate to simulate_non_bm()
# -------------------------------------------------------------------------

sim_bm <- function(tree, K, seed) {
  set.seed(seed)
  n <- length(tree$tip.label)
  L <- t(chol(ape::vcv(tree) + diag(1e-8, n)))
  out <- matrix(0, n, K)
  for (j in seq_len(K)) out[, j] <- as.numeric(L %*% stats::rnorm(n))
  rownames(out) <- tree$tip.label
  colnames(out) <- paste0("c", seq_len(K))
  as.data.frame(out)
}

sim_trait_by_model <- function(tree, model, K, seed) {
  if (identical(model, "BM")) return(sim_bm(tree, K, seed))
  # simulate_non_bm is internal; use `:::` to reach it without an export.
  df <- pigauto:::simulate_non_bm(tree, n_traits = K, scenario = model,
                                    seed = seed)
  # Z-score so metrics are comparable across models
  for (v in names(df)) df[[v]] <- as.numeric(scale(df[[v]]))
  df
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
      metric = c("rmse", "pearson_r"), value = c(rmse, pear))
  }
  do.call(rbind, rows)
}

cat(sprintf("Phase 8.2 evo-model sweep: %d models × %d reps × %d methods\n",
            length(CONFIG$models), CONFIG$n_reps, length(CONFIG$methods)))
results <- list()
t_start <- proc.time()[["elapsed"]]

for (rep_i in seq_len(CONFIG$n_reps)) {
  set.seed(rep_i)
  tree <- ape::rcoal(CONFIG$n_species,
                      tip.label = paste0("sp", seq_len(CONFIG$n_species)))
  for (model in CONFIG$models) {
    sim_seed <- 1000L * rep_i + match(model, CONFIG$models)
    truth <- sim_trait_by_model(tree, model, CONFIG$n_traits, sim_seed)
    miss  <- inject_mcar(truth, CONFIG$miss_frac, sim_seed + 1)
    for (meth in CONFIG$methods) {
      t0 <- proc.time()[["elapsed"]]
      completed <- tryCatch(switch(meth,
        mean_baseline   = method_mean(miss$data, tree, miss$mask),
        pigauto_default = method_pigauto(miss$data, tree, 0L, sim_seed + 2),
        pigauto_em5     = method_pigauto(miss$data, tree, 5L, sim_seed + 2)
      ), error = function(e) { message("  ", meth, " @ ", model, ": ",
                                          conditionMessage(e)); NULL })
      wall <- proc.time()[["elapsed"]] - t0
      if (is.null(completed)) next

      ev <- eval_cell(truth, completed, miss$mask)
      ev$method <- meth; ev$model <- model; ev$rep <- rep_i; ev$wall_s <- wall
      results[[length(results) + 1L]] <- ev
      cat(sprintf("  rep=%d model=%-14s %-18s %.1fs\n",
                  rep_i, model, meth, wall))
      saveRDS(list(results = do.call(rbind, results),
                    config = CONFIG,
                    script_wall = proc.time()[["elapsed"]] - t_start), out_rds)
    }
  }
}

all_res <- do.call(rbind, results)
saveRDS(list(results = all_res, config = CONFIG,
              script_wall = proc.time()[["elapsed"]] - t_start), out_rds)

summary_tbl <- aggregate(value ~ method + model + metric,
                          data = all_res, FUN = mean)
md <- c(
  "# Phase 8.2: evolutionary-model sweep",
  "",
  sprintf("n=%d, K=%d continuous traits, miss_frac=%.2f, n_reps=%d",
          CONFIG$n_species, CONFIG$n_traits, CONFIG$miss_frac, CONFIG$n_reps),
  sprintf("Total wall: %.1f min",
          (proc.time()[["elapsed"]] - t_start) / 60),
  "",
  "## Means per (method, model, metric), averaged over traits and reps",
  "",
  "```",
  capture.output(print(summary_tbl, row.names = FALSE)),
  "```"
)
writeLines(md, out_md)
cat("\n=== DONE === ", out_rds, " ", out_md, "\n")
