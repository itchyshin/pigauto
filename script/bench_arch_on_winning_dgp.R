#!/usr/bin/env Rscript
# script/bench_arch_on_winning_dgp.R
#
# DECISIVE architecture ablation on the DGP where pigauto wins.
#
# We've shown pigauto beats lm_nonlinear / phylolm-λ for λ ≥ 0.15 on
# the multi-obs nonlinear DGP.  But we have NOT yet decomposed that
# win: is it the graph transformer's multi-head attention + FFN, or
# is it the BM/GLS analytical baseline + obs-level aggregation, with
# the GNN delta gated to ~0?
#
# Earlier ablations (bench_transformer_ablation, bench_transformer_ablation_nosf)
# tested architecture toggles on cells where pigauto did not clearly win
# anyway -- the within-1.21% spread there is uninformative because no
# architecture was earning.  This script repeats the toggle on the
# WINNING cells, where the spread (if any) actually means something.
#
# DGP: multi-obs nonlinear / interactive at lambda=0.2/0.3/0.5,
#      n=500, ncov=10, beta=1.0, 3 reps.
#
# Three configs (all with safety_floor=TRUE since that's the default
# users will run):
#   1. pigauto_transformer   : use_transformer_blocks=TRUE,  use_attention=TRUE
#   2. pigauto_legacy_attn   : use_transformer_blocks=FALSE, use_attention=TRUE  (GAT-style)
#   3. pigauto_no_attn       : use_transformer_blocks=FALSE, use_attention=FALSE (plain GCN)
#
# Total: 3 lambdas × 2 f_types × 3 reps × 3 configs = 54 fits, ~60 min wall.
#
# Decision rule:
#   - If pigauto_transformer / pigauto_no_attn ≤ 0.95 (≥5% better) on
#     nonlinear cells: attention / FFN earns its keep.  Default to
#     transformer + claim the architecture as a contribution.
#   - If ratio in [0.98, 1.02]: architectures equivalent.  Keep
#     transformer as default for safety/forward-compatibility but
#     drop the architectural-novelty claim.

suppressPackageStartupMessages({
  library(ape); library(phylolm)
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH",
                          unset = "/Users/z3437171/Dropbox/Github Local/pigauto")
  devtools::load_all(pkg_path, quiet = TRUE)
})
options(warn = -1L)

HERE <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(HERE, "script", "bench_arch_on_winning_dgp.rds")
out_log <- file.path(HERE, "script", "bench_arch_on_winning_dgp.log")

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) cat(sprintf("[%6.0fs] ", proc.time()[["elapsed"]] - t0), ..., "\n", sep = "")

N_SPECIES <- 500L; OBS_PER <- 5L; N_COVS <- 10L
LAMBDAS <- c(0.2, 0.3, 0.5)
F_TYPES <- c("nonlinear", "interactive")
BETA    <- 1.0
N_REPS  <- 3L
EPOCHS  <- 250L
SP_MISSING_FRAC  <- 0.5
WITHIN_MISS_FRAC <- 0.2
SIGMA_RES <- 0.5

log_line("ARCHITECTURE ABLATION ON WINNING DGP")
log_line("Sweep: n=", N_SPECIES, " ncov=", N_COVS,
          " lambda=", paste(LAMBDAS, collapse=","),
          " f=", paste(F_TYPES, collapse=","),
          " beta=", BETA, " reps=", N_REPS, " epochs=", EPOCHS)

# ---------- DGP (same as bench_lambda_sweep / bench_phase1_extended) ----
sim_dgp <- function(lambda, f_type, seed) {
  set.seed(seed)
  tree <- ape::rtree(N_SPECIES); sp <- tree$tip.label
  phylo_vals <- as.numeric(ape::rTraitCont(tree, model = "BM",
                                              sigma = sqrt(lambda), root.value = 0)[sp])
  names(phylo_vals) <- sp
  n_per_sp <- pmin(pmax(rpois(N_SPECIES, OBS_PER), 1L), 15L)
  species_vec <- rep(sp, n_per_sp); n_total <- length(species_vec)
  Z <- matrix(rnorm(n_total * N_COVS), nrow = n_total, ncol = N_COVS)
  colnames(Z) <- paste0("z", seq_len(N_COVS))
  contrib <- switch(f_type,
    nonlinear = sin(Z[,1])*exp(0.3*Z[,2]) + 0.5*tanh(Z[,3]) +
                  0.3*sin(2*Z[,4]) - 0.2*cos(Z[,5]) +
                  0.4*tanh(Z[,6]+Z[,7]) - 0.3*sin(Z[,8]),
    interactive = Z[,1]*Z[,2] + 0.5*Z[,3]^2 + 0.3*Z[,1]*Z[,4] -
                    0.2*Z[,2]*Z[,5] + 0.25*Z[,6]*Z[,7] - 0.15*Z[,8]*Z[,9]^2
  )
  contrib <- as.numeric(scale(contrib))
  y <- phylo_vals[species_vec] + BETA*contrib + rnorm(n_total, 0, SIGMA_RES)
  miss_sp <- sample(sp, floor(N_SPECIES * SP_MISSING_FRAC))
  y_obs <- y; y_obs[species_vec %in% miss_sp] <- NA
  obs_idx <- which(!is.na(y_obs))
  n_within <- floor(WITHIN_MISS_FRAC * length(obs_idx))
  if (n_within > 0L) y_obs[sample(obs_idx, n_within)] <- NA
  mask <- is.na(y_obs)
  df_obs <- data.frame(species = species_vec, y = y_obs, Z, stringsAsFactors = FALSE)
  df_full <- data.frame(species = species_vec, y = y, Z, stringsAsFactors = FALSE)
  list(tree = tree, df_complete = df_full, df_observed = df_obs,
        mask = mask, predictor_names = colnames(Z))
}

score <- function(pred, truth, mask) {
  ok <- mask & is.finite(truth) & is.finite(pred)
  if (sum(ok) < 5L) return(NA_real_)
  sqrt(mean((pred[ok] - truth[ok])^2))
}

fit_pigauto_arch <- function(d, use_transformer_blocks, use_attention, seed) {
  ddf <- d$df_observed; pred_n <- d$predictor_names
  cov_df <- as.data.frame(lapply(pred_n, function(p) as.numeric(ddf[[p]])))
  names(cov_df) <- pred_n
  res <- pigauto::impute(traits = ddf[, c("species", "y"), drop = FALSE],
                          tree = d$tree, species_col = "species", covariates = cov_df,
                          missing_frac = 0.0, epochs = EPOCHS, verbose = FALSE,
                          seed = seed, safety_floor = TRUE,
                          use_transformer_blocks = use_transformer_blocks,
                          use_attention = use_attention)
  res$completed$y
}

m_transformer  <- function(d, s) fit_pigauto_arch(d, TRUE,  TRUE,  s)
m_legacy_attn  <- function(d, s) fit_pigauto_arch(d, FALSE, TRUE,  s)
m_no_attn      <- function(d, s) fit_pigauto_arch(d, FALSE, FALSE, s)

cells <- expand.grid(lambda = LAMBDAS, f_type = F_TYPES, rep_id = seq_len(N_REPS),
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
total <- nrow(cells); log_line("Total cells: ", total, " x 3 archs = ", total*3, " fits")

method_funs <- list(
  pigauto_transformer = m_transformer,
  pigauto_legacy_attn = m_legacy_attn,
  pigauto_no_attn     = m_no_attn
)

results <- list()
for (i in seq_len(total)) {
  row <- cells[i, ]; seed <- 31000L + i
  log_line(sprintf("Cell %d/%d: f=%s, lam=%.1f, rep=%d (seed=%d)",
                    i, total, row$f_type, row$lambda, row$rep_id, seed))
  d <- tryCatch(sim_dgp(row$lambda, row$f_type, seed),
                 error = function(e) { log_line("DGP ERR: ", e$message); NULL })
  if (is.null(d)) next
  truth <- d$df_complete$y
  for (m_name in names(method_funs)) {
    pred <- tryCatch(method_funs[[m_name]](d, seed),
                      error = function(e) { log_line("  ", m_name, " ERR: ", e$message); rep(NA_real_, nrow(d$df_observed)) })
    rmse_val <- score(pred, truth, d$mask)
    results[[length(results) + 1L]] <- data.frame(
      lambda = row$lambda, f_type = row$f_type, rep = row$rep_id,
      method = m_name, rmse = rmse_val, stringsAsFactors = FALSE)
  }
  if (i %% 2L == 0L || i == total) {
    df <- do.call(rbind, results); saveRDS(df, out_rds)
    log_line(sprintf("  [checkpoint] %d rows", nrow(df)))
  }
  invisible(gc(verbose = FALSE))
}

df <- do.call(rbind, results); saveRDS(df, out_rds)
log_line(sprintf("=== DONE: %d cells, %d rows ===", total, nrow(df)))
