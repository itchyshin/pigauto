#!/usr/bin/env Rscript
# script/bench_transformer_ablation.R
#
# Architecture ablation: graph transformer (current default) vs legacy
# attention-GNN vs no-attention GNN, on linear/nonlinear/interactive DGPs.
#
# Pigauto's `ResidualPhyloDAE` exposes two architecture toggles:
#   use_transformer_blocks (Phase 9 + B2 default = TRUE):
#     when TRUE  -> n_gnn_layers GraphTransformerBlock instances
#                   (multi-head attention + FFN + 2 residual skips +
#                   per-head Gaussian-bandwidth phylo bias on attention).
#     when FALSE -> legacy 2-layer single-head attention stack.
#   use_attention (default = TRUE):
#     when TRUE  -> attention with learnable log-adjacency bias
#                   (relevant for both transformer and legacy paths).
#     when FALSE -> plain message passing on the symmetric-normalised
#                   adjacency (no attention).
#
# This bench answers the question raised in FUTURE_WORK.md item 1:
# "is the rate-aware Gaussian-bandwidth attention bias doing real work,
#  or is plain self-attention enough?"  Concretely:
#
#   - on linear DGPs: all 3 configs should converge to "match phylolm-lambda"
#     (gate-closed safety) -> they should be ~tied.
#   - on nonlinear DGPs: the transformer's per-head bandwidths should
#     allow it to capture the function across different evolutionary
#     scales, so it should beat the legacy attention-GNN modestly and
#     beat the no-attention GNN clearly.
#   - on interactive DGPs: hard to predict; this is where multi-head
#     attention might or might not pay off.
#
# CONFIGS (3 pigauto + 2 baselines = 5 methods per cell)
#   1. column_mean              -- loss floor
#   2. phylolm_lambda_blup      -- linear smart baseline
#   3. pigauto_transformer      -- use_transformer_blocks=TRUE,  use_attention=TRUE
#   4. pigauto_legacy_attn      -- use_transformer_blocks=FALSE, use_attention=TRUE
#   5. pigauto_legacy_noattn    -- use_transformer_blocks=FALSE, use_attention=FALSE
#
# Tier (env var PIGAUTO_TIER): smoke / medium.
#
# Run:
#   PIGAUTO_TIER=smoke   PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_transformer_ablation.R
#   PIGAUTO_TIER=medium  PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_transformer_ablation.R

suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH",
                          unset = "/Users/z3437171/Dropbox/Github Local/pigauto")
  devtools::load_all(pkg_path, quiet = TRUE)
  source(file.path(pkg_path, "script", "sim_bace_dgp.R"))
})
options(warn = -1L)

TIER  <- toupper(Sys.getenv("PIGAUTO_TIER", "smoke"))
HERE  <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(HERE, "script",
                      sprintf("bench_transformer_ablation_%s.rds", tolower(TIER)))
out_log <- file.path(HERE, "script",
                      sprintf("bench_transformer_ablation_%s.log", tolower(TIER)))

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) cat(sprintf("[%6.0fs] ", proc.time()[["elapsed"]] - t0),
                                ..., "\n", sep = "")

if (TIER == "SMOKE") {
  # 24 cells (12 grid * 2 reps) * 5 methods = 120 fits, ~75-90 min wall
  N_SPECIES <- 200L
  MULTI_OBS <- c(1L, 4L)
  PHYLO_SIG <- 0.4
  BETA_STR  <- c(0.3, 0.6)
  F_TYPES   <- c("linear", "nonlinear", "interactive")
  N_REPS    <- 2L
} else if (TIER == "MEDIUM") {
  # 144 cells (72 grid * 2 reps) * 5 methods = 720 fits, ~9-12 hr wall
  N_SPECIES <- c(200L, 500L)
  MULTI_OBS <- c(1L, 4L)
  PHYLO_SIG <- c(0.2, 0.6)
  BETA_STR  <- c(0.3, 0.6)
  F_TYPES   <- c("linear", "nonlinear", "interactive")
  N_REPS    <- 3L
} else {
  stop("PIGAUTO_TIER must be one of: smoke, medium")
}

EPOCHS    <- as.integer(Sys.getenv("PIGAUTO_BENCH_EPOCHS", "150"))
MISS_FRAC <- 0.30

# Pull shared helpers (sim_nonlinear_dgp, score, method_*).
source(file.path(HERE, "script", "sim_nonlinear_helpers.R"))

log_line("=========================================================")
log_line("Tier: ", TIER, "   (TRANSFORMER ABLATION)")
log_line("Sweep: n_species=", paste(N_SPECIES, collapse=","),
          " | multi_obs=", paste(MULTI_OBS, collapse=","),
          " | phylo_sig=", paste(PHYLO_SIG, collapse=","),
          " | beta=", paste(BETA_STR, collapse=","),
          " | f_types=", paste(F_TYPES, collapse=","),
          " | reps=", N_REPS)
log_line("Output: ", out_rds)

# ---------------- score / baselines reused from the nonlinear bench --------
# (sourcing that file at the top exposes score(), method_column_mean(),
#  method_phylolm_blup(), and sim_nonlinear_dgp() in this scope.)

# pigauto with config knobs
fit_pigauto_with <- function(d, use_transformer_blocks, use_attention, seed) {
  ddf <- d$df_observed
  pred_n <- d$predictor_names
  cov_df <- as.data.frame(lapply(pred_n, function(p) {
    v <- ddf[[p]]; if (is.factor(v)) as.numeric(v) - 1 else as.numeric(v)
  }))
  names(cov_df) <- pred_n
  multi_obs <- d$meta$multi_obs_ratio > 1L
  res <- if (multi_obs) {
    pigauto::impute(traits = ddf[, c("species", "y"), drop = FALSE],
                     tree = d$tree, species_col = "species",
                     covariates = cov_df, missing_frac = 0.0,
                     epochs = EPOCHS, verbose = FALSE, seed = seed,
                     use_transformer_blocks = use_transformer_blocks,
                     use_attention = use_attention)
  } else {
    df_single <- ddf[, c("y", pred_n), drop = FALSE]
    rownames(df_single) <- ddf$species
    pigauto::impute(traits = df_single[, "y", drop = FALSE],
                     tree = d$tree,
                     covariates = cov_df, missing_frac = 0.0,
                     epochs = EPOCHS, verbose = FALSE, seed = seed,
                     use_transformer_blocks = use_transformer_blocks,
                     use_attention = use_attention)
  }
  if (multi_obs) res$completed$y else
    res$completed[as.character(ddf$species), "y"]
}

method_pigauto_transformer   <- function(d) fit_pigauto_with(d, TRUE, TRUE, sample.int(1e7L,1L))
method_pigauto_legacy_attn   <- function(d) fit_pigauto_with(d, FALSE, TRUE, sample.int(1e7L,1L))
method_pigauto_legacy_noattn <- function(d) fit_pigauto_with(d, FALSE, FALSE, sample.int(1e7L,1L))

# ---------------- main loop -----------------------------------------------

cells <- expand.grid(
  n_species  = N_SPECIES,
  multi_obs  = MULTI_OBS,
  phylo_sig  = PHYLO_SIG,
  beta_str   = BETA_STR,
  f_type     = F_TYPES,
  rep_id     = seq_len(N_REPS),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_cells <- nrow(cells)
log_line("Total cells: ", total_cells)

method_funs <- list(
  column_mean              = method_column_mean,
  phylolm_lambda_blup      = method_phylolm_blup,
  pigauto_transformer      = method_pigauto_transformer,
  pigauto_legacy_attn      = method_pigauto_legacy_attn,
  pigauto_legacy_noattn    = method_pigauto_legacy_noattn
)

results <- list()
for (i in seq_len(total_cells)) {
  row <- cells[i, ]
  seed <- 8000L + i
  log_line(sprintf("Cell %d/%d: f=%s, n=%d, multi=%d, beta=%.1f, rep=%d (seed=%d)",
                    i, total_cells, row$f_type, row$n_species, row$multi_obs,
                    row$beta_str, row$rep_id, seed))
  d <- tryCatch(
    sim_nonlinear_dgp(n_species = row$n_species,
                        multi_obs_ratio = row$multi_obs,
                        phylo_signal = row$phylo_sig,
                        beta_strength = row$beta_str,
                        f_type = row$f_type, seed = seed),
    error = function(e) { log_line("  DGP ERROR: ", e$message); NULL })
  if (is.null(d)) next
  truth <- d$df_complete$y

  for (m_name in names(method_funs)) {
    pred <- tryCatch(method_funs[[m_name]](d),
                      error = function(e) {
                        log_line("  ", m_name, " ERROR: ", e$message)
                        rep(NA_real_, nrow(d$df_observed))
                      })
    s <- score(pred, truth, d$mask)
    results[[length(results) + 1L]] <- data.frame(
      n_species = row$n_species, multi_obs = row$multi_obs,
      phylo_signal = row$phylo_sig, beta_strength = row$beta_str,
      f_type = row$f_type, rep = row$rep_id,
      method = m_name,
      rmse = s["rmse"], r = s["r"],
      stringsAsFactors = FALSE
    )
  }

  if (i %% 4L == 0L || i == total_cells) {
    df <- do.call(rbind, results)
    saveRDS(df, out_rds)
    log_line(sprintf("  [checkpoint] saved %d rows", nrow(df)))
  }
  invisible(gc(verbose = FALSE))
}

df <- do.call(rbind, results)
saveRDS(df, out_rds)
log_line(sprintf("=== DONE: %d cells, %d rows ===", total_cells, nrow(df)))
log_line(sprintf("  rds: %s", out_rds))
