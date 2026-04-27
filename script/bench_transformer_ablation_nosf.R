#!/usr/bin/env Rscript
# script/bench_transformer_ablation_nosf.R
#
# Architecture ablation with the SAFETY FLOOR DISABLED.
#
# The earlier bench_transformer_ablation.R (smoke tier) showed all 3
# pigauto configs (transformer / legacy_attn / legacy_noattn) collapse to
# within 1% of each other on every cell.  The most likely cause: the
# safety floor uniformly closes the gate when no GNN config can extract
# signal above noise, masking real architectural differences.  This
# script repeats the ablation with safety_floor = FALSE so we can see
# the raw GNN delta's contribution per architecture.
#
# Methods (3 pigauto + 2 baselines):
#   1. column_mean              -- loss floor
#   2. phylolm_lambda_blup      -- linear smart baseline
#   3. pigauto_tfm_nosf         -- use_transformer_blocks=TRUE,  safety_floor=FALSE
#   4. pigauto_attn_nosf        -- use_transformer_blocks=FALSE, use_attention=TRUE,  safety_floor=FALSE
#   5. pigauto_noattn_nosf      -- use_transformer_blocks=FALSE, use_attention=FALSE, safety_floor=FALSE
#
# Tier:
#   smoke:  24 cells, n=500, 300 epochs.  ~3-4 hr wall.
#   medium: 144 cells, ~14-18 hr.
#
# Run:
#   PIGAUTO_TIER=smoke PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_transformer_ablation_nosf.R

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
                      sprintf("bench_transformer_ablation_nosf_%s.rds", tolower(TIER)))
out_log <- file.path(HERE, "script",
                      sprintf("bench_transformer_ablation_nosf_%s.log", tolower(TIER)))

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) cat(sprintf("[%6.0fs] ", proc.time()[["elapsed"]] - t0),
                                ..., "\n", sep = "")

if (TIER == "SMOKE") {
  N_SPECIES <- 500L
  MULTI_OBS <- c(1L, 4L)
  PHYLO_SIG <- 0.4
  BETA_STR  <- c(0.5, 1.0)
  F_TYPES   <- c("linear", "nonlinear", "interactive")
  N_REPS    <- 2L
} else if (TIER == "MEDIUM") {
  N_SPECIES <- c(200L, 500L)
  MULTI_OBS <- c(1L, 4L)
  PHYLO_SIG <- c(0.2, 0.6)
  BETA_STR  <- c(0.5, 1.0)
  F_TYPES   <- c("linear", "nonlinear", "interactive")
  N_REPS    <- 3L
} else {
  stop("PIGAUTO_TIER must be one of: smoke, medium")
}

EPOCHS    <- as.integer(Sys.getenv("PIGAUTO_BENCH_EPOCHS", "300"))
MISS_FRAC <- 0.30

source(file.path(HERE, "script", "sim_nonlinear_helpers.R"))

log_line("=========================================================")
log_line("Tier: ", TIER, "   (TRANSFORMER ABLATION, NO SAFETY FLOOR)")
log_line("Sweep: n=", paste(N_SPECIES, collapse=","),
          " | multi=", paste(MULTI_OBS, collapse=","),
          " | beta=", paste(BETA_STR, collapse=","),
          " | f_types=", paste(F_TYPES, collapse=","),
          " | reps=", N_REPS,
          " | epochs=", EPOCHS,
          " | safety_floor=FALSE")
log_line("Output: ", out_rds)

fit_pigauto_arch <- function(d, use_transformer_blocks, use_attention, seed) {
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
                     safety_floor = FALSE,
                     use_transformer_blocks = use_transformer_blocks,
                     use_attention = use_attention)
  } else {
    df_single <- ddf[, c("y", pred_n), drop = FALSE]
    rownames(df_single) <- ddf$species
    pigauto::impute(traits = df_single[, "y", drop = FALSE],
                     tree = d$tree,
                     covariates = cov_df, missing_frac = 0.0,
                     epochs = EPOCHS, verbose = FALSE, seed = seed,
                     safety_floor = FALSE,
                     use_transformer_blocks = use_transformer_blocks,
                     use_attention = use_attention)
  }
  if (multi_obs) res$completed$y else
    res$completed[as.character(ddf$species), "y"]
}

method_pigauto_tfm_nosf    <- function(d) fit_pigauto_arch(d, TRUE,  TRUE,  sample.int(1e7L,1L))
method_pigauto_attn_nosf   <- function(d) fit_pigauto_arch(d, FALSE, TRUE,  sample.int(1e7L,1L))
method_pigauto_noattn_nosf <- function(d) fit_pigauto_arch(d, FALSE, FALSE, sample.int(1e7L,1L))

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
  column_mean         = method_column_mean,
  phylolm_lambda_blup = method_phylolm_blup,
  pigauto_tfm_nosf    = method_pigauto_tfm_nosf,
  pigauto_attn_nosf   = method_pigauto_attn_nosf,
  pigauto_noattn_nosf = method_pigauto_noattn_nosf
)

results <- list()
for (i in seq_len(total_cells)) {
  row <- cells[i, ]
  seed <- 11000L + i
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
  if (i %% 2L == 0L || i == total_cells) {
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
