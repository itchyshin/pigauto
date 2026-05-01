#!/usr/bin/env Rscript
# script/bench_gnn_capacity.R
#
# GNN-CAPACITY TEST -- the bench that should actually demonstrate the
# GNN's value-add over the analytical baselines, by removing the design
# choices that prevented it from doing so in earlier benches:
#
#   1. n_species = 500   (was 100): the GNN needs enough species to learn
#                                    nonlinear function structure across
#                                    phylogenetic neighborhoods.
#   2. beta in {0.5, 1.0}: stronger nonlinear contribution relative to
#                          BACE-sim's noise floor (was beta in {0, 0.5}).
#   3. epochs = 300       (was 150): more training for the GNN to converge
#                                     past the linear-baseline-matching phase.
#   4. ALSO test safety_floor = FALSE (in addition to TRUE): the safety
#                            floor uniformly closed the gate in the previous
#                            ablation, masking architectural differences.
#                            With it off, we see the raw GNN delta's
#                            contribution.
#
# Methods (6 per cell, run 6 fits per cell):
#   1. column_mean              -- loss floor
#   2. lm                       -- linear OLS (no phylo)
#   3. lm_nonlinear             -- poly(2) + bilinear interaction (no phylo)
#   4. phylolm_lambda_blup      -- linear + lambda-BM-BLUP (the linear-smart baseline)
#   5. pigauto_sfT              -- pigauto with safety_floor = TRUE  (the safe default)
#   6. pigauto_sfF              -- pigauto with safety_floor = FALSE (the GNN unleashed)
#
# DECISIVE PREDICTIONS (registered before harvest):
#
#   linear cells:        pigauto_sfT matches phylolm_lambda within 5%; sfF
#                         could overshoot (worse than phylolm, the cost of
#                         removing the safety floor on signal-free DGPs).
#   nonlinear cells, beta=1.0:
#                         pigauto_sfF beats phylolm_lambda by >= 15%.
#                         pigauto_sfT may also beat by 5-15% (gate stays
#                         partially open with strong nonlinear signal).
#                         lm_nonlinear is a meaningful bar -- pigauto needs
#                         to beat or match it to claim "phylo-aware GNN".
#   interactive cells, beta=1.0:
#                         lm_nonlinear should be near-optimal in no-phylo
#                         framing.  Pigauto wins only via phylo+nonlinear
#                         joint use.  Expect pigauto_sfF / lm_nonlinear ratio
#                         in [0.85, 0.95].
#
# Tier (env var PIGAUTO_TIER): smoke / medium.
#   smoke:  24 cells = 2 multi_obs * 2 beta * 3 f_type * 2 reps
#           Each cell: 6 methods, 2 of which are heavy pigauto fits at
#                     n=500, 300 epochs (~5 min each).
#           Total: 24 cells * (~10-12 min/cell) ≈ 4-5 hr wall.
#
#   medium: 72 cells = 2 multi_obs * 2 phylo_sig * 2 beta * 3 f_type * 3 reps
#           ~14-18 hr.
#
# Run:
#   PIGAUTO_TIER=smoke   PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_gnn_capacity.R

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
                      sprintf("bench_gnn_capacity_%s.rds", tolower(TIER)))
out_log <- file.path(HERE, "script",
                      sprintf("bench_gnn_capacity_%s.log", tolower(TIER)))

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
  N_SPECIES <- 500L
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

# Pull shared helpers (sim_nonlinear_dgp, score, method_*).  Includes the
# 2026-04-26 PM phylolm-singleobs fix.
source(file.path(HERE, "script", "sim_nonlinear_helpers.R"))

log_line("=========================================================")
log_line("Tier: ", TIER, "   (GNN CAPACITY)")
log_line("Sweep: n=", paste(N_SPECIES, collapse=","),
          " | multi=", paste(MULTI_OBS, collapse=","),
          " | ps=", paste(PHYLO_SIG, collapse=","),
          " | beta=", paste(BETA_STR, collapse=","),
          " | f_types=", paste(F_TYPES, collapse=","),
          " | reps=", N_REPS,
          " | epochs=", EPOCHS)
log_line("Output: ", out_rds)

# --- pigauto with safety_floor toggle --------------------------------------
fit_pigauto_sf <- function(d, safety_floor, seed) {
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
                     safety_floor = safety_floor)
  } else {
    df_single <- ddf[, c("y", pred_n), drop = FALSE]
    rownames(df_single) <- ddf$species
    pigauto::impute(traits = df_single[, "y", drop = FALSE],
                     tree = d$tree,
                     covariates = cov_df, missing_frac = 0.0,
                     epochs = EPOCHS, verbose = FALSE, seed = seed,
                     safety_floor = safety_floor)
  }
  if (multi_obs) res$completed$y else
    res$completed[as.character(ddf$species), "y"]
}
method_pigauto_sfT <- function(d) fit_pigauto_sf(d, TRUE,  sample.int(1e7L,1L))
method_pigauto_sfF <- function(d) fit_pigauto_sf(d, FALSE, sample.int(1e7L,1L))

# --- main loop -------------------------------------------------------------
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
  lm                  = method_lm,
  lm_nonlinear        = method_lm_nonlinear,
  phylolm_lambda_blup = method_phylolm_blup,
  pigauto_sfT         = method_pigauto_sfT,
  pigauto_sfF         = method_pigauto_sfF
)

results <- list()
for (i in seq_len(total_cells)) {
  row <- cells[i, ]
  seed <- 9000L + i
  log_line(sprintf("Cell %d/%d: f=%s, n=%d, multi=%d, ps=%.1f, beta=%.1f, rep=%d (seed=%d)",
                    i, total_cells, row$f_type, row$n_species, row$multi_obs,
                    row$phylo_sig, row$beta_str, row$rep_id, seed))
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
