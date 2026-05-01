#!/usr/bin/env Rscript
# script/bench_sim_bace_nonlinear.R
#
# GNN-ability sweep: BACE-simulated DGP + post-hoc nonlinear/interactive
# response injection x pigauto (Fix A-H + row-fix) vs analytical baselines.
#
# WHY THIS BENCH EXISTS
#
# bench_sim_bace_pigauto.R uses purely linear DGPs (y = X*beta + species + phylo
# + obs noise).  On linear DGPs, phylolm-lambda is mathematically optimal and
# pigauto's GNN should close its gate to MATCH (not beat) phylolm.  The
# safety_floor and gate-closure are pigauto's value there: robustness, not
# lift.
#
# To demonstrate the GNN's ABILITY, we need DGPs where the BM/GLS baseline
# mathematically can't capture the signal.  This script injects nonlinear /
# interactive structure on top of BACE-sim's linear gaussian output:
#
#   f_type = "linear"      -> y = (BACE-sim linear y) + 0  (control)
#   f_type = "nonlinear"   -> y = (BACE-sim base y) + beta * scale(sin(2*x1) * exp(0.3*x2))
#   f_type = "interactive" -> y = (BACE-sim base y) + beta * scale(x1*x2 + 0.5*x1^2)
#
# In linear cells, pigauto should match phylolm-lambda (gate-closed parity).
# In nonlinear / interactive cells, pigauto should clearly BEAT phylolm-lambda
# AND beat a nonlinear OLS baseline (lm with poly + interactions, no phylo).
#
# METHODS (5)
#   1. column_mean             -- loss floor
#   2. lm                      -- lm(y ~ x1 + x2 + x3) (linear, no phylo)
#   3. lm_nonlinear            -- lm(y ~ poly(x1, 2) + poly(x2, 2) + x3 + x1:x2)
#                                  (smart practitioner with poly + interaction features, no phylo)
#   4. phylolm_lambda_blup     -- GLS regression + lambda-fitted BM BLUP
#                                  (linear features + phylo, the "linear smart" baseline)
#   5. pigauto_cov_sfT         -- pigauto::impute() with covariates, safety_floor=TRUE
#                                  (nonlinear via GNN + phylo, what we are testing)
#
# DECISIVE PREDICTIONS
#
#   linear cells:        pigauto matches phylolm_lambda within 10 %; both clearly
#                         beat lm and lm_nonlinear.
#   nonlinear cells:     pigauto BEATS phylolm_lambda by >= 10 %; pigauto matches
#                         or beats lm_nonlinear; lm clearly worst.
#   interactive cells:   pigauto matches lm_nonlinear within 10 % (since that
#                         baseline has the bilinear feature) but BEATS
#                         phylolm_lambda by >= 10 %.
#
# Tier (env var PIGAUTO_TIER): smoke / medium / full.
#
# Run:
#   PIGAUTO_TIER=smoke   PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_sim_bace_nonlinear.R
#   PIGAUTO_TIER=medium  PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_sim_bace_nonlinear.R

suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH",
                          unset = "/Users/z3437171/Dropbox/Github Local/pigauto")
  devtools::load_all(pkg_path, quiet = TRUE)
  source(file.path(pkg_path, "script", "sim_bace_dgp.R"))
})
options(warn = -1L)

# MISS_FRAC and EPOCHS are referenced by sim_nonlinear_helpers.R; define
# them BEFORE sourcing so the helper picks up our values rather than the
# defaults.

TIER  <- toupper(Sys.getenv("PIGAUTO_TIER", "smoke"))
HERE  <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(HERE, "script",
                      sprintf("bench_sim_bace_nonlinear_%s.rds", tolower(TIER)))
out_log <- file.path(HERE, "script",
                      sprintf("bench_sim_bace_nonlinear_%s.log", tolower(TIER)))

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) cat(sprintf("[%6.0fs] ", proc.time()[["elapsed"]] - t0),
                                ..., "\n", sep = "")

# ---------------- Sweep grid (per tier) -----------------------------------

if (TIER == "SMOKE") {
  # 24 cells: tests basic feasibility, ~30-40 min wall.
  N_SPECIES <- c(100L)
  MULTI_OBS <- c(1L, 4L)
  PHYLO_SIG <- c(0.4)
  BETA_STR  <- c(0.0, 0.5)
  F_TYPES   <- c("linear", "nonlinear", "interactive")
  N_REPS    <- 2L
} else if (TIER == "MEDIUM") {
  # 216 cells: ~12 hr wall.
  N_SPECIES <- c(100L, 500L)
  MULTI_OBS <- c(1L, 4L)
  PHYLO_SIG <- c(0.2, 0.6)
  BETA_STR  <- c(0.0, 0.3, 0.6)
  F_TYPES   <- c("linear", "nonlinear", "interactive")
  N_REPS    <- 3L
} else if (TIER == "FULL") {
  # 1296 cells: not recommended without HPC.
  N_SPECIES <- c(100L, 500L, 1500L)
  MULTI_OBS <- c(1L, 4L)
  PHYLO_SIG <- c(0.2, 0.5, 0.8)
  BETA_STR  <- c(0.0, 0.3, 0.6)
  F_TYPES   <- c("linear", "nonlinear", "interactive")
  N_REPS    <- 3L
} else {
  stop("PIGAUTO_TIER must be one of: smoke, medium, full")
}

MISS_FRAC <- 0.30
EPOCHS    <- as.integer(Sys.getenv("PIGAUTO_BENCH_EPOCHS", "150"))

# Pull shared helpers (sim_nonlinear_dgp, score, method_*).
source(file.path(HERE, "script", "sim_nonlinear_helpers.R"))

log_line("=========================================================")
log_line("Tier: ", TIER, "   (NONLINEAR DGP)")
log_line("Sweep: n_species=", paste(N_SPECIES, collapse=","),
          " | multi_obs=", paste(MULTI_OBS, collapse=","),
          " | phylo_signal=", paste(PHYLO_SIG, collapse=","),
          " | beta=", paste(BETA_STR, collapse=","),
          " | f_types=", paste(F_TYPES, collapse=","),
          " | reps=", N_REPS)
log_line("Output: ", out_rds)

method_pigauto_cov_sfT <- function(d) {
  ddf <- d$df_observed
  pred_n <- d$predictor_names

  # Encode binary factor predictor as numeric
  cov_df <- as.data.frame(lapply(pred_n, function(p) {
    v <- ddf[[p]]; if (is.factor(v)) as.numeric(v) - 1 else as.numeric(v)
  }))
  names(cov_df) <- pred_n

  multi_obs <- d$meta$multi_obs_ratio > 1L

  res <- if (multi_obs) {
    pigauto::impute(traits = ddf[, c("species", "y"), drop = FALSE],
                     tree = d$tree, species_col = "species",
                     covariates = cov_df, missing_frac = 0.0,
                     epochs = EPOCHS, verbose = FALSE,
                     seed = sample.int(1e7L, 1L))
  } else {
    # single-obs: rownames must be species
    df_single <- ddf[, c("y", pred_n), drop = FALSE]
    rownames(df_single) <- ddf$species
    pigauto::impute(traits = df_single[, "y", drop = FALSE],
                     tree = d$tree,
                     covariates = cov_df, missing_frac = 0.0,
                     epochs = EPOCHS, verbose = FALSE,
                     seed = sample.int(1e7L, 1L))
  }
  if (multi_obs) res$completed$y else {
    # single-obs result has rownames = tree$tip.label.  Map to ddf order.
    res$completed[as.character(ddf$species), "y"]
  }
}

# ---------------- Main loop -----------------------------------------------

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

results <- list()
for (i in seq_len(total_cells)) {
  row <- cells[i, ]
  seed <- 7000L + i
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

  out <- list()
  method_funs <- list(
    column_mean         = method_column_mean,
    lm                  = method_lm,
    lm_nonlinear        = method_lm_nonlinear,
    phylolm_lambda_blup = method_phylolm_blup,
    pigauto_cov_sfT     = method_pigauto_cov_sfT
  )
  for (m_name in names(method_funs)) {
    fn <- method_funs[[m_name]]
    pred <- tryCatch(fn(d),
                      error = function(e) {
                        log_line("  ", m_name, " ERROR: ", e$message)
                        rep(NA_real_, nrow(d$df_observed))
                      })
    s <- score(pred, truth, d$mask)
    out[[m_name]] <- list(rmse = s["rmse"], r = s["r"])
  }

  for (m in names(out)) {
    results[[length(results) + 1L]] <- data.frame(
      n_species = row$n_species, multi_obs = row$multi_obs,
      phylo_signal = row$phylo_sig, beta_strength = row$beta_str,
      f_type = row$f_type, rep = row$rep_id,
      method = m,
      rmse = out[[m]]$rmse, r = out[[m]]$r,
      stringsAsFactors = FALSE
    )
  }

  # Incremental save
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
