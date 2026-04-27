#!/usr/bin/env Rscript
# script/bench_ou_regime.R
#
# OU and regime-shift DGP test -- where pigauto's GNN should beat
# phylolm-lambda BLUP because the analytical baselines assume BM.
#
# DGPs:
#   - OU:           rTraitCont(model="OU", alpha=2, theta=0)
#                   strong stabilising selection -- distantly-related
#                   species converge near theta, breaking BM-MVN structure.
#   - regime_shift: BM trait + clade-specific intercept shifts (one shift
#                   on each of two major root-descending clades).
#                   Bimodal trait distribution that BM (single Gaussian)
#                   cannot represent.
#
# Plus a covariate signal (linear or nonlinear) on top.
#
# Methods (5):
#   1. column_mean              -- loss floor
#   2. lm                       -- linear OLS, no phylo
#   3. phylolm_lambda_blup      -- linear + BM-BLUP (linear smart, but assumes BM)
#   4. pigauto_sfT              -- pigauto with safety_floor TRUE
#   5. pigauto_sfF              -- pigauto with safety_floor FALSE
#
# Decisive prediction: pigauto_sfF beats phylolm_lambda by >= 10% on OU
# cells and >= 15% on regime_shift cells.
#
# Tier:
#   smoke:  16 cells, n=300, 200 epochs.  ~2-3 hr wall.
#   medium: 64 cells, ~10 hr.

suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH",
                          unset = "/Users/z3437171/Dropbox/Github Local/pigauto")
  devtools::load_all(pkg_path, quiet = TRUE)
})
options(warn = -1L)

TIER  <- toupper(Sys.getenv("PIGAUTO_TIER", "smoke"))
HERE  <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(HERE, "script",
                      sprintf("bench_ou_regime_%s.rds", tolower(TIER)))
out_log <- file.path(HERE, "script",
                      sprintf("bench_ou_regime_%s.log", tolower(TIER)))

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) cat(sprintf("[%6.0fs] ", proc.time()[["elapsed"]] - t0),
                                ..., "\n", sep = "")

if (TIER == "SMOKE") {
  N_SPECIES <- 300L
  SCENARIOS <- c("OU", "regime_shift")
  BETA_STR  <- c(0.0, 0.5, 1.0)
  N_REPS    <- 3L
  EPOCHS    <- 200L
} else if (TIER == "MEDIUM") {
  N_SPECIES <- c(300L, 500L)
  SCENARIOS <- c("OU", "regime_shift")
  BETA_STR  <- c(0.0, 0.3, 0.6, 1.0)
  N_REPS    <- 4L
  EPOCHS    <- 300L
}

EPOCHS    <- as.integer(Sys.getenv("PIGAUTO_BENCH_EPOCHS", EPOCHS))
MISS_FRAC <- 0.30

log_line("=========================================================")
log_line("Tier: ", TIER, "   (OU + REGIME-SHIFT DGP)")
log_line("Sweep: n=", paste(N_SPECIES, collapse=","),
          " | scenarios=", paste(SCENARIOS, collapse=","),
          " | beta=", paste(BETA_STR, collapse=","),
          " | reps=", N_REPS,
          " | epochs=", EPOCHS)
log_line("Output: ", out_rds)

# --- DGP wrapper -----------------------------------------------------------
sim_ou_regime_dgp <- function(n_species, scenario, beta_strength, seed) {
  set.seed(seed)
  tree <- ape::rphylo(n = n_species, birth = 0.8, death = 0.4, T0 = 100)
  tree$tip.label <- paste0("sp", seq_len(n_species))

  # Generate the y-trait via simulate_non_bm()
  y_df <- simulate_non_bm(tree, n_traits = 1L,
                            scenario = scenario,
                            alpha = if (scenario == "OU") 2.0 else 0.5,
                            theta = 0,
                            sigma = 1.0,
                            shift_magnitude = 2.0,
                            seed = seed)
  y_full <- as.numeric(y_df[tree$tip.label, 1])

  # Two predictors: x1 BM, x2 noise (no signal)
  x1 <- ape::rTraitCont(tree, model = "BM", sigma = 1, root.value = 0)
  x1 <- as.numeric(x1[tree$tip.label])
  x2 <- rnorm(n_species)

  # Inject covariate signal (linear) at strength beta
  if (beta_strength > 0) {
    contrib <- 0.7 * x1 + 0.3 * x2
    contrib <- as.numeric(scale(contrib))
    y_full <- y_full + beta_strength * contrib
  }

  # Apply MCAR mask
  n <- length(y_full)
  mask <- logical(n)
  mask[sample.int(n, floor(n * MISS_FRAC))] <- TRUE
  y_obs <- y_full
  y_obs[mask] <- NA

  df_complete <- data.frame(y = y_full, x1 = x1, x2 = x2,
                              row.names = tree$tip.label)
  df_observed <- data.frame(y = y_obs, x1 = x1, x2 = x2,
                              row.names = tree$tip.label)
  list(tree = tree, df_complete = df_complete, df_observed = df_observed,
        mask = mask, predictor_names = c("x1", "x2"),
        scenario = scenario, beta_strength = beta_strength)
}

# --- methods (single-obs only here) ----------------------------------------
score <- function(pred, truth, mask) {
  ok <- mask & is.finite(truth) & is.finite(pred)
  if (sum(ok) < 5L) return(c(rmse = NA_real_, r = NA_real_))
  p <- pred[ok]; t <- truth[ok]
  rmse <- sqrt(mean((p - t)^2))
  r <- if (stats::sd(p) > 1e-10) suppressWarnings(stats::cor(p, t)) else NA_real_
  c(rmse = rmse, r = r)
}

method_column_mean <- function(d) {
  cm <- mean(d$df_observed$y, na.rm = TRUE)
  rep(cm, nrow(d$df_observed))
}

method_lm <- function(d) {
  fit_df <- d$df_observed[!is.na(d$df_observed$y), c("y", "x1", "x2")]
  fit <- tryCatch(lm(y ~ x1 + x2, data = fit_df), error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, nrow(d$df_observed)))
  predict(fit, newdata = d$df_observed)
}

method_phylolm_blup <- function(d) {
  ddf <- d$df_observed
  sp_df <- data.frame(species = rownames(ddf), y = ddf$y, x1 = ddf$x1, x2 = ddf$x2,
                       stringsAsFactors = FALSE)
  sp_df_obs <- sp_df[!is.na(sp_df$y), , drop = FALSE]
  if (nrow(sp_df_obs) < 8L) return(rep(NA_real_, nrow(ddf)))
  rownames(sp_df_obs) <- sp_df_obs$species
  tree_obs <- ape::keep.tip(d$tree, sp_df_obs$species)
  fit <- tryCatch(
    phylolm::phylolm(y ~ x1 + x2, data = sp_df_obs, phy = tree_obs, model = "lambda"),
    error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, nrow(ddf)))
  beta_hat <- coef(fit)
  Xmat <- model.matrix(~ x1 + x2, sp_df)
  fixed <- as.numeric(Xmat %*% beta_hat)

  R <- ape::vcv(d$tree); R <- stats::cov2cor(R); R <- R[sp_df$species, sp_df$species]
  obs_idx <- which(!is.na(sp_df$y))
  miss_idx <- which(is.na(sp_df$y))
  if (length(miss_idx) == 0L) return(sp_df$y[match(rownames(ddf), sp_df$species)])
  e_obs <- sp_df$y[obs_idx] - fixed[obs_idx]
  R_oo <- R[obs_idx, obs_idx, drop = FALSE]
  R_mo <- R[miss_idx, obs_idx, drop = FALSE]
  blup <- as.numeric(R_mo %*% solve(R_oo + diag(1e-6, nrow(R_oo)), e_obs))
  sp_pred <- numeric(nrow(sp_df))
  sp_pred[obs_idx] <- sp_df$y[obs_idx]
  sp_pred[miss_idx] <- fixed[miss_idx] + blup
  names(sp_pred) <- sp_df$species
  sp_pred[rownames(ddf)]
}

fit_pigauto_sf <- function(d, safety_floor, seed) {
  ddf <- d$df_observed
  cov_df <- data.frame(x1 = ddf$x1, x2 = ddf$x2)
  res <- pigauto::impute(traits = ddf[, "y", drop = FALSE],
                          tree = d$tree, covariates = cov_df,
                          missing_frac = 0.0,
                          epochs = EPOCHS, verbose = FALSE, seed = seed,
                          safety_floor = safety_floor)
  res$completed[rownames(ddf), "y"]
}

method_pigauto_sfT <- function(d) fit_pigauto_sf(d, TRUE,  sample.int(1e7L,1L))
method_pigauto_sfF <- function(d) fit_pigauto_sf(d, FALSE, sample.int(1e7L,1L))

# --- main loop -------------------------------------------------------------
cells <- expand.grid(
  n_species = N_SPECIES,
  scenario  = SCENARIOS,
  beta_str  = BETA_STR,
  rep_id    = seq_len(N_REPS),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_cells <- nrow(cells)
log_line("Total cells: ", total_cells)

method_funs <- list(
  column_mean         = method_column_mean,
  lm                  = method_lm,
  phylolm_lambda_blup = method_phylolm_blup,
  pigauto_sfT         = method_pigauto_sfT,
  pigauto_sfF         = method_pigauto_sfF
)

results <- list()
for (i in seq_len(total_cells)) {
  row <- cells[i, ]
  seed <- 13000L + i
  log_line(sprintf("Cell %d/%d: scenario=%s, n=%d, beta=%.1f, rep=%d (seed=%d)",
                    i, total_cells, row$scenario, row$n_species,
                    row$beta_str, row$rep_id, seed))
  d <- tryCatch(
    sim_ou_regime_dgp(row$n_species, row$scenario, row$beta_str, seed),
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
      n_species = row$n_species, scenario = row$scenario,
      beta_strength = row$beta_str, rep = row$rep_id,
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
