# =============================================================================
# script/sim_bench/pigauto_full_sweep.R - Wide parameter-sweep simulation
# =============================================================================
#
# Companion to Dan's BACE simulation (BACE/dev/02_benchmark_simulated_full.R).
# Same metrics, same DGP family, same crossed (signal x mechanism) design;
# adds the axes BACE can't afford to sweep because each MCMC chain costs
# ~30-60 min:
#
#                                 BACE-side          pigauto-side
#     Phylogenetic-signal scenario      4                  4
#     Missingness mechanism             3                  3
#     N_SPECIES                         {500}              {100, 500, 2000}
#     MISS_RATE                         {0.35}             {0.15, 0.35, 0.55}
#     beta_sparsity (DGP signal)        {0.3}              {0.1, 0.3, 0.6}
#     DGP family                        {sim_bace BM}      {BM, OU, EB}
#     N_SIMS per cell                   50                 20 (sweep, not depth)
#     N cells                           12                 12 x 3 x 3 x 3 x 3 = 972
#     Total replicates                  600                ~19,400
#
# pigauto's compute budget (post Hadfield-S^-1 acceleration): ~30 s per
# replicate at n=500, scaling near-linearly to ~3 min at n=2000. Full
# sweep is ~7 days of single-machine compute, sub-day on 24 cores. Real
# audience is "run once before a paper figure freeze, not on every
# commit".
#
# Output schema matches Dan's so the two studies can be merged into a
# single paper figure: same long-format columns + same metric vocabulary.
#
# Extra pigauto-only metrics (because BACE can't compute them):
#   - gate_mean / gate_active: GNN gate utilisation per trait
#   - n_imputations_used: confirms multi_impute n_final
#   - wall_sec: per-replicate wall time, for the speed-vs-rigour plot
#
# Pipeline order:
#   1) For each cell of the design, simulate complete data (BM/OU/EB)
#   2) Inject missingness under one of the three mechanisms
#   3) multi_impute() with pigauto's documented-default config
#   4) Evaluate per-cell accuracy + calibration + gate-utilisation
#   5) Pool to summary tables stratified by every axis
#
# Spec: TBD (write after exp 4/5 settle and we know which solver ships)
# Dan's design: BACE/dev/02_benchmark_simulated_full.R

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library(BACE)        # for sim_bace()
  library(ape)
  library(MASS)
  library(parallel)
})

set.seed(2026L)

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

# ---- Crossed-design axes ----------------------------------------------------
# Mirror Dan's signal + mechanism axes, add three pigauto-only sweep axes.
SCENARIOS    <- c("all_high", "all_moderate", "all_low", "mixed")
MECHANISMS   <- c("phylo_MAR", "trait_MAR", "trait_MNAR")
N_SPECIES_GRID <- c(100L, 500L, 2000L)
MISS_RATE_GRID <- c(0.15, 0.35, 0.55)
BETA_GRID    <- c(0.1, 0.3, 0.6)
DGP_GRID     <- c("BM", "OU", "EB")
N_SIMS       <- 20L            # per cell of the full crossed design
N_CORES      <- 4L

SIGNAL <- list(high = 0.90, moderate = 0.60, low = 0.20)
DEP_STRENGTH <- 1.5            # same as Dan's (matches Penone 2014 ICCs)

# ---- pigauto's documented-default CI config ---------------------------------
# (kept in sync with script/gha/_ci_config.R::PIGAUTO_CI_CONFIG)
PIGAUTO_SIM_CONFIG <- list(
  n_imputations    = 20L,       # matches Dan's n_final
  pool_method      = "median",
  clamp_outliers   = FALSE,
  phylo_signal_gate = TRUE,
  mc_cores         = 1L,
  seed             = 2026L
)

# ---- Output paths -----------------------------------------------------------
RESULTS_DIR <- file.path("dev", "simulation_results_pigauto")
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# =============================================================================
# 2. SIGNAL SCENARIO GENERATOR (verbatim from Dan's)
# =============================================================================

make_phylo_signal <- function(scenario) {
  switch(scenario,
    all_high     = rep(SIGNAL$high,     5),
    all_moderate = rep(SIGNAL$moderate, 5),
    all_low      = rep(SIGNAL$low,      5),
    mixed = {
      # Fixed canonical mixed configuration (one per release, not per
      # replicate; this differs from Dan's design — see simulation-check
      # review note 8). Easier to interpret aggregates.
      c(SIGNAL$high, SIGNAL$moderate, SIGNAL$low, SIGNAL$moderate, SIGNAL$high)
    },
    stop("Unknown scenario: ", scenario)
  )
}

# =============================================================================
# 3. DGP FAMILY (sim_bace BM + alternative-DGP wrappers)
# =============================================================================

# Wrap BACE::sim_bace under each evolutionary model. For OU and EB we
# perturb the simulated continuous traits by transforming the tree's
# branch lengths before generating values, then restoring the original
# tree so missingness mechanisms still work on the canonical phylogeny.
.sim_complete_data <- function(scenario, n_cases, n_species, beta_sparsity,
                                dgp) {
  phylo_signal <- make_phylo_signal(scenario)
  if (dgp == "BM") {
    return(sim_bace(
      response_type   = "gaussian",
      predictor_types = c("binary", "multinomial3", "poisson", "threshold3"),
      var_names       = c("y", "x1", "x2", "x3", "x4"),
      phylo_signal    = phylo_signal,
      n_cases         = n_cases,
      n_species       = n_species,
      beta_sparsity   = beta_sparsity,
      missingness     = c(0, 0, 0, 0, 0)
    ))
  }
  # OU / EB: transform tree branch lengths, simulate on the transformed
  # tree, then swap the canonical tree back into the returned object so
  # missingness mechanisms see the same tree across DGP cells.
  sim_BM <- sim_bace(
    response_type   = "gaussian",
    predictor_types = c("binary", "multinomial3", "poisson", "threshold3"),
    var_names       = c("y", "x1", "x2", "x3", "x4"),
    phylo_signal    = phylo_signal,
    n_cases         = n_cases,
    n_species       = n_species,
    beta_sparsity   = beta_sparsity,
    missingness     = c(0, 0, 0, 0, 0)
  )
  canon_tree <- sim_BM$tree
  alt_tree   <- switch(dgp,
    OU = .transform_tree_OU(canon_tree, alpha = 1.0),
    EB = .transform_tree_EB(canon_tree, r = -0.5)
  )
  # Re-simulate the continuous y on the alt-tree but keep discrete
  # predictors from the BM sim (BACE's threshold model coupling stays
  # consistent). Replace y in the complete_data with the alt-DGP y.
  R_alt <- ape::vcv(alt_tree, corr = TRUE)
  sigma_y <- sqrt(phylo_signal[1])
  z       <- MASS::mvrnorm(1, mu = rep(0, nrow(R_alt)),
                            Sigma = sigma_y^2 * R_alt)
  names(z) <- rownames(R_alt)
  sim_BM$complete_data$y <- z[as.character(sim_BM$complete_data$species)] +
                              sqrt(1 - phylo_signal[1]) * rnorm(n_cases)
  sim_BM
}

.transform_tree_OU <- function(tree, alpha = 1.0) {
  # Standard OU branch-length rescale (Hansen 1997).
  bl <- tree$edge.length
  tree$edge.length <- (1 - exp(-2 * alpha * bl)) / (2 * alpha)
  tree
}

.transform_tree_EB <- function(tree, r = -0.5) {
  # Early-burst (Blomberg & Garland 2003): branch lengths decay
  # exponentially with node depth.
  depths <- ape::node.depth.edgelength(tree)
  parent <- tree$edge[, 1]
  bl     <- tree$edge.length
  tree$edge.length <- bl * exp(r * depths[parent])
  tree
}

# =============================================================================
# 4. MISSINGNESS MECHANISMS (mirror Dan's verbatim)
# =============================================================================

.trait_to_numeric <- function(x) {
  if (is.factor(x)) as.integer(x) else as.numeric(x)
}

.calibrate_intercept <- function(linpred, rate) {
  f <- function(c) mean(plogis(c + linpred)) - rate
  if (sign(f(-20)) == sign(f(20))) {
    return(if (abs(f(-20)) < abs(f(20))) -20 else 20)
  }
  uniroot(f, c(-20, 20))$root
}

inject_missingness <- function(complete_data, tree, mechanism, rate,
                                vars = c("y", "x1", "x2", "x3", "x4")) {
  # ... identical to Dan's inject_missingness(), with species-coverage
  # safeguard included. Reproduced here so the simulation can run
  # standalone of BACE/dev/02_benchmark_simulated_full.R.

  n <- nrow(complete_data)
  miss_mask <- as.data.frame(matrix(FALSE, nrow = n, ncol = length(vars),
                                      dimnames = list(NULL, vars)))
  miss_data <- complete_data
  Sigma_phylo <- if (mechanism == "phylo_MAR")
    ape::vcv(tree, corr = TRUE) else NULL
  trait_mar_covariate <- c(y = "x3", x1 = "y", x2 = "y",
                            x3 = "x1", x4 = "y")

  for (v in vars) {
    linpred <- switch(mechanism,
      phylo_MAR = {
        z_sp <- MASS::mvrnorm(1, mu = rep(0, nrow(Sigma_phylo)),
                               Sigma = Sigma_phylo)
        names(z_sp) <- rownames(Sigma_phylo)
        -DEP_STRENGTH * z_sp[as.character(complete_data$species)]
      },
      trait_MAR = {
        cov_v <- trait_mar_covariate[[v]]
        val   <- .trait_to_numeric(complete_data[[cov_v]])
        -DEP_STRENGTH * as.numeric(scale(val))
      },
      trait_MNAR = {
        val <- .trait_to_numeric(complete_data[[v]])
        -DEP_STRENGTH * as.numeric(scale(val))
      }
    )
    c_hat <- .calibrate_intercept(linpred, rate)
    p     <- plogis(c_hat + linpred)
    miss  <- rbinom(n, 1, p) == 1
    miss_mask[[v]]  <- miss
    miss_data[[v]][miss] <- NA
  }

  # Species-coverage safeguard (pigauto doesn't have BACE's ginverse
  # requirement, but we keep this for cross-method parity).
  by_sp <- split(seq_len(n), complete_data$species)
  for (v in vars) {
    for (sp_rows in by_sp) {
      if (length(sp_rows) == 0L) next
      if (all(miss_mask[[v]][sp_rows])) {
        keep <- if (length(sp_rows) == 1L) sp_rows else sample(sp_rows, 1)
        miss_mask[[v]][keep] <- FALSE
        miss_data[[v]][keep] <- complete_data[[v]][keep]
      }
    }
  }
  list(miss_data = miss_data, miss_mask = miss_mask)
}

# =============================================================================
# 5. EVALUATOR (mirrors Dan's + adds gate-utilisation metric)
# =============================================================================

evaluate_imputation_ensemble <- function(true_data, imp_list, miss_mask,
                                          var_types) {
  # Per-variable: pool point-estimate metrics across imputations; treat
  # all imputations as a posterior sample for coverage95 / Brier.
  # Returns a long-format data.frame matching Dan's schema verbatim.
  out <- list()
  for (v in names(var_types)) {
    idx <- miss_mask[[v]]
    if (!any(idx)) next
    tv  <- true_data[[v]][idx]
    n_c <- sum(idx)
    if (var_types[v] %in% c("gaussian", "count")) {
      tv_num <- as.numeric(tv)
      iv_mat <- vapply(imp_list, function(d) as.numeric(d[[v]][idx]),
                        numeric(n_c))
      if (!is.matrix(iv_mat)) iv_mat <- matrix(iv_mat, nrow = n_c)
      sd_marg <- sd(as.numeric(true_data[[v]]), na.rm = TRUE)
      per_imp_nrmse <- apply(iv_mat, 2, function(iv) {
        if (!is.finite(sd_marg) || sd_marg <= 0) NA_real_
        else sqrt(mean((iv - tv_num)^2)) / sd_marg
      })
      per_imp_mae  <- apply(iv_mat, 2, function(iv) mean(abs(iv - tv_num)))
      per_imp_cor  <- apply(iv_mat, 2, function(iv) {
        if (length(unique(tv_num)) > 1 && length(unique(iv)) > 1)
          suppressWarnings(cor(iv, tv_num)) else NA_real_
      })
      lo <- apply(iv_mat, 1, stats::quantile, probs = 0.025, na.rm = TRUE,
                  names = FALSE)
      hi <- apply(iv_mat, 1, stats::quantile, probs = 0.975, na.rm = TRUE,
                  names = FALSE)
      coverage <- mean(tv_num >= lo & tv_num <= hi)
      metrics <- c("nrmse", "mae", "correlation", "coverage95")
      values  <- c(mean(per_imp_nrmse, na.rm = TRUE), mean(per_imp_mae),
                    mean(per_imp_cor, na.rm = TRUE), coverage)
      if (var_types[v] == "count") {
        per_imp_log_mae <- apply(iv_mat, 2, function(iv) {
          mean(abs(log1p(pmax(iv, 0)) - log1p(pmax(tv_num, 0))))
        })
        metrics <- c(metrics, "log_mae")
        values  <- c(values, mean(per_imp_log_mae))
      }
      out[[length(out) + 1]] <- data.frame(
        variable = v,
        type     = if (var_types[v] == "gaussian") "continuous" else "count",
        metric   = metrics, value = values,
        stringsAsFactors = FALSE
      )
    } else {
      tv_chr  <- as.character(tv)
      imp_chr <- vapply(imp_list, function(d) as.character(d[[v]][idx]),
                         character(n_c))
      if (!is.matrix(imp_chr)) imp_chr <- matrix(imp_chr, nrow = n_c)
      per_imp_acc <- apply(imp_chr, 2, function(iv) mean(iv == tv_chr))
      classes <- unique(tv_chr)
      per_imp_bal <- apply(imp_chr, 2, function(iv) {
        recalls <- vapply(classes, function(cl) {
          is_cl <- tv_chr == cl
          if (!any(is_cl)) NA_real_ else mean(iv[is_cl] == cl)
        }, numeric(1))
        mean(recalls, na.rm = TRUE)
      })
      all_levels <- if (is.factor(true_data[[v]]))
        levels(true_data[[v]]) else sort(unique(tv_chr))
      prob_mat <- matrix(0, nrow = n_c, ncol = length(all_levels),
                          dimnames = list(NULL, all_levels))
      for (k in all_levels) prob_mat[, k] <- rowMeans(imp_chr == k)
      y_indic <- matrix(0, nrow = n_c, ncol = length(all_levels),
                         dimnames = list(NULL, all_levels))
      for (i in seq_len(n_c)) y_indic[i, tv_chr[i]] <- 1
      brier <- mean(rowSums((prob_mat - y_indic)^2))
      metrics <- c("accuracy", "balanced_accuracy", "brier")
      values  <- c(mean(per_imp_acc), mean(per_imp_bal, na.rm = TRUE), brier)
      if (var_types[v] == "ordered") {
        lvls    <- levels(true_data[[v]])
        tv_rank <- match(tv_chr, lvls)
        per_imp_ord <- apply(imp_chr, 2, function(iv) {
          mean(abs(match(iv, lvls) - tv_rank), na.rm = TRUE)
        })
        metrics <- c(metrics, "ordinal_mae")
        values  <- c(values, mean(per_imp_ord))
      }
      out[[length(out) + 1]] <- data.frame(
        variable = v,
        type     = if (var_types[v] == "ordered") "ordered" else "categorical",
        metric   = metrics, value = values,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

# =============================================================================
# 6. SINGLE-REPLICATE DRIVER (pigauto-side)
# =============================================================================

run_one_sim <- function(sim_id, scenario, mechanism,
                         n_species, miss_rate, beta_sparsity, dgp) {
  sim <- .sim_complete_data(scenario,
                              n_cases       = n_species,
                              n_species     = n_species,
                              beta_sparsity = beta_sparsity,
                              dgp           = dgp)
  complete_data <- sim$complete_data
  tree          <- sim$tree

  miss   <- inject_missingness(complete_data, tree, mechanism, miss_rate)
  realised_rate <- colMeans(miss$miss_mask)

  # ---- pigauto fit -----------------------------------------------------
  t0 <- Sys.time()
  mi <- tryCatch(
    multi_impute(
      miss$miss_data, tree,
      m              = PIGAUTO_SIM_CONFIG$n_imputations,
      pool_method    = PIGAUTO_SIM_CONFIG$pool_method,
      clamp_outliers = PIGAUTO_SIM_CONFIG$clamp_outliers,
      seed           = PIGAUTO_SIM_CONFIG$seed + sim_id  # different per rep
    ),
    error = function(e) {
      message("  [", scenario, "/", mechanism, "/n=", n_species,
              "/miss=", miss_rate, "/beta=", beta_sparsity, "/dgp=", dgp,
              " sim ", sim_id, "] multi_impute error: ", e$message)
      NULL
    }
  )
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(mi)) return(NULL)

  var_types <- c(y  = "gaussian",  x1 = "categorical", x2 = "categorical",
                 x3 = "count",     x4 = "ordered")
  scores <- evaluate_imputation_ensemble(
    complete_data, mi$datasets, miss$miss_mask, var_types
  )

  # ---- pigauto-only extras: gate utilisation + wall time --------------
  gates <- if (!is.null(mi$fit) && !is.null(mi$fit$r_cal_gnn))
    mi$fit$r_cal_gnn else NA_real_
  gate_mean   <- if (!all(is.na(gates))) mean(gates, na.rm = TRUE) else NA_real_
  gate_active <- if (!all(is.na(gates))) mean(gates > 0, na.rm = TRUE) else NA_real_

  scores$sim_id        <- sim_id
  scores$scenario      <- scenario
  scores$mechanism     <- mechanism
  scores$n_species     <- n_species
  scores$miss_rate     <- miss_rate
  scores$beta_sparsity <- beta_sparsity
  scores$dgp           <- dgp
  scores$gate_mean     <- gate_mean
  scores$gate_active   <- gate_active
  scores$wall_sec      <- wall
  scores$realised_rate_mean <- mean(realised_rate)
  scores
}

# =============================================================================
# 7. MAIN SWEEP LOOP
# =============================================================================

design <- expand.grid(
  scenario      = SCENARIOS,
  mechanism     = MECHANISMS,
  n_species     = N_SPECIES_GRID,
  miss_rate     = MISS_RATE_GRID,
  beta_sparsity = BETA_GRID,
  dgp           = DGP_GRID,
  stringsAsFactors = FALSE
)

cat("=============================================================\n")
cat("  pigauto-side simulation sweep\n")
cat("  Design cells:", nrow(design), "; N_SIMS per cell:", N_SIMS, "\n")
cat("  Total replicates:", nrow(design) * N_SIMS, "\n")
cat("=============================================================\n\n")

all_results <- vector("list", nrow(design))
for (cell_i in seq_len(nrow(design))) {
  cell <- design[cell_i, ]
  cell_id <- paste(cell, collapse = "|")
  cat(sprintf("--- cell %d / %d: %s ---\n", cell_i, nrow(design), cell_id))
  t0 <- Sys.time()
  reps <- parallel::mclapply(seq_len(N_SIMS), function(i) {
    run_one_sim(sim_id        = i,
                 scenario      = cell$scenario,
                 mechanism     = cell$mechanism,
                 n_species     = cell$n_species,
                 miss_rate     = cell$miss_rate,
                 beta_sparsity = cell$beta_sparsity,
                 dgp           = cell$dgp)
  }, mc.cores = N_CORES)
  ok <- !vapply(reps, is.null, logical(1))
  cat(sprintf("  done %d/%d in %.1f min\n",
              sum(ok), N_SIMS,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  all_results[[cell_i]] <- do.call(rbind, reps[ok])
}

results_df <- do.call(rbind, all_results)
rownames(results_df) <- NULL
saveRDS(results_df, file.path(RESULTS_DIR, "pigauto_sweep_results.rds"))
write.csv(results_df, file.path(RESULTS_DIR, "pigauto_sweep_results.csv"),
          row.names = FALSE)

# =============================================================================
# 8. SUMMARY TABLES (stratified by every axis, with MC SE)
# =============================================================================

summary_table <- aggregate(
  value ~ variable + type + metric +
    scenario + mechanism + n_species + miss_rate + beta_sparsity + dgp,
  data = results_df,
  FUN  = function(x) c(mean = mean(x, na.rm = TRUE),
                        sd   = sd(x,   na.rm = TRUE),
                        mc_se = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))),
                        q025 = quantile(x, 0.025, na.rm = TRUE),
                        q975 = quantile(x, 0.975, na.rm = TRUE))
)
summary_flat <- cbind(
  summary_table[, c("variable", "type", "metric", "scenario", "mechanism",
                    "n_species", "miss_rate", "beta_sparsity", "dgp")],
  as.data.frame(summary_table$value)
)
names(summary_flat)[(ncol(summary_flat) - 4):ncol(summary_flat)] <-
  c("mean", "sd", "mc_se", "q025", "q975")
write.csv(summary_flat,
          file.path(RESULTS_DIR, "pigauto_sweep_summary.csv"),
          row.names = FALSE)
cat("\nSummary saved to:", RESULTS_DIR, "\n")

# =============================================================================
# 9. PLOTS (per-metric facet of n_species x miss_rate, faceted by signal)
# =============================================================================
# Deliberately deferred to a separate plotting script
# (script/sim_bench/plot_pigauto_sweep.R) to keep this script focused on
# computation. The summary CSV is the source of truth.
