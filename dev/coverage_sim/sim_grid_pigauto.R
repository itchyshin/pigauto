# =============================================================================
# Phase 2 scaffold: pigauto calibration grid (Dan's framework, extended)
# =============================================================================
#
# Grid:
#   4 phylogenetic-signal scenarios (all_high, all_moderate, all_low, mixed)
#   3 missingness mechanisms (phylo_MAR, trait_MAR, trait_MNAR)
#   2 conformal_method ("split", "bootstrap")
#   N_REPS replicates per cell (default 10, bump to 50 for final)
#
# Output per cell:
#   - coverage95 on each trait
#   - NRMSE / Brier / accuracy per trait
#   - calibrated gate and conformal score per latent col
#   - wall time
#
# Extensions vs Dan's BACE sim:
#   - 8 trait types instead of 5 (adds proportion, zi_count, multi_proportion
#     later; v1 keeps Dan's 5 to enable head-to-head comparison).
#   - 2 conformal methods crossed in (tests the Phase 3 fix).
#
# Expected runtime (single-obs n=150, epochs=500, ~3 min per fit):
#   4 signals x 3 mechanisms x 2 methods x 10 reps = 240 fits x 3 min = 12 hr
#
# For the "pre-Vulcan" local test: shrink to N_REPS=3 -> ~4 hr locally.
#
# NOT RUN YET. Waiting for Phase 3c bootstrap-vs-split probe to validate
# that the fix actually reduces conformal-score variance before investing
# 12+ hours of compute on the full grid.
# =============================================================================

library(BACE)
suppressMessages(devtools::load_all(".", quiet = TRUE))
library(ape); library(MASS)

# Helpers from Dan's sim
source_funcs <- function(path, fn_names) {
  src <- readLines(path)
  for (fn in fn_names) {
    s <- grep(paste0("^", fn, " <- function"), src)[1]
    d <- 0; e <- NA
    for (i in s:length(src)) {
      d <- d + sum(gregexpr("\\{", src[i])[[1]] > 0) -
               sum(gregexpr("\\}", src[i])[[1]] > 0)
      if (i > s && d == 0) { e <- i; break }
    }
    eval(parse(text = src[s:e]), envir = globalenv())
  }
}
SIGNAL <- list(high=0.90, moderate=0.60, low=0.20); DEP_STRENGTH <- 1.5
source_funcs("dev/coverage_sim/sim_bace_source.R",
             c(".trait_to_numeric", ".calibrate_intercept",
               "evaluate_imputation_ensemble", "make_phylo_signal"))

# Single-obs-safe injector
inject_singleobs <- function(cd, tree, mech, rate, vars=c("y","x1","x2","x3","x4")) {
  n <- nrow(cd)
  mm <- as.data.frame(matrix(FALSE, n, length(vars), dimnames=list(NULL, vars)))
  md <- cd
  Sp <- if (mech == "phylo_MAR") ape::vcv(tree, corr=TRUE) else NULL
  tmc <- c(y="x3", x1="y", x2="y", x3="x1", x4="y")
  for (v in vars) {
    lp <- switch(mech,
      phylo_MAR = {
        z <- MASS::mvrnorm(1, mu=rep(0, nrow(Sp)), Sigma=Sp)
        names(z) <- rownames(Sp); -DEP_STRENGTH * z[rownames(cd)]
      },
      trait_MAR = -DEP_STRENGTH * as.numeric(scale(.trait_to_numeric(cd[[tmc[[v]]]]))),
      trait_MNAR = -DEP_STRENGTH * as.numeric(scale(.trait_to_numeric(cd[[v]])))
    )
    c_hat <- .calibrate_intercept(lp, rate)
    mm[[v]] <- rbinom(n, 1, plogis(c_hat + lp)) == 1
    md[[v]][mm[[v]]] <- NA
  }
  list(miss_data=md, miss_mask=mm)
}

# ---- Configuration ---------------------------------------------------------
CONFIG <- list(
  scenarios   = c("all_high", "all_moderate", "all_low", "mixed"),
  mechanisms  = c("phylo_MAR", "trait_MAR", "trait_MNAR"),
  methods     = c("split", "bootstrap"),
  n_reps      = 10,          # bump to 50 for final publication grid
  n_species   = 150,         # single-obs; tree size matches Dan's setup
  miss_rate   = 0.35,
  epochs      = 500L,
  out_path    = "dev/coverage_sim/results/grid_v1.rds"
)

run_one <- function(scenario, mechanism, method, rep_id) {
  seed_i <- rep_id * 101L + which(CONFIG$scenarios == scenario) * 10L +
            which(CONFIG$mechanisms == mechanism)
  set.seed(seed_i)

  sim <- tryCatch(
    sim_bace(response_type="gaussian",
             predictor_types=c("binary","multinomial3","poisson","threshold3"),
             var_names=c("y","x1","x2","x3","x4"),
             phylo_signal=make_phylo_signal(scenario),
             n_cases=CONFIG$n_species, n_species=CONFIG$n_species,
             beta_sparsity=0.3, missingness=c(0,0,0,0,0)),
    error = function(e) NULL)
  if (is.null(sim)) return(NULL)

  cd <- sim$data
  tree <- sim$tree
  rownames(cd) <- cd$species; cd$species <- NULL
  cd$x1 <- factor(as.character(cd$x1), levels=c("0","1"))

  miss <- tryCatch(inject_singleobs(cd, tree, mechanism, CONFIG$miss_rate),
                   error=function(e) NULL)
  if (is.null(miss)) return(NULL)

  t0 <- Sys.time()
  mi <- tryCatch(
    multi_impute(traits=miss$miss_data, tree=tree, m=20L,
                 species_col=NULL, draws_method="conformal",
                 log_transform=FALSE, missing_frac=0.25,
                 epochs=CONFIG$epochs, verbose=FALSE, seed=seed_i,
                 conformal_method=method),
    error = function(e) e)
  if (inherits(mi, "error")) return(NULL)
  wall_s <- as.numeric(difftime(Sys.time(), t0, units="secs"))

  var_types <- c(y="gaussian", x1="categorical", x2="categorical",
                 x3="count",    x4="ordered")
  scores <- tryCatch(
    evaluate_imputation_ensemble(cd, mi$datasets, miss$miss_mask, var_types),
    error = function(e) NULL)
  if (is.null(scores)) return(NULL)

  scores$scenario  <- scenario
  scores$mechanism <- mechanism
  scores$method    <- method
  scores$rep       <- rep_id
  scores$wall_s    <- wall_s
  scores
}

# ---- Loop ------------------------------------------------------------------
cat("=== pigauto calibration grid v1 ===\n")
cat(sprintf("scenarios=%s, mechanisms=%s, methods=%s, n_reps=%d\n",
            paste(CONFIG$scenarios, collapse=","),
            paste(CONFIG$mechanisms, collapse=","),
            paste(CONFIG$methods, collapse=","),
            CONFIG$n_reps))
cat(sprintf("total fits = %d x %d x %d x %d = %d (~%.0f min at 3 min/fit)\n\n",
            length(CONFIG$scenarios), length(CONFIG$mechanisms),
            length(CONFIG$methods), CONFIG$n_reps,
            length(CONFIG$scenarios) * length(CONFIG$mechanisms) *
              length(CONFIG$methods) * CONFIG$n_reps,
            length(CONFIG$scenarios) * length(CONFIG$mechanisms) *
              length(CONFIG$methods) * CONFIG$n_reps * 3))

results <- list()
cell_count <- 0L
for (sc in CONFIG$scenarios) {
  for (mech in CONFIG$mechanisms) {
    for (meth in CONFIG$methods) {
      for (rep_i in seq_len(CONFIG$n_reps)) {
        cell_count <- cell_count + 1L
        cat(sprintf("[%3d] %s / %s / %s / rep %d ...",
                    cell_count, sc, mech, meth, rep_i))
        flush.console()

        r <- run_one(sc, mech, meth, rep_i)
        if (!is.null(r)) {
          results[[length(results) + 1]] <- r
          # save partial progress every cell
          saveRDS(do.call(rbind, results), CONFIG$out_path)
          cov_y <- r$value[r$variable == "y" & r$metric == "coverage95"]
          cat(sprintf(" cov_y=%.2f (%.0fs)\n",
                      if (length(cov_y)) cov_y else NA,
                      r$wall_s[1]))
        } else {
          cat(" FAILED\n")
        }
      }
    }
  }
}

cat("\n=== DONE ===\n")
cat(sprintf("Saved %d cells to %s\n", length(results), CONFIG$out_path))
