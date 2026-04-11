#!/usr/bin/env Rscript
#
# script/bench_continuous.R
#
# Per-type benchmark: CONTINUOUS traits under four evolutionary models
# (BM, OU, regime shift, nonlinear) at three missingness levels.
#
# Purpose
#   Isolate pigauto's performance on continuous traits across evolutionary
#   models that range from easy (BM — where the baseline is exact) to hard
#   (nonlinear inter-trait relationships that BM's linear covariance misses).
#
# Design
#   - Tree:      ape::rtree(300) — matches the bundled avonet300 scale.
#   - Traits:    4 continuous traits per scenario.
#   - Models:    BM, OU (alpha = 2), regime_shift, nonlinear.
#   - Missingness: primary sweep at missing_frac = 0.25;
#                  secondary sweep at {0.15, 0.30, 0.50} for BM and OU.
#   - Methods:   mean, BM baseline, pigauto (full pipeline).
#   - Replicates: 5 per cell.
#   - Seed:      rep * 100 + scenario_index.
#   - Parallelism: parallel::mclapply with mc.cores = 16
#
# Output
#   script/bench_continuous.rds    tidy results + metadata
#   script/bench_continuous.md     human-readable summary
#
# Run with
#   cd pigauto && /usr/local/bin/Rscript script/bench_continuous.R

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  library(parallel)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto",
    quiet = TRUE
  )
})

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "bench_continuous.rds")
out_md  <- file.path(here, "script", "bench_continuous.md")
MC_CORES <- 16L

script_start <- proc.time()[["elapsed"]]

log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start), ..., "\n",
      sep = "")
  flush.console()
}

# -------------------------------------------------------------------------
# Constants
# -------------------------------------------------------------------------

n_species   <- 300L
n_traits    <- 4L
n_reps      <- 5L
epochs      <- 500L

# Scenario definitions: name -> index mapping for seed derivation
scenarios_primary   <- c("BM", "OU", "regime_shift", "nonlinear")
scenario_index      <- setNames(seq_along(scenarios_primary), scenarios_primary)

# Primary sweep: single missingness level
primary_frac <- 0.25

# Secondary sweep: multiple missingness levels for BM and OU
secondary_scenarios <- c("BM", "OU")
secondary_fracs     <- c(0.15, 0.30, 0.50)

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

mean_mode_impute <- function(X_latent, splits) {
  X_train <- X_latent
  X_train[c(splits$val_idx, splits$test_idx)] <- NA
  col_means <- colMeans(X_train, na.rm = TRUE)
  col_means[!is.finite(col_means)] <- 0
  matrix(col_means,
         nrow = nrow(X_latent), ncol = ncol(X_latent),
         byrow = TRUE,
         dimnames = dimnames(X_latent))
}

tag_rows <- function(ev_df, method, scenario, frac, rep_id) {
  if (is.null(ev_df) || nrow(ev_df) == 0L) return(NULL)
  ev_df$method   <- method
  ev_df$scenario <- scenario
  ev_df$missing_frac <- frac
  ev_df$rep      <- rep_id
  ev_df
}

run_one_cell <- function(scenario, frac, rep_id) {
  rep_seed <- rep_id * 100L + scenario_index[scenario]

  # Generate tree
  set.seed(rep_seed)
  tree <- ape::rtree(n_species)

  # Generate traits
  if (scenario == "BM") {
    df <- simulate_bm_traits(tree, n_traits, seed = rep_seed)
  } else {
    df <- simulate_non_bm(tree, n_traits = n_traits,
                          scenario = scenario, seed = rep_seed)
  }

  # Pipeline
  result <- tryCatch({
    pd  <- preprocess_traits(df, tree, log_transform = FALSE)
    spl <- make_missing_splits(pd$X_scaled, missing_frac = frac,
                               seed = rep_seed, trait_map = pd$trait_map)
    graph <- build_phylo_graph(tree, k_eigen = "auto")

    # Method 1: mean
    mean_pred <- mean_mode_impute(pd$X_scaled, spl)
    ev_mean <- evaluate_imputation(mean_pred, pd$X_scaled, spl,
                                   trait_map = pd$trait_map)

    # Method 2: BM baseline
    bl <- fit_baseline(pd, tree, splits = spl, graph = graph)
    ev_bm <- evaluate_imputation(bl$mu, pd$X_scaled, spl,
                                 trait_map = pd$trait_map)

    # Method 3: pigauto
    fit <- fit_pigauto(pd, tree, splits = spl, baseline = bl,
                       graph = graph, epochs = epochs,
                       verbose = FALSE, seed = rep_seed)
    pred <- stats::predict(fit, return_se = TRUE, n_imputations = 1L)
    ev_pg <- evaluate_imputation(pred, pd$X_scaled, spl)

    rbind(
      tag_rows(ev_mean, "mean",    scenario, frac, rep_id),
      tag_rows(ev_bm,   "BM",     scenario, frac, rep_id),
      tag_rows(ev_pg,   "pigauto", scenario, frac, rep_id)
    )
  }, error = function(e) {
    data.frame(error = conditionMessage(e), scenario = scenario,
               missing_frac = frac, rep = rep_id, stringsAsFactors = FALSE)
  })

  result
}

# -------------------------------------------------------------------------
# Build the full grid of cells to run
# -------------------------------------------------------------------------

cells <- rbind(
  expand.grid(scenario = scenarios_primary,
              frac     = primary_frac,
              rep_id   = seq_len(n_reps),
              stringsAsFactors = FALSE),
  expand.grid(scenario = secondary_scenarios,
              frac     = secondary_fracs,
              rep_id   = seq_len(n_reps),
              stringsAsFactors = FALSE)
)
cells$cell_key <- sprintf("%s/%.2f/%d", cells$scenario, cells$frac, cells$rep_id)

# -------------------------------------------------------------------------
# Check for partial results and skip already-completed cells
# -------------------------------------------------------------------------

prior_results <- NULL
if (file.exists(out_rds)) {
  prior <- readRDS(out_rds)
  if (!is.null(prior$results) && nrow(prior$results) > 0L) {
    prior_results <- prior$results
    done_keys <- unique(sprintf("%s/%.2f/%d",
                                prior_results$scenario,
                                prior_results$missing_frac,
                                prior_results$rep))
    skip <- cells$cell_key %in% done_keys
    log_line(sprintf("Resuming: %d / %d cells already done, %d remaining",
                     sum(skip), nrow(cells), sum(!skip)))
    cells <- cells[!skip, , drop = FALSE]
  }
}

n_remaining <- nrow(cells)
log_line(sprintf("Running %d cells across %d cores (PSOCK cluster)", n_remaining, MC_CORES))

# -------------------------------------------------------------------------
# Parallel execution via PSOCK cluster (fork-safe for torch)
# -------------------------------------------------------------------------

if (n_remaining > 0L) {
  cell_list <- split(cells, seq_len(n_remaining))

  log_line("Starting PSOCK cluster...")
  cl <- parallel::makeCluster(min(MC_CORES, n_remaining))

  # Export helper functions and constants to workers
  parallel::clusterExport(cl, c("run_one_cell", "mean_mode_impute", "tag_rows",
                                 "n_species", "n_traits", "epochs",
                                 "scenario_index"),
                          envir = environment())

  # Load pigauto in each worker
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(ape)
      devtools::load_all(
        "/Users/z3437171/Dropbox/Github Local/pigauto",
        quiet = TRUE
      )
    })
  })
  log_line("Cluster ready. Dispatching cells...")

  par_results <- parallel::parLapply(cl, cell_list, function(row) {
    run_one_cell(row$scenario, row$frac, row$rep_id)
  })

  parallel::stopCluster(cl)
  log_line("Cluster stopped.")

  # Separate successes from errors
  good <- list()
  errs <- list()
  for (i in seq_along(par_results)) {
    r <- par_results[[i]]
    if (is.null(r)) {
      errs[[length(errs) + 1L]] <- cell_list[[i]]
    } else if ("error" %in% names(r)) {
      log_line(sprintf("  ERROR [%s]: %s", cell_list[[i]]$cell_key, r$error))
      errs[[length(errs) + 1L]] <- cell_list[[i]]
    } else {
      good[[length(good) + 1L]] <- r
    }
  }

  new_results <- if (length(good) > 0L) do.call(rbind, good) else NULL
  all_results <- rbind(prior_results, new_results)
  rownames(all_results) <- NULL
  log_line(sprintf("Completed: %d good, %d errors", length(good), length(errs)))
} else {
  all_results <- prior_results
  log_line("All cells already done. Skipping to summary.")
}

# -------------------------------------------------------------------------
# Finalise
# -------------------------------------------------------------------------

total_wall <- proc.time()[["elapsed"]] - script_start

saveRDS(list(
  results     = all_results,
  total_wall  = total_wall,
  n_species   = n_species,
  n_traits    = n_traits,
  n_reps      = n_reps,
  epochs      = epochs,
  primary_frac         = primary_frac,
  secondary_fracs      = secondary_fracs,
  scenarios_primary    = scenarios_primary,
  secondary_scenarios  = secondary_scenarios,
  commit = tryCatch(system("git rev-parse HEAD", intern = TRUE),
                    error = function(e) "unknown")
), out_rds)

log_line(sprintf("Total wall: %.1fs (%.1f min)", total_wall, total_wall / 60))

# ---- Markdown summary ----------------------------------------------------

machine <- tryCatch(
  sprintf("%s %s (%s), R %s",
          Sys.info()[["sysname"]], Sys.info()[["release"]],
          Sys.info()[["machine"]],
          paste(R.version$major, R.version$minor, sep = ".")),
  error = function(e) "machine info unavailable"
)

md <- c(
  "# Continuous-trait benchmark",
  "",
  sprintf("Run on: %s", format(Sys.time())),
  sprintf("Machine: %s", machine),
  sprintf("Species: %d, traits: %d, reps: %d", n_species, n_traits, n_reps),
  sprintf("Total wall: %.1f min", total_wall / 60),
  "",
  "## Methods",
  "",
  "- **mean**: column mean imputation (no phylogeny).",
  "- **BM**: Brownian-motion baseline (Rphylopars).",
  "- **pigauto**: full pipeline (BM baseline + calibrated GNN).",
  "",
  "## Primary sweep (missing_frac = 0.25)",
  ""
)

test_df <- all_results[all_results$split == "test", ]
show_cols <- c("method", "scenario", "trait", "rmse", "pearson_r", "mae")

for (scen in scenarios_primary) {
  md <- c(md, sprintf("### %s", scen), "")
  sub <- test_df[test_df$scenario == scen &
                   test_df$missing_frac == primary_frac, show_cols, drop = FALSE]
  if (nrow(sub) == 0L) {
    md <- c(md, "(no data)", "")
    next
  }
  avg <- aggregate(cbind(rmse, pearson_r, mae) ~ method + trait,
                   data = sub, FUN = mean, na.rm = TRUE)
  avg <- avg[order(avg$trait, match(avg$method, c("mean", "BM", "pigauto"))), ]
  md <- c(md, "```", capture.output(print(avg, row.names = FALSE)), "```", "")
}

md <- c(md, "## Secondary sweep (BM + OU, varying missingness)", "")
for (scen in secondary_scenarios) {
  for (frac in secondary_fracs) {
    md <- c(md, sprintf("### %s, missing_frac = %.2f", scen, frac), "")
    sub <- test_df[test_df$scenario == scen &
                     test_df$missing_frac == frac, show_cols, drop = FALSE]
    if (nrow(sub) == 0L) {
      md <- c(md, "(no data)", "")
      next
    }
    avg <- aggregate(cbind(rmse, pearson_r, mae) ~ method + trait,
                     data = sub, FUN = mean, na.rm = TRUE)
    avg <- avg[order(avg$trait, match(avg$method, c("mean", "BM", "pigauto"))), ]
    md <- c(md, "```", capture.output(print(avg, row.names = FALSE)), "```", "")
  }
}

writeLines(md, out_md)
log_line(sprintf("Wrote %s", out_md))
log_line("done")
