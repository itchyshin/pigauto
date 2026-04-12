#!/usr/bin/env Rscript
#
# script/bench_categorical.R
#
# Per-type benchmark: CATEGORICAL traits across numbers of categories (K)
# and phylogenetic signal levels.
#
# Purpose
#   Isolate pigauto's performance on categorical traits.  The primary sweep
#   varies K (number of categories) at fixed phylogenetic signal; the
#   secondary sweep varies signal at fixed K = 5.
#
# Design
#   - Tree:      ape::rtree(300) — matches the bundled avonet300 scale.
#   - Traits:    2 categorical traits per scenario.
#   - Primary sweep (x-axis = K):
#       K in {3, 5, 8, 12}, signal = 0.8, missing_frac = 0.25
#       Scenarios: "K_3", "K_5", "K_8", "K_12"
#   - Secondary sweep (x-axis = signal):
#       signal in {0.3, 0.6, 1.0} at K = 5, missing_frac = 0.25
#       Scenarios: "signal_0.3", "signal_0.6", "signal_1.0"
#   - Methods:   mode (column-mean imputation), baseline (label prop), pigauto.
#   - Replicates: 5 per cell.
#   - Seed:      rep * 100 + scenario_index.
#   - Parallelism: PSOCK cluster (torch is not fork-safe).
#
# Output
#   script/bench_categorical.rds    tidy results + metadata
#   script/bench_categorical.md     human-readable summary
#
# Run with
#   cd pigauto && /usr/local/bin/Rscript script/bench_categorical.R

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
out_rds <- file.path(here, "script", "bench_categorical.rds")
out_md  <- file.path(here, "script", "bench_categorical.md")
MC_CORES <- as.integer(Sys.getenv("MC_CORES", unset = "16"))

script_start <- proc.time()[["elapsed"]]

log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start), ..., "\n",
      sep = "")
  flush.console()
}

# -------------------------------------------------------------------------
# Constants
# -------------------------------------------------------------------------

n_species     <- 300L
n_traits      <- 2L
n_reps        <- 5L
epochs        <- 500L
missing_frac  <- 0.25

# Scenario definitions: name -> index mapping for seed derivation
all_scenarios  <- c("K_3", "K_5", "K_8", "K_12",
                    "signal_0.3", "signal_0.6", "signal_1.0")
scenario_index <- setNames(seq_along(all_scenarios), all_scenarios)

scenarios_primary   <- c("K_3", "K_5", "K_8", "K_12")
scenarios_secondary <- c("signal_0.3", "signal_0.6", "signal_1.0")

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

parse_scenario <- function(scenario) {
  # "K_8"      -> n_levels = 8,  signal = 0.8
  # "signal_0.6" -> n_levels = 5, signal = 0.6
  if (grepl("^K_", scenario)) {
    n_levels <- as.integer(sub("^K_", "", scenario))
    signal   <- 0.8
  } else {
    n_levels <- 5L
    signal   <- as.numeric(sub("^signal_", "", scenario))
  }
  list(n_levels = n_levels, signal = signal)
}

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
  ev_df$method       <- method
  ev_df$scenario     <- scenario
  ev_df$missing_frac <- frac
  ev_df$rep          <- rep_id
  ev_df
}

run_one_cell <- function(scenario, frac, rep_id) {
  rep_seed <- rep_id * 100L + scenario_index[scenario]
  params   <- parse_scenario(scenario)

  # Generate tree
  set.seed(rep_seed)
  tree <- ape::rtree(n_species)

  # Generate categorical traits
  df <- simulate_categorical_traits(tree, n_traits = n_traits,
                                    n_levels = params$n_levels,
                                    signal = params$signal,
                                    seed = rep_seed)

  # Pipeline
  result <- tryCatch({
    pd  <- preprocess_traits(df, tree, log_transform = FALSE)
    spl <- make_missing_splits(pd$X_scaled, missing_frac = frac,
                               seed = rep_seed, trait_map = pd$trait_map)
    graph <- build_phylo_graph(tree, k_eigen = "auto")

    # Method 1: mode (column-mean imputation in latent space)
    mode_pred <- mean_mode_impute(pd$X_scaled, spl)
    ev_mode   <- evaluate_imputation(mode_pred, pd$X_scaled, spl,
                                     trait_map = pd$trait_map)

    # Method 2: baseline (phylogenetic label propagation)
    bl    <- fit_baseline(pd, tree, splits = spl, graph = graph)
    ev_bl <- evaluate_imputation(bl$mu, pd$X_scaled, spl,
                                 trait_map = pd$trait_map)

    # Method 3: pigauto (full pipeline)
    fit <- fit_pigauto(pd, tree, splits = spl, baseline = bl,
                       graph = graph, epochs = epochs,
                       verbose = FALSE, seed = rep_seed)
    pred  <- stats::predict(fit, return_se = TRUE, n_imputations = 1L)
    ev_pg <- evaluate_imputation(pred, pd$X_scaled, spl)

    rbind(
      tag_rows(ev_mode, "mode",     scenario, frac, rep_id),
      tag_rows(ev_bl,   "baseline", scenario, frac, rep_id),
      tag_rows(ev_pg,   "pigauto",  scenario, frac, rep_id)
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

cells <- expand.grid(
  scenario = all_scenarios,
  frac     = missing_frac,
  rep_id   = seq_len(n_reps),
  stringsAsFactors = FALSE
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
log_line(sprintf("Running %d cells across %d cores (PSOCK cluster)",
                 n_remaining, MC_CORES))

# -------------------------------------------------------------------------
# Parallel execution via PSOCK cluster (fork-safe for torch)
# -------------------------------------------------------------------------

if (n_remaining > 0L) {
  cell_list <- split(cells, seq_len(n_remaining))

  log_line("Starting PSOCK cluster...")
  cl <- parallel::makeCluster(min(MC_CORES, n_remaining))

  # Export helper functions and constants to workers
  parallel::clusterExport(cl, c("run_one_cell", "mean_mode_impute", "tag_rows",
                                 "parse_scenario",
                                 "n_species", "n_traits", "epochs",
                                 "missing_frac", "scenario_index"),
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
  results              = all_results,
  total_wall           = total_wall,
  n_species            = n_species,
  n_traits             = n_traits,
  n_reps               = n_reps,
  epochs               = epochs,
  missing_frac         = missing_frac,
  scenarios_primary    = scenarios_primary,
  scenarios_secondary  = scenarios_secondary,
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
  "# Categorical-trait benchmark",
  "",
  sprintf("Run on: %s", format(Sys.time())),
  sprintf("Machine: %s", machine),
  sprintf("Species: %d, traits: %d, reps: %d, missing_frac: %.2f",
          n_species, n_traits, n_reps, missing_frac),
  sprintf("Total wall: %.1f min", total_wall / 60),
  "",
  "## Methods",
  "",
  "- **mode**: column-mean imputation in latent (one-hot) space; argmax determines predicted category.",
  "- **baseline**: phylogenetic label propagation.",
  "- **pigauto**: full pipeline (label-propagation baseline + calibrated GNN).",
  "",
  "## Primary sweep (varying K, signal = 0.8)",
  ""
)

test_df   <- all_results[all_results$split == "test", ]
show_cols <- c("method", "scenario", "trait", "accuracy")

for (scen in scenarios_primary) {
  md <- c(md, sprintf("### %s", scen), "")
  sub <- test_df[test_df$scenario == scen, show_cols, drop = FALSE]
  if (nrow(sub) == 0L) {
    md <- c(md, "(no data)", "")
    next
  }
  avg <- aggregate(accuracy ~ method + trait, data = sub, FUN = mean,
                   na.rm = TRUE)
  avg <- avg[order(avg$trait, match(avg$method,
                                    c("mode", "baseline", "pigauto"))), ]
  md <- c(md, "```", capture.output(print(avg, row.names = FALSE)), "```", "")
}

md <- c(md, "## Secondary sweep (K = 5, varying signal)", "")

for (scen in scenarios_secondary) {
  md <- c(md, sprintf("### %s", scen), "")
  sub <- test_df[test_df$scenario == scen, show_cols, drop = FALSE]
  if (nrow(sub) == 0L) {
    md <- c(md, "(no data)", "")
    next
  }
  avg <- aggregate(accuracy ~ method + trait, data = sub, FUN = mean,
                   na.rm = TRUE)
  avg <- avg[order(avg$trait, match(avg$method,
                                    c("mode", "baseline", "pigauto"))), ]
  md <- c(md, "```", capture.output(print(avg, row.names = FALSE)), "```", "")
}

# Grand summary across all scenarios
md <- c(md, "## Grand summary (accuracy averaged across traits and reps)", "")
if (nrow(test_df) > 0L) {
  grand <- aggregate(accuracy ~ method + scenario, data = test_df, FUN = mean,
                     na.rm = TRUE)
  grand <- grand[order(grand$scenario,
                       match(grand$method, c("mode", "baseline", "pigauto"))), ]
  md <- c(md, "```", capture.output(print(grand, row.names = FALSE)), "```", "")
}

writeLines(md, out_md)
log_line(sprintf("Wrote %s", out_md))
log_line("done")
