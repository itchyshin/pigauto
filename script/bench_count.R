#!/usr/bin/env Rscript
#
# script/bench_count.R
#
# Per-type benchmark: COUNT traits across mean-count levels and
# overdispersion settings.
#
# Purpose
#   Isolate pigauto's performance on count traits.  The primary sweep
#   varies mean count (5, 20, 100, 500) to test whether the log1p-z
#   pipeline handles sparse vs. dense counts.  The secondary sweep
#   compares Poisson vs. negative-binomial noise at mean_count = 20.
#
# Design
#   - Tree:      ape::rtree(300).
#   - Traits:    3 count traits per scenario.
#   - Primary sweep (x-axis = mean count):
#       mean_count in {5, 20, 100, 500} at missing_frac = 0.25.
#       Scenarios: mean_5, mean_20, mean_100, mean_500.
#   - Secondary sweep (x-axis = overdispersion):
#       Poisson (overdispersion = NULL) vs NegBin (overdispersion = 2)
#       at mean_count = 20, missing_frac = 0.25.
#       Scenarios: poisson, negbin.
#   - Methods:   mean (column-mean on log1p-z), baseline (BM via
#                Rphylopars), pigauto (full pipeline).
#   - Replicates: 5 per cell.
#   - Seed:      rep * 100 + scenario_index.
#   - Parallelism: PSOCK cluster (torch is not fork-safe).
#
# Output
#   script/bench_count.rds    tidy results + metadata
#   script/bench_count.md     human-readable summary
#
# Run with
#   cd pigauto && /usr/local/bin/Rscript script/bench_count.R

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
out_rds <- file.path(here, "script", "bench_count.rds")
out_md  <- file.path(here, "script", "bench_count.md")
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

n_species     <- 300L
n_traits      <- 3L
n_reps        <- 5L
epochs        <- 500L
missing_frac  <- 0.25

# All scenario names and their unique integer indices for seed derivation
all_scenarios  <- c("mean_5", "mean_20", "mean_100", "mean_500",
                    "poisson", "negbin")
scenario_index <- setNames(seq_along(all_scenarios), all_scenarios)

# Primary sweep: varying mean count
scenarios_primary <- c("mean_5", "mean_20", "mean_100", "mean_500")

# Secondary sweep: Poisson vs NegBin at mean_count = 20
scenarios_secondary <- c("poisson", "negbin")

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
  ev_df$method       <- method
  ev_df$scenario     <- scenario
  ev_df$missing_frac <- frac
  ev_df$rep          <- rep_id
  ev_df
}

parse_scenario <- function(scenario) {
  if (grepl("^mean_", scenario)) {
    mc <- as.integer(sub("^mean_", "", scenario))
    list(mean_count = mc, overdispersion = NULL)
  } else if (scenario == "poisson") {
    list(mean_count = 20L, overdispersion = NULL)
  } else if (scenario == "negbin") {
    list(mean_count = 20L, overdispersion = 2)
  } else {
    stop("Unknown scenario: ", scenario)
  }
}

run_one_cell <- function(scenario, frac, rep_id) {
  rep_seed <- rep_id * 100L + scenario_index[scenario]
  pars     <- parse_scenario(scenario)

  # Generate tree
  set.seed(rep_seed)
  tree <- ape::rtree(n_species)

  # Generate count traits
  df <- simulate_count_traits(tree, n_traits,
                              mean_count     = pars$mean_count,
                              overdispersion = pars$overdispersion,
                              seed           = rep_seed)

  # Pipeline
  result <- tryCatch({
    pd  <- preprocess_traits(df, tree, log_transform = FALSE)
    spl <- make_missing_splits(pd$X_scaled, missing_frac = frac,
                               seed = rep_seed, trait_map = pd$trait_map)
    graph <- build_phylo_graph(tree, k_eigen = "auto")

    # Method 1: mean imputation (column means on log1p-z scale)
    mean_pred <- mean_mode_impute(pd$X_scaled, spl)
    ev_mean <- evaluate_imputation(mean_pred, pd$X_scaled, spl,
                                   trait_map = pd$trait_map)

    # Method 2: BM baseline
    bl <- fit_baseline(pd, tree, splits = spl, graph = graph)
    ev_bl <- evaluate_imputation(bl$mu, pd$X_scaled, spl,
                                 trait_map = pd$trait_map)

    # Method 3: pigauto (full pipeline)
    fit <- fit_pigauto(pd, tree, splits = spl, baseline = bl,
                       graph = graph, epochs = epochs,
                       verbose = FALSE, seed = rep_seed)
    pred <- stats::predict(fit, return_se = TRUE, n_imputations = 1L)
    ev_pg <- evaluate_imputation(pred, pd$X_scaled, spl)

    rbind(
      tag_rows(ev_mean, "mean",     scenario, frac, rep_id),
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

cells <- rbind(
  expand.grid(scenario = scenarios_primary,
              frac     = missing_frac,
              rep_id   = seq_len(n_reps),
              stringsAsFactors = FALSE),
  expand.grid(scenario = scenarios_secondary,
              frac     = missing_frac,
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
# Parallel execution via PSOCK cluster (torch is not fork-safe)
# -------------------------------------------------------------------------

if (n_remaining > 0L) {
  cell_list <- split(cells, seq_len(n_remaining))

  log_line("Starting PSOCK cluster...")
  cl <- parallel::makeCluster(min(MC_CORES, n_remaining))

  # Export helper functions and constants to workers
  parallel::clusterExport(cl, c("run_one_cell", "mean_mode_impute", "tag_rows",
                                 "parse_scenario",
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
  "# Count-trait benchmark",
  "",
  sprintf("Run on: %s", format(Sys.time())),
  sprintf("Machine: %s", machine),
  sprintf("Species: %d, traits: %d, reps: %d, missing_frac: %.2f",
          n_species, n_traits, n_reps, missing_frac),
  sprintf("Total wall: %.1f min", total_wall / 60),
  "",
  "## Methods",
  "",
  "- **mean**: column-mean imputation on log1p-z scale (no phylogeny).",
  "- **baseline**: Brownian-motion baseline on log1p-z scale (Rphylopars).",
  "- **pigauto**: full pipeline (BM baseline + calibrated GNN).",
  "",
  "## Primary sweep (x-axis = mean count)",
  ""
)

test_df   <- all_results[all_results$split == "test", ]
show_cols <- c("method", "scenario", "trait", "rmse", "mae", "pearson_r")

for (scen in scenarios_primary) {
  md <- c(md, sprintf("### %s", scen), "")
  sub <- test_df[test_df$scenario == scen, show_cols, drop = FALSE]
  if (nrow(sub) == 0L) {
    md <- c(md, "(no data)", "")
    next
  }
  avg <- aggregate(cbind(rmse, mae, pearson_r) ~ method + trait,
                   data = sub, FUN = mean, na.rm = TRUE)
  avg <- avg[order(avg$trait, match(avg$method, c("mean", "baseline", "pigauto"))), ]
  md <- c(md, "```", capture.output(print(avg, row.names = FALSE)), "```", "")
}

md <- c(md, "## Secondary sweep (Poisson vs NegBin, mean_count = 20)", "")
for (scen in scenarios_secondary) {
  md <- c(md, sprintf("### %s", scen), "")
  sub <- test_df[test_df$scenario == scen, show_cols, drop = FALSE]
  if (nrow(sub) == 0L) {
    md <- c(md, "(no data)", "")
    next
  }
  avg <- aggregate(cbind(rmse, mae, pearson_r) ~ method + trait,
                   data = sub, FUN = mean, na.rm = TRUE)
  avg <- avg[order(avg$trait, match(avg$method, c("mean", "baseline", "pigauto"))), ]
  md <- c(md, "```", capture.output(print(avg, row.names = FALSE)), "```", "")
}

writeLines(md, out_md)
log_line(sprintf("Wrote %s", out_md))
log_line("done")
