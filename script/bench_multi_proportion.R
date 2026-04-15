#!/usr/bin/env Rscript
#
# script/bench_multi_proportion.R
#
# Per-type benchmark: MULTI_PROPORTION traits (K-category compositions,
# rows sum to 1). Encoded via CLR + z-score, Brownian-motion baseline
# on CLR space, softmax decode.
#
# Primary sweep: phylogenetic signal
# Secondary sweep: K (number of categories: 3, 5, 8, 12)
#
# Run with
#   cd pigauto && Rscript script/bench_multi_proportion.R

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
out_rds <- file.path(here, "script", "bench_multi_proportion.rds")
out_md  <- file.path(here, "script", "bench_multi_proportion.md")
MC_CORES <- as.integer(Sys.getenv("MC_CORES", unset = "4"))

script_start <- proc.time()[["elapsed"]]

log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start), ..., "\n",
      sep = "")
  flush.console()
}

# -------------------------------------------------------------------------
# Constants
# -------------------------------------------------------------------------

n_species <- 300L
n_reps    <- 5L
epochs    <- 500L

# Primary sweep: signal strength (fixed K = 5)
scenarios_primary <- c("signal_0.2", "signal_0.4", "signal_0.6",
                       "signal_0.8", "signal_1.0")

# Secondary sweep: number of components (fixed signal = 0.6)
secondary_scenarios <- c("K_3", "K_5", "K_8", "K_12")

all_scenarios <- c(scenarios_primary, secondary_scenarios)
scenario_index <- setNames(seq_along(all_scenarios), all_scenarios)

primary_frac <- 0.25

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

tag_rows <- function(ev_df, method, scenario, frac, rep_id) {
  if (is.null(ev_df) || !is.data.frame(ev_df)) return(NULL)
  ev_df$method       <- method
  ev_df$scenario     <- scenario
  ev_df$missing_frac <- frac
  ev_df$rep          <- rep_id
  ev_df
}

# Naive mean-imputation baseline for compositions: use the training-row
# geometric-mean composition (= row-wise geometric mean then renormalise).
mean_composition_impute <- function(pd, splits) {
  X <- pd$X_scaled
  mask <- splits$mask   # TRUE = observed + training
  X_out <- X
  # Fill held-out cells with the column mean of TRAINING cells only
  for (j in seq_len(ncol(X))) {
    train_j <- X[mask[, j], j]
    train_j <- train_j[is.finite(train_j)]
    mu_j <- if (length(train_j) > 0L) mean(train_j) else 0
    hold <- which(!mask[, j])
    X_out[hold, j] <- mu_j
  }
  X_out
}

run_one_cell <- function(scenario, frac, rep_id) {
  rep_seed <- rep_id * 100L + scenario_index[scenario]

  # Parse scenario
  if (grepl("^signal_", scenario)) {
    signal <- as.numeric(sub("signal_", "", scenario))
    K <- 5L
  } else if (grepl("^K_", scenario)) {
    signal <- 0.6
    K <- as.integer(sub("K_", "", scenario))
  } else {
    stop("Unknown scenario: ", scenario)
  }

  set.seed(rep_seed)
  tree <- ape::rtree(n_species)
  df   <- simulate_multi_proportion_traits(tree, K = K, signal = signal,
                                           seed = rep_seed)
  group_cols <- names(df)
  group_name <- "comp"   # distinct from any column name

  result <- tryCatch({
    pd  <- preprocess_traits(df, tree,
                             multi_proportion_groups =
                               setNames(list(group_cols), group_name))
    spl <- make_missing_splits(pd$X_scaled, missing_frac = frac,
                               seed = rep_seed, trait_map = pd$trait_map)
    graph <- build_phylo_graph(tree, k_eigen = "auto")

    # Method 1: naive column-mean on CLR (no phylogeny)
    mean_pred <- mean_composition_impute(pd, spl)
    ev_mean   <- evaluate_imputation(mean_pred, pd$X_scaled, spl,
                                     trait_map = pd$trait_map)

    # Method 2: BM baseline (K independent BMs on CLR)
    bl <- fit_baseline(pd, tree, splits = spl, graph = graph)
    ev_bl <- evaluate_imputation(bl$mu, pd$X_scaled, spl,
                                 trait_map = pd$trait_map)

    # Method 3: pigauto (BM + GNN + calibrated gate)
    fit <- fit_pigauto(pd, tree, splits = spl, baseline = bl,
                       graph = graph, epochs = epochs,
                       verbose = FALSE, seed = rep_seed)
    pred <- stats::predict(fit, return_se = FALSE, n_imputations = 1L)
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
# Build full grid
# -------------------------------------------------------------------------

cells <- rbind(
  expand.grid(scenario = scenarios_primary, frac = primary_frac,
              rep_id = seq_len(n_reps), stringsAsFactors = FALSE),
  expand.grid(scenario = secondary_scenarios, frac = primary_frac,
              rep_id = seq_len(n_reps), stringsAsFactors = FALSE)
)
cells$cell_key <- sprintf("%s/%.2f/%d", cells$scenario, cells$frac, cells$rep_id)

# -------------------------------------------------------------------------
# Resume from partial RDS
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
# Parallel execution
# -------------------------------------------------------------------------

if (n_remaining > 0L) {
  cell_list <- split(cells, seq_len(n_remaining))

  log_line("Starting PSOCK cluster...")
  cl <- parallel::makeCluster(min(MC_CORES, n_remaining))

  parallel::clusterExport(cl, c("run_one_cell", "mean_composition_impute",
                                 "tag_rows", "n_species", "epochs",
                                 "scenario_index"),
                          envir = environment())
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

  good <- list(); errs <- list()
  for (i in seq_along(par_results)) {
    r_i <- par_results[[i]]
    if (is.null(r_i)) {
      errs[[length(errs) + 1L]] <- cell_list[[i]]
    } else if ("error" %in% names(r_i)) {
      log_line(sprintf("  ERROR [%s]: %s", cell_list[[i]]$cell_key, r_i$error))
      errs[[length(errs) + 1L]] <- cell_list[[i]]
    } else {
      good[[length(good) + 1L]] <- r_i
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
  results             = all_results,
  total_wall          = total_wall,
  n_species           = n_species,
  n_reps              = n_reps,
  epochs              = epochs,
  primary_frac        = primary_frac,
  scenarios_primary   = scenarios_primary,
  secondary_scenarios = secondary_scenarios,
  commit = tryCatch(system("git rev-parse HEAD", intern = TRUE),
                    error = function(e) "unknown")
), out_rds)

log_line(sprintf("Total wall: %.1fs (%.1f min)", total_wall, total_wall / 60))

# ---- Markdown summary ----------------------------------------------------

md <- c(
  "# Multi-proportion trait benchmark (compositional data)",
  "",
  sprintf("Run on: %s", format(Sys.time())),
  sprintf("Species: %d, reps: %d", n_species, n_reps),
  sprintf("Total wall: %.1f min", total_wall / 60),
  "",
  "## Methods",
  "",
  "- **mean**: column-mean imputation on CLR-z latent scale.",
  "- **baseline**: per-component Brownian-motion imputation on CLR-z latent scale.",
  "- **pigauto**: full pipeline (BM baseline + calibrated GNN) on CLR-z latent scale.",
  "",
  "Metrics:",
  "",
  "- `aitchison`: average Aitchison distance (Euclidean in CLR space) across held-out rows.",
  "- `rmse_clr`: RMSE on z-scored CLR values.",
  "- `simplex_mae`: mean absolute error on the simplex after softmax decode.",
  ""
)

test_df <- all_results[all_results$split == "test", ]
for (scen in c(scenarios_primary, secondary_scenarios)) {
  md <- c(md, sprintf("### %s", scen), "")
  sub <- test_df[test_df$scenario == scen, , drop = FALSE]
  if (nrow(sub) == 0L) { md <- c(md, "(no data)", ""); next }
  metrics <- c("aitchison", "rmse_clr", "simplex_mae")
  agg_list <- lapply(metrics, function(m) {
    sub_m <- sub[sub$metric == m, , drop = FALSE]
    if (nrow(sub_m) == 0L) return(NULL)
    ag <- aggregate(value ~ method + trait, data = sub_m, FUN = mean, na.rm = TRUE)
    names(ag)[3] <- m
    ag
  })
  agg_list <- Filter(Negate(is.null), agg_list)
  if (length(agg_list) == 0L) { md <- c(md, "(no metrics)", ""); next }
  avg <- Reduce(function(a, b) merge(a, b, by = c("method", "trait"), all = TRUE),
                agg_list)
  avg <- avg[order(avg$trait, match(avg$method, c("mean", "baseline", "pigauto"))), ]
  md <- c(md, "```", capture.output(print(avg, row.names = FALSE)), "```", "")
}

writeLines(md, out_md)
log_line(sprintf("Wrote %s", out_md))
log_line("done")
