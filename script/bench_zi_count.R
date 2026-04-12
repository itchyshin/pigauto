#!/usr/bin/env Rscript
#
# script/bench_zi_count.R
#
# Per-type benchmark: ZERO-INFLATED COUNT traits
# Primary sweep: zero fraction
# Secondary sweep: mean of non-zero counts
#
# Run with
#   cd pigauto && /usr/local/bin/Rscript script/bench_zi_count.R

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
out_rds <- file.path(here, "script", "bench_zi_count.rds")
out_md  <- file.path(here, "script", "bench_zi_count.md")
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

n_species   <- 300L
n_traits    <- 3L
n_reps      <- 5L
epochs      <- 500L

# Primary sweep: zero fraction
scenarios_primary <- c("zf_0.2", "zf_0.4", "zf_0.6", "zf_0.8")

# Secondary sweep: mean of non-zero counts at zero_frac = 0.5
secondary_scenarios <- c("mean_nz_5", "mean_nz_20", "mean_nz_100")

all_scenarios <- c(scenarios_primary, secondary_scenarios)
scenario_index <- setNames(seq_along(all_scenarios), all_scenarios)

primary_frac <- 0.25

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
         byrow = TRUE, dimnames = dimnames(X_latent))
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

  # Parse scenario
  if (grepl("^zf_", scenario)) {
    zero_frac <- as.numeric(sub("zf_", "", scenario))
    mean_nz <- 20
  } else if (grepl("^mean_nz_", scenario)) {
    zero_frac <- 0.5
    mean_nz <- as.numeric(sub("mean_nz_", "", scenario))
  } else {
    stop("Unknown scenario: ", scenario)
  }

  set.seed(rep_seed)
  tree <- ape::rtree(n_species)

  df <- simulate_zi_count_traits(tree, n_traits,
                                 zero_frac = zero_frac,
                                 mean_nz = mean_nz,
                                 seed = rep_seed)

  result <- tryCatch({
    pd  <- preprocess_traits(df, tree,
                             trait_types = setNames(
                               rep("zi_count", n_traits),
                               colnames(df)))
    spl <- make_missing_splits(pd$X_scaled, missing_frac = frac,
                               seed = rep_seed, trait_map = pd$trait_map)
    graph <- build_phylo_graph(tree, k_eigen = "auto")

    # Method 1: mean
    mean_pred <- mean_mode_impute(pd$X_scaled, spl)
    ev_mean <- evaluate_imputation(mean_pred, pd$X_scaled, spl,
                                   trait_map = pd$trait_map)

    # Method 2: baseline (LP for gate + BM for magnitude)
    bl <- fit_baseline(pd, tree, splits = spl, graph = graph)
    ev_bl <- evaluate_imputation(bl$mu, pd$X_scaled, spl,
                                 trait_map = pd$trait_map)

    # Method 3: pigauto
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
# Build full grid
# -------------------------------------------------------------------------

cells <- rbind(
  expand.grid(scenario = scenarios_primary,
              frac     = primary_frac,
              rep_id   = seq_len(n_reps),
              stringsAsFactors = FALSE),
  expand.grid(scenario = secondary_scenarios,
              frac     = primary_frac,
              rep_id   = seq_len(n_reps),
              stringsAsFactors = FALSE)
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

  parallel::clusterExport(cl, c("run_one_cell", "mean_mode_impute", "tag_rows",
                                 "n_species", "n_traits", "epochs",
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

  good <- list()
  errs <- list()
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
  n_traits            = n_traits,
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
  "# Zero-inflated count benchmark",
  "",
  sprintf("Run on: %s", format(Sys.time())),
  sprintf("Species: %d, traits: %d, reps: %d", n_species, n_traits, n_reps),
  sprintf("Total wall: %.1f min", total_wall / 60),
  ""
)

test_df <- all_results[all_results$split == "test", ]
show_cols <- c("method", "scenario", "trait", "rmse", "mae", "pearson_r",
               "zero_accuracy", "brier")

for (scen in c(scenarios_primary, secondary_scenarios)) {
  md <- c(md, sprintf("### %s", scen), "")
  sub <- test_df[test_df$scenario == scen, intersect(show_cols, names(test_df)),
                 drop = FALSE]
  if (nrow(sub) == 0L) {
    md <- c(md, "(no data)", "")
    next
  }
  avail_metrics <- intersect(c("rmse", "mae", "pearson_r", "zero_accuracy", "brier"),
                             names(sub)[vapply(sub, is.numeric, logical(1))])
  if (length(avail_metrics) == 0L) {
    md <- c(md, "(no numeric metrics)", "")
    next
  }
  fm <- as.formula(paste("cbind(",
                         paste(avail_metrics, collapse = ", "),
                         ") ~ method + trait"))
  avg <- aggregate(fm, data = sub, FUN = mean, na.rm = TRUE)
  avg <- avg[order(avg$trait, match(avg$method, c("mean", "baseline", "pigauto"))), ]
  md <- c(md, "```", capture.output(print(avg, row.names = FALSE)), "```", "")
}

writeLines(md, out_md)
log_line(sprintf("Wrote %s", out_md))
log_line("done")
