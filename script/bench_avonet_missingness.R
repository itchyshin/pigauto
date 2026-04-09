#!/usr/bin/env Rscript
#
# script/bench_avonet_missingness.R
#
# Missingness-sweep benchmark on the full AVONET3 + BirdTree dataset
# (avonet_full / tree_full bundled in pigauto >= 0.3.2).
#
# Purpose
#   Compare three imputation methods at three MCAR missingness levels on
#   a real 9,993-species mixed-type dataset, showing how each method
#   degrades as the held-out fraction grows.
#
# Design
#   - Data:     avonet_full (9,993 x 7 traits) + tree_full.
#   - Methods:  mean/mode, BM baseline, pigauto (full pipeline).
#   - Levels:   missing_frac in {0.20, 0.50, 0.80}, val_frac = 0.5.
#   - Split:    trait-level (categorical traits masked as whole groups).
#   - Metrics:  per-trait RMSE + r + accuracy via evaluate_imputation().
#   - Seed:     single seed per cell (2026L).
#   - Hyperparams for pigauto: copied verbatim from
#     script/validate_avonet_full.R for comparability.
#
# Output
#   script/bench_avonet_missingness.rds    tidy results + timings
#   script/bench_avonet_missingness.md     human-readable summary
#   script/bench_avonet_missingness.log    stdout/stderr
#
# Run with
#   cd pigauto && /usr/local/bin/Rscript script/bench_avonet_missingness.R

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto",
    quiet = TRUE
  )
})

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "bench_avonet_missingness.rds")
out_md  <- file.path(here, "script", "bench_avonet_missingness.md")

script_start <- proc.time()[["elapsed"]]

log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start), ..., "\n",
      sep = "")
  flush.console()
}

time_stage <- function(label, fun) {
  gc(reset = TRUE, full = TRUE)
  t0 <- proc.time()[["elapsed"]]
  value <- fun()
  wall <- proc.time()[["elapsed"]] - t0
  gc_after <- gc(full = TRUE)
  mb <- tryCatch(as.numeric(gc_after["Vcells", 7L]),
                 error = function(e) NA_real_)
  log_line(sprintf("[%-25s] wall = %7.1fs  maxR = %6.0f MB", label, wall, mb))
  list(value = value, wall = wall, mb = mb)
}

# -------------------------------------------------------------------------
# Load bundled data
# -------------------------------------------------------------------------

log_line("Loading avonet_full + tree_full from installed package")

# These are already xz-compressed in the pigauto package; devtools::load_all
# makes them available via data().
e <- new.env(parent = emptyenv())
utils::data("avonet_full", package = "pigauto", envir = e)
utils::data("tree_full",   package = "pigauto", envir = e)
avonet_full <- e$avonet_full
tree_full   <- e$tree_full

# preprocess_traits() expects rownames to be species keys matching tip labels.
df <- avonet_full
rownames(df) <- df$Species_Key
df$Species_Key <- NULL
stopifnot(all(rownames(df) == tree_full$tip.label))

log_line(sprintf("Aligned dataset: %d species x %d traits",
                 nrow(df), ncol(df)))

# -------------------------------------------------------------------------
# Preprocess + graph ONCE (shared across all missingness levels)
# -------------------------------------------------------------------------

pp <- time_stage("preprocess", function() {
  preprocess_traits(df, tree_full, log_transform = TRUE)
})
pd <- pp$value

gr <- time_stage("graph", function() {
  build_phylo_graph(tree_full, k_eigen = "auto")
})
graph <- gr$value

# Keep a pristine copy of graph$D so we can restore it at the start of
# each iteration after Fix B drops it before pigauto training.
D_cache <- graph$D

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

# Mean / mode imputation in latent space. Computes column means of the
# TRAINING-only copy of pd$X_scaled (held-out cells set to NA first), then
# predicts that column mean for every cell. For z-scored continuous /
# ordinal columns the prediction collapses to ~0; for one-hot categorical
# columns it collapses to the class frequency vector (argmax = mode).
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

# Force predict() onto CPU by monkey-patching get_device() for one call.
# (Same pattern as script/validate_avonet_full.R: at n ~ 10k the MPS
# attention tensor during predict can OOM on a 64 GB unified-memory
# Mac. Training still happens on MPS via the normal fit_pigauto() path.)
predict_on_cpu <- function(fit) {
  dev <- torch::torch_device("cpu")
  orig_get_device <- get("get_device", envir = asNamespace("pigauto"))
  on.exit(
    assignInNamespace("get_device", orig_get_device, ns = "pigauto"),
    add = TRUE
  )
  assignInNamespace("get_device", function() dev, ns = "pigauto")
  stats::predict(fit, return_se = TRUE, n_imputations = 1L)
}

# Uniform result-row shape for all three methods, so rbind works cleanly.
tag_rows <- function(ev_df, method, frac, wall_sec) {
  if (is.null(ev_df) || nrow(ev_df) == 0L) return(NULL)
  ev_df$method       <- method
  ev_df$missing_frac <- frac
  ev_df$wall_sec     <- wall_sec
  ev_df
}

# -------------------------------------------------------------------------
# Sweep
# -------------------------------------------------------------------------

results <- list()          # list of per-cell data.frames; rbound at end
stages  <- list()          # per-cell wall-clock / heap logs

missingness_levels <- c(0.20, 0.50, 0.80)
base_seed <- 2026L

for (frac in missingness_levels) {
  frac_tag <- sprintf("%02dpct", as.integer(round(100 * frac)))
  log_line(sprintf("=== missing_frac = %.2f ===", frac))

  # Restore cached cophenetic matrix for this iteration's baseline call
  # (we drop it before pigauto training to avoid the per-epoch regression
  # noted in Fix B, then restore it at the next iteration's top).
  if (is.null(graph$D)) graph$D <- D_cache

  # ---- Splits ---------------------------------------------------------
  sp <- time_stage(sprintf("splits_%s", frac_tag), function() {
    make_missing_splits(pd$X_scaled, missing_frac = frac,
                        val_frac = 0.5, seed = base_seed,
                        trait_map = pd$trait_map)
  })
  splits <- sp$value

  # ---- Method 1: mean / mode ------------------------------------------
  mm <- time_stage(sprintf("mean_%s", frac_tag), function() {
    mean_mode_impute(pd$X_scaled, splits)
  })
  mean_pred <- mm$value
  ev_mean <- evaluate_imputation(mean_pred, pd$X_scaled, splits,
                                 trait_map = pd$trait_map)
  results[[length(results) + 1L]] <-
    tag_rows(ev_mean, "mean", frac, mm$wall)

  # ---- Method 2: BM baseline ------------------------------------------
  bl <- time_stage(sprintf("baseline_%s", frac_tag), function() {
    fit_baseline(pd, tree_full, splits = splits,
                 model = "BM", graph = graph)
  })
  baseline <- bl$value
  ev_bm <- evaluate_imputation(baseline$mu, pd$X_scaled, splits,
                               trait_map = pd$trait_map)
  results[[length(results) + 1L]] <-
    tag_rows(ev_bm, "BM", frac, bl$wall)

  stages[[frac_tag]] <- list(
    splits   = list(wall = sp$wall, mb = sp$mb),
    mean     = list(wall = mm$wall, mb = mm$mb),
    baseline = list(wall = bl$wall, mb = bl$mb)
  )
  saveRDS(list(results = do.call(rbind, results), stages = stages), out_rds)

  # Drop graph$D before pigauto training (Fix B regression note)
  graph$D <- NULL
  invisible(gc(full = TRUE, verbose = FALSE))

  # ---- Method 3: pigauto full pipeline --------------------------------
  # Hyperparameters copied verbatim from script/validate_avonet_full.R
  # for comparability with the scaling benchmark.
  fr <- time_stage(sprintf("pigauto_train_%s", frac_tag), function() {
    fit_pigauto(
      data            = pd,
      tree            = tree_full,
      splits          = splits,
      graph           = graph,
      baseline        = baseline,
      hidden_dim      = 64L,
      dropout         = 0.10,
      lr              = 3e-3,
      epochs          = 500L,
      corruption_rate = 0.55,
      lambda_shrink   = 0.03,
      eval_every      = 25L,
      patience        = 4L,
      verbose         = TRUE,
      seed            = 1L
    )
  })
  fit <- fr$value

  pr <- time_stage(sprintf("pigauto_pred_%s", frac_tag), function() {
    predict_on_cpu(fit)
  })
  pred <- pr$value

  ev_pg <- evaluate_imputation(pred, pd$X_scaled, splits)
  results[[length(results) + 1L]] <-
    tag_rows(ev_pg, "pigauto", frac, fr$wall + pr$wall)

  stages[[frac_tag]]$pigauto_train <- list(wall = fr$wall, mb = fr$mb)
  stages[[frac_tag]]$pigauto_pred  <- list(wall = pr$wall, mb = pr$mb)

  saveRDS(list(results = do.call(rbind, results), stages = stages), out_rds)

  rm(baseline, fit, pred, mean_pred, splits)
  invisible(gc(full = TRUE, verbose = FALSE))
}

# -------------------------------------------------------------------------
# Finalise + markdown summary
# -------------------------------------------------------------------------

total_wall <- proc.time()[["elapsed"]] - script_start
all_results <- do.call(rbind, results)
rownames(all_results) <- NULL

# Drop the cached D reference so the saved RDS stays small.
rm(D_cache, graph); invisible(gc(full = TRUE, verbose = FALSE))

saveRDS(list(
  results    = all_results,
  stages     = stages,
  total_wall = total_wall,
  n_species  = nrow(df),
  n_traits   = ncol(df),
  trait_names = colnames(df),
  missingness_levels = missingness_levels,
  seed       = base_seed,
  commit     = tryCatch(system("git rev-parse HEAD", intern = TRUE),
                        error = function(e) "unknown")
), out_rds)

log_line(sprintf("Total wall: %.1fs (%.1f min)", total_wall, total_wall / 60))

# Human-readable markdown: a compact test-set-only table per missingness.
machine <- tryCatch(
  sprintf("%s %s (%s), R %s",
          Sys.info()[["sysname"]], Sys.info()[["release"]],
          Sys.info()[["machine"]],
          paste(R.version$major, R.version$minor, sep = ".")),
  error = function(e) "machine info unavailable"
)

md <- c(
  "# AVONET full-scale missingness sweep",
  "",
  sprintf("Run on: %s", format(Sys.time())),
  sprintf("Machine: %s", machine),
  sprintf("Species: %d, traits: %d", nrow(df), ncol(df)),
  sprintf("Total wall: %.1f min", total_wall / 60),
  "",
  "## Methods",
  "",
  "- **mean**: column mean / mode imputation (no phylogeny).",
  "- **BM**: Brownian-motion baseline (Rphylopars + phylogenetic label",
  "  propagation for discrete traits).",
  "- **pigauto**: full pipeline (BM baseline + calibrated GNN delta).",
  "",
  "## Test-set metrics (per trait)",
  ""
)

test_df <- all_results[all_results$split == "test", ]
show_cols <- c("method", "trait", "type", "n",
               "rmse", "pearson_r", "spearman_rho", "accuracy")
for (frac in missingness_levels) {
  md <- c(md, sprintf("### missing_frac = %.2f", frac), "")
  sub <- test_df[test_df$missing_frac == frac, show_cols, drop = FALSE]
  sub <- sub[order(sub$trait, sub$method), ]
  md <- c(md, "```", capture.output(print(sub, row.names = FALSE)), "```", "")
}

writeLines(md, out_md)
log_line(sprintf("Wrote %s", out_md))
log_line("done")
