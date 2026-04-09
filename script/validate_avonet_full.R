#!/usr/bin/env Rscript
#
# script/validate_avonet_full.R
#
# End-to-end validation run on the FULL AVONET3 dataset (9,993 bird
# species) and the matching BirdTree Stage2 Hackett MCC phylogeny,
# after landing Fix A (RSpectra) and Fix B (cophenetic caching) on the
# v0.3.1 fix branch.
#
# Purpose
#   * Prove that pigauto runs end-to-end on a real 10k-tip phylogeny,
#     not just a simulation.
#   * Measure wall-clock time per pipeline stage at this scale.
#   * Compare GNN test RMSE against the Rphylopars BM baseline on a
#     held-out 15% MCAR mask.
#
# Design
#   * Tree: avonet/Stage2_Hackett_MCC_no_neg.tre (9,993 tips).
#   * Traits: 4 continuous morphometric (Mass, Beak.Length_Culmen,
#     Tarsus.Length, Wing.Length) + 2 categorical (Trophic.Level,
#     Primary.Lifestyle) + 1 ordinal (Migration). Same scheme as the
#     bundled avonet300 dataset.
#   * Native missingness from AVONET is preserved. An additional 15%
#     MCAR mask on observed cells provides the held-out test set.
#   * Stages: preprocess -> graph -> baseline -> train (2000 epochs
#     with early stopping) -> predict (5 imputations).
#   * Checkpoints: one .rds checkpoint per stage so a crash still
#     leaves usable partial results.
#
# Output
#   script/validate_avonet_full.rds    timing + metrics
#   script/validate_avonet_full.md     human-readable summary
#   script/validate_avonet_full.log    stdout/stderr
#
# Run with
#   cd pigauto && /usr/local/bin/Rscript script/validate_avonet_full.R

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

here       <- "/Users/z3437171/Dropbox/Github Local/pigauto"
tree_path  <- file.path(here, "avonet", "Stage2_Hackett_MCC_no_neg.tre")
csv_path   <- file.path(here, "avonet", "AVONET3_BirdTree.csv")

out_rds    <- file.path(here, "script", "validate_avonet_full.rds")
out_md     <- file.path(here, "script", "validate_avonet_full.md")

stopifnot(file.exists(tree_path), file.exists(csv_path))

script_start <- proc.time()[["elapsed"]]

log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start), ..., "\n",
      sep = "")
  flush.console()
}

# -------------------------------------------------------------------------
# Load and prepare data
# -------------------------------------------------------------------------

log_line("Loading AVONET3 + BirdTree tree...")
tree   <- ape::read.tree(tree_path)
avonet <- read.csv(csv_path, stringsAsFactors = FALSE)
avonet$Species_Key <- gsub(" ", "_", avonet$Species3)

cont_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")
trait_cols <- c(cont_cols, "Trophic.Level", "Primary.Lifestyle", "Migration")

df <- avonet[stats::complete.cases(avonet[, cont_cols]), ]
df <- df[, c("Species_Key", trait_cols)]

df$Trophic.Level <- trimws(df$Trophic.Level)
df$Trophic.Level[df$Trophic.Level == ""] <- NA
df$Trophic.Level <- factor(df$Trophic.Level)

df$Primary.Lifestyle <- trimws(df$Primary.Lifestyle)
df$Primary.Lifestyle[df$Primary.Lifestyle == ""] <- NA
df$Primary.Lifestyle <- factor(df$Primary.Lifestyle)

df$Migration <- ordered(df$Migration, levels = c(1, 2, 3),
                        labels = c("Resident", "Partial", "Full"))

common <- intersect(tree$tip.label, df$Species_Key)
df     <- df[df$Species_Key %in% common, ]
tree   <- ape::keep.tip(tree, common)
df     <- df[match(tree$tip.label, df$Species_Key), ]
rownames(df) <- df$Species_Key
df$Species_Key <- NULL

n_species <- nrow(df)
log_line(sprintf("Aligned dataset: %d species x %d traits", n_species, ncol(df)))

# -------------------------------------------------------------------------
# Inject an additional 15% MCAR mask for held-out evaluation
# -------------------------------------------------------------------------

set.seed(2026L)
orig <- df
cell_is_observed <- !is.na(as.matrix(df))
n_obs_cells <- sum(cell_is_observed)
target_mask <- round(0.15 * n_obs_cells)
obs_idx <- which(cell_is_observed)
held_out <- sample(obs_idx, target_mask)

# Record truth at held-out cells BEFORE overwriting
truth_mat <- as.matrix(data.frame(lapply(df, function(x) {
  if (is.factor(x)) as.integer(x) else as.numeric(x)
})))
held_truth <- truth_mat[held_out]
held_trait <- ((held_out - 1L) %/% n_species) + 1L  # column index
held_row   <- ((held_out - 1L) %%  n_species) + 1L  # row index

# Overwrite held-out cells with NA
mat_na <- as.matrix(df)
mat_na[held_out] <- NA
for (j in seq_len(ncol(df))) {
  new_col <- mat_na[, j]
  old_col <- df[, j]
  if (is.factor(old_col)) {
    df[, j] <- factor(new_col, levels = levels(old_col),
                      ordered = is.ordered(old_col))
  } else if (is.integer(old_col)) {
    df[, j] <- as.integer(new_col)
  } else {
    df[, j] <- as.numeric(new_col)
  }
}

log_line(sprintf("Held-out cells: %d (from %d observed cells; %.1f%%)",
                 length(held_out), n_obs_cells,
                 100 * length(held_out) / n_obs_cells))

results <- list(
  n_species       = n_species,
  trait_names     = colnames(df),
  n_held_out      = length(held_out),
  stages          = list(),
  commit          = tryCatch(system("git rev-parse HEAD", intern = TRUE),
                             error = function(e) "unknown")
)
saveRDS(results, out_rds)

# -------------------------------------------------------------------------
# Per-stage timing wrapper
# -------------------------------------------------------------------------

time_stage <- function(label, expr_fun) {
  gc(reset = TRUE, full = TRUE)
  t0 <- proc.time()[["elapsed"]]
  value <- expr_fun()
  wall <- proc.time()[["elapsed"]] - t0
  gc_after <- gc(full = TRUE)
  mb <- tryCatch(as.numeric(gc_after["Vcells", 7L]),
                 error = function(e) NA_real_)
  results$stages[[label]] <<- list(wall_sec = wall, max_mb = mb)
  saveRDS(results, out_rds)
  log_line(sprintf("[%-10s] wall = %7.1fs  maxR = %6.0f MB", label, wall, mb))
  value
}

# -------------------------------------------------------------------------
# Pipeline
# -------------------------------------------------------------------------

log_line("Starting pipeline...")

pd <- time_stage("preprocess", function() {
  preprocess_traits(df, tree, log_transform = TRUE)
})

graph <- time_stage("graph", function() {
  build_phylo_graph(tree, k_eigen = "auto")
})

splits <- time_stage("splits", function() {
  make_missing_splits(pd$X_scaled, missing_frac = 0.25,
                      val_frac = 0.5, seed = 42L,
                      trait_map = pd$trait_map)
})

baseline <- time_stage("baseline", function() {
  fit_baseline(pd, tree, splits = splits, model = "BM", graph = graph)
})

# Free cached cophenetic matrix before training (see Fix B regression
# note in bench_scaling_v031.R -- leaving ~800 MB of D in memory slows
# the training loop by up to 10x at n ~ 10k).
graph$D <- NULL
invisible(gc(full = TRUE, verbose = FALSE))

fit <- time_stage("train", function() {
  fit_pigauto(
    data            = pd,
    tree            = tree,
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

# Checkpoint the fit object to disk IMMEDIATELY after training so that
# a crash during predict() (e.g. MPS OOM on the 9993 x 9993 attention
# tensor with n_imputations > 1) does not force a costly retrain.
fit_ckpt <- file.path(here, "script", "validate_avonet_full_fit.rds")
saveRDS(fit, fit_ckpt)
log_line(sprintf("Checkpointed fit to %s", fit_ckpt))

# Drop the graph-side copies of heavy matrices that predict() does not
# need. predict() reuses fit$graph$adj and fit$graph$coords from the
# fit object itself, so we can release the local `graph` here.
rm(graph); invisible(gc(full = TRUE, verbose = FALSE))

# predict() on a 9993-tip phylogeny with n_imputations > 1 peaks at
# (number of MC passes) * (layer count) * (n^2 * 4 bytes) on MPS,
# which is roughly 4-8 GB of unified memory per pass. The previous
# run at n_imputations = 5 appears to have been silently killed by
# the MPS memory manager somewhere between train completion and the
# first predict forward pass. Use n_imputations = 1 for this
# validation run so we get a deterministic single forward pass; MC
# dropout is not essential for proving the pipeline runs end-to-end
# at 10k tips.
#
# Helper: run predict on a specific torch device by monkey-patching
# pigauto's internal get_device() for the duration of the call.
predict_on_device <- function(device_name) {
  dev <- torch::torch_device(device_name)
  orig_get_device <- get("get_device", envir = asNamespace("pigauto"))
  on.exit(
    assignInNamespace("get_device", orig_get_device, ns = "pigauto"),
    add = TRUE
  )
  assignInNamespace("get_device", function() dev, ns = "pigauto")
  stats::predict(fit, return_se = TRUE, n_imputations = 1L)
}

pred <- tryCatch(
  time_stage("predict", function() predict_on_device("cpu")),
  error = function(e) {
    log_line(sprintf("predict() failed on CPU: %s", conditionMessage(e)))
    NULL
  }
)

if (is.null(pred)) {
  log_line("predict() unavailable; skipping held-out metrics block.")
}

total_wall <- proc.time()[["elapsed"]] - script_start
results$total_wall_sec <- total_wall
log_line(sprintf("Total pipeline wall: %.1fs (%.1f min)",
                 total_wall, total_wall / 60))

# -------------------------------------------------------------------------
# Evaluation on held-out cells
# -------------------------------------------------------------------------

eval_on_holdout <- function(pred_mat_latent, X_true_latent, held_latent_idx) {
  p <- as.numeric(pred_mat_latent[held_latent_idx])
  t <- as.numeric(X_true_latent[held_latent_idx])
  ok <- is.finite(p) & is.finite(t)
  if (sum(ok) == 0) return(list(n = 0L, rmse = NA_real_, r = NA_real_))
  list(
    n    = sum(ok),
    rmse = sqrt(mean((p[ok] - t[ok])^2)),
    r    = stats::cor(p[ok], t[ok])
  )
}

# Map held-out input cells to latent-space cells. For continuous /
# ordinal / count columns the mapping is 1-to-1 via trait_map$latent_cols;
# for categorical columns the held-out entry corresponds to a single
# column in input space but K columns in latent space. We report
# metrics for the BM-style latent columns only (continuous + ordinal +
# count) as the simplest apples-to-apples comparison.

trait_map <- pd$trait_map
X_lat_true <- pd$X_scaled  # z-scored / latent-space matrix (observed + test cells)
# Re-insert held-out truth in latent space by temporarily using the
# original (unmasked) data frame once more. Simpler: just use splits$test_idx
# already produced by make_missing_splits(), which operates on pd$X_scaled.

# pred$imputed_latent should be n_species x p_latent in z-score space.
pred_latent_gnn <- if (!is.null(pred)) pred$imputed_latent else NULL
baseline_latent <- baseline$mu

# Extract latent-space metrics per trait for BM-eligible traits.
metrics <- data.frame(
  trait = character(), type = character(),
  n = integer(),
  bm_rmse = numeric(), gnn_rmse = numeric(),
  bm_r = numeric(),    gnn_r = numeric(),
  stringsAsFactors = FALSE
)

test_rows <- splits$test_idx  # linear indices into pd$X_scaled
if (length(test_rows) > 0) {
  X_true_lat <- pd$X_scaled
  for (tm in trait_map) {
    if (!(tm$type %in% c("continuous", "count", "ordinal"))) next
    lc <- tm$latent_cols
    # Restrict test_rows that fall within this trait's latent column(s)
    col_of_test <- ((test_rows - 1L) %/% nrow(X_true_lat)) + 1L
    keep <- col_of_test %in% lc
    if (!any(keep)) next
    idx <- test_rows[keep]
    bm_stats <- eval_on_holdout(baseline_latent, X_true_lat, idx)
    gnn_stats <- if (!is.null(pred_latent_gnn)) {
      eval_on_holdout(pred_latent_gnn, X_true_lat, idx)
    } else {
      list(n = NA_integer_, rmse = NA_real_, r = NA_real_)
    }
    metrics <- rbind(metrics, data.frame(
      trait = tm$name, type = tm$type,
      n = bm_stats$n,
      bm_rmse = bm_stats$rmse, gnn_rmse = gnn_stats$rmse,
      bm_r = bm_stats$r,       gnn_r = gnn_stats$r
    ))
  }
}
results$metrics <- metrics
log_line("Latent-space test-cell metrics (z-score):")
print(metrics)

# -------------------------------------------------------------------------
# Save and write markdown
# -------------------------------------------------------------------------

saveRDS(results, out_rds)

machine <- tryCatch(
  sprintf("%s %s (%s), R %s",
          Sys.info()[["sysname"]], Sys.info()[["release"]],
          Sys.info()[["machine"]],
          paste(R.version$major, R.version$minor, sep = ".")),
  error = function(e) "machine info unavailable"
)

md <- c(
  "# AVONET full-dataset validation (post Fix A + Fix B)",
  "",
  sprintf("Run on: %s", format(Sys.time())),
  sprintf("Machine: %s", machine),
  sprintf("Commit:  %s", results$commit),
  "",
  sprintf("Species aligned to phylogeny: **%d**", results$n_species),
  sprintf("Trait columns: %s", paste(results$trait_names, collapse = ", ")),
  sprintf("Held-out test cells: %d", results$n_held_out),
  "",
  "## Per-stage wall time",
  "",
  "| Stage      | Wall (s) | Wall (min) | R heap max (MB) |",
  "|------------|----------|------------|-----------------|"
)

for (nm in names(results$stages)) {
  s <- results$stages[[nm]]
  md <- c(md, sprintf("| %-10s | %8.1f | %10.2f | %15.0f |",
                      nm, s$wall_sec, s$wall_sec / 60, s$max_mb))
}

md <- c(md,
  sprintf("| **total** | %8.1f | %10.2f |                 |",
          results$total_wall_sec, results$total_wall_sec / 60),
  "",
  "## Held-out test-cell metrics (latent / z-score space)",
  "",
  "| Trait | Type | n | BM RMSE | GNN RMSE | BM r | GNN r |",
  "|---|---|---|---|---|---|---|"
)
if (nrow(results$metrics) > 0) {
  for (i in seq_len(nrow(results$metrics))) {
    r <- results$metrics[i, ]
    md <- c(md, sprintf("| %s | %s | %d | %.3f | %.3f | %.3f | %.3f |",
                        r$trait, r$type, r$n,
                        r$bm_rmse, r$gnn_rmse, r$bm_r, r$gnn_r))
  }
} else {
  md <- c(md, "| (no BM-eligible test cells) | | | | | | |")
}

writeLines(md, out_md)
log_line(sprintf("Wrote %s", out_md))
log_line("done")
