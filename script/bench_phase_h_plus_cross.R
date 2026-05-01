#!/usr/bin/env Rscript
#
# script/bench_phase_h_plus_cross.R
#
# Phase H+: cross-K ordinal mode-vs-median pooling bench on simulated
# trees, to gather evidence for whether to flip pool_method default to
# "mode" for ordinal traits in v0.9.2.
#
# Pre-registered framing (from useful/MEMO_2026-05-01_phase_h_results.md):
#   Phase H showed mode > median on AVONET Migration K=3 at N_IMP=20
#   (+6.6 pp).  Single-trait single-dataset evidence.  Phase H+ checks
#   whether mode wins, ties, or loses for K in {3, 5, 7, 10} on
#   simulated phylogenetic ordinal data.
#
# Pre-registered acceptance criterion to flip default to "mode":
#   - Mode acc >= median acc on at least 80% of (K, lambda, rep) cells
#     at N_IMP = 20.
#   - Mode never regresses by > 2 pp on any (K, lambda) combo.
#
# Bench design:
#   K in {3, 5, 7}     (drop K=10 for compute; K=3 is the AVONET regime)
#   lambda in {0.6}    (moderate phylo signal; K-sweep at fixed signal)
#   reps = 3 per cell
#   pool_method in {"median", "mode"}
#   n_imputations = 20
#   n_species = 200 per simulated tree
#   = 3 K x 1 lambda x 3 reps x 2 pools x 1 N_IMP = 18 cells
#   Wall: ~3 min per cell, ~54 min total.
#
# Output:
#   script/bench_phase_h_plus_cross.{rds,md}

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH", unset = "")
  if (nzchar(pkg_path) && dir.exists(pkg_path) &&
      file.exists(file.path(pkg_path, "DESCRIPTION"))) {
    devtools::load_all(pkg_path, quiet = TRUE)
  } else {
    library(pigauto)
  }
})

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "bench_phase_h_plus_cross.rds")
out_md  <- file.path(here, "script", "bench_phase_h_plus_cross.md")

K_LIST    <- c(3L, 5L, 7L)
LAMBDA    <- 0.6
N_REPS    <- 3L
N_IMP     <- 20L
N_SPECIES <- 200L
MISS_FRAC <- 0.30
POOLS     <- c("median", "mode")

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - t0),
      ..., "\n", sep = "")
  flush.console()
}

# ---- Simulate one ordinal trait via thresholded BM -------------------------

sim_ordinal <- function(n, K, lambda, seed) {
  set.seed(seed)
  tree <- ape::rtree(n)
  z    <- ape::rTraitCont(tree, sigma = lambda)
  # Equal-prob class boundaries (so each class is roughly 1/K of the tips)
  q    <- stats::quantile(z, probs = seq_len(K - 1L) / K)
  cls  <- as.integer(cut(z, breaks = c(-Inf, q, Inf), labels = FALSE))
  levs <- paste0("c", seq_len(K))
  list(
    tree   = tree,
    factor = ordered(levs[cls], levels = levs),
    z      = z
  )
}

# ---- One bench cell --------------------------------------------------------

run_one <- function(K, lambda, rep_id) {
  seed <- as.integer(K * 1e4 + lambda * 1e3 + rep_id)
  log_line(sprintf("K=%d lambda=%.1f rep=%d: simulating ...",
                    K, lambda, rep_id))
  sim <- sim_ordinal(N_SPECIES, K, lambda, seed)

  df <- data.frame(o = sim$factor, row.names = sim$tree$tip.label)
  set.seed(seed + 1L)
  miss_idx <- sample(N_SPECIES, ceiling(MISS_FRAC * N_SPECIES))
  truth <- df$o
  df$o[miss_idx] <- NA

  out <- list()
  for (pm in POOLS) {
    log_line(sprintf("  K=%d rep=%d pool=%s: impute() N_IMP=%d ...",
                      K, rep_id, pm, N_IMP))
    res <- pigauto::impute(df, sim$tree,
                            epochs = 200L,
                            n_imputations = N_IMP,
                            pool_method = pm,
                            verbose = FALSE,
                            seed = seed + 2L)
    pred <- res$completed$o[miss_idx]
    acc  <- mean(as.character(pred) == as.character(truth[miss_idx]),
                 na.rm = TRUE)

    # Per-class-mode baseline (most frequent observed class)
    obs_table   <- sort(table(df$o), decreasing = TRUE)
    mode_class  <- names(obs_table)[1L]
    base_acc    <- mean(as.character(truth[miss_idx]) == mode_class,
                        na.rm = TRUE)

    out[[pm]] <- data.frame(
      K            = K,
      lambda       = lambda,
      rep          = rep_id,
      pool_method  = pm,
      pigauto_acc  = acc,
      baseline_acc = base_acc,
      lift_pp      = 100 * (acc - base_acc),
      n_test       = length(miss_idx),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

# ---- Sweep -----------------------------------------------------------------

results <- list()
i <- 0L
total <- length(K_LIST) * N_REPS
for (K in K_LIST) {
  for (rp in seq_len(N_REPS)) {
    i <- i + 1L
    log_line(sprintf("==== cell %d/%d: K=%d rep=%d ====", i, total, K, rp))
    results[[length(results) + 1L]] <- run_one(K, LAMBDA, rp)
  }
}
all_rows <- do.call(rbind, results)
saveRDS(all_rows, out_rds)

# ---- Summary -------------------------------------------------------------

# Per-cell mode > median?
cells <- unique(all_rows[, c("K", "rep")])
cell_mode_wins <- sapply(seq_len(nrow(cells)), function(i) {
  c_med <- subset(all_rows, K == cells$K[i] & rep == cells$rep[i] &
                  pool_method == "median")$pigauto_acc
  c_mod <- subset(all_rows, K == cells$K[i] & rep == cells$rep[i] &
                  pool_method == "mode")$pigauto_acc
  c_mod >= c_med
})
n_cells_mode_wins <- sum(cell_mode_wins, na.rm = TRUE)
n_cells_total     <- length(cell_mode_wins)
pct_wins          <- 100 * n_cells_mode_wins / n_cells_total

# Worst regression: where median > mode by the most
deltas <- sapply(seq_len(nrow(cells)), function(i) {
  c_med <- subset(all_rows, K == cells$K[i] & rep == cells$rep[i] &
                  pool_method == "median")$pigauto_acc
  c_mod <- subset(all_rows, K == cells$K[i] & rep == cells$rep[i] &
                  pool_method == "mode")$pigauto_acc
  100 * (c_mod - c_med)   # positive = mode wins
})
worst_regression <- if (length(deltas) > 0L) min(deltas) else NA_real_

# Per-K aggregate
agg <- stats::aggregate(
  cbind(pigauto_acc, lift_pp) ~ K + pool_method,
  data = all_rows,
  FUN = function(x) c(mean = mean(x), sd = stats::sd(x))
)

# Acceptance verdict
crit_a <- (pct_wins >= 80)
crit_b <- !is.na(worst_regression) && worst_regression >= -2
verdict <- if (crit_a && crit_b) {
  sprintf("**PHASE H+ PASSES**: mode wins on %d / %d cells (%.0f%%); worst per-cell regression %.1f pp (>= -2 pp).  Strong evidence to flip default to 'mode' in v0.9.2.",
          n_cells_mode_wins, n_cells_total, pct_wins, worst_regression)
} else if (crit_a && !crit_b) {
  sprintf("**PHASE H+ PARTIAL**: mode wins on %d / %d cells (%.0f%%) but worst per-cell regression %.1f pp violates the -2 pp guard.  Investigate that cell before flipping default.",
          n_cells_mode_wins, n_cells_total, pct_wins, worst_regression)
} else {
  sprintf("**PHASE H+ INCONCLUSIVE**: mode wins on only %d / %d cells (%.0f%%, target >= 80%%).  Default should remain 'median'.",
          n_cells_mode_wins, n_cells_total, pct_wins)
}

# ---- Markdown report -----------------------------------------------------

md <- c(
  "# Phase H+ smoke bench: cross-K ordinal pooling (median vs mode)",
  "",
  sprintf("Run: %s", format(Sys.time())),
  sprintf("n_species=%d, miss_frac=%.2f, N_IMP=%d, K=%s, lambda=%g, reps=%d",
          N_SPECIES, MISS_FRAC, N_IMP,
          paste(K_LIST, collapse = ","), LAMBDA, N_REPS),
  "",
  "## Per-cell results",
  "",
  "```",
  capture.output(print(all_rows[order(all_rows$K, all_rows$rep,
                                       all_rows$pool_method),
                                 c("K", "rep", "pool_method",
                                   "pigauto_acc", "baseline_acc", "lift_pp")],
                       row.names = FALSE)),
  "```",
  "",
  "## Aggregated by K and pool_method",
  "",
  "```",
  capture.output(print(agg, row.names = FALSE)),
  "```",
  "",
  sprintf("Mode wins (>= median acc) on %d / %d cells (%.0f%%).",
          n_cells_mode_wins, n_cells_total, pct_wins),
  sprintf("Worst per-cell mode-vs-median delta: %+.1f pp.",
          worst_regression),
  "",
  "## Verdict",
  "",
  verdict,
  ""
)
writeLines(md, out_md)
log_line("DONE -- ", out_rds)
log_line("DONE -- ", out_md)
