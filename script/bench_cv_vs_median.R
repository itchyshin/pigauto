#!/usr/bin/env Rscript
#
# script/bench_cv_vs_median.R
#
# Focused 3-way comparison: gate_method in
#   {"single_split", "median_splits", "cv_folds"}
# on synthetic discrete + continuous traits with controlled signal.
# Goal: measure whether cv_folds closes more val->test drift than
# median_splits on the regime where median_splits already closed
# ~50% of the discrete bench regression
# (useful/MEMO_2026-04-29_discrete_bench_reruns.md).
#
# DGP
#   * Tree: ape::rtree(300)
#   * Continuous trait: BM with phylo signal lambda = 1
#   * Binary trait: threshold-binarised BM (50/50 split)
#   * Categorical trait (K=4): independent BM per class, argmax decode
#   * 30% MCAR mask
#
# For each gate_method:
#   * Fit pigauto on training+val (val_frac = 0.30 of mask)
#   * Score on test set (the other 0.70 of mask)
#   * Per-trait metrics:
#       continuous: RMSE, Pearson r
#       binary / categorical: accuracy
#
# Replicates: 5 seeds.
#
# Output:
#   script/bench_cv_vs_median.rds
#   script/bench_cv_vs_median.md
#
# Wall: ~10-15 min per scenario; 3 scenarios; 5 reps -> ~3-5 min total
# at moderate parallelism (set MC_CORES env var; default 4).

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  here_path <- "/Users/z3437171/Dropbox/Github Local/pigauto"
  devtools::load_all(here_path, quiet = TRUE)
})

out_rds <- file.path(here_path, "script", "bench_cv_vs_median.rds")
out_md  <- file.path(here_path, "script", "bench_cv_vs_median.md")

script_start <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start),
      ..., "\n", sep = "")
  flush.console()
}

simulate_dataset <- function(n = 300L, signal = 1.0, seed = 1L) {
  set.seed(seed)
  tree <- ape::rtree(n)
  # Continuous BM
  cont <- ape::rTraitCont(tree, sigma = signal)
  # Binary: threshold BM
  bm_b <- ape::rTraitCont(tree, sigma = signal)
  bin <- factor(ifelse(bm_b > stats::median(bm_b), "1", "0"),
                levels = c("0", "1"))
  # Categorical K=4: argmax of 4 independent BMs
  K <- 4L
  bm_c <- vapply(seq_len(K), function(k)
                   ape::rTraitCont(tree, sigma = signal),
                 numeric(n))
  cat3 <- factor(LETTERS[max.col(bm_c, ties.method = "first")],
                 levels = LETTERS[seq_len(K)])
  df <- data.frame(cont = cont, bin = bin, cat3 = cat3,
                   row.names = tree$tip.label)
  list(df = df, tree = tree)
}

apply_mask <- function(df, mask_frac = 0.30, seed = 1L) {
  set.seed(seed)
  out <- df
  mask <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df),
                 dimnames = list(NULL, names(df)))
  for (v in names(df)) {
    n_mask <- round(mask_frac * nrow(df))
    idx    <- sample.int(nrow(df), n_mask)
    mask[idx, v] <- TRUE
    if (is.factor(df[[v]])) {
      out[[v]][idx] <- NA
    } else {
      out[[v]][idx] <- NA
    }
  }
  list(masked = out, mask = mask)
}

eval_one_method <- function(gate_method, df_truth, df_masked, mask, tree,
                              seed, epochs) {
  res <- pigauto::impute(df_masked, tree,
                           safety_floor = TRUE,
                           epochs = epochs, n_imputations = 1L,
                           gate_method = gate_method,
                           gate_cv_folds = 5L, gate_splits_B = 31L,
                           verbose = FALSE, seed = seed)
  comp <- res$completed
  out <- list()
  for (v in names(df_truth)) {
    if (!any(mask[, v])) next
    truth_v <- df_truth[[v]][mask[, v]]
    pred_v  <- comp[[v]][mask[, v]]
    if (is.factor(df_truth[[v]])) {
      acc <- mean(as.character(pred_v) == as.character(truth_v),
                  na.rm = TRUE)
      out[[v]] <- data.frame(trait = v, type = "discrete",
                              metric = "accuracy", value = acc,
                              stringsAsFactors = FALSE)
    } else {
      diffs <- as.numeric(truth_v) - as.numeric(pred_v)
      rmse  <- sqrt(mean(diffs^2, na.rm = TRUE))
      out[[v]] <- data.frame(trait = v, type = "continuous",
                              metric = "rmse", value = rmse,
                              stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

run_one_cell <- function(rep_id, gate_method, signal = 1.0, n = 300L,
                          mask_frac = 0.30, epochs = 60L) {
  sim <- simulate_dataset(n = n, signal = signal, seed = rep_id * 100L + 1L)
  msk <- apply_mask(sim$df, mask_frac = mask_frac,
                     seed = rep_id * 100L + 2L)
  ev <- tryCatch(
    eval_one_method(gate_method, sim$df, msk$masked, msk$mask,
                     sim$tree, seed = rep_id * 100L + 3L, epochs = epochs),
    error = function(e) {
      data.frame(trait = NA_character_, type = NA_character_,
                  metric = NA_character_, value = NA_real_,
                  error = conditionMessage(e),
                  stringsAsFactors = FALSE)
    })
  ev$gate_method <- gate_method
  ev$rep_id      <- rep_id
  ev$signal      <- signal
  ev
}

# -------------------------------------------------------------------------
# Sweep
# -------------------------------------------------------------------------

methods <- c("single_split", "median_splits", "cv_folds")
n_reps  <- 5L

log_line("Starting bench_cv_vs_median: ", n_reps, " reps x ", length(methods),
          " methods = ", n_reps * length(methods), " cells.")

cells <- expand.grid(rep_id = seq_len(n_reps),
                     gate_method = methods,
                     stringsAsFactors = FALSE)
results <- list()
for (i in seq_len(nrow(cells))) {
  log_line(sprintf("[%d/%d] rep=%d method=%s",
                    i, nrow(cells), cells$rep_id[i], cells$gate_method[i]))
  results[[i]] <- run_one_cell(rep_id = cells$rep_id[i],
                                 gate_method = cells$gate_method[i])
}
all_rows <- do.call(rbind, results)

# -------------------------------------------------------------------------
# Summarise
# -------------------------------------------------------------------------

log_line("Summarising across reps...")
agg <- stats::aggregate(value ~ gate_method + trait + metric,
                         data = all_rows,
                         FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                              sd   = stats::sd(x, na.rm = TRUE)))
agg_df <- data.frame(
  gate_method = agg$gate_method,
  trait       = agg$trait,
  metric      = agg$metric,
  mean        = as.numeric(agg$value[, "mean"]),
  sd          = as.numeric(agg$value[, "sd"]),
  stringsAsFactors = FALSE
)
agg_df <- agg_df[order(agg_df$trait, agg_df$gate_method), ]

saveRDS(list(per_rep = all_rows, summary = agg_df, n_reps = n_reps,
              methods = methods),
         out_rds)

# -------------------------------------------------------------------------
# Markdown
# -------------------------------------------------------------------------

md <- c(
  "# CV-folds vs median_splits vs single_split bench",
  "",
  sprintf("Run on: %s", format(Sys.time())),
  sprintf("n=%d species, %d reps, %d gate methods, mask_frac=0.30",
          300L, n_reps, length(methods)),
  "",
  "## Per-trait mean +/- SD across reps",
  "",
  "```",
  capture.output(print(agg_df, row.names = FALSE)),
  "```",
  "",
  "## Comparison: cv_folds vs median_splits per trait",
  ""
)
for (tr in unique(agg_df$trait)) {
  if (is.na(tr)) next
  rows_t <- agg_df[agg_df$trait == tr, ]
  if (nrow(rows_t) < 2L) next
  metric <- rows_t$metric[1]
  ss <- rows_t[rows_t$gate_method == "single_split", "mean"]
  ms <- rows_t[rows_t$gate_method == "median_splits", "mean"]
  cv <- rows_t[rows_t$gate_method == "cv_folds", "mean"]
  md <- c(md,
          sprintf("- **%s** (%s): single=%.4g, median=%.4g, cv=%.4g",
                  tr, metric, ss, ms, cv))
}
writeLines(md, out_md)

log_line("DONE -- ", out_rds)
log_line("DONE -- ", out_md)
