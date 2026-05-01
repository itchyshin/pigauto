#!/usr/bin/env Rscript
#
# script/bench_cv_default_flip.R
#
# Expanded cv_folds vs single_split comparison to inform whether to
# flip the gate_method default for v0.9.2.
#
# Sweeps 3 regimes (low/moderate/high phylo signal) x 3 trait types
# (continuous, binary, categorical) x 10 reps = 90 cells per method.
# Compares single_split (current default) vs cv_folds (proposed
# default).
#
# Decision rule (pre-registered):
# Flip default to cv_folds if:
#   - mean continuous RMSE lift >= 1 % across regimes
#   - mean discrete accuracy gap (cv - single) >= -0.5 pp
#     across regimes (cv_folds may not improve discrete, but
#     should not regress)
# Otherwise keep single_split as default and document cv_folds as
# recommended-but-opt-in.
#
# Output:
#   script/bench_cv_default_flip.{rds,md}
#
# Wall: ~30-45 min sequential.

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  here_path <- "/Users/z3437171/Dropbox/Github Local/pigauto"
  devtools::load_all(here_path, quiet = TRUE)
})

out_rds <- file.path(here_path, "script", "bench_cv_default_flip.rds")
out_md  <- file.path(here_path, "script", "bench_cv_default_flip.md")

script_start <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start),
      ..., "\n", sep = "")
  flush.console()
}

simulate_dataset <- function(n = 300L, lambda = 1.0, seed = 1L) {
  set.seed(seed)
  tree <- ape::rtree(n)
  cont <- ape::rTraitCont(tree, sigma = lambda)
  bm_b <- ape::rTraitCont(tree, sigma = lambda)
  bin  <- factor(ifelse(bm_b > stats::median(bm_b), "1", "0"),
                 levels = c("0", "1"))
  K <- 4L
  bm_c <- vapply(seq_len(K), function(k) ape::rTraitCont(tree, sigma = lambda),
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
    out[[v]][idx] <- if (is.factor(df[[v]])) NA else NA_real_
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

run_one_cell <- function(rep_id, gate_method, lambda,
                          n = 300L, mask_frac = 0.30, epochs = 60L) {
  sim <- simulate_dataset(n = n, lambda = lambda,
                           seed = rep_id * 1000L + round(lambda * 10L) + 1L)
  msk <- apply_mask(sim$df, mask_frac = mask_frac,
                     seed = rep_id * 1000L + round(lambda * 10L) + 2L)
  ev <- tryCatch(
    eval_one_method(gate_method, sim$df, msk$masked, msk$mask,
                     sim$tree,
                     seed = rep_id * 1000L + round(lambda * 10L) + 3L,
                     epochs = epochs),
    error = function(e) {
      data.frame(trait = NA_character_, type = NA_character_,
                  metric = NA_character_, value = NA_real_,
                  stringsAsFactors = FALSE)
    })
  ev$gate_method <- gate_method
  ev$rep_id      <- rep_id
  ev$lambda      <- lambda
  ev
}

# Sweep ----
methods <- c("single_split", "cv_folds")
lambdas <- c(0.4, 0.7, 1.0)
n_reps  <- 10L

cells <- expand.grid(rep_id = seq_len(n_reps),
                     gate_method = methods,
                     lambda = lambdas,
                     stringsAsFactors = FALSE)

log_line("bench_cv_default_flip starting: ", nrow(cells),
          " cells (", n_reps, " reps x ", length(methods),
          " methods x ", length(lambdas), " lambdas)")

results <- list()
for (i in seq_len(nrow(cells))) {
  log_line(sprintf("[%d/%d] rep=%d method=%s lambda=%.1f",
                    i, nrow(cells), cells$rep_id[i],
                    cells$gate_method[i], cells$lambda[i]))
  results[[i]] <- run_one_cell(rep_id = cells$rep_id[i],
                                 gate_method = cells$gate_method[i],
                                 lambda = cells$lambda[i])
}
all_rows <- do.call(rbind, results)

# Summarise ----
log_line("Aggregating ...")
agg <- stats::aggregate(value ~ gate_method + trait + metric + lambda,
                         data = all_rows,
                         FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                              sd   = stats::sd(x, na.rm = TRUE),
                                              n    = sum(!is.na(x))))
agg_df <- data.frame(
  gate_method = agg$gate_method,
  trait       = agg$trait,
  metric      = agg$metric,
  lambda      = agg$lambda,
  mean        = as.numeric(agg$value[, "mean"]),
  sd          = as.numeric(agg$value[, "sd"]),
  n           = as.numeric(agg$value[, "n"]),
  stringsAsFactors = FALSE
)
agg_df <- agg_df[order(agg_df$trait, agg_df$lambda, agg_df$gate_method), ]

saveRDS(list(per_rep = all_rows, summary = agg_df,
              n_reps = n_reps, lambdas = lambdas, methods = methods),
         out_rds)

# Decision rule ----
# Compute paired (cv_folds - single_split) per cell
wide <- stats::reshape(all_rows[, c("gate_method", "trait", "value",
                                      "rep_id", "lambda")],
                        direction = "wide",
                        idvar = c("trait", "rep_id", "lambda"),
                        timevar = "gate_method")
wide$delta_cv_minus_single <- wide$value.cv_folds - wide$value.single_split
delta_summary <- stats::aggregate(delta_cv_minus_single ~ trait + lambda,
                                    data = wide,
                                    FUN = function(x)
                                      c(mean = mean(x, na.rm = TRUE),
                                        sd   = stats::sd(x, na.rm = TRUE),
                                        n_pos = sum(x > 0, na.rm = TRUE),
                                        n_neg = sum(x < 0, na.rm = TRUE)))
delta_df <- data.frame(
  trait = delta_summary$trait,
  lambda = delta_summary$lambda,
  delta_mean = as.numeric(delta_summary$delta_cv_minus_single[, "mean"]),
  delta_sd = as.numeric(delta_summary$delta_cv_minus_single[, "sd"]),
  n_pos = as.numeric(delta_summary$delta_cv_minus_single[, "n_pos"]),
  n_neg = as.numeric(delta_summary$delta_cv_minus_single[, "n_neg"]),
  stringsAsFactors = FALSE
)
delta_df <- delta_df[order(delta_df$trait, delta_df$lambda), ]

# Markdown ----
md <- c(
  "# CV-folds vs single_split: default-flip evidence",
  "",
  sprintf("Run: %s", format(Sys.time())),
  sprintf("n=%d species, %d reps x %d gate methods x %d lambdas = %d cells",
          300L, n_reps, length(methods), length(lambdas), nrow(cells)),
  "",
  "## Per-trait per-lambda mean +/- SD across reps",
  "",
  "```",
  capture.output(print(agg_df, row.names = FALSE)),
  "```",
  "",
  "## Paired delta (cv_folds - single_split) per cell, summarised",
  "",
  "Lower delta is better for RMSE; higher is better for accuracy.",
  "n_pos / n_neg = number of reps where cv_folds was bigger / smaller",
  "than single_split.",
  "",
  "```",
  capture.output(print(delta_df, row.names = FALSE)),
  "```",
  "",
  "## Decision",
  ""
)

# Apply pre-registered decision rule
cont_deltas <- delta_df[delta_df$trait == "cont", "delta_mean"]
disc_deltas <- delta_df[delta_df$trait %in% c("bin", "cat3"), "delta_mean"]
mean_cont_lift <- mean(cont_deltas) * -1   # cv RMSE - single RMSE; negate so positive = cv better
mean_disc_lift <- mean(disc_deltas)        # cv acc - single acc; positive = cv better

verdict_cont <- mean_cont_lift >= 0.01 * mean(agg_df$mean[agg_df$trait == "cont"])
verdict_disc <- mean_disc_lift >= -0.005

md <- c(md,
        sprintf("- Mean continuous RMSE lift (cv - single, negated): %+.4g",
                mean_cont_lift),
        sprintf("- Mean discrete accuracy gap (cv - single): %+.4g",
                mean_disc_lift),
        sprintf("- Continuous lift >= 1%% of baseline RMSE? %s",
                if (verdict_cont) "YES" else "NO"),
        sprintf("- Discrete gap >= -0.5 pp? %s",
                if (verdict_disc) "YES" else "NO"),
        "",
        if (verdict_cont && verdict_disc)
          "**RECOMMENDATION: flip default to cv_folds in v0.9.2.**"
        else
          "**RECOMMENDATION: keep single_split as default; cv_folds remains opt-in.**")

writeLines(md, out_md)
log_line("DONE -- ", out_rds)
log_line("DONE -- ", out_md)
