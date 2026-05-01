#!/usr/bin/env Rscript
#
# script/bench_avonet_full_local.R
#
# AVONET full n=9,993 bird bench, pigauto vs mean baseline, Apple
# MPS local (GPU failed on Vulcan L40S due to predict-stage OOM at
# n>=5000; MPS unified memory absorbs the n^2 attention at full
# scale).
#
# Mirrors bench_fishbase.R eval pattern: RMSE / pearson_r per trait
# plus conformal + MC-dropout 95% coverage.
#
# Run with:
#   PIGAUTO_PKG_PATH=$(pwd) Rscript script/bench_avonet_full_local.R
# Optional env vars:
#   PIGAUTO_N_SUBSET=NNNN  random subset from 9,993
#   PIGAUTO_N_IMPUTATIONS=K  defaults to 20
#
# Output:
#   script/bench_avonet_full_local.rds
#   script/bench_avonet_full_local.md

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
out_rds <- file.path(here, "script", "bench_avonet_full_local.rds")
out_md  <- file.path(here, "script", "bench_avonet_full_local.md")

SEED      <- 2026L
MISS_FRAC <- 0.30
N_SUB <- {
  v <- Sys.getenv("PIGAUTO_N_SUBSET", unset = "")
  if (nzchar(v)) as.integer(v) else NA_integer_
}
N_IMP <- {
  v <- Sys.getenv("PIGAUTO_N_IMPUTATIONS", unset = "")
  if (nzchar(v)) as.integer(v) else 20L
}

script_start <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start),
      ..., "\n", sep = "")
  flush.console()
}

log_line("Loading avonet_full + tree_full from pigauto ...")
e <- new.env(parent = emptyenv())
utils::data("avonet_full", package = "pigauto", envir = e)
utils::data("tree_full",   package = "pigauto", envir = e)
df   <- e$avonet_full
tree <- e$tree_full
rownames(df) <- df$Species_Key
df$Species_Key <- NULL

if (!is.na(N_SUB) && N_SUB < nrow(df)) {
  set.seed(SEED)
  keep <- sample(rownames(df), N_SUB)
  df   <- df[keep, , drop = FALSE]
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
  df   <- df[tree$tip.label, , drop = FALSE]
  log_line(sprintf("N_SUBSET = %d: random subset -> %d species",
                    N_SUB, nrow(df)))
}
stopifnot(all(rownames(df) == tree$tip.label))
log_line(sprintf("Aligned: %d species x %d traits", nrow(df), ncol(df)))

set.seed(SEED)
df_truth <- df
mask_test <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df),
                     dimnames = list(rownames(df), names(df)))
for (v in names(df)) {
  obs_idx <- which(!is.na(df[[v]]))
  to_hide <- sample(obs_idx, ceiling(length(obs_idx) * MISS_FRAC))
  mask_test[to_hide, v] <- TRUE
}
df_miss <- df
for (v in names(df)) df_miss[[v]][mask_test[, v]] <- NA

log_line("Held-out test cells per trait:")
print(colSums(mask_test))

safe_cor <- function(x, y) {
  idx <- which(is.finite(x) & is.finite(y))
  if (length(idx) < 2L) return(NA_real_)
  if (stats::sd(x[idx]) == 0 || stats::sd(y[idx]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[idx], y[idx]))
}

eval_completed <- function(completed, truth, mask, method, wall_s,
                            res_obj = NULL) {
  if (is.null(completed)) return(NULL)
  lo      <- if (!is.null(res_obj)) res_obj$prediction$conformal_lower else NULL
  hi      <- if (!is.null(res_obj)) res_obj$prediction$conformal_upper else NULL
  mi_list <- if (!is.null(res_obj)) res_obj$prediction$imputed_datasets else NULL
  rows <- list()
  for (v in colnames(mask)) {
    idx <- which(mask[, v])
    if (!length(idx)) next
    t_v <- truth[[v]][idx]; c_v <- completed[[v]][idx]
    if (is.factor(t_v) || is.ordered(t_v) || is.character(t_v)) {
      acc <- mean(as.character(c_v) == as.character(t_v), na.rm = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(
        method = method, trait = v, metric = "accuracy",
        value = acc, n_cells = length(idx), wall_s = wall_s)
    } else {
      rmse <- sqrt(mean((as.numeric(t_v) - as.numeric(c_v))^2, na.rm = TRUE))
      pear <- safe_cor(as.numeric(t_v), as.numeric(c_v))
      rows[[length(rows) + 1L]] <- data.frame(
        method = method, trait = v, metric = "rmse",
        value = rmse, n_cells = length(idx), wall_s = wall_s)
      rows[[length(rows) + 1L]] <- data.frame(
        method = method, trait = v, metric = "pearson_r",
        value = pear, n_cells = length(idx), wall_s = wall_s)
      t_num <- as.numeric(t_v)
      if (!is.null(lo) && !is.null(hi) && v %in% colnames(lo)) {
        lo_v <- lo[idx, v]; hi_v <- hi[idx, v]
        valid <- is.finite(lo_v) & is.finite(hi_v) & is.finite(t_num)
        if (any(valid)) {
          hits <- t_num[valid] >= lo_v[valid] & t_num[valid] <= hi_v[valid]
          rows[[length(rows) + 1L]] <- data.frame(
            method = method, trait = v, metric = "coverage95_conformal",
            value = mean(hits), n_cells = sum(valid), wall_s = wall_s)
        }
      }
      if (!is.null(mi_list) && length(mi_list) > 1L &&
          v %in% names(mi_list[[1]])) {
        draws_mat <- vapply(mi_list, function(d) as.numeric(d[idx, v]),
                             numeric(length(idx)))
        if (!is.matrix(draws_mat))
          draws_mat <- matrix(draws_mat, ncol = length(mi_list))
        if (nrow(draws_mat) == length(idx) && ncol(draws_mat) > 1L) {
          q_lo <- apply(draws_mat, 1L, stats::quantile, probs = 0.025,
                         na.rm = TRUE)
          q_hi <- apply(draws_mat, 1L, stats::quantile, probs = 0.975,
                         na.rm = TRUE)
          valid <- is.finite(q_lo) & is.finite(q_hi) & is.finite(t_num)
          if (any(valid)) {
            hits <- t_num[valid] >= q_lo[valid] & t_num[valid] <= q_hi[valid]
            rows[[length(rows) + 1L]] <- data.frame(
              method = method, trait = v, metric = "coverage95_mcdropout",
              value = mean(hits), n_cells = sum(valid), wall_s = wall_s)
          }
        }
      }
    }
  }
  if (!length(rows)) NULL else do.call(rbind, rows)
}

run_mean <- function() {
  out <- df_miss
  for (v in names(out)) {
    if (is.factor(out[[v]])) {
      mode <- names(sort(table(out[[v]], useNA = "no"), decreasing = TRUE))[1]
      out[[v]][mask_test[, v]] <- factor(mode, levels = levels(out[[v]]),
                                           ordered = is.ordered(out[[v]]))
    } else {
      out[[v]][mask_test[, v]] <- mean(out[[v]], na.rm = TRUE)
    }
  }
  list(completed = out, res = NULL)
}

run_pigauto <- function() {
  res <- pigauto::impute(df_miss, tree,
                           log_transform = TRUE,
                           missing_frac  = 0.20,
                           n_imputations = N_IMP,
                           epochs        = 500L,
                           verbose       = TRUE,
                           seed          = SEED)
  list(completed = res$completed, res = res)
}

timed <- function(tag, fn) {
  log_line("[STAGE] ", tag, " starting")
  t0 <- proc.time()[["elapsed"]]
  val <- fn()
  wall <- proc.time()[["elapsed"]] - t0
  log_line("[STAGE] ", tag, " done in ", round(wall, 1), " s")
  list(val = val, wall = wall)
}

r_mean <- timed("mean_baseline",   run_mean)
r_pig  <- timed("pigauto_default", run_pigauto)

ev_mean <- eval_completed(r_mean$val$completed, df_truth, mask_test,
                            "mean_baseline",   r_mean$wall, NULL)
ev_pig  <- eval_completed(r_pig$val$completed,  df_truth, mask_test,
                            "pigauto_default", r_pig$wall, r_pig$val$res)

all_rows <- do.call(rbind, Filter(Negate(is.null),
                                     list(ev_mean, ev_pig)))

saveRDS(list(
  results   = all_rows,
  seed      = SEED,
  miss_frac = MISS_FRAC,
  n_species = nrow(df),
  n_subset  = N_SUB,
  n_imputations = N_IMP
), out_rds)

md <- c(
  "# AVONET full n=9,993 x pigauto (local Mac MPS)",
  "",
  sprintf("n = %d species x %d traits (bundled avonet_full + tree_full).",
          nrow(df), ncol(df)),
  sprintf("Seed = %d, miss_frac = %.2f, n_imputations = %d.",
          SEED, MISS_FRAC, N_IMP),
  "",
  "## Per-trait metrics (pigauto_default vs mean_baseline)",
  "",
  "```",
  capture.output(print(all_rows, row.names = FALSE, max = 500)),
  "```"
)
writeLines(md, out_md)

log_line("=== DONE ===")
log_line("  rds: ", out_rds)
log_line("  md : ", out_md)
