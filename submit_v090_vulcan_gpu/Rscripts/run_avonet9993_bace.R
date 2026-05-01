#!/usr/bin/env Rscript
#
# Rscripts/run_avonet9993_bace.R
#
# AVONET 9,993 × BACE head-to-head on Vulcan GPU.
#
# Data: pigauto's bundled `avonet_full` + `tree_full` (9,993 bird species,
#   7 traits: Mass, Beak.Length_Culmen, Tarsus.Length, Wing.Length [continuous]
#   + Trophic.Level, Primary.Lifestyle [categorical] + Migration [ordinal]).
#
# Methods compared:
#   mean_baseline   — column mean / mode
#   pigauto_default — em_iterations = 0, n_imputations = 20
#   pigauto_em5     — em_iterations = 5, n_imputations = 20
#   bace_default    — BACE::bace() with OVR, MCMCglmm 2000 iter (CPU)
#
# Smoke mode:
#   Set PIGAUTO_SMOKE=1 to subset to the first 500 species for a
#   ~2–5 min end-to-end plumbing check locally before submitting to
#   Vulcan. Everything else is identical.
#
# Output:
#   bench_avonet9993_bace.rds — tidy data.frame; one row per (method, trait,
#                                metric) with value, n_cells, wall_s
#   bench_avonet9993_bace.md  — human-readable summary
#
# Run:
#   Rscript Rscripts/run_avonet9993_bace.R

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  library(torch)
  # Try library(pigauto) first (the Vulcan install path); fall back to
  # devtools::load_all() on a local dev machine.
  if (!requireNamespace("pigauto", quietly = TRUE)) {
    pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH",
                             "/Users/z3437171/Dropbox/Github Local/pigauto")
    if (!dir.exists(pkg_path)) {
      stop("pigauto not installed and no local source at ", pkg_path,
           ". Install via pak::pak('local::~/pigauto_src') or set ",
           "PIGAUTO_PKG_PATH to the checkout path.", call. = FALSE)
    }
    devtools::load_all(pkg_path, quiet = TRUE)
  } else {
    library(pigauto)
  }
})

SMOKE <- identical(Sys.getenv("PIGAUTO_SMOKE"), "1")
SEED  <- 2026L
MISS_FRAC <- 0.30

out_rds <- "bench_avonet9993_bace.rds"
out_md  <- "bench_avonet9993_bace.md"

log_line <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
      ..., "\n", sep = "")
  flush.console()
}

# -------------------------------------------------------------------------
# Load bundled AVONET 9,993 + BirdTree full
# -------------------------------------------------------------------------

log_line("Loading avonet_full + tree_full from pigauto package ...")
e <- new.env(parent = emptyenv())
utils::data("avonet_full", package = "pigauto", envir = e)
utils::data("tree_full",   package = "pigauto", envir = e)
df   <- e$avonet_full
tree <- e$tree_full

rownames(df) <- df$Species_Key
df$Species_Key <- NULL

if (SMOKE) {
  keep <- rownames(df)[seq_len(min(500L, nrow(df)))]
  df   <- df[keep, , drop = FALSE]
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
  df   <- df[match(tree$tip.label, rownames(df)), , drop = FALSE]
  log_line("SMOKE mode: subsetting to ", nrow(df), " species")
}

stopifnot(all(rownames(df) == tree$tip.label))
log_line(sprintf("Aligned: %d species × %d traits", nrow(df), ncol(df)))

# -------------------------------------------------------------------------
# Inject MCAR 30% hold-out mask
# -------------------------------------------------------------------------

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

# -------------------------------------------------------------------------
# Eval helper: captures RMSE/accuracy + both coverage types
# -------------------------------------------------------------------------

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
      if (!is.null(mi_list) && length(mi_list) > 1L && v %in% names(mi_list[[1]])) {
        draws_mat <- vapply(mi_list, function(d) as.numeric(d[idx, v]),
                             numeric(length(idx)))
        if (!is.matrix(draws_mat))
          draws_mat <- matrix(draws_mat, ncol = length(mi_list))
        if (nrow(draws_mat) == length(idx) && ncol(draws_mat) > 1L) {
          q_lo <- apply(draws_mat, 1L, stats::quantile, probs = 0.025, na.rm = TRUE)
          q_hi <- apply(draws_mat, 1L, stats::quantile, probs = 0.975, na.rm = TRUE)
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

# -------------------------------------------------------------------------
# Method: mean / mode baseline
# -------------------------------------------------------------------------

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

# -------------------------------------------------------------------------
# Method: pigauto (default / em5)
# -------------------------------------------------------------------------

run_pigauto <- function(em_iter) {
  res <- pigauto::impute(df_miss, tree,
                           log_transform  = TRUE,
                           missing_frac   = 0.20,
                           n_imputations  = 20L,
                           epochs         = 500L,
                           verbose        = TRUE,
                           seed           = SEED,
                           em_iterations  = em_iter)
  list(completed = res$completed, res = res)
}

# -------------------------------------------------------------------------
# Method: BACE (CPU; may skip on install miss or error)
# -------------------------------------------------------------------------

run_bace <- function() {
  if (!requireNamespace("BACE", quietly = TRUE)) {
    log_line("BACE not installed — skipping")
    return(list(completed = NULL, res = NULL))
  }
  out <- tryCatch({
    BACE::bace(
      data    = df_miss,
      tree    = tree,
      n_iter  = 2000L,
      burnin  = 500L,
      thin    = 5L,
      ovr     = TRUE,
      verbose = FALSE
    )
  }, error = function(e) {
    log_line("BACE run failed: ", conditionMessage(e))
    NULL
  })
  if (is.null(out)) return(list(completed = NULL, res = NULL))
  # BACE shape varies by version — try common slots
  completed <- tryCatch({
    if ("completed"    %in% names(out)) out$completed
    else if ("data"    %in% names(out)) out$data
    else if ("imputed_data" %in% names(out)) out$imputed_data
    else NULL
  }, error = function(e) NULL)
  if (is.null(completed)) log_line("BACE output shape not recognised")
  list(completed = completed, res = NULL)
}

# -------------------------------------------------------------------------
# Run all four methods
# -------------------------------------------------------------------------

timed <- function(tag, fn) {
  log_line("[STAGE] ", tag, " starting")
  t0 <- proc.time()[["elapsed"]]
  val <- fn()
  wall <- proc.time()[["elapsed"]] - t0
  log_line("[STAGE] ", tag, " done in ", round(wall, 1), " s")
  list(val = val, wall = wall)
}

r_mean <- timed("mean_baseline",   run_mean)
r_def  <- timed("pigauto_default", function() run_pigauto(0L))
r_em5  <- timed("pigauto_em5",     function() run_pigauto(5L))
r_bace <- timed("bace_default",    run_bace)

ev_mean <- eval_completed(r_mean$val$completed, df_truth, mask_test,
                            "mean_baseline",   r_mean$wall, NULL)
ev_def  <- eval_completed(r_def$val$completed,  df_truth, mask_test,
                            "pigauto_default", r_def$wall, r_def$val$res)
ev_em5  <- eval_completed(r_em5$val$completed,  df_truth, mask_test,
                            "pigauto_em5",     r_em5$wall, r_em5$val$res)
ev_bace <- eval_completed(r_bace$val$completed, df_truth, mask_test,
                            "bace_default",    r_bace$wall, NULL)

all_rows <- do.call(rbind, Filter(Negate(is.null),
                                     list(ev_mean, ev_def, ev_em5, ev_bace)))

saveRDS(list(
  results   = all_rows,
  seed      = SEED,
  miss_frac = MISS_FRAC,
  n_species = nrow(df),
  smoke     = SMOKE,
  bace_ran  = !is.null(ev_bace)
), out_rds)

# -------------------------------------------------------------------------
# Markdown summary
# -------------------------------------------------------------------------

md <- c(
  "# AVONET 9,993 × pigauto + BACE head-to-head",
  "",
  sprintf("Mode: %s. n = %d species × %d traits.",
          if (SMOKE) "SMOKE (subset)" else "full", nrow(df), ncol(df)),
  sprintf("Seed = %d, miss_frac = %.2f.", SEED, MISS_FRAC),
  if (is.null(ev_bace)) "**BACE skipped** (not installed or failed)." else "",
  "",
  "## Per-trait metrics",
  "",
  "```",
  capture.output(print(all_rows, row.names = FALSE, max = 500)),
  "```"
)
writeLines(md, out_md)

log_line("=== DONE ===")
log_line("  rds: ", out_rds)
log_line("  md : ", out_md)
