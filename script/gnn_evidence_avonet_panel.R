#!/usr/bin/env Rscript
#
# script/gnn_evidence_avonet_panel.R
#
# AVONET300 mixed-type multi-seed panel (real-data corroboration).
# 5 seeds × 3 method arms = 15 jobs; 30% MCAR held-out per trait.
#
# Methods:
#   pigauto_fixed1  — primary arm (lambda_mode = fixed_1); also records
#                     baseline (= Rphylopars joint path when available)
#   pigauto_bayes   — sensitivity arm (lambda_mode = bayes; continuous focus)
#   bace            — BACE::bace chained equations when installed
#
# Pin: candidate SHA 6fddd79, driver git SHA, host, seeds.
#
# Env:
#   PIGAUTO_PKG_PATH      package root (devtools::load_all)
#   PIGAUTO_OUT_DIR       results directory (default: cwd/results_avonet_panel)
#   PIGAUTO_JOB_ID        0-based job index
#   PIGAUTO_SKIP_EXISTING skip if output RDS exists (default: 0)
#   PIGAUTO_SMOKE         if 1, run only job 0
#
# Output: avonet_panel_job_<id>.rds

options(warn = 1, stringsAsFactors = FALSE)

CANDIDATE_SHA <- "6fddd79"
MISS_FRAC     <- 0.30
N_IMP         <- 20L
EPOCHS        <- 500L
SEEDS         <- c(2026L, 2027L, 2028L, 2029L, 2030L)
METHODS       <- c("pigauto_fixed1", "pigauto_bayes", "bace")

pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH", unset = "")
if (!nzchar(pkg_path)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    this_file <- sub("^--file=", "", file_arg[[1]])
    pkg_path <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  } else {
    pkg_path <- getwd()
  }
}
stopifnot(file.exists(file.path(pkg_path, "DESCRIPTION")))

suppressPackageStartupMessages({
  library(ape)
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(pkg_path, quiet = TRUE)
  } else {
    library(pigauto)
  }
})

out_dir <- Sys.getenv("PIGAUTO_OUT_DIR", unset = file.path(getwd(), "results_avonet_panel"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

skip_existing <- identical(Sys.getenv("PIGAUTO_SKIP_EXISTING", "0"), "1")

driver_sha <- tryCatch(
  system2("git", c("-C", pkg_path, "rev-parse", "HEAD"),
          stdout = TRUE, stderr = FALSE),
  error = function(e) character()
)
if (length(driver_sha) == 0L) driver_sha <- NA_character_
driver_sha <- driver_sha[[1L]]
hostname <- Sys.info()[["nodename"]]

build_jobs <- function() {
  grid <- expand.grid(
    seed   = SEEDS,
    method = METHODS,
    stringsAsFactors = FALSE
  )
  grid$job_id <- seq_len(nrow(grid)) - 1L
  grid[, c("job_id", "seed", "method")]
}

jobs <- build_jobs()
n_jobs <- nrow(jobs)

if (identical(Sys.getenv("PIGAUTO_SMOKE", "0"), "1")) {
  jobs <- jobs[1L, , drop = FALSE]
  n_jobs <- 1L
}

job_env <- Sys.getenv("PIGAUTO_JOB_ID", unset = "")
run_all <- !nzchar(job_env)
job_ids <- if (run_all) jobs$job_id else as.integer(strsplit(job_env, ",", fixed = TRUE)[[1L]])
if (any(is.na(job_ids)) || any(job_ids < 0L) || any(job_ids >= nrow(build_jobs()))) {
  stop("PIGAUTO_JOB_ID must be in 0..", nrow(build_jobs()) - 1L)
}

script_start <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start), ..., "\n", sep = "")
  flush.console()
}

log_line("=== AVONET300 multi-seed panel ===")
log_line(sprintf(
  "candidate=%s  driver=%s  host=%s  jobs=%s/%d",
  CANDIDATE_SHA, driver_sha, hostname,
  paste(job_ids, collapse = ","), nrow(build_jobs())
))
log_line(sprintf("miss_frac=%.2f  seeds=%d  methods=%s",
                 MISS_FRAC, length(SEEDS), paste(METHODS, collapse = ",")))
log_line(sprintf("Rphylopars=%s  BACE=%s  MCMCglmm=%s",
                 requireNamespace("Rphylopars", quietly = TRUE),
                 requireNamespace("BACE", quietly = TRUE),
                 requireNamespace("MCMCglmm", quietly = TRUE)))

# ---- data load (AVONET300) -------------------------------------------------

load_avonet300 <- function() {
  e <- new.env(parent = emptyenv())
  utils::data("avonet300", package = "pigauto", envir = e)
  utils::data("tree300",   package = "pigauto", envir = e)
  df <- e$avonet300
  tree <- e$tree300
  rownames(df) <- df$Species_Key
  df$Species_Key <- NULL
  stopifnot(all(rownames(df) == tree$tip.label))
  list(df = df, tree = tree)
}

make_mcar_mask <- function(df, seed, miss_frac) {
  set.seed(seed)
  mask <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df),
                 dimnames = list(rownames(df), names(df)))
  for (v in names(df)) {
    obs_idx <- which(!is.na(df[[v]]))
    if (!length(obs_idx)) next
    n_hide <- max(1L, ceiling(length(obs_idx) * miss_frac))
    mask[sample(obs_idx, n_hide), v] <- TRUE
  }
  mask
}

apply_mask <- function(df, mask) {
  out <- df
  for (v in names(df)) out[[v]][mask[, v]] <- NA
  out
}

safe_cor <- function(x, y) {
  idx <- which(is.finite(x) & is.finite(y))
  if (length(idx) < 2L) return(NA_real_)
  if (stats::sd(x[idx]) == 0 || stats::sd(y[idx]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[idx], y[idx]))
}

eval_completed <- function(completed, truth, mask, method, wall_s) {
  if (is.null(completed)) return(NULL)
  rows <- list()
  for (v in colnames(mask)) {
    idx <- which(mask[, v])
    if (!length(idx)) next
    t_v <- truth[[v]][idx]
    c_v <- completed[[v]][idx]
    tp <- if (is.ordered(t_v)) "ordinal"
          else if (is.factor(t_v)) "categorical"
          else "continuous"
    if (tp %in% c("categorical", "ordinal")) {
      acc <- mean(as.character(c_v) == as.character(t_v), na.rm = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(
        method = method, trait = v, type = tp, metric = "accuracy",
        value = acc, n_cells = length(idx), wall_s = wall_s,
        stringsAsFactors = FALSE)
    } else {
      rmse <- sqrt(mean((as.numeric(t_v) - as.numeric(c_v))^2, na.rm = TRUE))
      pear <- safe_cor(as.numeric(t_v), as.numeric(c_v))
      rows[[length(rows) + 1L]] <- data.frame(
        method = method, trait = v, type = tp, metric = "rmse",
        value = rmse, n_cells = length(idx), wall_s = wall_s,
        stringsAsFactors = FALSE)
      rows[[length(rows) + 1L]] <- data.frame(
        method = method, trait = v, type = tp, metric = "pearson_r",
        value = pear, n_cells = length(idx), wall_s = wall_s,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) NULL else do.call(rbind, rows)
}

eval_impute_result <- function(res, method_label, wall_s) {
  eval_df <- evaluate(res$fit, data = res$data, splits = res$splits)
  rows <- list()
  for (m in unique(eval_df$method)) {
    sub <- eval_df[eval_df$method == m, , drop = FALSE]
    tag <- if (m == "pigauto") method_label else paste0("baseline_", method_label)
    for (i in seq_len(nrow(sub))) {
      rows[[length(rows) + 1L]] <- data.frame(
        method = tag,
        trait = sub$trait[i],
        type = sub$type[i],
        metric = sub$metric[i],
        value = sub$value[i],
        n_cells = sub$n_test[i],
        wall_s = wall_s,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

run_pigauto <- function(df_miss, tree, seed, lambda_mode) {
  impute(
    traits        = df_miss,
    tree          = tree,
    seed          = seed,
    epochs        = EPOCHS,
    missing_frac  = 0.20,
    n_imputations = N_IMP,
    verbose       = FALSE,
    safety_floor  = TRUE,
    lambda_mode   = lambda_mode,
    log_transform = TRUE
  )
}

run_bace <- function(df_miss, tree, mask_test) {
  if (!requireNamespace("BACE", quietly = TRUE)) {
    return(list(completed = NULL, error = "BACE not installed"))
  }
  tree_b <- tree
  if (any(tree_b$edge.length == 0, na.rm = TRUE)) {
    tree_b$edge.length[tree_b$edge.length == 0] <- 1e-8
  }
  df_b <- df_miss
  df_b$Species <- rownames(df_miss)
  all_traits <- setdiff(names(df_b), "Species")
  fixformula <- lapply(all_traits, function(v) {
    others <- setdiff(all_traits, v)
    paste0(v, " ~ ", paste(others, collapse = " + "))
  })
  out <- tryCatch({
    BACE::bace(
      fixformula     = fixformula,
      ran_phylo_form = "~ 1 |Species",
      phylo          = tree_b,
      data           = df_b,
      nitt           = 2000L,
      burnin         = 500L,
      thin           = 5L,
      runs           = 2L,
      n_final        = 2L,
      verbose        = FALSE,
      skip_conv      = TRUE
    )
  }, error = function(e) {
    list(.error = conditionMessage(e))
  })
  if (!is.null(out$.error)) {
    return(list(completed = NULL, error = out$.error))
  }
  imputed_sets <- if ("imputed_datasets" %in% names(out)) out$imputed_datasets
                 else if ("imputed_data" %in% names(out)) out$imputed_data
                 else if ("data" %in% names(out)) list(out$data)
                 else NULL
  if (is.null(imputed_sets) || !length(imputed_sets)) {
    return(list(completed = NULL, error = "BACE output shape not recognised"))
  }
  M <- length(imputed_sets)
  completed <- df_miss
  for (v in names(completed)) {
    if (!any(mask_test[, v])) next
    draws <- sapply(imputed_sets, function(d) d[[v]])
    if (!is.matrix(draws)) draws <- matrix(draws, ncol = M)
    if (is.factor(completed[[v]])) {
      for (i in which(mask_test[, v])) {
        vals <- as.character(draws[i, ])
        vals <- vals[!is.na(vals)]
        if (!length(vals)) next
        mode_val <- names(sort(table(vals), decreasing = TRUE))[1]
        completed[[v]][i] <- factor(mode_val, levels = levels(completed[[v]]),
                                    ordered = is.ordered(completed[[v]]))
      }
    } else {
      for (i in which(mask_test[, v])) {
        vals <- as.numeric(draws[i, ])
        completed[[v]][i] <- stats::median(vals, na.rm = TRUE)
      }
    }
  }
  list(completed = completed, error = NA_character_)
}

result_stub <- function(row, fit_ok, error = NA_character_, fit_sec = NA_real_,
                        metrics = NULL, skipped = FALSE, extra = list()) {
  c(list(
    campaign      = "avonet_panel",
    job_id        = row$job_id,
    seed          = row$seed,
    method        = row$method,
    miss_frac     = MISS_FRAC,
    n_species     = 300L,
    n_traits      = 7L,
    candidate_sha = CANDIDATE_SHA,
    driver_sha    = driver_sha,
    hostname      = hostname,
    fit_ok        = fit_ok,
    skipped       = skipped,
    error         = error,
    fit_sec       = fit_sec,
    metrics       = metrics,
    rphylopars_available = requireNamespace("Rphylopars", quietly = TRUE),
    bace_available       = requireNamespace("BACE", quietly = TRUE)
  ), extra)
}

run_one_job <- function(job_id) {
  all_jobs <- build_jobs()
  row <- all_jobs[all_jobs$job_id == job_id, , drop = FALSE]
  stopifnot(nrow(row) == 1L)

  out_path <- file.path(out_dir, sprintf("avonet_panel_job_%04d.rds", job_id))
  if (skip_existing && file.exists(out_path)) {
    prev <- readRDS(out_path)
    if (isTRUE(prev$fit_ok)) {
      log_line(sprintf("--- job %d  SKIP (ok exists) ---", job_id))
      prev$skipped <- TRUE
      return(prev)
    }
    log_line(sprintf("--- job %d  RE-RUN (prior fit_ok=FALSE) ---", job_id))
  }

  seed <- as.integer(row$seed)
  method <- row$method
  log_line(sprintf("job %d  seed=%d  method=%s", job_id, seed, method))

  av <- load_avonet300()
  df_truth <- av$df
  tree <- av$tree
  mask_test <- make_mcar_mask(df_truth, seed, MISS_FRAC)
  df_miss <- apply_mask(df_truth, mask_test)

  t0 <- proc.time()[["elapsed"]]
  metrics <- NULL
  err <- NA_character_

  if (method %in% c("pigauto_fixed1", "pigauto_bayes")) {
    lambda_mode <- if (method == "pigauto_fixed1") "fixed_1" else "bayes"
    res <- tryCatch({
      run_pigauto(df_miss, tree, seed, lambda_mode)
    }, error = function(e) {
      list(.error = conditionMessage(e))
    })
    fit_sec <- proc.time()[["elapsed"]] - t0
    if (!is.null(res$.error)) {
      log_line("  FIT FAILED: ", res$.error)
      return(result_stub(row, fit_ok = FALSE, error = res$.error, fit_sec = fit_sec))
    }
    metrics <- eval_impute_result(res, method, fit_sec)
    gates <- as.numeric(res$fit$r_cal_gnn)
    log_line(sprintf("  ok  fit_sec=%.1f  gate_mean=%.3f  rows=%d",
                     fit_sec, mean(gates, na.rm = TRUE), nrow(metrics)))
    out <- result_stub(row, fit_ok = TRUE, fit_sec = fit_sec, metrics = metrics,
                       extra = list(
                         lambda_mode = lambda_mode,
                         r_cal_gnn_mean = mean(gates, na.rm = TRUE),
                         gate_open_frac = mean(gates > 0, na.rm = TRUE)
                       ))
    saveRDS(out, out_path)
    return(out)
  }

  if (method == "bace") {
    bace_res <- run_bace(df_miss, tree, mask_test)
    fit_sec <- proc.time()[["elapsed"]] - t0
    if (is.null(bace_res$completed)) {
      log_line("  BACE FAILED: ", bace_res$error)
      return(result_stub(row, fit_ok = FALSE, error = bace_res$error, fit_sec = fit_sec))
    }
    metrics <- eval_completed(bace_res$completed, df_truth, mask_test,
                              "bace", fit_sec)
    log_line(sprintf("  ok  fit_sec=%.1f  rows=%d", fit_sec, nrow(metrics)))
    out <- result_stub(row, fit_ok = TRUE, fit_sec = fit_sec, metrics = metrics)
    saveRDS(out, out_path)
    return(out)
  }

  stop("unknown method: ", method)
}

results <- lapply(job_ids, run_one_job)

for (i in seq_along(job_ids)) {
  jid <- job_ids[[i]]
  out_path <- file.path(out_dir, sprintf("avonet_panel_job_%04d.rds", jid))
  if (!isTRUE(results[[i]]$skipped) && !file.exists(out_path)) {
    saveRDS(results[[i]], out_path)
    log_line("wrote ", out_path)
  }
}

log_line("done")
