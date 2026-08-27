#!/usr/bin/env Rscript
#
# script/gnn_evidence_campaign_phase_b.R
#
# GNN evidence Phase B campaign (G0-approved 2026-08-27).
# 243 cells × 30 paired seeds = 7,290 fits (phylo_MAR / covariate_MAR / MNAR).
#
# Parallel lanes (interleaved cell order from merge(cells, reps)):
#   lane 3a phylo+cov: 4,860 fits — PIGAUTO_LIST_JOBS=phylo_cov
#   lane 3b MNAR:      2,430 fits — PIGAUTO_LIST_JOBS=mnar
#
# Pin: candidate SHA 6fddd79, driver git SHA, host, seeds.
# Threads: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 (set by launcher).
#
# Env:
#   PIGAUTO_PKG_PATH        package root (devtools::load_all)
#   PIGAUTO_OUT_DIR         results directory (default: cwd/results_phase_b)
#   PIGAUTO_JOB_ID          0-based job index within the active arm (default: all)
#   PIGAUTO_LAMBDA_MODE     impute lambda_mode (default: fixed_1)
#   PIGAUTO_SKIP_EXISTING   skip job if output RDS exists (default: 0)
#   PIGAUTO_SMOKE           if 1, run only job 0
#   PIGAUTO_MECHANISM_ARM   all | MAR | MNAR (default: all)
#                           MAR = phylo_MAR + covariate_MAR (lane 3a)
#                           MNAR = MNAR only (lane 3b)
#   PIGAUTO_JOB_PREFIX      RDS prefix (default: gnn_phase_b)
#
# Output per job: <PIGAUTO_JOB_PREFIX>_job_<id>.rds
#
# G6 labeling: only covariate_MAR is "MAR"; phylo_MAR and MNAR are not MAR.

options(warn = 1, stringsAsFactors = FALSE)

CANDIDATE_SHA <- "6fddd79"
EPOCHS        <- 500L
N_REPS        <- 30L
DEP_STRENGTH  <- 1.5
ALL_MECHANISMS <- c("phylo_MAR", "covariate_MAR", "MNAR")
DRIVER_COL     <- 1L  # first trait column — fully observed under covariate_MAR

mechanism_arm <- Sys.getenv("PIGAUTO_MECHANISM_ARM", unset = "all")
MECHANISMS <- switch(
  mechanism_arm,
  all = ALL_MECHANISMS,
  MAR = c("phylo_MAR", "covariate_MAR"),
  MNAR = "MNAR",
  stop("PIGAUTO_MECHANISM_ARM must be all, MAR, or MNAR; got: ", mechanism_arm)
)
job_prefix <- Sys.getenv("PIGAUTO_JOB_PREFIX", unset = "gnn_phase_b")

pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH", unset = "")
if (!nzchar(pkg_path)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    this_file <- sub("^--file=", "", file_arg[[1]])
    pkg_path <- normalizePath(file.path(dirname(this_file), ".."),
                              mustWork = FALSE)
  } else {
    pkg_path <- getwd()
  }
}
stopifnot(file.exists(file.path(pkg_path, "DESCRIPTION")))

suppressPackageStartupMessages({
  library(ape)
  if (requireNamespace("MASS", quietly = TRUE)) NULL
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(pkg_path, quiet = TRUE)
  } else {
    library(pigauto)
  }
})

out_dir <- Sys.getenv("PIGAUTO_OUT_DIR", unset = file.path(getwd(), "results_phase_b"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lambda_mode <- Sys.getenv("PIGAUTO_LAMBDA_MODE", unset = "fixed_1")
skip_existing <- identical(Sys.getenv("PIGAUTO_SKIP_EXISTING", "0"), "1")

driver_sha <- tryCatch(
  system2("git", c("-C", pkg_path, "rev-parse", "HEAD"),
          stdout = TRUE, stderr = FALSE),
  error = function(e) character()
)
if (length(driver_sha) == 0L) driver_sha <- NA_character_
driver_sha <- driver_sha[[1L]]
hostname <- Sys.info()[["nodename"]]

build_all_cells <- function() {
  base <- expand.grid(
    family       = c("F1", "F2", "F3"),
    n_species    = c(100L, 300L, 1000L),
    phylo_lambda = c(0.2, 0.5, 1.0),
    miss_frac    = c(0.10, 0.30, 0.50),
    stringsAsFactors = FALSE
  )
  cells <- do.call(rbind, lapply(ALL_MECHANISMS, function(mech) {
    cbind(
      data.frame(missing_mechanism = mech, stringsAsFactors = FALSE),
      base
    )
  }))
  rownames(cells) <- NULL
  cells$cell_id <- seq_len(nrow(cells)) - 1L
  cells
}

build_jobs <- function() {
  cells <- build_all_cells()
  cells <- cells[cells$missing_mechanism %in% MECHANISMS, , drop = FALSE]
  reps <- data.frame(rep = seq_len(N_REPS))
  jobs <- merge(cells, reps, by = NULL)
  jobs$job_id <- seq_len(nrow(jobs)) - 1L
  jobs$seed <- 30000L + jobs$cell_id * 100L + jobs$rep
  jobs[, c("job_id", "cell_id", "missing_mechanism", "family", "n_species",
           "phylo_lambda", "miss_frac", "rep", "seed")]
}

build_jobs_full <- function() {
  cells <- build_all_cells()
  reps <- data.frame(rep = seq_len(N_REPS))
  jobs <- merge(cells, reps, by = NULL)
  jobs$job_id <- seq_len(nrow(jobs)) - 1L
  jobs$seed <- 30000L + jobs$cell_id * 100L + jobs$rep
  jobs[, c("job_id", "cell_id", "missing_mechanism", "family", "n_species",
           "phylo_lambda", "miss_frac", "rep", "seed")]
}

job_rds_path <- function(job_id) {
  file.path(out_dir, sprintf("%s_job_%04d.rds", job_prefix, job_id))
}

jobs <- build_jobs()
n_jobs <- nrow(jobs)

list_jobs <- Sys.getenv("PIGAUTO_LIST_JOBS", unset = "")
if (nzchar(list_jobs)) {
  full <- build_jobs_full()
  filt <- switch(list_jobs,
                 phylo_cov = full$missing_mechanism != "MNAR",
                 mnar      = full$missing_mechanism == "MNAR",
                 stop("PIGAUTO_LIST_JOBS must be phylo_cov or mnar"))
  cat(paste(full$job_id[filt], collapse = "\n"), "\n", sep = "")
  quit(save = "no")
}

if (identical(Sys.getenv("PIGAUTO_SMOKE", "0"), "1")) {
  jobs <- jobs[1L, , drop = FALSE]
  n_jobs <- 1L
}

job_env <- Sys.getenv("PIGAUTO_JOB_ID", unset = "")
run_all <- !nzchar(job_env)
job_ids <- if (run_all) jobs$job_id else as.integer(strsplit(job_env, ",", fixed = TRUE)[[1L]])
n_full <- nrow(build_jobs_full())
if (any(is.na(job_ids)) || any(job_ids < 0L) || any(job_ids >= n_full)) {
  stop("PIGAUTO_JOB_ID must be in 0..", n_full - 1L)
}

script_start <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start), ..., "\n",
      sep = "")
  flush.console()
}

log_line("=== GNN evidence Phase B campaign ===")
log_line(sprintf(
  "candidate=%s  driver=%s  host=%s  jobs=%s/%d  lambda_mode=%s",
  CANDIDATE_SHA, driver_sha, hostname,
  paste(job_ids, collapse = ","), nrow(build_jobs()), lambda_mode
))
log_line(sprintf(
  "arm=%s  mechanisms=%s  epochs=%d  prefix=%s  pkg=%s  out=%s",
  mechanism_arm, paste(MECHANISMS, collapse = ","), EPOCHS, job_prefix,
  pkg_path, out_dir
))

# ---- DGP helpers (regime-map design, issue #135) ---------------------------

transform_tree_pagel_local <- function(tree, lambda) {
  n_tip <- ape::Ntip(tree)
  node_depths <- ape::node.depth.edgelength(tree)
  parents <- tree$edge[, 1L]
  children <- tree$edge[, 2L]
  is_terminal <- children <= n_tip
  out <- tree
  out$edge.length[!is_terminal] <- lambda * tree$edge.length[!is_terminal]
  parent_depth_term <- node_depths[parents[is_terminal]]
  out$edge.length[is_terminal] <-
    tree$edge.length[is_terminal] + (1 - lambda) * parent_depth_term
  out
}

simulate_latents <- function(tree, K = 4L, seed, cross_r = 0.7) {
  set.seed(seed)
  n <- ape::Ntip(tree)
  indep <- matrix(NA_real_, n, K)
  for (k in seq_len(K)) {
    tr <- ape::rTraitCont(tree, model = "BM", sigma = 1.0, root.value = 0)
    indep[, k] <- as.numeric(tr[tree$tip.label])
  }
  Sigma <- matrix(cross_r, K, K)
  diag(Sigma) <- 1
  L <- chol(Sigma)
  Z <- indep %*% t(L)
  colnames(Z) <- paste0("Z", seq_len(K))
  rownames(Z) <- tree$tip.label
  Z
}

family_to_traits <- function(Z, family) {
  spp <- rownames(Z)
  n <- nrow(Z)
  if (family == "F1") {
    data.frame(
      trait1 = Z[, 1L],
      trait2 = Z[, 2L],
      trait3 = Z[, 3L],
      trait4 = Z[, 4L],
      row.names = spp
    )
  } else if (family == "F2") {
    z2 <- Z[, 2L]
    data.frame(
      trait1 = Z[, 1L],
      trait2 = z2,
      trait3 = sin(2 * Z[, 1L]) + 0.5 * Z[, 3L],
      trait4 = (z2^2) * sign(z2) + 0.5 * Z[, 4L],
      row.names = spp
    )
  } else if (family == "F3") {
    lat_cat <- sin(2 * Z[, 1L]) + 0.5 * Z[, 3L]
    q <- stats::quantile(lat_cat, probs = c(1/3, 2/3), type = 7)
    cat_vals <- ifelse(lat_cat <= q[[1L]], "A",
                ifelse(lat_cat <= q[[2L]], "B", "C"))
    rate <- exp(0.5 * Z[, 4L] + 0.3 * Z[, 1L])
    data.frame(
      cont1 = Z[, 1L],
      binary_trait = factor(ifelse(Z[, 2L] > 0, "yes", "no"),
                            levels = c("no", "yes")),
      cat_trait = factor(cat_vals, levels = c("A", "B", "C")),
      count_trait = as.integer(rpois(n, rate)),
      row.names = spp
    )
  } else {
    stop("unknown family: ", family)
  }
}

scale_vec <- function(x) {
  sx <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx < 1e-8) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / sx)
}

calibrate_plogis <- function(linpred, miss_frac) {
  obj <- function(c) (mean(stats::plogis(c + linpred)) - miss_frac)^2
  c_hat <- stats::optimize(obj, interval = c(-10, 10))$minimum
  stats::plogis(c_hat + linpred)
}

sample_col_miss <- function(p_row, n, frac, seed) {
  set.seed(seed)
  n_miss <- max(1L, round(n * frac))
  if (n_miss >= n) return(seq_len(n))
  sample.int(n, size = n_miss, prob = p_row)
}

punch_phylo_mar <- function(df, tree, frac, seed) {
  n <- nrow(df)
  R <- stats::cov2cor(ape::vcv(tree))
  R <- R[tree$tip.label, tree$tip.label, drop = FALSE]
  set.seed(seed)
  z <- as.numeric(MASS::mvrnorm(1L, mu = rep(0, n), Sigma = R))
  linpred <- -DEP_STRENGTH * scale_vec(z)
  p_row <- calibrate_plogis(linpred, frac)
  out <- df
  for (col in names(df)) {
    idx <- sample_col_miss(p_row, n, frac, seed + as.integer(charToRaw(col)[1L]))
    out[[col]][idx] <- NA
  }
  out
}

punch_covariate_mar <- function(df, frac, seed) {
  n <- nrow(df)
  driver <- df[[DRIVER_COL]]
  linpred <- DEP_STRENGTH * scale_vec(driver)
  p_row <- calibrate_plogis(linpred, frac)
  out <- df
  for (j in seq_along(names(df))) {
    if (j == DRIVER_COL) next  # driver fully observed — genuine MAR (G6)
    col <- names(df)[j]
    idx <- sample_col_miss(p_row, n, frac, seed + j)
    out[[col]][idx] <- NA
  }
  out
}

punch_mnar <- function(df, frac, seed) {
  n <- nrow(df)
  out <- df
  for (j in seq_along(names(df))) {
    col <- names(df)[j]
    y <- df[[j]]
    if (is.factor(y)) y <- as.numeric(y)
    linpred <- DEP_STRENGTH * scale_vec(y)
    p_row <- calibrate_plogis(linpred, frac)
    idx <- sample_col_miss(p_row, n, frac, seed + j)
    out[[col]][idx] <- NA
  }
  out
}

punch_mechanism <- function(df, tree, mechanism, frac, seed) {
  switch(mechanism,
         phylo_MAR = punch_phylo_mar(df, tree, frac, seed),
         covariate_MAR = punch_covariate_mar(df, frac, seed),
         MNAR = punch_mnar(df, frac, seed),
         stop("unknown mechanism: ", mechanism))
}

trait_test_loss <- function(eval_df, method) {
  sub <- eval_df[eval_df$method == method, , drop = FALSE]
  if (nrow(sub) == 0L) return(NA_real_)
  traits <- unique(sub$trait)
  losses <- vapply(traits, function(tr) {
    d <- sub[sub$trait == tr, , drop = FALSE]
    tp <- d$type[1L]
    if (tp %in% c("continuous", "count", "ordinal", "proportion")) {
      v <- d$value[d$metric == "rmse"]
      if (length(v) == 0L) NA_real_ else v[[1L]]
    } else if (tp %in% c("binary", "categorical")) {
      v <- d$value[d$metric == "accuracy"]
      if (length(v) == 0L) NA_real_ else 1 - v[[1L]]
    } else {
      NA_real_
    }
  }, numeric(1))
  mean(losses, na.rm = TRUE)
}

gate_summary <- function(fit) {
  rg <- as.numeric(fit$r_cal_gnn)
  rb <- as.numeric(fit$r_cal_bm)
  rm <- as.numeric(fit$r_cal_mean)
  nm <- names(fit$r_cal_gnn)
  if (length(nm) == 0L) nm <- paste0("col", seq_along(rg))
  data.frame(
    latent_col = nm,
    r_cal_gnn  = rg,
    r_cal_bm   = rb,
    r_cal_mean = rm,
    floor_mean_corner = rm > 0.99,
    stringsAsFactors = FALSE
  )
}

result_stub <- function(row, fit_ok, error = NA_character_, fit_sec = NA_real_,
                        blend_loss = NA_real_, baseline_loss = NA_real_,
                        paired_delta = NA_real_, gates = NULL,
                        floor_fired = NA, skipped = FALSE) {
  list(
    job_id = row$job_id,
    cell_id = row$cell_id,
    rep = row$rep,
    family = row$family,
    n_species = row$n_species,
    seed = row$seed,
    tree_seed = 20260827L + row$job_id,
    candidate_sha = CANDIDATE_SHA,
    driver_sha = driver_sha,
    hostname = hostname,
    phylo_lambda = row$phylo_lambda,
    miss_frac = row$miss_frac,
    missing_mechanism = row$missing_mechanism,
    g6_mar_label = identical(row$missing_mechanism, "covariate_MAR"),
    lambda_mode = lambda_mode,
    fit_ok = fit_ok,
    skipped = skipped,
    error = error,
    fit_sec = fit_sec,
    blend_loss = blend_loss,
    baseline_loss = baseline_loss,
    paired_delta = paired_delta,
    r_cal_gnn_mean = if (is.null(gates)) NA_real_ else mean(gates$r_cal_gnn, na.rm = TRUE),
    r_cal_bm_mean  = if (is.null(gates)) NA_real_ else mean(gates$r_cal_bm, na.rm = TRUE),
    r_cal_mean_mean = if (is.null(gates)) NA_real_ else mean(gates$r_cal_mean, na.rm = TRUE),
    gate_open_frac = if (is.null(gates)) NA_real_ else mean(gates$r_cal_gnn > 0, na.rm = TRUE),
    floor_fired = floor_fired,
    gates = gates
  )
}

run_one_job <- function(job_id) {
  all_jobs <- build_jobs_full()
  row <- all_jobs[all_jobs$job_id == job_id, , drop = FALSE]
  stopifnot(nrow(row) == 1L)
  if (!row$missing_mechanism %in% MECHANISMS) {
    stop("job ", job_id, " mechanism ", row$missing_mechanism,
         " not in active arm (", paste(MECHANISMS, collapse = ","), ")")
  }

  out_path <- job_rds_path(job_id)
  if (skip_existing && file.exists(out_path)) {
    log_line(sprintf("--- job %d  SKIP (exists) ---", job_id))
    prev <- readRDS(out_path)
    prev$skipped <- TRUE
    return(prev)
  }

  fam <- row$family
  mech <- row$missing_mechanism
  n_sp <- as.integer(row$n_species)
  seed <- as.integer(row$seed)
  phylo_lam <- row$phylo_lambda
  miss_frac <- row$miss_frac
  tree_seed <- 20260827L + job_id

  log_line(sprintf(
    "job %d  cell=%d  %s  mech=%s  n=%d  lambda=%.1f  miss=%.0f%%  rep=%d  seed=%d",
    job_id, row$cell_id, fam, mech, n_sp, phylo_lam, miss_frac * 100, row$rep, seed
  ))

  set.seed(tree_seed)
  tree_raw <- ape::rtree(n_sp)
  tree <- transform_tree_pagel_local(tree_raw, phylo_lam)

  Z <- simulate_latents(tree, seed = seed + 1000L)
  df_full <- family_to_traits(Z, fam)
  df_miss <- punch_mechanism(df_full, tree, mech, miss_frac, seed = seed + 2000L)

  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch({
    impute(
      traits        = df_miss,
      tree          = tree,
      seed          = seed,
      epochs        = EPOCHS,
      missing_frac  = 0.25,
      verbose       = FALSE,
      safety_floor  = TRUE,
      lambda_mode   = lambda_mode,
      log_transform = FALSE,
      trait_types   = if (fam == "F3") c(count_trait = "count") else NULL
    )
  }, error = function(e) {
    list(error = conditionMessage(e))
  })
  fit_sec <- proc.time()[["elapsed"]] - t0

  if (!is.null(res$error)) {
    log_line("  FIT FAILED: ", res$error)
    return(result_stub(row, fit_ok = FALSE, error = res$error, fit_sec = fit_sec))
  }

  fit  <- res$fit
  pd   <- res$data
  spl  <- res$splits
  eval_df <- evaluate(fit, data = pd, splits = spl)

  blend_loss    <- trait_test_loss(eval_df, "pigauto")
  baseline_loss <- trait_test_loss(eval_df, "baseline")
  paired_delta  <- blend_loss - baseline_loss

  gates <- gate_summary(fit)
  floor_any <- any(gates$floor_mean_corner, na.rm = TRUE)

  log_line(sprintf(
    "  blend=%.4f  baseline=%.4f  delta=%+.4f  fit_sec=%.1f  floor=%s",
    blend_loss, baseline_loss, paired_delta, fit_sec,
    if (floor_any) "yes" else "no"
  ))

  out <- result_stub(
    row,
    fit_ok = TRUE,
    fit_sec = fit_sec,
    blend_loss = blend_loss,
    baseline_loss = baseline_loss,
    paired_delta = paired_delta,
    gates = gates,
    floor_fired = floor_any
  )
  out$evaluation <- eval_df
  out
}

results <- lapply(job_ids, run_one_job)

for (i in seq_along(job_ids)) {
  jid <- job_ids[[i]]
  out_path <- job_rds_path(jid)
  if (!isTRUE(results[[i]]$skipped)) {
    saveRDS(results[[i]], out_path)
    log_line("wrote ", out_path)
  }
}

if (length(job_ids) == nrow(build_jobs()) && !identical(Sys.getenv("PIGAUTO_SMOKE", "0"), "1")) {
  fs <- sort(Sys.glob(file.path(out_dir, paste0(job_prefix, "_job_*.rds"))))
  rows <- lapply(fs, readRDS)
  summary_df <- do.call(rbind, lapply(rows, function(r) {
    data.frame(
      job_id = r$job_id,
      cell_id = r$cell_id,
      rep = r$rep,
      missing_mechanism = r$missing_mechanism,
      g6_mar_label = isTRUE(r$g6_mar_label),
      family = r$family,
      n_species = r$n_species,
      phylo_lambda = r$phylo_lambda,
      miss_frac = r$miss_frac,
      seed = r$seed,
      lambda_mode = r$lambda_mode,
      fit_ok = r$fit_ok,
      fit_sec = r$fit_sec,
      blend_loss = r$blend_loss,
      baseline_loss = r$baseline_loss,
      paired_delta = r$paired_delta,
      r_cal_gnn_mean = r$r_cal_gnn_mean,
      r_cal_bm_mean = r$r_cal_bm_mean,
      r_cal_mean_mean = r$r_cal_mean_mean,
      floor_fired = r$floor_fired,
      error = if (is.null(r$error)) NA_character_ else r$error,
      stringsAsFactors = FALSE
    )
  }))
  summary_stub <- paste0(job_prefix, "_summary")
  saveRDS(summary_df, file.path(out_dir, paste0(summary_stub, ".rds")))
  write.csv(summary_df, file.path(out_dir, paste0(summary_stub, ".csv")),
            row.names = FALSE)
  log_line("summary written (", nrow(summary_df), " rows)")
}

log_line("done")
