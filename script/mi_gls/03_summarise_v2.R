#!/usr/bin/env Rscript
# script/mi_gls/03_summarise_v2.R
#
# Summarise script/mi_gls/01_cell_v2.R output into:
#   (a) a regime x method x downstream table, read by 04_acceptance.R (G6).
#       Methods: "complete" (the reference, D2), posterior_full, and
#       posterior_none in both-missing regimes only (D6). Columns: paired
#       bias vs complete and its MCSE, SE ratio (mean SE / empirical SD),
#       coverage and MCSE, the complete-data reference values, fit-failure
#       rate, and per-regime x method convergence (n_fits, n_converged,
#       median max R-hat, median min ESS).
#   (b) a regime x method x trait x mechanism per-cell coverage table
#       (n, covered_sum, mean_width, convergence), read by
#       05_cell_coverage.R (G7). Methods: posterior_full, posterior_none
#       (both-missing regimes) and conformal (descriptive comparator, D7).
#
# Expected grid, fail-closed (design.md section 5c, D3): rows are built
# from script/mi_gls/regimes.R, never from the files present. Env
# MI_N_REPS (default 200) sets the expected reps 1..MI_N_REPS; env
# MI_REGIMES (e.g. 1,21 or 17-24; default all) restricts the regimes for
# staged checks, and the gates then withhold their pass token. A missing
# rep file counts as a fit failure (fit_failure_rate = 1 - finite
# estimates / expected reps) and as non-converged (n_fits = expected reps).
# Every expected row is written even when no rep file exists for it, so
# the gates see the gap. The default grid is regimes 1-40 (25-40 are the
# in-model twins of 1-16, CP2 follow-up 2026-09-24). Both CSVs carry the
# regime's `dgp` ("tree_raw", "kronecker", "twin") and `gated` (TRUE for
# 17-40; 1-16 are the reported stress test) from regimes.R.
#
# Downstream-coverage truth (D1): regimes 1-16 and their twins 25-40 use
# regimes$true_beta_pop (rho = 0.7, exact for any analysis model under
# the proportional DGP). The Kronecker regimes 17-24 use a pseudo-truth:
# the mean complete-data slope over the expected reps of that (regime,
# downstream). `covered` is recomputed here for complete and every method
# from the saved estimate/se/df (t quantile with the saved df; normal
# quantile if df is missing). The saved true_beta_pop coverage is kept as
# the descriptive column coverage_pop.
#
# Design-review item B3: paired bias, SE ratio and coverage for the
# posterior methods use converged reps only (posterior_none shares the
# chains and diagnostics of the same seed). "complete" has no convergence
# concept and uses every expected rep with a finite estimate.
#
# Usage:
#   Rscript script/mi_gls/03_summarise_v2.R <indir> <out_summary.csv> [out.md] [out_cell_coverage.csv]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("expected: indir out_summary.csv [out.md] [out_cell_coverage.csv]", call. = FALSE)
indir       <- args[[1L]]
out_csv     <- args[[2L]]
out_md      <- if (length(args) >= 3L && nzchar(args[[3L]])) args[[3L]] else NA_character_
out_cell_csv <- if (length(args) >= 4L && nzchar(args[[4L]])) args[[4L]] else NA_character_

source(file.path("script", "mi_gls", "regimes.R"))

ex <- mi_gls_v2_expected()
exp_regimes <- regimes[regimes$regime_id %in% ex$regime_ids, ]
n_reps <- ex$n_reps
cat(sprintf("EXPECTED_GRID %s%s\n", ex$label, if (ex$full) "" else " (NOT the planned grid)"))

files <- list.files(indir, pattern = "^regime_[0-9]+_rep_[0-9]+\\.rds$", full.names = TRUE)
if (length(files) == 0L) stop("no .rds files found in ", indir, call. = FALSE)

cells <- list(); diags <- list(); cell_dets <- list(); prov <- list(); conf_err <- list()
n_ignored <- 0L
for (f in files) {
  r <- readRDS(f)
  if (!(r$regime_id %in% ex$regime_ids) || r$rep < 1L || r$rep > n_reps) {
    n_ignored <- n_ignored + 1L
    next
  }
  d <- r$results
  # Harness before 2026-09-24 saved the true_beta_pop coverage as `covered`.
  if (!("covered_pop" %in% names(d)) && ("covered" %in% names(d))) d$covered_pop <- d$covered
  d$covered <- NULL
  if (!("n_ok" %in% names(d))) d$n_ok <- NA_integer_
  d$regime_id <- r$regime_id; d$rep <- r$rep
  cells[[length(cells) + 1L]] <- d
  if (!is.null(r$diagnostics) && nrow(r$diagnostics)) {
    dg <- r$diagnostics; dg$regime_id <- r$regime_id; dg$rep <- r$rep
    diags[[length(diags) + 1L]] <- dg
  }
  if (!is.null(r$cell_detail) && nrow(r$cell_detail)) {
    cd <- r$cell_detail; cd$regime_id <- r$regime_id; cd$rep <- r$rep
    cell_dets[[length(cell_dets) + 1L]] <- cd
  }
  prov[[length(prov) + 1L]] <- data.frame(
    regime_id = r$regime_id, rep = r$rep,
    code_sha = if (is.null(r$code_sha)) NA_character_ else as.character(r$code_sha),
    pigauto_version = if (is.null(r$pigauto_version)) NA_character_ else as.character(r$pigauto_version))
  if (!is.null(r$conformal_error) && !is.na(r$conformal_error)) {
    conf_err[[length(conf_err) + 1L]] <- sprintf("regime %d rep %d: %s", r$regime_id, r$rep, r$conformal_error)
  }
}
if (n_ignored > 0L) cat(sprintf("IGNORED %d rep file(s) outside the expected grid\n", n_ignored))
if (!length(cells)) stop("no rep files inside the expected grid (", ex$label, ") in ", indir, call. = FALSE)

all_cells <- do.call(rbind, cells)
all_diag  <- if (length(diags))     do.call(rbind, diags)     else NULL
all_cellD <- if (length(cell_dets)) do.call(rbind, cell_dets) else NULL
prov_df   <- do.call(rbind, prov)

# ---- presence and provenance ------------------------------------------------
present <- unique(prov_df[, c("regime_id", "rep")])
for (rid in exp_regimes$regime_id) {
  k <- sum(present$regime_id == rid)
  if (k < n_reps) cat(sprintf("INCOMPLETE regime %d: %d/%d rep files\n", rid, k, n_reps))
}
shas <- sort(unique(ifelse(is.na(prov_df$code_sha), "unrecorded", prov_df$code_sha)))
vers <- sort(unique(ifelse(is.na(prov_df$pigauto_version), "unrecorded", prov_df$pigauto_version)))
code_sha_label <- if (length(shas) == 1L) shas else paste0("MIXED:", paste(shas, collapse = "|"))
cat(sprintf("PROVENANCE code_sha=%s pigauto_version=%s\n", code_sha_label, paste(vers, collapse = "|")))
if (length(shas) > 1L) cat("MIXED_CODE_SHA: rep files come from more than one code version\n")
if (length(conf_err)) {
  cat(sprintf("CONFORMAL_FAILED in %d rep file(s); first: %s\n", length(conf_err), conf_err[[1L]]))
}

# ---- downstream-coverage truth (D1) -----------------------------------------
covers_t <- function(est, se, dfr, truth, conf = 0.95) {
  ok <- is.finite(est) & is.finite(se) & is.finite(truth)
  crit <- ifelse(is.finite(dfr) & dfr > 0,
                 stats::qt(1 - (1 - conf) / 2, pmax(dfr, 1)),
                 stats::qnorm(1 - (1 - conf) / 2))
  out <- abs(est - truth) <= crit * se
  out[!ok] <- NA
  out
}

truth_rows <- do.call(rbind, lapply(exp_regimes$regime_id, function(rid) {
  do.call(rbind, lapply(c("gls", "phylolm"), function(ds) {
    if (regimes$dgp[regimes$regime_id == rid] != "kronecker") {   # tree_raw 1-16, twin 25-40
      data.frame(regime_id = rid, downstream = ds,
                 truth = regimes$true_beta_pop[regimes$regime_id == rid],
                 truth_source = "rho", n_truth = NA_integer_)
    } else {
      ce <- all_cells$estimate[all_cells$regime_id == rid & all_cells$method == "complete" &
                               all_cells$downstream == ds]
      ce <- ce[is.finite(ce)]
      data.frame(regime_id = rid, downstream = ds,
                 truth = if (length(ce)) mean(ce) else NA_real_,
                 truth_source = "pseudo_complete_mean", n_truth = length(ce))
    }
  }))
}))
all_cells <- merge(all_cells, truth_rows[, c("regime_id", "downstream", "truth")],
                   by = c("regime_id", "downstream"), all.x = TRUE)
all_cells$covered <- covers_t(all_cells$estimate, all_cells$se, all_cells$df, all_cells$truth)

converged_reps <- function(regime_id, method) {
  if (is.null(all_diag)) return(integer(0))
  sub <- all_diag[all_diag$regime_id == regime_id & all_diag$method == method &
                  (all_diag$converged %in% TRUE), ]
  unique(sub$rep)
}

# ---- (a) regime x method x downstream summary --------------------------------
ratio_stats <- function(est, se) {
  fe <- is.finite(est)
  emp_sd  <- if (sum(fe) > 1) stats::sd(est[fe]) else NA_real_
  mean_se <- if (any(is.finite(se))) mean(se[is.finite(se)]) else NA_real_
  se_ratio <- if (is.finite(mean_se) && is.finite(emp_sd) && emp_sd > 0) mean_se / emp_sd else NA_real_
  list(emp_sd = emp_sd, mean_se = mean_se, se_ratio = se_ratio)
}
cov_stats <- function(covered) {
  ok <- !is.na(covered)
  cv <- if (any(ok)) mean(covered[ok]) else NA_real_
  list(coverage = cv,
       coverage_mcse = if (is.finite(cv)) sqrt(cv * (1 - cv) / sum(ok)) else NA_real_)
}
mean_or_na <- function(x) if (any(!is.na(x))) mean(x[!is.na(x)]) else NA_real_

summarise_complete <- function(sub) {
  fin <- is.finite(sub$estimate)
  rs <- ratio_stats(sub$estimate[fin], sub$se[fin])
  cs <- cov_stats(sub$covered)
  data.frame(R = sum(fin), paired_bias = NA_real_, paired_bias_mcse = NA_real_,
             emp_sd = rs$emp_sd, mean_se = rs$mean_se, se_ratio = rs$se_ratio,
             coverage = cs$coverage, coverage_mcse = cs$coverage_mcse,
             coverage_pop = mean_or_na(sub$covered_pop))
}

summarise_mi <- function(sub, complete_sub) {
  m <- merge(sub, complete_sub[, c("rep", "estimate", "se")],
             by = "rep", suffixes = c("", "_complete"))
  ok <- is.finite(m$estimate) & is.finite(m$estimate_complete)
  R  <- sum(ok)
  paired <- m$estimate[ok] - m$estimate_complete[ok]
  rs <- ratio_stats(m$estimate, m$se)
  cs <- cov_stats(m$covered)
  data.frame(R = R,
             paired_bias = if (R > 0) mean(paired) else NA_real_,
             paired_bias_mcse = if (R > 1) stats::sd(paired) / sqrt(R) else NA_real_,
             emp_sd = rs$emp_sd, mean_se = rs$mean_se, se_ratio = rs$se_ratio,
             coverage = cs$coverage, coverage_mcse = cs$coverage_mcse,
             coverage_pop = mean_or_na(m$covered_pop))
}

conv_summary <- function(regime_id, method) {
  sub <- if (is.null(all_diag)) NULL else
    all_diag[all_diag$regime_id == regime_id & all_diag$method == method, ]
  # Denominator is the EXPECTED rep count: a missing rep file is a
  # non-converged fit (D3).
  data.frame(n_fits = n_reps,
             n_converged = if (is.null(sub)) 0L else length(unique(sub$rep[sub$converged %in% TRUE])),
             median_max_rhat = if (is.null(sub) || !nrow(sub)) NA_real_ else stats::median(sub$max_rhat, na.rm = TRUE),
             median_min_ess  = if (is.null(sub) || !nrow(sub)) NA_real_ else stats::median(sub$min_ess, na.rm = TRUE))
}
no_conv <- data.frame(n_fits = NA_integer_, n_converged = NA_integer_,
                      median_max_rhat = NA_real_, median_min_ess = NA_real_)

keys <- do.call(rbind, lapply(seq_len(nrow(exp_regimes)), function(i) {
  expand.grid(regime_id = exp_regimes$regime_id[i],
              method = c("complete", mi_gls_v2_methods(exp_regimes$missing[i])),
              downstream = c("gls", "phylolm"),
              KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
}))

rows <- lapply(seq_len(nrow(keys)), function(i) {
  k <- keys[i, ]
  sub_all <- all_cells[all_cells$regime_id == k$regime_id & all_cells$method == k$method &
                       all_cells$downstream == k$downstream, ]
  n_present <- length(unique(sub_all$rep))
  n_finite  <- sum(is.finite(sub_all$estimate))
  fit_failure_rate <- 1 - n_finite / n_reps   # missing files count as failures (D3)
  complete_sub <- all_cells[all_cells$regime_id == k$regime_id & all_cells$method == "complete" &
                            all_cells$downstream == k$downstream, ]
  tr <- truth_rows[truth_rows$regime_id == k$regime_id & truth_rows$downstream == k$downstream, ]
  if (k$method == "complete") {
    s <- summarise_complete(sub_all)
    cv <- no_conv
  } else {
    keep <- converged_reps(k$regime_id, k$method)
    s <- summarise_mi(sub_all[sub_all$rep %in% keep, , drop = FALSE], complete_sub)
    cv <- conv_summary(k$regime_id, k$method)
  }
  cbind(k, n_expected = n_reps, n_present = n_present, n_finite = n_finite,
        fit_failure_rate = fit_failure_rate, truth = tr$truth, truth_source = tr$truth_source,
        s, cv)
})
summary_df <- do.call(rbind, rows)

# Complete-data reference beside each row (D2); 04 reads the complete rows.
cref <- summary_df[summary_df$method == "complete",
                   c("regime_id", "downstream", "se_ratio", "coverage", "coverage_pop")]
names(cref) <- c("regime_id", "downstream", "complete_se_ratio", "complete_coverage",
                 "complete_coverage_pop")
summary_df <- merge(summary_df, cref, by = c("regime_id", "downstream"), all.x = TRUE)
summary_df$code_sha <- code_sha_label
summary_df$dgp   <- regimes$dgp[match(summary_df$regime_id, regimes$regime_id)]
summary_df$gated <- regimes$gated[match(summary_df$regime_id, regimes$regime_id)]
method_order <- c("complete", "posterior_full", "posterior_none")
summary_df <- summary_df[order(summary_df$regime_id, match(summary_df$method, method_order),
                               summary_df$downstream), ]
lead <- c("regime_id", "method", "downstream", "dgp", "gated")
summary_df <- summary_df[, c(lead, setdiff(names(summary_df), lead))]
rownames(summary_df) <- NULL

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(summary_df, out_csv, row.names = FALSE)
cat("wrote", out_csv, "(", nrow(summary_df), "rows )\n")

if (!is.na(out_md)) {
  fmt <- function(x, digits = 4) ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
  md <- c("# Posterior MI simulation summary", "",
         sprintf("Generated from %d rep file(s) in `%s`; expected grid %s.",
                 nrow(present), indir, ex$label),
         sprintf("Code SHA: %s. pigauto version: %s.", code_sha_label, paste(vers, collapse = ", ")),
         "",
         paste0("| ", paste(names(summary_df), collapse = " | "), " |"),
         paste0("|", paste(rep("---", ncol(summary_df)), collapse = "|"), "|"))
  for (i in seq_len(nrow(summary_df))) {
    r <- summary_df[i, ]
    md <- c(md, paste0("| ", paste(vapply(r, function(v)
      if (is.numeric(v)) fmt(v) else as.character(v), ""), collapse = " | "), " |"))
  }
  dir.create(dirname(out_md), recursive = TRUE, showWarnings = FALSE)
  writeLines(md, out_md)
  cat("wrote", out_md, "\n")
}

# ---- (b) per-cell coverage table (regime x method x trait x mechanism) -------
if (!is.na(out_cell_csv)) {
  ck <- do.call(rbind, lapply(seq_len(nrow(exp_regimes)), function(i) {
    expand.grid(regime_id = exp_regimes$regime_id[i],
                method = c(mi_gls_v2_methods(exp_regimes$missing[i]), "conformal"),
                trait = mi_gls_v2_traits(exp_regimes$missing[i]),
                KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  }))
  cell_rows <- lapply(seq_len(nrow(ck)), function(i) {
    k <- ck[i, ]
    is_post <- k$method %in% c("posterior_full", "posterior_none")
    # Posterior methods: converged reps only (B3). Conformal: every present rep.
    keep <- if (is_post) converged_reps(k$regime_id, k$method) else
      present$rep[present$regime_id == k$regime_id]
    sub <- if (is.null(all_cellD)) NULL else
      all_cellD[all_cellD$regime_id == k$regime_id & all_cellD$method == k$method &
                all_cellD$trait == k$trait & all_cellD$rep %in% keep, ]
    n <- if (is.null(sub)) 0L else nrow(sub)
    reg_row <- regimes[regimes$regime_id == k$regime_id, ]
    cbind(k, mechanism = reg_row$mechanism, dgp = reg_row$dgp, gated = reg_row$gated, n = n,
          covered_sum = if (n) sum(sub$covered) else NA_integer_,   # NA propagates (fail-closed)
          mean_width = if (n) mean(sub$width) else NA_real_,
          n_expected = n_reps,
          n_reps_scored = if (n) length(unique(sub$rep)) else 0L,
          if (is_post) conv_summary(k$regime_id, k$method) else no_conv)
  })
  cell_cov_df <- do.call(rbind, cell_rows)
  cell_cov_df <- cell_cov_df[order(cell_cov_df$regime_id,
                                   match(cell_cov_df$method, c("posterior_full", "posterior_none", "conformal")),
                                   cell_cov_df$trait), ]
  rownames(cell_cov_df) <- NULL
  dir.create(dirname(out_cell_csv), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(cell_cov_df, out_cell_csv, row.names = FALSE)
  cat("wrote", out_cell_csv, "(", nrow(cell_cov_df), "rows )\n")
}
