#!/usr/bin/env Rscript
# Usage: Rscript 06b_summarise_mi_se.R <indir> [out.md]
#
# Aggregates the rep_<i>.rds files written by 06_mi_se_sim.R into the
# cross-rep metrics pre-registered in
# useful/mondrian-mi-se-justification.md: SE ratio (with a delta-method
# Monte Carlo SE), coverage, bias, mean FMI, and per-cell draw PIT /
# coverage by near/far stratum, for the split and Mondrian arms.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("expected: indir [out.md]", call. = FALSE)
indir  <- args[[1L]]
out_md <- if (length(args) >= 2L) args[[2L]] else file.path(indir, "summary.md")

files <- list.files(indir, pattern = "^rep_\\d+\\.rds$", full.names = TRUE)
if (!length(files)) stop("no rep_*.rds files found in ", indir, call. = FALSE)

reps <- lapply(files, readRDS)
true_beta <- reps[[1]]$true_beta

# ---- SE ratio: mean pooled total SE / empirical SD of the pooled estimate,
# across reps, for one arm. Delta-method MCSE on the ratio via the standard
# first-order approximation Var(A/B) ~= (1/B^2) Var(A) + (A^2/B^4) Var(B)
# for A = mean(SE), B = sd(estimate), treating reps as iid.
se_ratio_summary <- function(arm_name) {
  ok <- vapply(reps, function(r) !isTRUE(r[[arm_name]]$failed), logical(1))
  est <- vapply(reps[ok], function(r) r[[arm_name]]$estimate, numeric(1))
  se  <- vapply(reps[ok], function(r) r[[arm_name]]$std.error, numeric(1))
  fmi <- vapply(reps[ok], function(r) r[[arm_name]]$fmi %||% NA_real_, numeric(1))

  n_ok <- sum(ok)
  mean_se <- mean(se, na.rm = TRUE)
  sd_est  <- stats::sd(est, na.rm = TRUE)
  se_ratio <- mean_se / sd_est

  # delta-method MCSE
  var_mean_se <- stats::var(se, na.rm = TRUE) / n_ok
  var_sd_est  <- stats::var(est, na.rm = TRUE)^2 * 2 / (n_ok - 1) / (4 * sd_est^2)
  mcse <- sqrt(var_mean_se / sd_est^2 +
              mean_se^2 / sd_est^4 * var_sd_est)

  bias <- mean(est, na.rm = TRUE) - true_beta
  # nominal 95% CI coverage of the true slope, using the pooled t reference
  df  <- vapply(reps[ok], function(r) r[[arm_name]]$df %||% NA_real_, numeric(1))
  tcrit <- stats::qt(0.975, df)
  covered <- abs(est - true_beta) <= tcrit * se
  coverage <- mean(covered, na.rm = TRUE)

  list(n_ok = n_ok, n_failed = length(reps) - n_ok,
       mean_se = mean_se, sd_est = sd_est, se_ratio = se_ratio, mcse = mcse,
       bias = bias, coverage95 = coverage, mean_fmi = mean(fmi, na.rm = TRUE))
}
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- per-cell draw PIT / coverage by stratum, for one arm -----------------
stratum_summary <- function(arm_name) {
  cells <- do.call(rbind, lapply(reps, function(r) {
    a <- r[[arm_name]]
    if (isTRUE(a$failed) || is.null(a$cell)) return(NULL)
    a$cell
  }))
  if (is.null(cells) || !nrow(cells)) {
    return(data.frame(stratum = character(0), n = integer(0),
                      mean_pit = numeric(0), coverage = numeric(0)))
  }
  cells$stratum[is.na(cells$stratum)] <- "unknown"
  agg <- stats::aggregate(cbind(pit, covered) ~ stratum, data = cells,
                          FUN = mean)
  n_tab <- table(cells$stratum)
  agg$n <- as.integer(n_tab[agg$stratum])
  agg[, c("stratum", "n", "pit", "covered")]
}

split_se    <- se_ratio_summary("split")
mondrian_se <- se_ratio_summary("mondrian")
split_strat    <- stratum_summary("split")
mondrian_strat <- stratum_summary("mondrian")

fmt <- function(x, d = 4) formatC(x, digits = d, format = "f")

lines <- c(
  "# Mondrian vs split multi_impute() Rubin-SE campaign summary",
  "",
  sprintf("Reps found: %d. True slope (rho): %s.", length(reps), fmt(true_beta)),
  "",
  "## SE ratio (mean pooled SE / empirical SD of pooled estimate)",
  "",
  "| arm | n_ok | n_failed | mean SE | empirical SD | SE ratio | MCSE(ratio) | bias | 95% CI coverage | mean FMI |",
  "|---|---|---|---|---|---|---|---|---|---|",
  sprintf("| split | %d | %d | %s | %s | %s | %s | %s | %s | %s |",
          split_se$n_ok, split_se$n_failed, fmt(split_se$mean_se),
          fmt(split_se$sd_est), fmt(split_se$se_ratio), fmt(split_se$mcse),
          fmt(split_se$bias), fmt(split_se$coverage95), fmt(split_se$mean_fmi)),
  sprintf("| mondrian | %d | %d | %s | %s | %s | %s | %s | %s | %s |",
          mondrian_se$n_ok, mondrian_se$n_failed, fmt(mondrian_se$mean_se),
          fmt(mondrian_se$sd_est), fmt(mondrian_se$se_ratio), fmt(mondrian_se$mcse),
          fmt(mondrian_se$bias), fmt(mondrian_se$coverage95), fmt(mondrian_se$mean_fmi)),
  "",
  sprintf("Pre-registered claim: |SE ratio - 1| smaller for mondrian than split by >= 0.03. Observed: split=%s, mondrian=%s, delta=%s.",
          fmt(abs(split_se$se_ratio - 1)), fmt(abs(mondrian_se$se_ratio - 1)),
          fmt(abs(split_se$se_ratio - 1) - abs(mondrian_se$se_ratio - 1))),
  "",
  "## Per-cell draw PIT / coverage by stratum -- split arm",
  "",
  "| stratum | n | mean PIT | coverage |",
  "|---|---|---|---|",
  if (nrow(split_strat))
    sprintf("| %s | %d | %s | %s |", split_strat$stratum, split_strat$n,
            fmt(split_strat$pit), fmt(split_strat$covered))
  else "| (no cells) | | | |",
  "",
  "## Per-cell draw PIT / coverage by stratum -- mondrian arm",
  "",
  "| stratum | n | mean PIT | coverage |",
  "|---|---|---|---|",
  if (nrow(mondrian_strat))
    sprintf("| %s | %d | %s | %s |", mondrian_strat$stratum, mondrian_strat$n,
            fmt(mondrian_strat$pit), fmt(mondrian_strat$covered))
  else "| (no cells) | | | |",
  "",
  sprintf("Pre-registered claim: far-stratum coverage in [0.93, 0.97]. Observed (mondrian arm, far): %s.",
          {
            far <- mondrian_strat[mondrian_strat$stratum == "far", "covered"]
            if (length(far)) fmt(far) else "NA (no far-stratum cells)"
          })
)

writeLines(lines, out_md)
cat(paste(lines, collapse = "\n"), "\n")
cat("\nWrote", out_md, "\n")
