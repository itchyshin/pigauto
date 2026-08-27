#!/usr/bin/env Rscript
#
# Collect AVONET300 multi-seed panel results with MCSE across seeds.
#
# Usage:
#   Rscript script/collect_gnn_evidence_avonet_panel.R
#
# Env:
#   PIGAUTO_AVONET_DIR   default: script/returned_gnn_campaign/results_avonet_panel
#   PIGAUTO_OUT_MD       default: script/returned_gnn_campaign/AVONET_PANEL_SUMMARY.md

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

out_dir <- Sys.getenv(
  "PIGAUTO_AVONET_DIR",
  unset = file.path("script", "returned_gnn_campaign", "results_avonet_panel")
)
out_md <- Sys.getenv(
  "PIGAUTO_OUT_MD",
  unset = file.path("script", "returned_gnn_campaign", "AVONET_PANEL_SUMMARY.md")
)

stopifnot(dir.exists(out_dir))
rds_files <- sort(Sys.glob(file.path(out_dir, "avonet_panel_job_*.rds")))
if (!length(rds_files)) stop("no panel RDS in ", out_dir)

rows <- lapply(rds_files, readRDS)

job_df <- do.call(rbind, lapply(rows, function(r) {
  data.frame(
    job_id = r$job_id,
    seed = r$seed,
    method = r$method,
    fit_ok = isTRUE(r$fit_ok),
    fit_sec = r$fit_sec,
    lambda_mode = r$lambda_mode %||% NA_character_,
    r_cal_gnn_mean = r$r_cal_gnn_mean %||% NA_real_,
    gate_open_frac = r$gate_open_frac %||% NA_real_,
    error = if (is.null(r$error)) NA_character_ else r$error,
    stringsAsFactors = FALSE
  )
}))

eval_df <- do.call(rbind, lapply(rows, function(r) {
  if (is.null(r$metrics) || !nrow(r$metrics)) return(NULL)
  d <- r$metrics
  d$job_id <- r$job_id
  d$seed <- r$seed
  d$panel_method <- r$method
  d
}))
if (is.null(eval_df)) stop("no metrics rows in panel RDS")
rownames(eval_df) <- NULL

fmt_num <- function(x, digits = 4L) {
  if (length(x) != 1L || is.na(x)) "NA" else format(round(x, digits), nsmall = digits)
}

summarize_metric <- function(d, method, trait, metric) {
  sub <- d[d$method == method & d$trait == trait & d$metric == metric, , drop = FALSE]
  if (!nrow(sub)) return(NULL)
  vals <- sub$value
  n <- length(vals)
  m <- mean(vals, na.rm = TRUE)
  sd_v <- if (n > 1L) stats::sd(vals, na.rm = TRUE) else NA_real_
  mcse <- if (n > 1L) sd_v / sqrt(n) else NA_real_
  data.frame(
    method = method,
    trait = trait,
    type = sub$type[1L],
    metric = metric,
    n_seeds = n,
    mean = m,
    sd = sd_v,
    mcse = mcse,
    stringsAsFactors = FALSE
  )
}

metric_keys <- unique(eval_df[, c("method", "trait", "metric", "type")])
metric_summ <- do.call(rbind, lapply(seq_len(nrow(metric_keys)), function(i) {
  summarize_metric(
    eval_df,
    metric_keys$method[i],
    metric_keys$trait[i],
    metric_keys$metric[i]
  )
}))
rownames(metric_summ) <- NULL

write.csv(job_df, file.path(out_dir, "avonet_panel_job_summary.csv"), row.names = FALSE)
write.csv(eval_df, file.path(out_dir, "avonet_panel_eval_long.csv"), row.names = FALSE)
write.csv(metric_summ, file.path(out_dir, "avonet_panel_metric_summary.csv"), row.names = FALSE)

# Continuous bayes sensitivity: pigauto_bayes vs pigauto_fixed1 RMSE
cont_traits <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")
cont_rmse <- metric_summ[
  metric_summ$trait %in% cont_traits & metric_summ$metric == "rmse",
, drop = FALSE]
bayes_sens <- NULL
for (tr in cont_traits) {
  fix <- cont_rmse[cont_rmse$trait == tr & cont_rmse$method == "pigauto_fixed1", , drop = FALSE]
  bay <- cont_rmse[cont_rmse$trait == tr & cont_rmse$method == "pigauto_bayes", , drop = FALSE]
  if (nrow(fix) && nrow(bay)) {
    bayes_sens <- rbind(bayes_sens, data.frame(
      trait = tr,
      rmse_fixed_1 = fix$mean,
      rmse_bayes = bay$mean,
      delta_bayes_minus_fixed = bay$mean - fix$mean,
      stringsAsFactors = FALSE
    ))
  }
}

# Head-to-head vs baseline_pigauto_fixed1 (Rphylopars joint path when available)
paired_compare <- function(method_a, method_b, metric_name = "rmse") {
  traits <- unique(eval_df$trait)
  out <- list()
  for (tr in traits) {
    tp <- eval_df$type[eval_df$trait == tr][1L]
    met <- if (tp == "continuous") metric_name else "accuracy"
    a <- metric_summ[metric_summ$method == method_a & metric_summ$trait == tr &
                       metric_summ$metric == met, , drop = FALSE]
    b <- metric_summ[metric_summ$method == method_b & metric_summ$trait == tr &
                       metric_summ$metric == met, , drop = FALSE]
    if (nrow(a) && nrow(b)) {
      out[[length(out) + 1L]] <- data.frame(
        trait = tr,
        type = tp,
        metric = met,
        method_a = method_a,
        method_b = method_b,
        mean_a = a$mean,
        mean_b = b$mean,
        delta_b_minus_a = b$mean - a$mean,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(out)) do.call(rbind, out) else NULL
}

h2h_blend_vs_base <- paired_compare("pigauto_fixed1", "baseline_pigauto_fixed1")
h2h_bace_vs_blend <- paired_compare("pigauto_fixed1", "bace")

meta <- rows[[1L]]
n_fail <- sum(!job_df$fit_ok)
md <- c(
  "# AVONET300 multi-seed external validation panel",
  "",
  sprintf("Generated: %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  sprintf("Candidate SHA: `%s`", meta$candidate_sha %||% "unknown"),
  sprintf("Driver SHA: `%s`", meta$driver_sha %||% "unknown"),
  sprintf("Dataset: avonet300 (%d species, %d traits), miss_frac = %.0f%% MCAR (paired per seed).",
          meta$n_species %||% 300L, meta$n_traits %||% 7L, 100 * (meta$miss_frac %||% 0.30)),
  sprintf("Seeds: %s", paste(sort(unique(job_df$seed)), collapse = ", ")),
  sprintf("Jobs: %d total, %d failures", nrow(job_df), n_fail),
  "",
  "**Framing:** Real-data corroboration only — not pre-registered in Phase A.",
  "BACE may win on continuous morphometrics; pigauto's claim is unified mixed-type",
  "imputation with calibrated gates. `baseline_pigauto_fixed1` is pigauto's internal",
  "phylogenetic baseline (Rphylopars joint path when installed).",
  "",
  "## Per-trait metrics (mean ± MCSE across seeds)",
  "",
  "| Method | Trait | Type | Metric | Mean | MCSE | n |",
  "|---|---|---|---|---:|---:|---:|"
)
for (i in seq_len(nrow(metric_summ))) {
  r <- metric_summ[i, ]
  md <- c(md, sprintf(
    "| %s | %s | %s | %s | %s | %s | %d |",
    r$method, r$trait, r$type, r$metric,
    fmt_num(r$mean), fmt_num(r$mcse), r$n_seeds
  ))
}

if (!is.null(bayes_sens) && nrow(bayes_sens)) {
  md <- c(md, "",
    "## Bayes λ sensitivity (continuous traits: pigauto_bayes vs pigauto_fixed1)",
    "",
    "| Trait | RMSE fixed_1 | RMSE bayes | Δ (bayes − fixed) |",
    "|---|---:|---:|---:|"
  )
  for (i in seq_len(nrow(bayes_sens))) {
    r <- bayes_sens[i, ]
    md <- c(md, sprintf(
      "| %s | %s | %s | %s |",
      r$trait, fmt_num(r$rmse_fixed_1), fmt_num(r$rmse_bayes),
      fmt_num(r$delta_bayes_minus_fixed)
    ))
  }
}

if (!is.null(h2h_blend_vs_base) && nrow(h2h_blend_vs_base)) {
  md <- c(md, "",
    "## pigauto blend vs internal baseline (fixed_1)",
    "",
    "| Trait | Type | Metric | Blend | Baseline | Δ (blend − base) |",
    "|---|---|---|---:|---:|---:|"
  )
  for (i in seq_len(nrow(h2h_blend_vs_base))) {
    r <- h2h_blend_vs_base[i, ]
    md <- c(md, sprintf(
      "| %s | %s | %s | %s | %s | %s |",
      r$trait, r$type, r$metric,
      fmt_num(r$mean_a), fmt_num(r$mean_b), fmt_num(r$delta_b_minus_a)
    ))
  }
}

if (!is.null(h2h_bace_vs_blend) && nrow(h2h_bace_vs_blend)) {
  md <- c(md, "",
    "## BACE vs pigauto_fixed1 (held-out cells)",
    "",
    "| Trait | Type | Metric | pigauto | BACE | Δ (BACE − pigauto) |",
    "|---|---|---|---:|---:|---:|"
  )
  for (i in seq_len(nrow(h2h_bace_vs_blend))) {
    r <- h2h_bace_vs_blend[i, ]
    md <- c(md, sprintf(
      "| %s | %s | %s | %s | %s | %s |",
      r$trait, r$type, r$metric,
      fmt_num(r$mean_a), fmt_num(r$mean_b), fmt_num(r$delta_b_minus_a)
    ))
  }
}

md <- c(md, "",
  "## Wall time by method (mean seconds per job)",
  "",
  "```",
  capture.output(print(aggregate(fit_sec ~ method, job_df[job_df$fit_ok, ], mean)),
                 row.names = FALSE),
  "```",
  "",
  sprintf("Artifacts: `%s`", out_dir)
)

dir.create(dirname(out_md), recursive = TRUE, showWarnings = FALSE)
writeLines(md, out_md)
cat("Wrote ", out_md, "\n", sep = "")
cat("Jobs:", nrow(job_df), " failures:", n_fail, " metric rows:", nrow(metric_summ), "\n")
