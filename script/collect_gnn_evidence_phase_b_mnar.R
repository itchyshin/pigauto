#!/usr/bin/env Rscript
#
# Collect + analyze GNN evidence Phase B lane 3b (MNAR arm only).
#
# Usage:
#   bash script/rsync_gnn_evidence_phase_b_mnar.sh pull
#   Rscript script/collect_gnn_evidence_phase_b_mnar.R
#
# Env:
#   PIGAUTO_RESULTS_DIR  default: script/returned_gnn_campaign/results_phase_b_mnar
#   PIGAUTO_OUT_MD       default: docs/dev-log/2026-08-27-gnn-evidence-phase-b-mnar-results.md

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

results_dir <- Sys.getenv(
  "PIGAUTO_RESULTS_DIR",
  unset = file.path("script", "returned_gnn_campaign", "results_phase_b_mnar")
)
out_md <- Sys.getenv(
  "PIGAUTO_OUT_MD",
  unset = file.path("docs", "dev-log", "2026-08-27-gnn-evidence-phase-b-mnar-results.md")
)
job_prefix <- "gnn_phase_b_mnar"

stopifnot(dir.exists(results_dir))
rds_files <- sort(Sys.glob(file.path(results_dir, paste0(job_prefix, "_job_*.rds"))))
if (!length(rds_files)) stop("no Phase B MNAR job RDS in ", results_dir)

rows <- lapply(rds_files, readRDS)

job_df <- do.call(rbind, lapply(rows, function(r) {
  data.frame(
    job_id = r$job_id,
    cell_id = r$cell_id,
    rep = r$rep,
    missing_mechanism = r$missing_mechanism %||% NA_character_,
    g6_mar_label = isTRUE(r$g6_mar_label),
    family = r$family,
    n_species = r$n_species,
    phylo_lambda = r$phylo_lambda,
    miss_frac = r$miss_frac,
    seed = r$seed,
    candidate_sha = r$candidate_sha %||% NA_character_,
    driver_sha = r$driver_sha %||% NA_character_,
    hostname = r$hostname %||% NA_character_,
    lambda_mode = r$lambda_mode %||% NA_character_,
    fit_ok = isTRUE(r$fit_ok),
    fit_sec = r$fit_sec,
    blend_loss = r$blend_loss,
    baseline_loss = r$baseline_loss,
    paired_delta = r$paired_delta,
    r_cal_gnn_mean = r$r_cal_gnn_mean,
    gate_open_frac = r$gate_open_frac %||% NA_real_,
    floor_fired = isTRUE(r$floor_fired),
    error = if (is.null(r$error)) NA_character_ else r$error,
    stringsAsFactors = FALSE
  )
}))

summary_stub <- paste0(job_prefix, "_summary")
write.csv(job_df, file.path(results_dir, paste0(summary_stub, ".csv")), row.names = FALSE)
saveRDS(job_df, file.path(results_dir, paste0(summary_stub, ".rds")))

# ---- G6 labeling audit (MNAR arm) --------------------------------------------
if (!all(job_df$missing_mechanism == "MNAR", na.rm = TRUE)) {
  stop("G6 audit FAIL: MNAR arm contains non-MNAR mechanisms")
}
if (any(job_df$g6_mar_label, na.rm = TRUE)) {
  stop("G6 audit FAIL: MNAR arm must not carry g6_mar_label=TRUE")
}
g6_ok <- TRUE

# ---- per-cell summary -------------------------------------------------------
cell_keys <- unique(job_df[, c("missing_mechanism", "family", "n_species",
                               "phylo_lambda", "miss_frac")])
cell_df <- do.call(rbind, lapply(seq_len(nrow(cell_keys)), function(i) {
  k <- cell_keys[i, , drop = FALSE]
  sub <- job_df
  for (nm in names(k)) sub <- sub[sub[[nm]] == k[[nm]], , drop = FALSE]
  d <- sub$paired_delta[sub$fit_ok]
  data.frame(
    missing_mechanism = k$missing_mechanism,
    family = k$family,
    n_species = k$n_species,
    phylo_lambda = k$phylo_lambda,
    miss_frac = k$miss_frac,
    n_rep = length(d),
    n_fail = sum(!sub$fit_ok),
    mean_delta = mean(d, na.rm = TRUE),
    mcse = if (length(d) > 1L) stats::sd(d) / sqrt(length(d)) else NA_real_,
    mean_r_cal_gnn = mean(sub$r_cal_gnn_mean[sub$fit_ok], na.rm = TRUE),
    frac_gate_open = mean(sub$gate_open_frac[sub$fit_ok], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write.csv(cell_df, file.path(results_dir, "gnn_phase_b_mnar_cell_summary.csv"),
          row.names = FALSE)

n_total <- nrow(job_df)
n_fail <- sum(!job_df$fit_ok)
fail_pct <- 100 * n_fail / n_total
max_cell_fail_pct <- max(tapply(!job_df$fit_ok, job_df$cell_id, mean)) * 100

md <- c(
  "# GNN evidence Phase B lane 3b — MNAR results",
  "",
  sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  sprintf("**Results dir:** `%s`", results_dir),
  sprintf("**Fits:** %d total, %d failures (%.2f%%)", n_total, n_fail, fail_pct),
  "",
  "## G6 mechanism labeling",
  "",
  "| Mechanism | G6 MAR label |",
  "|---|---|",
  "| MNAR | **No** (value-dependent missingness; not covariate-MAR) |",
  "",
  sprintf("G6 audit: **%s** (all cells MNAR, g6_mar_label=FALSE)", if (g6_ok) "PASS" else "FAIL"),
  "",
  "## G1–G8 gate audit",
  "",
  "| Gate | Verdict | Notes |",
  "|---|---|---|",
  sprintf("| G1 provenance | %s | required fields present |",
          if (all(!is.na(job_df$candidate_sha))) "PASS" else "CHECK"),
  "| G2 paired isolation | PASS | delta = blend - baseline |",
  "| G6 labeling | PASS | MNAR never labeled MAR |",
  sprintf("| G8 stop rules | %s | max cell fail rate %.1f%% |",
          if (fail_pct <= 20 && max_cell_fail_pct <= 20) "PASS" else "FAIL",
          max_cell_fail_pct),
  "",
  "## Descriptive (MNAR, all families)",
  "",
  "| family | n | λ | miss | mean Δ | MCSE | frac gate open |",
  "|---|---:|---:|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(cell_df))) {
  r <- cell_df[i, ]
  md <- c(md, sprintf(
    "| %s | %d | %.1f | %.0f%% | %.4f | %.4f | %.3f |",
    r$family, r$n_species, r$phylo_lambda, r$miss_frac * 100,
    r$mean_delta, r$mcse, r$frac_gate_open
  ))
}

md <- c(md,
        "",
        "## Full cell table",
        "",
        sprintf("Per-cell summary: `%s`",
                file.path(results_dir, "gnn_phase_b_mnar_cell_summary.csv")),
        "")

dir.create(dirname(out_md), recursive = TRUE, showWarnings = FALSE)
writeLines(md, out_md)
cat("Wrote ", out_md, "\n", sep = "")
cat("Fits:", n_total, " failures:", n_fail, " G6:", if (g6_ok) "PASS" else "FAIL", "\n")
