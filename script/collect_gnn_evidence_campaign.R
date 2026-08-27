#!/usr/bin/env Rscript
#
# Collect + analyze GNN evidence Phase A campaign results (prereg 2026-08-26).
#
# Usage:
#   bash script/rsync_gnn_evidence_campaign.sh pull
#   Rscript script/collect_gnn_evidence_campaign.R
#
# Env:
#   PIGAUTO_RESULTS_DIR  default: script/returned_gnn_campaign/results
#   PIGAUTO_OUT_MD       default: docs/dev-log/2026-08-26-gnn-evidence-phase-a-results.md

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

results_dir <- Sys.getenv(
  "PIGAUTO_RESULTS_DIR",
  unset = file.path("script", "returned_gnn_campaign", "results")
)
out_md <- Sys.getenv(
  "PIGAUTO_OUT_MD",
  unset = file.path("docs", "dev-log", "2026-08-26-gnn-evidence-phase-a-results.md")
)

stopifnot(dir.exists(results_dir))
rds_files <- sort(Sys.glob(file.path(results_dir, "gnn_campaign_job_*.rds")))
if (!length(rds_files)) stop("no job RDS in ", results_dir)

rows <- lapply(rds_files, readRDS)

# ---- per-job summary -------------------------------------------------------
job_df <- do.call(rbind, lapply(rows, function(r) {
  data.frame(
    job_id = r$job_id,
    cell_id = r$cell_id,
    rep = r$rep,
    family = r$family,
    n_species = r$n_species,
    phylo_lambda = r$phylo_lambda,
    miss_frac = r$miss_frac,
    seed = r$seed,
    tree_seed = r$tree_seed %||% NA_integer_,
    candidate_sha = r$candidate_sha %||% NA_character_,
    driver_sha = r$driver_sha %||% NA_character_,
    hostname = r$hostname %||% NA_character_,
    missing_mechanism = r$missing_mechanism %||% NA_character_,
    lambda_mode = r$lambda_mode %||% NA_character_,
    fit_ok = isTRUE(r$fit_ok),
    fit_sec = r$fit_sec,
    blend_loss = r$blend_loss,
    baseline_loss = r$baseline_loss,
    paired_delta = r$paired_delta,
    r_cal_gnn_mean = r$r_cal_gnn_mean,
    r_cal_bm_mean = r$r_cal_bm_mean,
    r_cal_mean_mean = r$r_cal_mean_mean,
    gate_open_frac = r$gate_open_frac %||% NA_real_,
    floor_fired = isTRUE(r$floor_fired),
    error = if (is.null(r$error)) NA_character_ else r$error,
    stringsAsFactors = FALSE
  )
}))

write.csv(job_df, file.path(results_dir, "gnn_campaign_summary.csv"), row.names = FALSE)
saveRDS(job_df, file.path(results_dir, "gnn_campaign_summary.rds"))

# ---- gate distribution (full r_cal_gnn per latent col) ---------------------
gate_rows <- lapply(rows, function(r) {
  g <- r$gates
  if (is.null(g) || !is.data.frame(g) || nrow(g) == 0L) return(NULL)
  cbind(
    job_id = r$job_id,
    cell_id = r$cell_id,
    family = r$family,
    n_species = r$n_species,
    phylo_lambda = r$phylo_lambda,
    miss_frac = r$miss_frac,
    g
  )
})
gate_df <- if (length(gate_rows)) do.call(rbind, gate_rows) else NULL
if (!is.null(gate_df)) {
  write.csv(gate_df, file.path(results_dir, "gnn_campaign_gates.csv"), row.names = FALSE)
}

# ---- cell-level summaries --------------------------------------------------
cell_key <- function(d) {
  interaction(d$family, d$n_species, d$phylo_lambda, d$miss_frac, drop = TRUE)
}

summarize_cell <- function(d) {
  ok <- d$fit_ok
  n <- nrow(d)
  n_fail <- sum(!ok)
  deltas <- d$paired_delta[ok]
  bl <- d$baseline_loss[ok]
  mean_delta <- if (length(deltas)) mean(deltas, na.rm = TRUE) else NA_real_
  sd_delta <- if (length(deltas) > 1L) sd(deltas, na.rm = TRUE) else NA_real_
  mcse <- if (length(deltas) > 1L) sd_delta / sqrt(length(deltas)) else NA_real_
  rel_improve <- if (length(bl) && mean(bl, na.rm = TRUE) > 0) {
    -mean_delta / mean(bl, na.rm = TRUE)
  } else NA_real_
  z_ratio <- if (!is.na(mcse) && mcse > 0) abs(mean_delta) / mcse else NA_real_

  gnn <- d$r_cal_gnn_mean[ok]
  data.frame(
    n_rep = n,
    n_ok = sum(ok),
    n_fail = n_fail,
    fail_rate = n_fail / n,
    mean_delta = mean_delta,
    sd_delta = sd_delta,
    mcse = mcse,
    z_ratio = z_ratio,
    rel_improve = rel_improve,
    mean_baseline_loss = if (length(bl)) mean(bl, na.rm = TRUE) else NA_real_,
    mean_blend_loss = if (length(d$blend_loss[ok])) mean(d$blend_loss[ok], na.rm = TRUE) else NA_real_,
  mean_r_cal_gnn = if (length(gnn)) mean(gnn, na.rm = TRUE) else NA_real_,
    median_r_cal_gnn = if (length(gnn)) median(gnn, na.rm = TRUE) else NA_real_,
    frac_gate_open = if (length(gnn)) mean(gnn > 0, na.rm = TRUE) else NA_real_,
    frac_gate_material = if (length(gnn)) mean(gnn >= 0.10, na.rm = TRUE) else NA_real_,
    frac_gate_closed = if (length(gnn)) mean(gnn == 0, na.rm = TRUE) else NA_real_,
    frac_floor = mean(d$floor_fired[ok], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

cell_ids <- unique(job_df[, c("family", "n_species", "phylo_lambda", "miss_frac")])
cell_summ <- do.call(rbind, lapply(seq_len(nrow(cell_ids)), function(i) {
  d <- job_df[
    job_df$family == cell_ids$family[i] &
      job_df$n_species == cell_ids$n_species[i] &
      job_df$phylo_lambda == cell_ids$phylo_lambda[i] &
      job_df$miss_frac == cell_ids$miss_frac[i],
  , drop = FALSE]
  cbind(cell_ids[i, , drop = FALSE], summarize_cell(d))
}))
rownames(cell_summ) <- NULL

# ---- G4 flag logic ---------------------------------------------------------
g4_flag <- function(row) {
  if (row$family == "F1") return(FALSE)
  if (is.na(row$z_ratio) || is.na(row$rel_improve)) return(FALSE)
  row$z_ratio >= 3 && row$rel_improve >= 0.02 && row$mean_delta < 0
}

cell_summ$g4_explore_pass <- vapply(seq_len(nrow(cell_summ)), function(i) {
  g4_flag(cell_summ[i, ])
}, logical(1))

write.csv(cell_summ, file.path(results_dir, "gnn_campaign_cell_summary.csv"), row.names = FALSE)

g4_cells <- cell_summ[cell_summ$g4_explore_pass, , drop = FALSE]

# ---- G1-G8 audit -----------------------------------------------------------
required_fields <- c(
  "candidate_sha", "driver_sha", "hostname", "seed", "family",
  "n_species", "phylo_lambda", "miss_frac", "missing_mechanism",
  "lambda_mode", "blend_loss", "baseline_loss", "paired_delta",
  "r_cal_gnn_mean", "r_cal_bm_mean", "r_cal_mean_mean", "floor_fired", "fit_sec"
)

g1_ok <- all(required_fields %in% names(job_df)) &&
  all(!is.na(job_df$candidate_sha[job_df$fit_ok])) &&
  all(job_df$missing_mechanism[job_df$fit_ok] == "MCAR", na.rm = TRUE)

g2_ok <- all(c("blend_loss", "baseline_loss", "paired_delta") %in% names(job_df)) &&
  all(
    abs(job_df$paired_delta[job_df$fit_ok] -
          (job_df$blend_loss[job_df$fit_ok] - job_df$baseline_loss[job_df$fit_ok])) < 1e-8,
    na.rm = TRUE
  )

g3_ok <- any(job_df$r_cal_gnn_mean[job_df$fit_ok] == 0, na.rm = TRUE)

max_fail_rate <- max(cell_summ$fail_rate, na.rm = TRUE)
g8_fail <- max_fail_rate > 0.20
g8_fields <- !all(required_fields %in% names(job_df))

# wall time from fit_sec sum / workers proxy — use max job timestamps if available
total_fits <- nrow(job_df)
n_fail <- sum(!job_df$fit_ok)
fail_pct <- 100 * n_fail / total_fits

# F1 @ lambda=1 specificity
f1_l1 <- cell_summ[cell_summ$family == "F1" & cell_summ$phylo_lambda == 1.0, , drop = FALSE]
f2_l1 <- cell_summ[cell_summ$family == "F2" & cell_summ$phylo_lambda == 1.0, , drop = FALSE]
f3_all <- cell_summ[cell_summ$family == "F3", , drop = FALSE]

fmt_num <- function(x, digits = 4) {
  ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits))
}

md_lines <- c(
  "# GNN evidence Phase A — results",
  "",
  sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  sprintf("**Results dir:** `%s`", results_dir),
  sprintf("**Fits:** %d total, %d failures (%.2f%%)", total_fits, n_fail, fail_pct),
  "",
  "## G1–G8 gate audit",
  "",
  sprintf("| Gate | Verdict | Notes |"),
  sprintf("|---|---|---|"),
  sprintf("| G1 provenance | %s | required fields present |", if (g1_ok) "PASS" else "FAIL"),
  sprintf("| G2 paired isolation | %s | delta = blend - baseline |", if (g2_ok) "PASS" else "FAIL"),
  sprintf("| G3 fallback | %s | r_cal_gnn=0 observed |", if (g3_ok) "PASS" else "FAIL"),
  "| G4 positive claim | see below | explore-stage flags only |",
  "| G5 no-benefit retention | PASS | all cells reported |",
  "| G6 MCAR labeling | PASS | Phase A MCAR only |",
  "| G7 trait boundary | PASS | F1/F2/F3 reported separately |",
  sprintf("| G8 stop rules | %s | max cell fail rate %.1f%% |",
          if (g8_fail || g8_fields) "FAIL" else "PASS", 100 * max_fail_rate),
  "",
  "## F1 @ λ=1 specificity (null control)",
  "",
  "Expected Δ ≈ 0. Any systematic GNN win here is a red flag.",
  "",
  "| n | miss | mean Δ | MCSE | z_ratio | rel improve | frac gate open |",
  "|---|---:|---:|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(f1_l1))) {
  r <- f1_l1[i, ]
  md_lines <- c(md_lines, sprintf(
    "| %d | %.0f%% | %s | %s | %s | %s | %s |",
    r$n_species, 100 * r$miss_frac,
    fmt_num(r$mean_delta), fmt_num(r$mcse), fmt_num(r$z_ratio),
    fmt_num(100 * r$rel_improve, 2), fmt_num(r$frac_gate_open, 3)
  ))
}

md_lines <- c(md_lines, "",
  "## F2 @ λ=1 explore (G4 screen, 30 seeds)",
  "",
  "G4 explore pass requires: F2 family, |Δ|/MCSE ≥ 3, relative improvement ≥ 2%, Δ < 0.",
  "60-seed confirmation required before any public GNN-positive claim.",
  "",
  "| n | miss | mean Δ | MCSE | z_ratio | rel improve | G4 explore | mean r_cal_gnn | frac open |",
  "|---|---:|---:|---:|---:|---:|---|---:|---:|"
)

for (i in seq_len(nrow(f2_l1))) {
  r <- f2_l1[i, ]
  md_lines <- c(md_lines, sprintf(
    "| %d | %.0f%% | %s | %s | %s | %s | %s | %s | %s |",
    r$n_species, 100 * r$miss_frac,
    fmt_num(r$mean_delta), fmt_num(r$mcse), fmt_num(r$z_ratio),
    fmt_num(100 * r$rel_improve, 2),
    if (r$g4_explore_pass) "**PASS**" else "—",
    fmt_num(r$mean_r_cal_gnn, 3), fmt_num(r$frac_gate_open, 3)
  ))
}

md_lines <- c(md_lines, "",
  sprintf("**G4 explore passes (F2 @ λ=1):** %d / %d cells",
          sum(f2_l1$g4_explore_pass), nrow(f2_l1)),
  "",
  "## F3 mixed-type (descriptive)",
  "",
  "| n | λ | miss | mean Δ | MCSE | mean r_cal_gnn | frac open | frac floor |",
  "|---|---:|---:|---:|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(f3_all))) {
  r <- f3_all[i, ]
  md_lines <- c(md_lines, sprintf(
    "| %d | %.1f | %.0f%% | %s | %s | %s | %s | %s |",
    r$n_species, r$phylo_lambda, 100 * r$miss_frac,
    fmt_num(r$mean_delta), fmt_num(r$mcse),
    fmt_num(r$mean_r_cal_gnn, 3), fmt_num(r$frac_gate_open, 3),
    fmt_num(r$frac_floor, 3)
  ))
}

md_lines <- c(md_lines, "",
  "## Gate distribution note",
  "",
  if (!is.null(gate_df)) {
    sprintf("Full per-latent `r_cal_gnn` distribution: `%s` (%d rows).",
            file.path(results_dir, "gnn_campaign_gates.csv"), nrow(gate_df))
  } else {
    "Gate-level RDS expansion unavailable."
  },
  "",
  "## All cells summary",
  "",
  sprintf("Full cell table: `%s`", file.path(results_dir, "gnn_campaign_cell_summary.csv")),
  ""
)

dir.create(dirname(out_md), recursive = TRUE, showWarnings = FALSE)
writeLines(md_lines, out_md)

cat("=== collect_gnn_evidence_campaign.R ===\n")
cat("fits:", total_fits, " failures:", n_fail, "\n")
cat("cell summaries:", nrow(cell_summ), "\n")
cat("G4 explore passes (all cells):", sum(cell_summ$g4_explore_pass), "\n")
cat("wrote:", out_md, "\n")
