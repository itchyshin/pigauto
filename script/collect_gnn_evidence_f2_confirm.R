#!/usr/bin/env Rscript
#
# G4 confirmatory audit — F2 @ λ=1, 60 seeds (prereg §3.4).
#
# Usage:
#   Rscript script/collect_gnn_evidence_f2_confirm.R
#
# Env:
#   PIGAUTO_CONFIRM_DIR   default: script/returned_gnn_campaign/results_confirm
#   PIGAUTO_PRIMARY_DIR   default: script/returned_gnn_campaign/results
#   PIGAUTO_OUT_MD        default: append section to phase-a-results md

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

confirm_dir <- Sys.getenv(
  "PIGAUTO_CONFIRM_DIR",
  unset = file.path("script", "returned_gnn_campaign", "results_confirm")
)
primary_dir <- Sys.getenv(
  "PIGAUTO_PRIMARY_DIR",
  unset = file.path("script", "returned_gnn_campaign", "results")
)
out_csv <- file.path(confirm_dir, "gnn_confirm_cell_summary.csv")

stopifnot(dir.exists(confirm_dir))
rds_files <- sort(Sys.glob(file.path(confirm_dir, "gnn_confirm_job_*.rds")))
if (!length(rds_files)) stop("no confirm RDS in ", confirm_dir)

rows <- lapply(rds_files, readRDS)

job_df <- do.call(rbind, lapply(rows, function(r) {
  data.frame(
    job_id = r$job_id,
    confirm_cell_id = r$confirm_cell_id,
    explore_cell_id = r$explore_cell_id,
    rep = r$rep,
    family = r$family,
    n_species = r$n_species,
    phylo_lambda = r$phylo_lambda,
    miss_frac = r$miss_frac,
    seed = r$seed,
    lambda_mode = r$lambda_mode,
    fit_ok = isTRUE(r$fit_ok),
    fit_sec = r$fit_sec,
    blend_loss = r$blend_loss,
    baseline_loss = r$baseline_loss,
    paired_delta = r$paired_delta,
    r_cal_gnn_mean = r$r_cal_gnn_mean,
    floor_fired = isTRUE(r$floor_fired),
    error = if (is.null(r$error)) NA_character_ else r$error,
    stringsAsFactors = FALSE
  )
}))

write.csv(job_df, file.path(confirm_dir, "gnn_confirm_job_summary.csv"), row.names = FALSE)

summarize_arm <- function(d) {
  ok <- d$fit_ok
  deltas <- d$paired_delta[ok]
  bl <- d$baseline_loss[ok]
  n_ok <- sum(ok)
  mean_delta <- if (length(deltas)) mean(deltas, na.rm = TRUE) else NA_real_
  sd_delta <- if (length(deltas) > 1L) sd(deltas, na.rm = TRUE) else NA_real_
  mcse <- if (length(deltas) > 1L) sd_delta / sqrt(length(deltas)) else NA_real_
  rel_improve <- if (length(bl) && mean(bl, na.rm = TRUE) > 0) {
    -mean_delta / mean(bl, na.rm = TRUE)
  } else NA_real_
  z_ratio <- if (!is.na(mcse) && mcse > 0) abs(mean_delta) / mcse else NA_real_
  gnn <- d$r_cal_gnn_mean[ok]
  data.frame(
    n_rep = nrow(d),
    n_ok = n_ok,
    n_fail = sum(!ok),
    mean_delta = mean_delta,
    sd_delta = sd_delta,
    mcse = mcse,
    z_ratio = z_ratio,
    rel_improve = rel_improve,
    mean_r_cal_gnn = if (length(gnn)) mean(gnn, na.rm = TRUE) else NA_real_,
    frac_gate_open = if (length(gnn)) mean(gnn > 0, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}

g4_pass <- function(mean_delta, z_ratio, rel_improve) {
  !is.na(z_ratio) && !is.na(rel_improve) &&
    z_ratio >= 3 && rel_improve >= 0.02 && mean_delta < 0
}

cell_ids <- unique(job_df[, c(
  "confirm_cell_id", "explore_cell_id", "family",
  "n_species", "phylo_lambda", "miss_frac"
)])

confirm_summ <- do.call(rbind, lapply(seq_len(nrow(cell_ids)), function(i) {
  d <- job_df[
    job_df$confirm_cell_id == cell_ids$confirm_cell_id[i],
  , drop = FALSE]
  s <- summarize_arm(d)
  cbind(cell_ids[i, , drop = FALSE], s)
}))
rownames(confirm_summ) <- NULL

confirm_summ$g4_confirm_pass <- mapply(
  g4_pass,
  confirm_summ$mean_delta,
  confirm_summ$z_ratio,
  confirm_summ$rel_improve
)

# ---- explore (30-seed primary) comparison ----------------------------------
primary_summ <- NULL
if (file.exists(file.path(primary_dir, "gnn_campaign_cell_summary.csv"))) {
  primary_summ <- read.csv(
    file.path(primary_dir, "gnn_campaign_cell_summary.csv"),
    stringsAsFactors = FALSE
  )
} else if (dir.exists(primary_dir)) {
  source_lines <- readLines(file.path("script", "collect_gnn_evidence_campaign.R"))
  invisible(source_lines) # not sourcing; re-read primary if needed
}

explore_cols <- c(
  "mean_delta", "mcse", "z_ratio", "rel_improve", "g4_explore_pass"
)
if (!is.null(primary_summ)) {
  for (i in seq_len(nrow(confirm_summ))) {
    cid <- confirm_summ$explore_cell_id[i]
    # explore cells are F2 @ lambda=1 matching n/miss
    pr <- primary_summ[
      primary_summ$family == "F2" &
        primary_summ$phylo_lambda == 1.0 &
        primary_summ$n_species == confirm_summ$n_species[i] &
        abs(primary_summ$miss_frac - confirm_summ$miss_frac[i]) < 1e-9,
    , drop = FALSE]
    if (nrow(pr) == 1L) {
      for (col in explore_cols) {
        confirm_summ[i, paste0("explore_", col)] <- pr[[col]]
      }
    }
  }
}

confirm_summ$delta_shift <- confirm_summ$mean_delta - confirm_summ$explore_mean_delta
confirm_summ$z_shift <- confirm_summ$z_ratio - confirm_summ$explore_z_ratio
confirm_summ$stable_sign <- sign(confirm_summ$mean_delta) == sign(confirm_summ$explore_mean_delta)

write.csv(confirm_summ, out_csv, row.names = FALSE)

fmt_num <- function(x, digits = 4) {
  ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits))
}

md_lines <- c(
  "",
  "## F2 @ λ=1 confirm (G4, 60 seeds)",
  "",
  sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  sprintf("**Results dir:** `%s`", confirm_dir),
  sprintf("**Fits:** %d total, %d failures", nrow(job_df), sum(!job_df$fit_ok)),
  "",
  "G4 confirm pass requires: F2 family, |Δ|/MCSE ≥ 3, relative improvement ≥ 2%, Δ < 0.",
  "These cells passed G4 explore (30 seeds) and were re-run at 60 seeds before any manuscript claim.",
  "",
  "| cell | n | miss | mean Δ | MCSE | z_ratio | rel improve | G4 confirm | explore z | explore rel | stable |",
  "|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---|"
)

for (i in seq_len(nrow(confirm_summ))) {
  r <- confirm_summ[i, ]
  md_lines <- c(md_lines, sprintf(
    "| %d | %d | %.0f%% | %s | %s | %s | %s | %s | %s | %s | %s |",
    r$confirm_cell_id, r$n_species, 100 * r$miss_frac,
    fmt_num(r$mean_delta), fmt_num(r$mcse), fmt_num(r$z_ratio),
    fmt_num(100 * r$rel_improve, 2),
    if (r$g4_confirm_pass) "**PASS**" else "FAIL",
    fmt_num(r$explore_z_ratio), fmt_num(100 * r$explore_rel_improve, 2),
    if (isTRUE(r$stable_sign)) "yes" else "no"
  ))
}

n_pass <- sum(confirm_summ$g4_confirm_pass)
md_lines <- c(md_lines, "",
  sprintf("**G4 confirm passes (F2 @ λ=1):** %d / %d cells", n_pass, nrow(confirm_summ)),
  "",
  sprintf("Per-job table: `%s`", file.path(confirm_dir, "gnn_confirm_job_summary.csv")),
  sprintf("Cell summary: `%s`", out_csv),
  ""
)

out_md <- Sys.getenv(
  "PIGAUTO_OUT_MD",
  unset = file.path("docs", "dev-log", "2026-08-26-gnn-evidence-phase-a-results.md")
)
if (file.exists(out_md)) {
  existing <- readLines(out_md, warn = FALSE)
  marker <- "## F2 @ λ=1 confirm (G4, 60 seeds)"
  if (any(grepl(marker, existing, fixed = TRUE))) {
    cut <- which(grepl(marker, existing, fixed = TRUE))[1L]
    existing <- existing[seq_len(cut - 1L)]
  }
  writeLines(c(existing, md_lines), out_md)
} else {
  writeLines(c("# GNN evidence Phase A — results", md_lines), out_md)
}

cat("=== collect_gnn_evidence_f2_confirm.R ===\n")
cat("confirm fits:", nrow(job_df), " failures:", sum(!job_df$fit_ok), "\n")
cat("G4 confirm passes:", n_pass, "/", nrow(confirm_summ), "\n")
cat("wrote:", out_csv, "\n")
cat("appended:", out_md, "\n")
