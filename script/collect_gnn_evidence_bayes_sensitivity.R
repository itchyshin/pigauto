#!/usr/bin/env Rscript
#
# Bayes sensitivity arm analysis — low-λ cells only (prereg §4).
#
# Compares lambda_mode = "bayes" vs primary "fixed_1" on λ_DGP ∈ {0.2, 0.5}.
# gnn_res = paired_delta under bayes arm (blend − λ-fitted baseline).
# Closure = fixed_1 arm passed G4 explore but bayes arm does not.
#
# Usage:
#   Rscript script/collect_gnn_evidence_bayes_sensitivity.R

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

bayes_dir <- Sys.getenv(
  "PIGAUTO_BAYES_DIR",
  unset = file.path("script", "returned_gnn_campaign", "results_bayes")
)
primary_dir <- Sys.getenv(
  "PIGAUTO_PRIMARY_DIR",
  unset = file.path("script", "returned_gnn_campaign", "results")
)

stopifnot(dir.exists(bayes_dir))
rds_files <- sort(Sys.glob(file.path(bayes_dir, "gnn_campaign_job_*.rds")))
if (!length(rds_files)) stop("no bayes RDS in ", bayes_dir)

rows <- lapply(rds_files, readRDS)

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
    lambda_mode = r$lambda_mode,
    fit_ok = isTRUE(r$fit_ok),
    blend_loss = r$blend_loss,
    baseline_loss = r$baseline_loss,
    paired_delta = r$paired_delta,
    r_cal_gnn_mean = r$r_cal_gnn_mean,
    floor_fired = isTRUE(r$floor_fired),
    stringsAsFactors = FALSE
  )
}))

write.csv(job_df, file.path(bayes_dir, "gnn_bayes_job_summary.csv"), row.names = FALSE)

summarize_cell <- function(d) {
  ok <- d$fit_ok
  deltas <- d$paired_delta[ok]
  bl <- d$baseline_loss[ok]
  mean_delta <- if (length(deltas)) mean(deltas, na.rm = TRUE) else NA_real_
  sd_delta <- if (length(deltas) > 1L) sd(deltas, na.rm = TRUE) else NA_real_
  mcse <- if (length(deltas) > 1L) sd_delta / sqrt(length(deltas)) else NA_real_
  rel_improve <- if (length(bl) && mean(bl, na.rm = TRUE) > 0) {
    -mean_delta / mean(bl, na.rm = TRUE)
  } else NA_real_
  z_ratio <- if (!is.na(mcse) && mcse > 0) abs(mean_delta) / mcse else NA_real_
  data.frame(
    n_rep = nrow(d),
    n_ok = sum(ok),
    mean_delta = mean_delta,
    mcse = mcse,
    z_ratio = z_ratio,
    rel_improve = rel_improve,
    mean_r_cal_gnn = mean(d$r_cal_gnn_mean[ok], na.rm = TRUE),
    frac_gate_open = mean(d$r_cal_gnn_mean[ok] > 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

g4_residual_pass <- function(family, mean_delta, z_ratio, rel_improve) {
  family != "F1" &&
    !is.na(z_ratio) && !is.na(rel_improve) &&
    z_ratio >= 3 && rel_improve >= 0.02 && mean_delta < 0
}

cell_ids <- unique(job_df[, c("family", "n_species", "phylo_lambda", "miss_frac")])
bayes_summ <- do.call(rbind, lapply(seq_len(nrow(cell_ids)), function(i) {
  d <- job_df[
    job_df$family == cell_ids$family[i] &
      job_df$n_species == cell_ids$n_species[i] &
      job_df$phylo_lambda == cell_ids$phylo_lambda[i] &
      abs(job_df$miss_frac - cell_ids$miss_frac[i]) < 1e-9,
  , drop = FALSE]
  cbind(cell_ids[i, , drop = FALSE], summarize_cell(d))
}))
rownames(bayes_summ) <- NULL

bayes_summ$gnn_res_pass <- mapply(
  g4_residual_pass,
  bayes_summ$family,
  bayes_summ$mean_delta,
  bayes_summ$z_ratio,
  bayes_summ$rel_improve
)

write.csv(bayes_summ, file.path(bayes_dir, "gnn_bayes_cell_summary.csv"), row.names = FALSE)

# ---- merge with primary fixed_1 --------------------------------------------
primary_summ <- read.csv(
  file.path(primary_dir, "gnn_campaign_cell_summary.csv"),
  stringsAsFactors = FALSE
)

merge_key <- function(d) {
  paste(d$family, d$n_species, d$phylo_lambda, d$miss_frac, sep = "|")
}

primary_low <- primary_summ[primary_summ$phylo_lambda %in% c(0.2, 0.5), , drop = FALSE]
primary_low$merge_key <- merge_key(primary_low)
bayes_summ$merge_key <- merge_key(bayes_summ)

cmp <- merge(
  primary_low[, c(
    "merge_key", "family", "n_species", "phylo_lambda", "miss_frac",
    "mean_delta", "mcse", "z_ratio", "rel_improve", "g4_explore_pass"
  )],
  bayes_summ[, c(
    "merge_key", "mean_delta", "mcse", "z_ratio", "rel_improve",
    "gnn_res_pass", "mean_r_cal_gnn", "frac_gate_open"
  )],
  by = "merge_key",
  suffixes = c("_fixed1", "_bayes")
)

cmp$closed <- cmp$g4_explore_pass & !cmp$gnn_res_pass
cmp$retained <- cmp$g4_explore_pass & cmp$gnn_res_pass
cmp$delta_attenuation <- cmp$mean_delta_bayes - cmp$mean_delta_fixed1

write.csv(cmp, file.path(bayes_dir, "gnn_bayes_closure_comparison.csv"), row.names = FALSE)

closure_by_family_lambda <- do.call(rbind, lapply(
  split(cmp, list(cmp$family, cmp$phylo_lambda), drop = TRUE),
  function(d) {
    g4_cells <- sum(d$g4_explore_pass, na.rm = TRUE)
    data.frame(
      family = d$family[1L],
      phylo_lambda = d$phylo_lambda[1L],
      n_cells = nrow(d),
      n_g4_fixed1 = g4_cells,
      n_gnn_res_bayes = sum(d$gnn_res_pass, na.rm = TRUE),
      n_closed = sum(d$closed, na.rm = TRUE),
      pct_closed_of_g4 = if (g4_cells > 0) 100 * sum(d$closed) / g4_cells else NA_real_,
      pct_gnn_res_of_cells = 100 * sum(d$gnn_res_pass) / nrow(d),
      stringsAsFactors = FALSE
    )
  }
))
rownames(closure_by_family_lambda) <- NULL
write.csv(
  closure_by_family_lambda,
  file.path(bayes_dir, "gnn_bayes_closure_by_family_lambda.csv"),
  row.names = FALSE
)

fmt_num <- function(x, digits = 4) {
  ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits))
}

md_lines <- c(
  "",
  "## Bayes sensitivity (λ_DGP ∈ {0.2, 0.5}, lambda_mode = bayes)",
  "",
  sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  sprintf("**Results dir:** `%s`", bayes_dir),
  sprintf("**Fits:** %d total, %d failures", nrow(job_df), sum(!job_df$fit_ok)),
  "",
  "gnn_res = blend − λ-fitted baseline under `lambda_mode = \"bayes\"`.",
  "Closure = fixed_1 arm passed G4 explore but bayes arm does not (|Δ|/MCSE < 3 or rel improve < 2%).",
  "",
  "### Closure summary by family × λ",
  "",
  "| family | λ | cells | G4 fixed_1 | gnn_res bayes | closed | % closed of G4 |",
  "|---|---:|---:|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(closure_by_family_lambda))) {
  r <- closure_by_family_lambda[i, ]
  md_lines <- c(md_lines, sprintf(
    "| %s | %.1f | %d | %d | %d | %d | %s |",
    r$family, r$phylo_lambda, r$n_cells,
    r$n_g4_fixed1, r$n_gnn_res_bayes, r$n_closed,
    fmt_num(r$pct_closed_of_g4, 1)
  ))
}

md_lines <- c(md_lines, "",
  "### F1 vs F2 detail (low-λ cells)",
  "",
  "| family | λ | n | miss | z fixed_1 | z bayes | G4 fixed_1 | gnn_res bayes | closed |",
  "|---|---:|---:|---:|---:|---:|---|---:|---|"
)

for (i in seq_len(nrow(cmp))) {
  r <- cmp[i, ]
  md_lines <- c(md_lines, sprintf(
    "| %s | %.1f | %d | %.0f%% | %s | %s | %s | %s | %s |",
    r$family, r$phylo_lambda, r$n_species, 100 * r$miss_frac,
    fmt_num(r$z_ratio_fixed1), fmt_num(r$z_ratio_bayes),
    if (isTRUE(r$g4_explore_pass)) "PASS" else "—",
    if (isTRUE(r$gnn_res_pass)) "PASS" else "—",
    if (isTRUE(r$closed)) "yes" else "no"
  ))
}

out_md <- Sys.getenv(
  "PIGAUTO_OUT_MD",
  unset = file.path("docs", "dev-log", "2026-08-26-gnn-evidence-phase-a-results.md")
)
if (file.exists(out_md)) {
  existing <- readLines(out_md, warn = FALSE)
  marker <- "## Bayes sensitivity"
  if (any(grepl(marker, existing, fixed = TRUE))) {
    cut <- which(grepl(marker, existing, fixed = TRUE))[1L]
    existing <- existing[seq_len(cut - 1L)]
  }
  writeLines(c(existing, md_lines), out_md)
} else {
  writeLines(c("# GNN evidence Phase A — results", md_lines), out_md)
}

cat("=== collect_gnn_evidence_bayes_sensitivity.R ===\n")
cat("bayes fits:", nrow(job_df), " failures:", sum(!job_df$fit_ok), "\n")
print(closure_by_family_lambda)
cat("wrote:", file.path(bayes_dir, "gnn_bayes_closure_comparison.csv"), "\n")
cat("appended:", out_md, "\n")
