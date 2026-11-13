#!/usr/bin/env Rscript

source("script/validation_20261018/roster.R")

args <- commandArgs(trailingOnly = TRUE)
receipt_dir <- if (length(args) >= 1L) args[[1L]] else
  "docs/dev-log/validation-20261018/receipts"
out_dir <- if (length(args) >= 2L) args[[2L]] else receipt_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bien_roster_path <- file.path(receipt_dir, "production-bien-source-roster.csv")
bien_verify_path <- file.path(receipt_dir, "production-bien-roster-verification.txt")
bien_roster <- utils::read.csv(bien_roster_path, stringsAsFactors = FALSE)
bien_verify_lines <- if (file.exists(bien_verify_path)) readLines(bien_verify_path) else character()
marker <- which(grepl("^BIEN_TOTORO_ROSTER_VERIFIED$", bien_verify_lines))
if (length(marker) != 1L || marker[[1L]] != length(bien_verify_lines))
  stop("fresh read-only Totoro BIEN roster verification receipt is missing")
verified_bien_roster <- utils::read.csv(text = paste(bien_verify_lines[-marker], collapse = "\n"),
  stringsAsFactors = FALSE, check.names = FALSE)
if (!isTRUE(all.equal(verified_bien_roster, bien_roster, check.attributes = FALSE)))
  stop("Totoro BIEN verification payload does not match the retained source roster")
source_identity <- readRDS(file.path(receipt_dir, "proposed-source.rds"))
expected_cache_hashes <- source_identity$inputs[c("VALIDATION_BIEN_TRAITS", "VALIDATION_BIEN_TREE")]
if (!identical(unname(as.character(unique(bien_roster$traits_cache_md5))),
               unname(as.character(expected_cache_hashes[["VALIDATION_BIEN_TRAITS"]]))) ||
    !identical(unname(as.character(unique(bien_roster$tree_cache_md5))),
               unname(as.character(expected_cache_hashes[["VALIDATION_BIEN_TREE"]]))) ||
    any(bien_roster$matched_species != 6425L) || any(bien_roster$subset_n != 300L) ||
    any(bien_roster$subset_seed != 20261018L) || anyDuplicated(bien_roster$trait))
  stop("BIEN source-derived trait roster does not match the approved input identity")

denominator_path <- file.path(receipt_dir, "production-posterior-denominator-audit.csv.gz")
denominator_audit <- utils::read.csv(gzfile(denominator_path), stringsAsFactors = FALSE)
if (!all(c("job_id", "trait", "n_expected_truth", "n_expected_receipt", "status") %in%
         names(denominator_audit)) || nrow(denominator_audit) != 8000L ||
    any(denominator_audit$status != "match"))
  stop("independently regenerated posterior denominator audit is incomplete or mismatched")

audit <- validation_roster_audit_tables(
  file.path(receipt_dir, "jobs.csv.gz"),
  file.path(receipt_dir, "production-result_rows.csv.gz"),
  file.path(receipt_dir, "production-calibration_bins_by_replicate.csv.gz"),
  file.path(receipt_dir, "production-posterior_intervals_by_replicate.csv.gz"),
  file.path(receipt_dir, "production-attempts.csv.gz"),
  root = ".", bien_trait_names = as.character(bien_roster$trait))

write_gzip <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.csv(x, con, row.names = FALSE, na = "")
}

write_gzip(audit$audit, file.path(out_dir, "production-roster-audit.csv.gz"))
write_gzip(audit$metric_audit, file.path(out_dir, "production-metric-roster-audit.csv.gz"))
write_gzip(audit$unexpected_rows, file.path(out_dir, "production-unexpected-rows.csv.gz"))
write_gzip(audit$unrostered_rows, file.path(out_dir, "production-unrostered-rows.csv.gz"))
write_gzip(audit$roster_issues, file.path(out_dir, "production-roster-issues.csv.gz"))

summary <- data.frame(
  jobs = 32240L,
  roster_groups = nrow(audit$audit),
  planned_method_estimand = sum(audit$audit$n_planned),
  present_method_estimand = sum(audit$audit$n_present),
  absent_method_estimand = sum(audit$audit$n_absent),
  failed_method_estimand = sum(audit$audit$n_failed),
  unavailable_method_estimand = sum(audit$audit$n_unavailable),
  absent_all_groups = sum(audit$audit$absent_all_replicates),
  partial_groups = sum(audit$audit$audit_status == "partially_absent"),
  unexpected_rows = nrow(audit$unexpected_rows),
  unrostered_rows = nrow(audit$unrostered_rows),
  roster_issues = nrow(audit$roster_issues),
  metric_groups = nrow(audit$metric_audit),
  metric_groups_missing = sum(audit$metric_audit$n_missing),
  metric_groups_failed = sum(audit$metric_audit$n_failed),
  metric_groups_unavailable = sum(audit$metric_audit$n_unavailable),
  posterior_denominator_rows = nrow(denominator_audit),
  posterior_denominator_mismatches = sum(denominator_audit$status != "match"),
  stringsAsFactors = FALSE)
write_gzip(summary, file.path(out_dir, "production-roster-summary.csv.gz"))
print(summary)
cat("ROSTER_RECEIPTS_AUDITED\n")
