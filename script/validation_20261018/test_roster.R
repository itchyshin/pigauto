#!/usr/bin/env Rscript

source("script/validation_20261018/roster.R")

.test_job <- function(id, adapter, scenario, rep = 1L, ...) {
  x <- data.frame(job_id = id, scenario_id = sub("__.*$", "", id), slice = "test",
    adapter = adapter, scenario = scenario, rep = rep, seed = 100L + rep,
    claim_tier = "descriptive", model = NA_character_, regime_id = NA_integer_,
    n = NA_integer_, mechanism = NA_character_, frac = NA_real_,
    stringsAsFactors = FALSE)
  extras <- list(...)
  for (nm in names(extras)) x[[nm]] <- extras[[nm]]
  x
}

.test_result <- function(job, method, estimand, status = "ok") {
  data.frame(job_id = job$job_id, method = method, estimand = estimand,
    status = status, truth = NA_real_, estimate = if (status == "ok") 1 else NA_real_,
    se = NA_real_, lower = NA_real_, upper = NA_real_, complete_reference = NA_real_,
    failure_reason = if (status == "ok") NA_character_ else "planned failure",
    stringsAsFactors = FALSE)
}

jobs <- rbind(
  .test_job("analysis_a__0001", "analysis", "lm_n100_MCAR", 1L,
            model = "lm", n = 100L, mechanism = "MCAR"),
  .test_job("analysis_a__0002", "analysis", "lm_n100_MCAR", 2L,
            model = "lm", n = 100L, mechanism = "MCAR"),
  .test_job("posterior_b__0001", "posterior", "regime_33", 1L,
            model = "gls", regime_id = 33L),
  .test_job("prediction_c__0001", "prediction", "types_mixed_BM_l1_r0_mcar0.3_n100", 1L,
            mechanism = "mcar", frac = 0.3)
)
records <- list(
  list(job_id = jobs$job_id[[1]], results = rbind(
    .test_result(jobs[1, ], "complete_data", "beta_x"),
    .test_result(jobs[1, ], "analysis_lm", "beta_x", "failed")), metrics = list()),
  list(job_id = jobs$job_id[[2]], results = .test_result(jobs[2, ], "complete_data", "beta_x"), metrics = list()),
  list(job_id = jobs$job_id[[3]], results = rbind(
    .test_result(jobs[3, ], "complete_data", "beta_x"),
    .test_result(jobs[3, ], "posterior_full", "beta_x")),
    metrics = list(posterior_cell_interval = data.frame(
      trait = c("x", "y"), n = c(2L, 2L), n_expected = c(2L, 2L),
      n_finite_scored = c(2L, 1L), n_missing_or_nonfinite = c(0L, 1L),
      coverage = c(1, .5), mean_width = c(1, 2), status = c("ok", "failed")))),
  list(job_id = jobs$job_id[[4]], results = data.frame(
    job_id = character(), method = character(), estimand = character(), status = character(),
    truth = numeric(), estimate = numeric(), se = numeric(), lower = numeric(), upper = numeric(),
    complete_reference = numeric(), failure_reason = character()),
    metrics = list(calibration_bins = data.frame()))
)

audit <- validation_roster_audit(jobs, records)
stopifnot(any(audit$method == "analysis_lm" & audit$estimand == "beta_x" &
              audit$n_failed == 1L & audit$n_absent == 1L & !audit$absent_all_replicates))
stopifnot(any(audit$method == "posterior_none" & audit$estimand == "beta_x" &
              audit$n_absent == 1L & audit$absent_all_replicates))
stopifnot(any(audit$method == "gnn_off" & audit$estimand == "cat3_ece" &
              audit$n_absent == 1L & audit$absent_all_replicates))
stopifnot(any(audit$method == "floor" & audit$estimand == "ord_brier" &
              audit$n_absent == 1L & audit$absent_all_replicates))

metric_audit <- validation_metric_roster_audit(jobs, records)
stopifnot(any(metric_audit$metric_family == "calibration_bins" &
              metric_audit$status == "missing"))
stopifnot(any(metric_audit$metric_family == "posterior_cell_interval" &
              metric_audit$trait == "y" & metric_audit$status == "failed" &
              metric_audit$n_missing_or_nonfinite == 1L))

flat_results <- do.call(rbind, lapply(records, function(e) {
  if (nrow(e$results)) e$results else NULL
}))
flat_posterior <- data.frame(job_id = jobs$job_id[[3]],
  scenario_id = jobs$scenario_id[[3]], rep = jobs$rep[[3]], seed = jobs$seed[[3]],
  status = "ok", interval_status = c("ok", "failed"), trait = c("x", "y"),
  n_expected = c(2L, 2L), n_finite_scored = c(2L, 1L),
  n_missing_or_nonfinite = c(0L, 1L), n_cells = c(2L, 1L),
  coverage = c(1, .5), mean_width = c(1, 2), scope = "descriptive")
flat_bins <- data.frame(job_id = character(), scenario_id = character(), rep = integer(),
  seed = integer(), method = character(), trait = character(), bin = integer(),
  n = integer(), status = character())
table_audit <- validation_roster_audit_tables(jobs, flat_results,
  calibration_bins = flat_bins, posterior_intervals = flat_posterior)
stopifnot(any(table_audit$audit$method == "posterior_none" &
              table_audit$audit$absent_all_replicates),
          any(table_audit$metric_audit$metric_family == "posterior_cell_interval" &
              table_audit$metric_audit$trait == "y" &
              table_audit$metric_audit$n_failed == 1L &
              table_audit$metric_audit$n_missing == 0L),
          any(table_audit$metric_audit$metric_family == "calibration_bins" &
              table_audit$metric_audit$n_missing == 1L &
              table_audit$metric_audit$absent_all_replicates))
bad_status_rows <- flat_results
bad_status_rows$status[[1L]] <- "bogus"
stopifnot(inherits(try(validation_roster_audit_tables(jobs, bad_status_rows,
  calibration_bins = flat_bins, posterior_intervals = flat_posterior), silent = TRUE),
  "try-error"))
bad_status_records <- records
bad_status_records[[1L]]$results$status[[1L]] <- "bogus"
stopifnot(inherits(try(validation_roster_audit(jobs, bad_status_records), silent = TRUE),
  "try-error"))
prediction_job <- jobs[4, , drop = FALSE]
partial_bins <- data.frame(job_id = prediction_job$job_id,
  method = "gnn_off", trait = "bin", bin = 1:9,
  bin_lower = (0:8) / 10, bin_upper = (1:9) / 10, n = 1L,
  mean_confidence = 0.5, mean_correct = 0.5, absolute_gap = 0,
  status = "ok")
failed_ece <- .test_result(prediction_job, "gnn_off", "bin_ece", "failed")
partial_table_audit <- validation_roster_audit_tables(jobs,
  rbind(flat_results, failed_ece), calibration_bins = partial_bins,
  posterior_intervals = flat_posterior)
stopifnot(any(partial_table_audit$metric_audit$metric_family == "calibration_bins" &
              partial_table_audit$metric_audit$method == "gnn_off" &
              partial_table_audit$metric_audit$trait == "bin" &
          partial_table_audit$metric_audit$n_failed == 1L))

# Malformed but complete metric payloads must fail closed.
valid_bins <- data.frame(bin = 1:10, bin_lower = (0:9) / 10,
  bin_upper = (1:10) / 10, n = 1L, mean_confidence = 0.5,
  mean_correct = 0.5, absolute_gap = 0, status = "ok")
valid_ece <- data.frame(status = "ok", estimate = 0)
stopifnot(identical(.vr_calibration_payload_status(valid_bins, 10L, valid_ece), "present"))
bad_bins <- valid_bins; bad_bins$n[[1L]] <- -1L
stopifnot(identical(.vr_calibration_payload_status(bad_bins, 10L, valid_ece), "failed"))
bad_bins <- valid_bins; bad_bins$mean_confidence[[1L]] <- 1.5
stopifnot(identical(.vr_calibration_payload_status(bad_bins, 10L, valid_ece), "failed"))
bad_bins <- valid_bins; bad_bins$status[[1L]] <- "bogus"
stopifnot(identical(.vr_calibration_payload_status(bad_bins, 10L, valid_ece), "failed"))
bad_bins <- valid_bins; bad_bins$status[[1L]] <- NA_character_
stopifnot(identical(.vr_calibration_payload_status(bad_bins, 10L, valid_ece), "failed"))
bad_interval <- data.frame(interval_status = "ok", n_expected = 4L,
  n_finite_scored = 4L, n_missing_or_nonfinite = 0L, coverage = 3,
  mean_width = -1)
stopifnot(identical(.vr_posterior_interval_status(bad_interval), "failed"))
denominator_fixture <- data.frame(job_id = c("p1", "p1"), trait = c("x", "y"),
  n_expected = c(3L, 2L))
expected_fixture <- data.frame(job_id = c("p1", "p1"), trait = c("x", "y"),
  n_expected_truth = c(3L, 2L))
stopifnot(all(validation_compare_posterior_denominators(
  denominator_fixture, expected_fixture)$status == "match"))
denominator_fixture$n_expected[[1L]] <- 4L
stopifnot(any(validation_compare_posterior_denominators(
  denominator_fixture, expected_fixture)$status == "mismatch"))

# Exercise strict payload checks through the production compact-table audit.
ece_ok <- .test_result(prediction_job, "gnn_off", "bin_ece")
ece_ok$estimate <- 0
complete_bins <- valid_bins
complete_bins$job_id <- prediction_job$job_id
complete_bins$method <- "gnn_off"
complete_bins$trait <- "bin"
complete_bins <- complete_bins[, c("job_id", "method", "trait", "bin",
  "bin_lower", "bin_upper", "n", "mean_confidence", "mean_correct",
  "absolute_gap", "status")]
malformed_bins <- complete_bins
malformed_bins$n[[1L]] <- -1L
integrated_bad_bin <- validation_roster_audit_tables(jobs,
  rbind(flat_results, ece_ok), calibration_bins = malformed_bins,
  posterior_intervals = flat_posterior)
stopifnot(any(integrated_bad_bin$metric_audit$metric_family == "calibration_bins" &
              integrated_bad_bin$metric_audit$method == "gnn_off" &
              integrated_bad_bin$metric_audit$trait == "bin" &
              integrated_bad_bin$metric_audit$n_failed == 1L))
bad_posterior_receipt <- flat_posterior
bad_posterior_receipt$coverage[[1L]] <- 3
bad_posterior_receipt$mean_width[[1L]] <- -1
integrated_bad_interval <- validation_roster_audit_tables(jobs, flat_results,
  calibration_bins = flat_bins, posterior_intervals = bad_posterior_receipt)
stopifnot(any(integrated_bad_interval$metric_audit$metric_family == "posterior_cell_interval" &
              integrated_bad_interval$metric_audit$trait == "x" &
              integrated_bad_interval$metric_audit$n_failed == 1L))

write_gzip_csv <- function(x) {
  path <- tempfile(fileext = ".csv.gz"); con <- gzfile(path, "wt")
  utils::write.csv(x, con, row.names = FALSE); close(con); path
}
receipt_paths <- lapply(list(jobs, flat_results, flat_bins, flat_posterior,
                             data.frame(job_id = jobs$job_id, status = "ok")), write_gzip_csv)
compressed_audit <- validation_roster_audit_tables(
  receipt_paths[[1]], receipt_paths[[2]], receipt_paths[[3]], receipt_paths[[4]],
  receipt_paths[[5]])
stopifnot(nrow(compressed_audit$expected_roster) == nrow(table_audit$expected_roster),
          sum(compressed_audit$audit$absent_all_replicates) ==
            sum(table_audit$audit$absent_all_replicates))
unlink(unlist(receipt_paths))

# Dynamic BIEN keys are derived from explicitly named source caches, never
# copied from an observed result envelope. This exercises the schema path only.
if (requireNamespace("ape", quietly = TRUE)) {
  tree <- ape::rtree(300L)
  traits <- setNames(lapply(c("height_m", "leaf_area"), function(nm)
    data.frame(species = tree$tip.label, mean_value = seq_along(tree$tip.label))),
    c("height_m", "leaf_area"))
  tree_path <- tempfile(fileext = ".rds"); traits_path <- tempfile(fileext = ".rds")
  saveRDS(list(scenario.3 = tree), tree_path); saveRDS(traits, traits_path)
  bien_job <- .test_job("real_bien__0001", "real", "bien", 1L)
  bien_roster <- validation_expected_method_roster(bien_job,
    bien_traits_path = traits_path, bien_tree_path = tree_path)
  stopifnot(setequal(bien_roster$keys$estimand,
                     c("height_m_rmse", "leaf_area_rmse")),
            !nrow(bien_roster$issues))
  unlink(c(tree_path, traits_path))
}
cat("ROSTER_PASS\n")
