#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[1L] != "--run-dir") {
  stop("Usage: Rscript audit_run.R --run-dir PATH", call. = FALSE)
}
run_dir <- args[2L]
source("script/gnn_incremental/campaign.R")
config_path <- file.path(run_dir, "run_config.rds")
tasks_path <- file.path(run_dir, "tasks.csv")
replicate_dir <- file.path(run_dir, "replicates")
if (!all(file.exists(c(config_path, tasks_path))) ||
    !dir.exists(replicate_dir)) {
  stop("Run directory is missing its config, tasks, or replicate folder.",
       call. = FALSE)
}
config <- readRDS(config_path)
tasks <- utils::read.csv(tasks_path, stringsAsFactors = FALSE)
verify_config_fingerprint(config)
verify_task_manifest(config, tasks)
manifest_live <- read.csv(
  "script/gnn_incremental/scenario_manifest.csv", stringsAsFactors = FALSE
)
expected_task_count <- nrow(manifest_live) * config$n_replicates
if (nrow(tasks) != expected_task_count || anyDuplicated(tasks$task_index) ||
    anyDuplicated(paste(tasks$scenario_id, tasks$replicate))) {
  stop("Task manifest is incomplete or has duplicate task keys.")
}
if (identical(config$run_type, "production") && config$n_replicates != 500L) {
  stop("Production run must allocate exactly 500 replicates per scenario.")
}
if (!identical(as.integer(config$master_seed), 20261019L)) {
  stop("Run config does not use the frozen master seed 20261019.")
}
for (i in seq_len(nrow(tasks))) {
  scenario_index <- tasks$scenario_index[i]
  if (is.na(scenario_index) || scenario_index < 1L ||
      scenario_index > nrow(config$scenario_manifest) ||
      !identical(as.character(tasks$scenario_id[i]),
        as.character(config$scenario_manifest$scenario_id[scenario_index]))) {
    stop("Task scenario_index does not map to its declared scenario_id.")
  }
}
for (i in seq_len(nrow(config$scenario_manifest))) {
  reps <- sort(as.integer(tasks$replicate[tasks$scenario_index == i]))
  if (!identical(reps, seq_len(config$n_replicates))) {
    stop("Scenario does not contain exactly replicates 1..N_REPS: ",
         config$scenario_manifest$scenario_id[i])
  }
}
if (!identical(source_identity(), config$source)) {
  stop("Current source/environment identity differs from campaign initialization.")
}
files <- list.files(replicate_dir, pattern = "[.]rds$", full.names = TRUE)
actual <- basename(files)
expected <- sprintf("%s__rep_%04d.rds", tasks$scenario_id, tasks$replicate)
if (anyDuplicated(actual)) stop("Duplicate replicate files are present.")
missing <- setdiff(expected, actual)
unexpected <- setdiff(actual, expected)

records <- lapply(files, function(path) tryCatch(readRDS(path),
  error = function(e) structure(list(error = conditionMessage(e)),
                                class = "audit_read_error")))
if (any(vapply(records, inherits, logical(1), "audit_read_error"))) {
  stop("At least one replicate file cannot be read.")
}
manifest_index <- vapply(records, function(x) {
  task <- x$identity$task
  if (!identical(as.character(task$scenario_id),
      as.character(config$scenario_manifest$scenario_id[task$scenario_index]))) NA_integer_
  else as.integer(task$scenario_index)
}, integer(1))
expected_fingerprints <- vapply(seq_along(records), function(i) {
  x <- records[[i]]
  task <- x$identity$task
  scenario <- config$scenario_manifest[task$scenario_index, , drop = FALSE]
  seed_for <- function(stream) as.integer((as.double(config$master_seed) +
    task$scenario_index * 1000003 + task$replicate * 7919 +
    stream * 104729) %% 2147483646) + 1L
  identity <- list(
    fingerprint = config$fingerprint, task = task,
    dgp_seed = seed_for(1L), fit_seed = seed_for(2L),
    f_type = scenario$f_type, variance_multiplier = scenario$variance_multiplier,
    beta_strength = scenario$beta_strength, epochs = config$epochs
  )
  digest::digest(identity, algo = "sha256", serialize = TRUE)
}, character(1))
if (anyNA(manifest_index) || any(vapply(seq_along(records), function(i) {
  !identical(records[[i]]$fingerprint, expected_fingerprints[i]) ||
    !identical(records[[i]]$source, config$source)
}, logical(1)))) {
  stop("A replicate has a stale, mismatched, or unrecognized fingerprint.")
}

statuses <- unlist(lapply(records, function(x) {
  vapply(x$gate_results, function(y) y$status, character(1))
}), use.names = FALSE)
fit_counts <- table(factor(statuses,
  levels = c("completed", "invalid_output", "failed", "nonconverged", "infrastructure_skip")))
gls_nonconverged <- sum(vapply(records, function(x) sum(vapply(
  x$gate_results, function(y) identical(y$fit_diagnostics$gls_status,
                                        "nonconverged"), logical(1))),
  integer(1)))
signal_status <- vapply(records, function(rec) {
  gate <- rec$gate_results[["phylo_signal_gate_TRUE"]]
  if (is.null(gate$fit_diagnostics$gate_signal_status)) "missing" else
    as.character(gate$fit_diagnostics$gate_signal_status)
}, character(1L))
expected_fits <- 2L * nrow(tasks)
observed_fits <- length(statuses)
scenario_counts <- table(factor(vapply(records, function(x) {
  x$identity$task$scenario_id
}, character(1)), levels = config$scenario_manifest$scenario_id))
complete_counts <- vapply(config$scenario_manifest$scenario_id, function(id) {
  sum(vapply(records, function(x) identical(x$identity$task$scenario_id, id) &&
                all(vapply(x$gate_results, function(y) y$status == "completed",
                           logical(1))), logical(1)))
}, integer(1))
requested_by_scenario <- vapply(config$scenario_manifest$scenario_id,
  function(id) sum(tasks$scenario_id == id), integer(1))

summary <- data.frame(
  requested_datasets = nrow(tasks),
  completed_dataset_records = length(records),
  missing_datasets = length(missing), unexpected_results = length(unexpected),
  requested_gate_fits = expected_fits, recorded_gate_fits = observed_fits,
  completed_gate_fits = unname(fit_counts["completed"]),
  invalid_output_gate_fits = unname(fit_counts["invalid_output"]),
  failed_gate_fits = unname(fit_counts["failed"]),
  nonconverged_gate_fits = unname(fit_counts["nonconverged"]),
  infrastructure_skips = unname(fit_counts["infrastructure_skip"]),
  gls_nonconverged_attempts = gls_nonconverged,
  gate_true_signal_estimated = sum(signal_status == "estimated"),
  gate_true_signal_unestimable = sum(signal_status == "unestimable"),
  gate_true_signal_missing = sum(signal_status == "missing"),
  stringsAsFactors = FALSE
)
cat("Scenario accounting:\n")
print(data.frame(scenario_id = config$scenario_manifest$scenario_id,
  requested = requested_by_scenario,
  dataset_records = as.integer(scenario_counts),
  complete_paired_gate_fits = complete_counts))
cat("Run accounting:\n")
print(summary)
method_names <- c("pigauto_full", "pigauto_neural_off", "mean", "lm_linear",
                  "lm_poly_interaction", "phylo_gls")
method_statuses <- setNames(lapply(method_names, function(method) {
  vapply(records, function(rec) {
    unlist(lapply(rec$gate_results, function(g) {
      status <- g$fit_diagnostics$method_status
      if (is.null(status)) "not_attempted" else unname(status[method])
    }), use.names = FALSE)
  }, character(2L))
}), method_names)
method_accounting <- do.call(rbind, lapply(method_names, function(method) {
  status <- unlist(method_statuses[[method]], use.names = FALSE)
  data.frame(method = method, requested = expected_fits,
    attempted = sum(status != "not_attempted"),
    successful = sum(status == "ok"), failed = sum(grepl("failed|invalid_output", status)),
    nonconverged = sum(status == "nonconverged"),
    stringsAsFactors = FALSE)
}))
cat("Method accounting:\n")
print(method_accounting)
if (length(missing) || length(unexpected) || observed_fits != expected_fits) {
  stop("RUN_ACCOUNTING_INCOMPLETE", call. = FALSE)
}
known_statuses <- c("completed", "invalid_output", "failed", "nonconverged",
                    "infrastructure_skip")
if (any(!statuses %in% known_statuses) || sum(fit_counts) != expected_fits) {
  stop("Gate-fit attempt statuses do not reconcile to the requested denominator.")
}
cat("RUN_ACCOUNTING_PASS\n")
