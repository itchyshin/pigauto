#!/usr/bin/env Rscript

source("script/gnn_incremental/campaign.R")
split_fixture <- list(
  X_scaled = matrix(c(10, 20, NA, 40, 50, 60), ncol = 1L),
  n_obs = 6L
)
split_contract <- list(mask = matrix(c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE),
                                     ncol = 1L))
stopifnot(identical(training_rows_from_splits(split_fixture, split_contract),
                    c(1L, 4L, 6L)))
run_dir <- tempfile("gnn-accounting-")
init_campaign(run_dir, n_replicates = 2L, master_seed = 20261019L, epochs = 2L)
config <- readRDS(file.path(run_dir, "run_config.rds"))
tasks <- utils::read.csv(file.path(run_dir, "tasks.csv"),
                         stringsAsFactors = FALSE)

for (i in seq_len(nrow(tasks))) {
  task <- tasks[i, , drop = FALSE]
  scenario <- config$scenario_manifest[task$scenario_index, , drop = FALSE]
  dgp_seed <- derive_seed(config$master_seed, task$scenario_index,
                          task$replicate, 1L)
  fit_seed <- derive_seed(config$master_seed, task$scenario_index,
                          task$replicate, 2L)
  identity <- list(
    fingerprint = config$fingerprint, task = task,
    dgp_seed = dgp_seed, fit_seed = fit_seed,
    f_type = scenario$f_type, variance_multiplier = scenario$variance_multiplier,
    beta_strength = scenario$beta_strength, epochs = config$epochs
  )
  replicate_record <- list(
    identity = identity, fingerprint = configuration_fingerprint(identity),
    source = config$source,
    gate_results = lapply(c(FALSE, TRUE), function(enabled) {
      p <- if (enabled) 0.8 else 0.7
      predictions <- data.frame(
        row_id = 1:4, species = paste0("sp", 1:4),
        stratum = rep(c("fully_unobserved_species",
                        "partly_observed_species"), each = 2L),
        truth = c(0, 1, 1, 2),
        pigauto_full = c(0, 0.8, 0.8, 1.8),
        pigauto_neural_off = c(0, 0.5, 0.5, 1.5),
        mean = rep(1, 4), lm_linear = c(0, 1, 1, 2),
        lm_poly_interaction = c(0, 1, 1, 2),
        phylo_gls = c(0, 1, 1, 2), stringsAsFactors = FALSE
      )
      list(status = "completed", failure_reason = NA_character_,
        elapsed_seconds = 1, predictions = predictions,
        fit_diagnostics = list(
          phylo_signal_per_trait = if (enabled) c(y = 0.3) else c(y = NA_real_),
          gate_signal_status = if (enabled) "estimated" else "not_requested",
          phylo_gate_triggered = c(y = NA), r_cal_bm = c(y = 0.3),
          r_cal_gnn = c(y = 0.5), r_cal_mean = c(y = 0.2),
          gls_status = "failed: synthetic comparator control",
          method_status = c(pigauto_full = "ok", pigauto_neural_off = "ok",
            mean = "ok", lm_linear = "ok", lm_poly_interaction = "ok",
            phylo_gls = "failed: synthetic comparator control")
        ))
    })
  )
  names(replicate_record$gate_results) <- c("phylo_signal_gate_FALSE",
                                             "phylo_signal_gate_TRUE")
  if (i == 1L) {
    replicate_record$gate_results[[1L]]$status <- "invalid_output"
    replicate_record$gate_results[[1L]]$failure_reason <-
      "synthetic non-finite neural prediction control"
    replicate_record$gate_results[[1L]]$predictions$pigauto_full[1L] <- NA_real_
    replicate_record$gate_results[[1L]]$fit_diagnostics$method_status["pigauto_full"] <-
      "invalid_output"
  }
  out <- file.path(run_dir, "replicates",
    sprintf("%s__rep_%04d.rds", task$scenario_id, task$replicate))
  atomic_save_rds(replicate_record, out)
}

drifted <- config
drifted$source$git_head <- "deliberately-stale-source"
drifted$fingerprint <- NULL
drifted$fingerprint <- configuration_fingerprint(drifted)
drift_error <- tryCatch({
  run_one_dataset(drifted, tasks[1L, , drop = FALSE])
  NULL
}, error = function(e) e)
if (is.null(drift_error) || !grepl("source/environment identity", conditionMessage(drift_error))) {
  stop("Source-drift negative control was not detected.")
}

audit_out <- system2("Rscript", c("--vanilla", "script/gnn_incremental/audit_run.R",
  "--run-dir", run_dir), stdout = TRUE, stderr = TRUE)
if (!any(grepl("RUN_ACCOUNTING_PASS", audit_out, fixed = TRUE))) {
  stop(paste(audit_out, collapse = "\n"))
}

audit_must_reject <- function(path, label) {
  out <- system2("Rscript", c("--vanilla", "script/gnn_incremental/audit_run.R",
    "--run-dir", path), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  if (status == 0L) stop("Audit accepted tampered ", label, ".")
}

tampered_config_parent <- tempfile("gnn-tampered-config-")
dir.create(tampered_config_parent)
file.copy(run_dir, tampered_config_parent, recursive = TRUE)
tampered_config_dir <- file.path(tampered_config_parent, basename(run_dir))
tampered_config <- readRDS(file.path(tampered_config_dir, "run_config.rds"))
tampered_config$fit_settings$dropout <- 0.9
saveRDS(tampered_config, file.path(tampered_config_dir, "run_config.rds"))
audit_must_reject(tampered_config_dir, "run configuration fingerprint")

tampered_tasks_parent <- tempfile("gnn-tampered-tasks-")
dir.create(tampered_tasks_parent)
file.copy(run_dir, tampered_tasks_parent, recursive = TRUE)
tampered_tasks_dir <- file.path(tampered_tasks_parent, basename(run_dir))
tampered_tasks <- utils::read.csv(file.path(tampered_tasks_dir, "tasks.csv"))
tampered_tasks$task_index[1L] <- max(tampered_tasks$task_index) + 1L
utils::write.csv(tampered_tasks, file.path(tampered_tasks_dir, "tasks.csv"),
                 row.names = FALSE)
audit_must_reject(tampered_tasks_dir, "noncanonical task_index")

summary_out <- system2("Rscript", c("--vanilla", "script/gnn_incremental/summarise.R",
  "--verify", "--run-dir", run_dir), stdout = TRUE, stderr = TRUE)
if (!any(grepl("SUMMARY_REGENERATION_PASS", summary_out, fixed = TRUE))) {
  stop(paste(summary_out, collapse = "\n"))
}
method_accounting <- read.csv(file.path(run_dir, "summaries", "method_accounting.csv"))
if (!any(method_accounting$method == "phylo_gls" & method_accounting$failed > 0L) ||
    !any(method_accounting$method == "pigauto_full" & method_accounting$failed > 0L)) {
  stop("Summary accounting did not retain synthetic method failures.")
}
dataset_metrics <- read.csv(file.path(run_dir, "summaries", "paired_dataset_metrics.csv"))
if (!any(dataset_metrics$phylo_signal_gate == "TRUE" &
         dataset_metrics$gate_signal_status == "estimated") ||
    any(dataset_metrics$gate_signal_status == "unestimable" &
        dataset_metrics$gate_triggered %in% FALSE)) {
  stop("Signal-estimation status or unknown gate trigger was not retained.")
}
paired_summary <- read.csv(file.path(run_dir, "summaries", "paired_summary.csv"))
pilot_confirmatory <- paired_summary$scenario_id %in%
  c("nonlinear_primary", "interactive_secondary") &
  paired_summary$stratum == "all_outer_holdouts"
if (any(grepl("^supported$|^failed$", paired_summary$classification[pilot_confirmatory]))) {
  stop("Pilot-size data were incorrectly assigned an efficacy classification.")
}
cat("ACCOUNTING_AND_SUMMARY_VALIDATION_PASS\n")
