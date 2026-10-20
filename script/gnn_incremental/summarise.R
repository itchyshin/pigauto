#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
verify <- "--verify" %in% args
run_pos <- match("--run-dir", args)
if (is.na(run_pos) || length(args) < run_pos + 1L) {
  stop("Usage: Rscript summarise.R [--verify] --run-dir PATH",
       call. = FALSE)
}
run_dir <- args[run_pos + 1L]
config <- readRDS(file.path(run_dir, "run_config.rds"))
tasks <- utils::read.csv(file.path(run_dir, "tasks.csv"),
                         stringsAsFactors = FALSE)
files <- list.files(file.path(run_dir, "replicates"), pattern = "[.]rds$",
                    full.names = TRUE)
records <- lapply(files, readRDS)
if (!length(records)) stop("No replicate outputs are available to summarize.",
                            call. = FALSE)
scenario_order <- config$scenario_manifest$scenario_id
gate_labels <- c("FALSE" = "default_gate_off", "TRUE" = "gate_on_sensitivity")
methods <- c("pigauto_full", "pigauto_neural_off", "mean", "lm_linear",
             "lm_poly_interaction", "phylo_gls")
or_default <- function(x, default) if (is.null(x)) default else x

dataset_rows <- list()
method_rows <- list()
fail_rows <- list()
method_accounting_rows <- list()
counter <- 0L
for (rec in records) {
  scenario_id <- rec$identity$task$scenario_id
  replicate <- rec$identity$task$replicate
  for (gate_name in names(rec$gate_results)) {
    gate_result <- rec$gate_results[[gate_name]]
    gate_setting <- sub("phylo_signal_gate_", "", gate_name, fixed = TRUE)
    gate_label <- unname(gate_labels[gate_setting])
    if (!identical(gate_result$status, "completed") ||
        is.null(gate_result$predictions)) {
      fail_rows[[length(fail_rows) + 1L]] <- data.frame(
        scenario_id = scenario_id, replicate = replicate,
        phylo_signal_gate = gate_setting, status = gate_result$status,
        failure_reason = or_default(gate_result$failure_reason, NA_character_),
        elapsed_seconds = or_default(gate_result$elapsed_seconds, NA_real_),
        stringsAsFactors = FALSE
      )
    }
    status_map <- gate_result$fit_diagnostics$method_status
    for (method in methods) {
      method_accounting_rows[[length(method_accounting_rows) + 1L]] <- data.frame(
        scenario_id = scenario_id, replicate = replicate,
        phylo_signal_gate = gate_setting, method = method,
        status = if (is.null(status_map)) "not_attempted" else unname(status_map[method]),
        stringsAsFactors = FALSE
      )
    }
    if (is.null(gate_result$predictions)) next
    pred <- gate_result$predictions
    diag <- gate_result$fit_diagnostics
    for (stratum in c("all_outer_holdouts", "fully_unobserved_species",
                      "partly_observed_species")) {
      dat <- if (stratum == "all_outer_holdouts") pred else
        pred[pred$stratum == stratum, , drop = FALSE]
      if (!nrow(dat)) next
      mse <- vapply(methods, function(method) {
        nm <- if (method == "pigauto_full") "pigauto_full" else method
        if (all(c("truth", nm) %in% names(dat)) &&
            all(is.finite(dat[[nm]]))) mean((dat$truth - dat[[nm]])^2) else NA_real_
      }, numeric(1))
      rmse <- sqrt(mse)
      counter <- counter + 1L
      dataset_rows[[counter]] <- data.frame(
        scenario_id = scenario_id, replicate = replicate,
        phylo_signal_gate = gate_setting, gate_policy = gate_label,
        stratum = stratum, n_cells = nrow(dat),
        mse_full = mse[["pigauto_full"]],
        mse_neural_off = mse[["pigauto_neural_off"]],
        rmse_full = rmse[["pigauto_full"]],
        rmse_neural_off = rmse[["pigauto_neural_off"]],
        mse_difference = mse[["pigauto_full"]] - mse[["pigauto_neural_off"]],
        relative_rmse_reduction = 1 - rmse[["pigauto_full"]] /
          rmse[["pigauto_neural_off"]],
        phylo_signal = if (length(diag$phylo_signal_per_trait))
          unname(diag$phylo_signal_per_trait[[1L]]) else NA_real_,
        gate_signal_status = or_default(diag$gate_signal_status, NA_character_),
        gate_triggered = if (length(diag$phylo_gate_triggered) &&
          !is.na(diag$phylo_gate_triggered[[1L]]))
            as.logical(diag$phylo_gate_triggered[[1L]]) else NA,
        r_cal_bm = if (length(diag$r_cal_bm)) unname(diag$r_cal_bm[[1L]]) else NA_real_,
        r_cal_gnn = if (length(diag$r_cal_gnn)) unname(diag$r_cal_gnn[[1L]]) else NA_real_,
        r_cal_mean = if (length(diag$r_cal_mean)) unname(diag$r_cal_mean[[1L]]) else NA_real_,
        gls_status = or_default(diag$gls_status, NA_character_),
        stringsAsFactors = FALSE
      )
      for (method in methods) {
        method_rows[[length(method_rows) + 1L]] <- data.frame(
          scenario_id = scenario_id, replicate = replicate,
          phylo_signal_gate = gate_setting, gate_policy = gate_label,
          stratum = stratum, method = method,
          status = if (is.null(status_map)) "not_attempted" else unname(status_map[method]),
          mse = mse[[method]], rmse = rmse[[method]], n_cells = nrow(dat),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
dataset_metrics <- if (length(dataset_rows)) do.call(rbind, dataset_rows) else data.frame()
method_metrics <- if (length(method_rows)) do.call(rbind, method_rows) else data.frame()
failures <- if (length(fail_rows)) do.call(rbind, fail_rows) else data.frame()
method_accounting <- if (length(method_accounting_rows))
  do.call(rbind, method_accounting_rows) else data.frame()
method_accounting_summary <- if (nrow(method_accounting)) {
  groups <- split(method_accounting,
    interaction(method_accounting$scenario_id,
      method_accounting$phylo_signal_gate, method_accounting$method, drop = TRUE))
  do.call(rbind, lapply(groups, function(d) data.frame(
    scenario_id = d$scenario_id[1L], phylo_signal_gate = d$phylo_signal_gate[1L],
    method = d$method[1L], requested = nrow(d),
    attempted = sum(d$status != "not_attempted"), successful = sum(d$status == "ok"),
    failed = sum(grepl("failed|invalid_output", d$status)),
    nonconverged = sum(d$status == "nonconverged"), stringsAsFactors = FALSE
  )))
} else data.frame()

set.seed(config$master_seed + 900001L)
bootstrap_reduction <- function(d, B = 2000L) {
  ok <- is.finite(d$mse_full) & is.finite(d$mse_neural_off) &
    d$mse_neural_off > 0
  d <- d[ok, , drop = FALSE]
  if (nrow(d) < 2L) return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
  calc <- function(idx) 1 - sqrt(mean(d$mse_full[idx]) /
                                   mean(d$mse_neural_off[idx]))
  est <- calc(seq_len(nrow(d)))
  boot <- replicate(B, calc(sample.int(nrow(d), nrow(d), replace = TRUE)))
  c(estimate = est, lower = unname(stats::quantile(boot, 0.025)),
    upper = unname(stats::quantile(boot, 0.975)))
}

all_groups <- split(dataset_metrics,
  interaction(dataset_metrics$scenario_id, dataset_metrics$phylo_signal_gate,
              dataset_metrics$stratum, drop = TRUE))
paired_summary <- lapply(all_groups, function(d) {
  d <- d[is.finite(d$mse_difference), , drop = FALSE]
  n <- nrow(d)
  if (n < 2L) return(data.frame(
    scenario_id = d$scenario_id[1L], phylo_signal_gate = d$phylo_signal_gate[1L],
    gate_policy = d$gate_policy[1L], stratum = d$stratum[1L],
    n_paired = n, mean_mse_difference = NA_real_, mcse_mse_difference = NA_real_,
    ci_bonf_lower = NA_real_, ci_bonf_upper = NA_real_, p_one_sided = NA_real_,
    relative_rmse_reduction = NA_real_, reduction_boot_lower = NA_real_,
    reduction_boot_upper = NA_real_, classification = "indeterminate"))
  delta <- d$mse_difference
  se <- stats::sd(delta) / sqrt(n)
  critical <- stats::qt(1 - 0.05 / (2 * 2), df = n - 1L)
  t_stat <- mean(delta) / se
  p_less <- stats::pt(t_stat, df = n - 1L)
  rel <- bootstrap_reduction(d)
  data.frame(
    scenario_id = d$scenario_id[1L], phylo_signal_gate = d$phylo_signal_gate[1L],
    gate_policy = d$gate_policy[1L], stratum = d$stratum[1L], n_paired = n,
    mean_mse_difference = mean(delta), mcse_mse_difference = se,
    ci_bonf_lower = mean(delta) - critical * se,
    ci_bonf_upper = mean(delta) + critical * se,
    p_one_sided = p_less,
    relative_rmse_reduction = rel[["estimate"]],
    reduction_boot_lower = rel[["lower"]],
    reduction_boot_upper = rel[["upper"]],
    classification = "descriptive until multiplicity adjustment",
    stringsAsFactors = FALSE
  )
})
paired_summary <- if (length(paired_summary)) do.call(rbind, paired_summary) else data.frame()
if (nrow(paired_summary)) {
  paired_summary$n_requested <- vapply(paired_summary$scenario_id, function(id) {
    sum(tasks$scenario_id == id)
  }, integer(1L))
  paired_summary$n_unpaired_or_invalid <- paired_summary$n_requested -
    paired_summary$n_paired
  paired_summary$paired_fraction <- paired_summary$n_paired /
    paired_summary$n_requested
  estimated_by_key <- vapply(seq_len(nrow(paired_summary)), function(i) {
    d <- dataset_metrics[dataset_metrics$scenario_id == paired_summary$scenario_id[i] &
      dataset_metrics$phylo_signal_gate == paired_summary$phylo_signal_gate[i] &
      dataset_metrics$stratum == paired_summary$stratum[i], , drop = FALSE]
    if (!nrow(d)) return(0L)
    sum(d$gate_signal_status == "estimated", na.rm = TRUE)
  }, integer(1L))
  paired_summary$n_signal_estimated <- estimated_by_key
  paired_summary$signal_estimated_fraction <- estimated_by_key /
    paired_summary$n_requested
}

confirmatory <- c("nonlinear_primary", "interactive_secondary")
paired_summary$p_holm <- NA_real_
paired_summary$classification <- ifelse(paired_summary$scenario_id %in% confirmatory &
  paired_summary$stratum == "all_outer_holdouts", "indeterminate", "descriptive")
paired_summary$classification[paired_summary$phylo_signal_gate == "TRUE" &
  paired_summary$scenario_id %in% confirmatory &
  paired_summary$stratum == "all_outer_holdouts"] <- "policy_sensitivity"
paired_summary$classification[paired_summary$phylo_signal_gate == "TRUE" &
  paired_summary$scenario_id %in% confirmatory &
  paired_summary$stratum == "all_outer_holdouts" &
  paired_summary$signal_estimated_fraction < 0.95] <-
  "policy_sensitivity_not_exercised: signal estimated below 95%"
for (gate_setting in "FALSE") {
  idx_all <- which(paired_summary$phylo_signal_gate == gate_setting &
    paired_summary$scenario_id %in% confirmatory &
    paired_summary$stratum == "all_outer_holdouts")
  if (length(idx_all)) paired_summary$p_holm[idx_all] <- stats::p.adjust(
    paired_summary$p_one_sided[idx_all], method = "holm")
  idx <- idx_all[is.finite(paired_summary$p_one_sided[idx_all])]
  for (i in idx) {
    row <- paired_summary[i, ]
    if (row$paired_fraction < 0.95) {
      paired_summary$classification[i] <- "indeterminate: paired completion below 95%"
    } else if (is.finite(row$p_holm) && row$p_holm < 0.05 &&
        row$ci_bonf_upper < 0 && row$relative_rmse_reduction >= 0.02) {
      paired_summary$classification[i] <- "supported"
    } else if (is.finite(row$reduction_boot_upper) &&
               row$reduction_boot_upper < 0.02) {
      paired_summary$classification[i] <- "failed"
    } else {
      paired_summary$classification[i] <- "indeterminate"
    }
  }
}
pilot_confirmatory <- paired_summary$scenario_id %in% confirmatory &
  paired_summary$stratum == "all_outer_holdouts" &
  paired_summary$n_requested < 500L
paired_summary$classification[pilot_confirmatory] <-
  "indeterminate: pilot sample below 500"

method_summary <- if (nrow(method_metrics)) {
  groups <- split(method_metrics,
    interaction(method_metrics$scenario_id, method_metrics$phylo_signal_gate,
      method_metrics$stratum, method_metrics$method, drop = TRUE))
  do.call(rbind, lapply(groups, function(d) data.frame(
    scenario_id = d$scenario_id[1L], phylo_signal_gate = d$phylo_signal_gate[1L],
    gate_policy = d$gate_policy[1L], stratum = d$stratum[1L],
    method = d$method[1L], n_datasets = sum(is.finite(d$rmse)),
    mean_rmse = mean(d$rmse, na.rm = TRUE),
    mcse_rmse = stats::sd(d$rmse, na.rm = TRUE) /
      sqrt(sum(is.finite(d$rmse))),
    mean_mse = mean(d$mse, na.rm = TRUE),
    stringsAsFactors = FALSE
  )))
} else data.frame()

out_dir <- file.path(run_dir, "summaries")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(dataset_metrics, file.path(out_dir, "paired_dataset_metrics.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(method_summary, file.path(out_dir, "method_summary.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(paired_summary, file.path(out_dir, "paired_summary.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(failures, file.path(out_dir, "failures.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(method_accounting_summary,
  file.path(out_dir, "method_accounting.csv"), row.names = FALSE, na = "")
if (verify) {
  rerun <- read.csv(file.path(out_dir, "paired_dataset_metrics.csv"),
                    stringsAsFactors = FALSE)
  rerun$phylo_signal_gate <- as.character(rerun$phylo_signal_gate)
  same <- all.equal(dataset_metrics, rerun, check.attributes = FALSE)
  if (!isTRUE(same)) {
    stop("Summary regeneration did not reproduce retained metrics: ",
         paste(same, collapse = "; "))
  }
  cat("SUMMARY_REGENERATION_PASS\n")
}
