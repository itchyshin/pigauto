#!/usr/bin/env Rscript

# Summarise completed pigauto v0.10.0 MI validation tasks and apply the
# pre-specified release gates. This script never treats a pilot as release proof.

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else ""
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path)) else getwd()
source(file.path(script_dir, "0_prepare.R"))

mi_methods <- c(
  "conformal", "mc_dropout", "pmm", "oracle_conditional", "standard_smc"
)
release_methods <- c("conformal", "mc_dropout", "pmm")
control_methods <- c("oracle_conditional", "standard_smc")
calibration_only <- "--calibration-only" %in% args ||
  tolower(.env_value("PIGAUTO_MI_CALIBRATION_ONLY", "false")) %in%
    c("true", "1", "yes")

config <- mi_validation_config(args)
manifest_path <- file.path(config$output_dir, "manifest.rds")
if (!file.exists(manifest_path)) {
  stop("No manifest at ", manifest_path, ". Run 0_prepare.R first.", call. = FALSE)
}
manifest <- readRDS(manifest_path)
files <- list.files(config$output_dir, pattern = "^task-[0-9]+\\.rds$",
                    full.names = TRUE)
if (length(files) == 0L) {
  stop("No task RDS files found under ", config$output_dir, ".", call. = FALSE)
}
results <- lapply(files, readRDS)

.meta_value <- function(result, name) unlist(result$meta[[name]], use.names = FALSE)[1]

.bind_component <- function(name) {
  rows <- lapply(results, function(result) {
    value <- result[[name]]
    if (is.null(value) || !is.data.frame(value) || nrow(value) == 0L) return(NULL)
    value$task_id <- as.integer(.meta_value(result, "task_id"))
    value$dgp <- as.character(.meta_value(result, "dgp"))
    value$regime <- as.character(.meta_value(result, "regime"))
    value$replicate <- as.integer(.meta_value(result, "replicate"))
    value$m_requested <- as.integer(.meta_value(result, "m"))
    value$task_status <- result$status
    value
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) data.frame() else do.call(rbind, rows)
}

fixed <- .bind_component("fixed_effects")
variance <- .bind_component("variance_components")
gates <- .bind_component("gate")
training <- .bind_component("training")
processed <- data.frame(
  task_id = vapply(results, function(x) as.integer(.meta_value(x, "task_id")),
                   integer(1)),
  dgp = vapply(results, function(x) as.character(.meta_value(x, "dgp")),
               character(1)),
  regime = vapply(results, function(x) as.character(.meta_value(x, "regime")),
                  character(1)),
  status = vapply(results, `[[`, character(1), "status"),
  elapsed_seconds = vapply(results, function(x) x$elapsed_seconds %||% NA_real_,
                           numeric(1)),
  stringsAsFactors = FALSE
)

comparator_diagnostics <- do.call(rbind, lapply(results, function(result) {
  value <- result$comparator_diagnostics
  if (is.null(value) || length(value) == 0L) return(NULL)
  data.frame(
    task_id = as.integer(.meta_value(result, "task_id")),
    dgp = as.character(.meta_value(result, "dgp")),
    regime = as.character(.meta_value(result, "regime")),
    replicate = as.integer(.meta_value(result, "replicate")),
    engine = value$engine %||% NA_character_,
    elapsed_seconds = value$elapsed_seconds %||% NA_real_,
    warning_count = length(value$warnings %||% character(0)),
    warnings = paste(value$warnings %||% character(0), collapse = " | "),
    stringsAsFactors = FALSE
  )
}))
if (is.null(comparator_diagnostics)) comparator_diagnostics <- data.frame()

if (nrow(fixed) > 0L) {
  oracle_reference <- fixed[
    fixed$method == "oracle",
    c("task_id", "dgp", "regime", "replicate", "term", "estimate"),
    drop = FALSE
  ]
  names(oracle_reference)[names(oracle_reference) == "estimate"] <-
    "oracle_estimate"
  fixed <- merge(
    fixed, oracle_reference,
    by = c("task_id", "dgp", "regime", "replicate", "term"),
    all.x = TRUE, sort = FALSE
  )
}

.cell_counts <- function(dgp, regime) {
  n_planned <- sum(manifest$dgp == dgp & manifest$regime == regime)
  cell <- processed$dgp == dgp & processed$regime == regime
  c(n_planned = n_planned, n_processed = sum(cell))
}

.summarise_fixed_group <- function(data) {
  counts <- .cell_counts(data$dgp[[1]], data$regime[[1]])
  mi_method <- data$method[[1]] %in% mi_methods
  downstream_ok <- if (mi_method) {
    data$n_fits >= ceiling(0.95 * data$m_requested)
  } else {
    rep(TRUE, nrow(data))
  }
  qualified <- data[downstream_ok, , drop = FALSE]
  estimates <- qualified$estimate
  empirical_sd <- if (length(estimates) > 1L) {
    stats::sd(estimates, na.rm = TRUE)
  } else NA_real_
  truth <- unique(data$truth)
  truth <- truth[is.finite(truth)][1]
  finite <- is.finite(qualified$estimate) & is.finite(qualified$std.error) &
    is.finite(qualified$conf.low) & is.finite(qualified$conf.high) &
    (is.finite(qualified$df) | is.infinite(qualified$df)) & qualified$df > 0
  rubin_valid <- if (mi_method) {
    is.finite(qualified$fmi) & qualified$fmi >= 0 & qualified$fmi <= 1 &
      (is.finite(qualified$riv) | is.infinite(qualified$riv)) & qualified$riv >= 0
  } else {
    rep(TRUE, nrow(qualified))
  }
  valid_interval <- is.finite(qualified$conf.low) &
    is.finite(qualified$conf.high)
  coverage_denominator <- sum(valid_interval)
  coverage <- if (coverage_denominator > 0L) {
    mean(qualified$conf.low[valid_interval] <= truth &
           qualified$conf.high[valid_interval] >= truth)
  } else NA_real_
  n_success <- nrow(qualified)
  success_rate <- n_success / counts[["n_processed"]]
  bias <- if (n_success > 0L) mean(estimates, na.rm = TRUE) - truth else NA_real_
  standardized_bias <- abs(bias) / empirical_sd
  oracle_empirical_sd <- if (nrow(qualified) > 1L) {
    stats::sd(qualified$oracle_estimate, na.rm = TRUE)
  } else NA_real_
  added_bias <- if (n_success > 0L) {
    mean(qualified$estimate - qualified$oracle_estimate, na.rm = TRUE)
  } else NA_real_
  standardized_added_bias <- abs(added_bias) / oracle_empirical_sd
  bias_gate_value <- if (mi_method) standardized_added_bias else standardized_bias
  mean_se <- if (n_success > 0L) {
    mean(qualified$std.error, na.rm = TRUE)
  } else NA_real_
  se_ratio <- mean_se / empirical_sd
  coverage_mcse <- if (coverage_denominator > 0L && is.finite(coverage)) {
    sqrt(coverage * (1 - coverage) / coverage_denominator)
  } else NA_real_
  planned_coverage_mcse <- if (counts[["n_planned"]] > 0L) {
    sqrt(0.925 * 0.075 / counts[["n_planned"]])
  } else NA_real_
  coverage_successes <- if (coverage_denominator > 0L) {
    sum(qualified$conf.low[valid_interval] <= truth &
          qualified$conf.high[valid_interval] >= truth)
  } else 0L
  coverage_ci <- if (coverage_denominator > 0L) {
    stats::binom.test(coverage_successes, coverage_denominator)$conf.int
  } else c(NA_real_, NA_real_)
  finite_rate <- if (n_success > 0L) mean(finite & rubin_valid) else 0
  evidence_complete <- counts[["n_processed"]] >= 1000L &&
    counts[["n_processed"]] == counts[["n_planned"]]

  data.frame(
    dgp = data$dgp[[1]], regime = data$regime[[1]],
    method = data$method[[1]], term = data$term[[1]],
    n_planned = counts[["n_planned"]],
    n_processed = counts[["n_processed"]], n_success = n_success,
    success_rate = success_rate,
    mean_downstream_fit_rate = mean(data$n_fits / data$m_requested),
    bias = bias,
    empirical_sd = empirical_sd,
    standardized_bias = standardized_bias,
    added_bias = added_bias,
    oracle_empirical_sd = oracle_empirical_sd,
    standardized_added_bias = standardized_added_bias,
    mean_se = mean_se, se_ratio = se_ratio,
    coverage = coverage, coverage_mcse = coverage_mcse,
    coverage_denominator = coverage_denominator,
    planned_coverage_mcse = planned_coverage_mcse,
    coverage_ci_low = coverage_ci[[1]], coverage_ci_high = coverage_ci[[2]],
    finite_valid_rate = finite_rate,
    evidence_complete = evidence_complete,
    pilot_criteria_pass = success_rate >= 0.95 &&
      is.finite(bias_gate_value) && bias_gate_value <= 0.10 &&
      is.finite(se_ratio) && se_ratio >= 0.90 && se_ratio <= 1.10 &&
      is.finite(coverage) && coverage >= 0.925 && coverage <= 0.975 &&
      is.finite(finite_rate) && finite_rate >= 0.99,
    pass = evidence_complete && success_rate >= 0.95 &&
      is.finite(bias_gate_value) && bias_gate_value <= 0.10 &&
      is.finite(se_ratio) && se_ratio >= 0.90 && se_ratio <= 1.10 &&
      is.finite(coverage) && coverage >= 0.925 && coverage <= 0.975 &&
      is.finite(planned_coverage_mcse) && planned_coverage_mcse <= 0.01 &&
      is.finite(finite_rate) && finite_rate >= 0.99,
    stringsAsFactors = FALSE
  )
}

fixed_summary <- data.frame()
if (nrow(fixed) > 0L) {
  groups <- split(fixed, interaction(fixed$dgp, fixed$regime, fixed$method,
                                     fixed$term, drop = TRUE))
  fixed_summary <- do.call(rbind, lapply(groups, .summarise_fixed_group))
  rownames(fixed_summary) <- NULL
}

.paired_coverage_group <- function(data) {
  oracle <- data[data$method == "oracle_conditional",
                 c("replicate", "conf.low", "conf.high", "truth"),
                 drop = FALSE]
  standard <- data[data$method == "standard_smc",
                   c("replicate", "conf.low", "conf.high", "truth"),
                   drop = FALSE]
  names(oracle)[2:4] <- c("oracle_low", "oracle_high", "oracle_truth")
  names(standard)[2:4] <- c("standard_low", "standard_high", "standard_truth")
  paired <- merge(oracle, standard, by = "replicate", sort = FALSE)
  if (nrow(paired) == 0L) return(NULL)
  oracle_cover <- paired$oracle_low <= paired$oracle_truth &
    paired$oracle_high >= paired$oracle_truth
  standard_cover <- paired$standard_low <= paired$standard_truth &
    paired$standard_high >= paired$standard_truth
  oracle_only <- sum(oracle_cover & !standard_cover)
  standard_only <- sum(!oracle_cover & standard_cover)
  discordant <- oracle_only + standard_only
  data.frame(
    dgp = data$dgp[[1]], regime = data$regime[[1]], term = data$term[[1]],
    n_paired = nrow(paired),
    oracle_coverage = mean(oracle_cover),
    standard_coverage = mean(standard_cover),
    paired_difference = mean(standard_cover - oracle_cover),
    both_cover = sum(oracle_cover & standard_cover),
    oracle_only = oracle_only, standard_only = standard_only,
    neither_cover = sum(!oracle_cover & !standard_cover),
    mcnemar_exact_p = if (discordant > 0L) {
      stats::binom.test(standard_only, discordant, p = 0.5)$p.value
    } else 1,
    stringsAsFactors = FALSE
  )
}

paired_coverage <- data.frame()
if (nrow(fixed) > 0L) {
  paired_rows <- fixed[fixed$method %in% control_methods &
                         fixed$term %in% c("x", "z"), , drop = FALSE]
  if (nrow(paired_rows) > 0L) {
    paired_groups <- split(
      paired_rows,
      interaction(paired_rows$dgp, paired_rows$regime, paired_rows$term,
                  drop = TRUE)
    )
    paired_coverage_result <- do.call(
      rbind, lapply(paired_groups, .paired_coverage_group)
    )
    if (!is.null(paired_coverage_result)) {
      paired_coverage <- paired_coverage_result
      rownames(paired_coverage) <- NULL
    }
  }
}

.summarise_variance_group <- function(data) {
  truth <- unique(data$truth)
  truth <- truth[is.finite(truth)][1]
  estimate <- mean(data$estimate, na.rm = TRUE)
  data.frame(
    dgp = data$dgp[[1]], regime = data$regime[[1]],
    method = data$method[[1]], component = data$component[[1]],
    n_success = nrow(data), estimate = estimate, truth = truth,
    relative_bias = (estimate - truth) / truth,
    boundary_rate = mean(data$boundary_rate, na.rm = TRUE),
    singular_rate = mean(data$singular_rate, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

variance_summary <- data.frame()
if (nrow(variance) > 0L) {
  groups <- split(variance, interaction(variance$dgp, variance$regime,
                                        variance$method, variance$component,
                                        drop = TRUE))
  variance_summary <- do.call(rbind, lapply(groups, .summarise_variance_group))
  rownames(variance_summary) <- NULL
  oracle <- variance_summary[variance_summary$method == "oracle",
                             c("dgp", "regime", "component", "relative_bias",
                               "boundary_rate"), drop = FALSE]
  names(oracle)[4:5] <- c("oracle_relative_bias", "oracle_boundary_rate")
  variance_summary <- merge(variance_summary, oracle,
                            by = c("dgp", "regime", "component"), all.x = TRUE,
                            sort = FALSE)
  variance_summary$added_abs_relative_bias <-
    abs(variance_summary$relative_bias) - abs(variance_summary$oracle_relative_bias)
  variance_summary$boundary_rate_increase <-
    variance_summary$boundary_rate - variance_summary$oracle_boundary_rate
  variance_summary$diagnostic_flag <-
    variance_summary$method %in% mi_methods &
    (variance_summary$added_abs_relative_bias > 0.20 |
       variance_summary$boundary_rate_increase > 0.05)
}

gate_summary <- data.frame()
if (nrow(gates) > 0L) {
  groups <- split(gates, interaction(gates$dgp, gates$regime, drop = TRUE))
  gate_summary <- do.call(rbind, lapply(groups, function(data) {
    data.frame(
      dgp = data$dgp[[1]], regime = data$regime[[1]], n = nrow(data),
      mean_r_bm = mean(data$r_bm, na.rm = TRUE),
      mean_r_gnn = mean(data$r_gnn, na.rm = TRUE),
      mean_r_mean = mean(data$r_mean, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(gate_summary) <- NULL
}

training_summary <- data.frame()
if (nrow(training) > 0L) {
  usable <- training[is.finite(training$epochs_run) &
                       !is.na(training$suggests_more_epochs), , drop = FALSE]
  if (nrow(usable) > 0L) {
    elapsed <- processed$elapsed_seconds[
      match(usable$task_id, processed$task_id)
    ]
    median_elapsed <- stats::median(elapsed, na.rm = TRUE)
    training_summary <- data.frame(
      n_models = nrow(usable),
      median_epochs_run = stats::median(usable$epochs_run),
      median_elapsed_seconds = median_elapsed,
      projected_full_hours_96 = median_elapsed * 6000 * 1.20 / 96 / 3600,
      proportion_reaching_cap = mean(usable$reached_epoch_cap),
      proportion_suggesting_more = mean(usable$suggests_more_epochs),
      escalate_above_500 = mean(usable$suggests_more_epochs) > 0.05,
      stringsAsFactors = FALSE
    )
  }
}

core <- fixed_summary[
  fixed_summary$method %in% mi_methods &
    fixed_summary$term %in% c("x", "z"), , drop = FALSE
]
expected_cells <- length(unique(manifest$dgp)) * length(unique(manifest$regime)) * 2L
method_pass <- vapply(mi_methods, function(method) {
  cells <- core[core$method == method, , drop = FALSE]
  nrow(cells) == expected_cells && all(cells$pass)
}, logical(1))
git_shas <- vapply(results, function(x) {
  x$provenance$git_sha %||% NA_character_
}, character(1))
git_dirty <- vapply(results, function(x) {
  x$provenance$git_dirty %||% NA
}, logical(1))
sha_consistent <- !anyNA(git_shas) && length(unique(git_shas)) == 1L &&
  !anyNA(git_dirty) && !any(git_dirty)
setting_names <- c(
  "smcfcs_numit", "smcfcs_rjlimit", "jomo_nburn", "jomo_nbetween"
)
if (all(setting_names %in% names(manifest))) {
  settings <- list(
    smcfcs = unique(manifest[manifest$dgp != "lmer",
                             c("smcfcs_numit", "smcfcs_rjlimit"), drop = FALSE]),
    jomo = unique(manifest[manifest$dgp == "lmer",
                           c("jomo_nburn", "jomo_nbetween"), drop = FALSE])
  )
  settings_consistent <- nrow(settings$smcfcs) == 1L &&
    nrow(settings$jomo) == 1L
} else {
  settings <- list(
    smcfcs = data.frame(smcfcs_numit = NA_integer_,
                        smcfcs_rjlimit = NA_integer_),
    jomo = data.frame(jomo_nburn = NA_integer_, jomo_nbetween = NA_integer_)
  )
  settings_consistent <- FALSE
}
dependency_version <- function(package) {
  vapply(results, function(x) {
    versions <- x$provenance$dependency_versions %||% character(0)
    unname(versions[[package]] %||% NA_character_)
  }, character(1))
}
comparator_versions <- list(
  smcfcs = unique(dependency_version("smcfcs")),
  jomo = unique(dependency_version("jomo"))
)
versions_consistent <- all(vapply(comparator_versions, function(x) {
  length(x) == 1L && !is.na(x)
}, logical(1)))
warning_free <- nrow(comparator_diagnostics) == nrow(processed) &&
  sum(comparator_diagnostics$warning_count) == 0L
cell_sizes <- table(manifest$dgp, manifest$regime)
full_grid <- setequal(unique(manifest$dgp), c("lm", "glm", "lmer")) &&
  setequal(unique(manifest$regime), c("phylogeny", "auxiliary"))
campaign_complete <- nrow(processed) == nrow(manifest) &&
  length(unique(processed$task_id)) == nrow(manifest) &&
  all(cell_sizes >= 1000L) && sha_consistent && settings_consistent &&
  versions_consistent && warning_free && full_grid

oracle_attainable <- isTRUE(method_pass[["oracle_conditional"]])
standard_smc_ready <- isTRUE(method_pass[["standard_smc"]])
decision <- if (calibration_only || !full_grid) {
  "CALIBRATION ONLY: restricted evidence cannot support a release decision"
} else if (!campaign_complete) {
  "INCOMPLETE: pilot/dry-run evidence cannot support a release decision"
} else if (any(method_pass[release_methods])) {
  paste0(
    "PASS ", paste(release_methods[method_pass[release_methods]],
                    collapse = ", "),
    " ONLY: unsupported draw methods remain experimental"
  )
} else if (oracle_attainable && standard_smc_ready) {
  paste0(
    "CONTROL PASS ONLY: gate is attainable; block CRAN pending an ",
    "analysis-aware inferential backend"
  )
} else if (oracle_attainable && !standard_smc_ready) {
  "ORACLE PASS / STANDARD SMC FAIL: gate attainable; comparator not ready"
} else if (!oracle_attainable && standard_smc_ready) {
  "STANDARD SMC PASS / ORACLE FAIL: audit the oracle and harness"
} else {
  "FAIL: neither positive control passes the MI gate"
}

summary_object <- list(
  generated_at = Sys.time(), config = config, manifest = manifest,
  processed = processed, fixed_effects = fixed_summary,
  paired_coverage = paired_coverage,
  comparator_diagnostics = comparator_diagnostics,
  variance_components = variance_summary, gate_weights = gate_summary,
  training = training_summary,
  git_shas = unique(git_shas), git_dirty = unique(git_dirty),
  sha_consistent = sha_consistent, comparator_settings = settings,
  settings_consistent = settings_consistent,
  comparator_versions = comparator_versions,
  versions_consistent = versions_consistent,
  warning_free = warning_free,
  method_pass = method_pass, campaign_complete = campaign_complete,
  oracle_attainable = oracle_attainable,
  standard_smc_ready = standard_smc_ready,
  decision = decision
)
atomic_save_rds(summary_object, file.path(config$output_dir, "summary.rds"))

cat("pigauto v0.10.0 MI validation summary\n")
cat("Processed:", nrow(processed), "of", nrow(manifest), "tasks\n")
cat("Decision:", decision, "\n\n")
if (nrow(training_summary) > 0L) {
  cat(sprintf(
    "Epoch diagnostic: %.1f%% suggest >500 epochs; escalation %s\n\n",
    100 * training_summary$proportion_suggesting_more,
    if (training_summary$escalate_above_500) "triggered" else "not triggered"
  ))
}
if (nrow(core) > 0L) {
  print(core[, c("dgp", "regime", "method", "term", "n_processed",
                 "success_rate", "standardized_bias",
                 "standardized_added_bias", "se_ratio", "coverage",
                 "coverage_mcse", "planned_coverage_mcse",
                 "coverage_denominator", "coverage_ci_low", "coverage_ci_high",
                 "finite_valid_rate",
                 "pilot_criteria_pass",
                 "pass")],
        row.names = FALSE)
}
if (nrow(paired_coverage) > 0L) {
  cat("\nPaired standard-SMC minus oracle-conditional coverage:\n")
  print(paired_coverage, row.names = FALSE)
}
if (nrow(comparator_diagnostics) > 0L) {
  cat("\nComparator diagnostics:\n")
  cat("Warnings:", sum(comparator_diagnostics$warning_count), "across",
      nrow(comparator_diagnostics), "tasks\n")
  cat("Median engine seconds:",
      stats::median(comparator_diagnostics$elapsed_seconds, na.rm = TRUE),
      "\n")
}
if (nrow(variance_summary) > 0L && any(variance_summary$diagnostic_flag)) {
  cat("\nVariance-component diagnostic flags (descriptive; not Rubin-pooled):\n")
  print(variance_summary[variance_summary$diagnostic_flag, ], row.names = FALSE)
}
