#!/usr/bin/env Rscript

# Run one or more manifest tasks for the pigauto v0.10.0 MI validation.
# Each task trains pigauto exactly once, then derives paired conformal and
# Brownian-posterior/MC-dropout imputations from that same fitted object.

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else ""
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path)) else getwd()
source(file.path(script_dir, "0_prepare.R"))
source(file.path(script_dir, "comparator_helpers.R"))

config <- mi_validation_config(args)
manifest <- write_manifest(config)

.has_flag <- function(flag) paste0("--", flag) %in% args

.task_from_args <- function() {
  from_arg <- .arg_value(args, "task", NULL)
  from_env <- .env_value("PIGAUTO_MI_TASK_ID",
                         .env_value("SLURM_ARRAY_TASK_ID", NULL))
  value <- from_arg %||% from_env
  if (is.null(value)) return(NULL)
  .as_int(value, "task", 1L)
}

.load_package <- function() {
  if ("pigauto" %in% loadedNamespaces()) return(invisible(TRUE))
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(getwd(), quiet = TRUE)
  } else if (!requireNamespace("pigauto", quietly = TRUE)) {
    stop("Install pigauto or devtools before running the validation harness.",
         call. = FALSE)
  } else {
    loadNamespace("pigauto")
  }
  invisible(TRUE)
}

.provenance <- function() {
  git_sha <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1],
    error = function(e) NA_character_
  )
  git_dirty <- tryCatch(
    length(system2("git", c("status", "--porcelain"), stdout = TRUE,
                   stderr = FALSE)) > 0L,
    error = function(e) NA
  )
  description <- read.dcf("DESCRIPTION", fields = "Version")
  package_version <- unname(description[1, "Version"])
  dependencies <- c("ape", "lme4", "torch", "smcfcs", "jomo")
  versions <- vapply(dependencies, function(package) {
    if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
    as.character(utils::packageVersion(package))
  }, character(1))
  list(
    git_sha = git_sha,
    git_dirty = git_dirty,
    package_version = package_version,
    r_version = R.version.string,
    platform = R.version$platform,
    dependency_versions = versions
  )
}

.fixed_effect_table <- function(fit, method, truth) {
  estimates <- if (inherits(fit, "merMod")) lme4::fixef(fit) else stats::coef(fit)
  covariance <- as.matrix(stats::vcov(fit))
  terms <- intersect(names(truth), names(estimates))
  se <- sqrt(diag(covariance))[terms]
  if (inherits(fit, "lm") && !inherits(fit, "glm")) {
    df <- stats::df.residual(fit)
    critical <- stats::qt(0.975, df = df)
  } else {
    df <- Inf
    critical <- stats::qnorm(0.975)
  }
  data.frame(
    method = method,
    term = terms,
    estimate = unname(estimates[terms]),
    std.error = unname(se),
    df = df,
    conf.low = unname(estimates[terms] - critical * se),
    conf.high = unname(estimates[terms] + critical * se),
    fmi = NA_real_,
    riv = NA_real_,
    truth = unname(truth[terms]),
    n_fits = 1L,
    stringsAsFactors = FALSE
  )
}

.fit_many <- function(datasets, dgp) {
  fits <- lapply(datasets, function(data) {
    tryCatch(fit_downstream(data, dgp), error = function(e) e)
  })
  ok <- !vapply(fits, inherits, logical(1), "error")
  list(fits = fits[ok], n_attempted = length(fits),
       errors = vapply(fits[!ok], conditionMessage, character(1)))
}

.pool_fixed_effects <- function(fits, dgp, method, truth) {
  if (length(fits) < 2L) {
    stop(method, " produced fewer than two successful downstream fits.",
         call. = FALSE)
  }
  if (dgp == "lmer") {
    pooled <- pigauto::pool_mi(
      fits,
      coef_fun = lme4::fixef,
      vcov_fun = function(fit) as.matrix(stats::vcov(fit))
    )
  } else {
    pooled <- pigauto::pool_mi(
      fits,
      coef_fun = stats::coef,
      vcov_fun = stats::vcov,
      df_fun = stats::df.residual
    )
  }
  pooled <- as.data.frame(pooled)
  pooled <- pooled[pooled$term %in% names(truth), , drop = FALSE]
  pooled$method <- method
  pooled$truth <- unname(truth[pooled$term])
  pooled$n_fits <- length(fits)
  pooled[, c("method", "term", "estimate", "std.error", "df", "conf.low",
              "conf.high", "fmi", "riv", "truth", "n_fits")]
}

.analyse_method <- function(datasets, dgp, method, truth, truth_vc) {
  fitted <- .fit_many(datasets, dgp)
  pooled <- tryCatch(
    .pool_fixed_effects(fitted$fits, dgp, method, truth),
    error = function(e) e
  )
  list(
    fixed_effects = if (inherits(pooled, "error")) data.frame() else pooled,
    variance_components = if (dgp == "lmer" && length(fitted$fits) > 0L) {
      .variance_components(fitted$fits, method, truth_vc)
    } else data.frame(),
    fit_errors = fitted$errors,
    pool_error = if (inherits(pooled, "error")) conditionMessage(pooled) else NULL
  )
}

.variance_components <- function(fits, method, truth_vc) {
  if (length(fits) == 0L || length(truth_vc) == 0L) return(data.frame())
  one <- lapply(fits, function(fit) {
    vc <- as.data.frame(lme4::VarCorr(fit))
    c(
      tau2 = vc$vcov[vc$grp == "species" & is.na(vc$var2)][1],
      sigma2 = stats::sigma(fit)^2,
      singular = as.numeric(lme4::isSingular(fit, tol = 1e-5))
    )
  })
  values <- do.call(rbind, one)
  data.frame(
    method = method,
    component = names(truth_vc),
    estimate = colMeans(values[, names(truth_vc), drop = FALSE], na.rm = TRUE),
    between_sd = apply(values[, names(truth_vc), drop = FALSE], 2L,
                       stats::sd, na.rm = TRUE),
    boundary_rate = vapply(names(truth_vc), function(component) {
      mean(values[, component] <= 1e-8, na.rm = TRUE)
    }, numeric(1)),
    singular_rate = mean(values[, "singular"], na.rm = TRUE),
    truth = unname(truth_vc),
    n_fits = nrow(values),
    stringsAsFactors = FALSE
  )
}

.conformal_datasets_from_fit <- function(mi, original_traits, species_col,
                                         m, seed) {
  pred <- stats::predict(mi$fit, return_se = TRUE, n_imputations = 1L)
  sample_draw <- getFromNamespace(".sample_conformal_draw", "pigauto")
  build_completed <- getFromNamespace("build_completed", "pigauto")
  lapply(seq_len(m), function(i) {
    imputed <- sample_draw(
      pred = pred,
      imputed_mask = mi$imputed_mask,
      trait_map = mi$fit$trait_map,
      seed_i = as.integer(seed + 50000L + i),
      input_row_order = mi$data$input_row_order
    )
    build_completed(
      original_traits, imputed, species_col,
      input_row_order = mi$data$input_row_order
    )$completed
  })
}

.pmm_datasets_from_fit <- function(mi, original_traits, species_col,
                                    m, seed, k = 5L) {
  pred <- stats::predict(
    mi$fit, return_se = FALSE, n_imputations = m,
    match_observed = "none"
  )
  if (is.null(pred$imputed_datasets) ||
      length(pred$imputed_datasets) != m) {
    stop("PMM validation did not receive ", m, " stochastic predictions.",
         call. = FALSE)
  }

  input_row_order <- mi$data$input_row_order
  truth_input <- as.numeric(original_traits$x)
  truth_internal <- if (!is.null(input_row_order)) {
    truth_input[input_row_order]
  } else {
    truth_input
  }
  pmm_one <- getFromNamespace("pmm_impute_one_trait", "pigauto")
  build_completed <- getFromNamespace("build_completed", "pigauto")

  lapply(seq_len(m), function(i) {
    imputed <- pred$imputed_datasets[[i]]
    imputed$x <- pmm_one(
      predictions = as.numeric(imputed$x), truth = truth_internal,
      K = k, seed = as.integer(seed + 70000L + i)
    )
    build_completed(
      original_traits, imputed, species_col,
      input_row_order = input_row_order
    )$completed
  })
}

.x_gate <- function(fit) {
  tm_index <- which(vapply(fit$trait_map, function(x) identical(x$name, "x"),
                           logical(1)))
  if (length(tm_index) != 1L) {
    return(data.frame(r_bm = NA_real_, r_gnn = NA_real_, r_mean = NA_real_))
  }
  latent_col <- fit$trait_map[[tm_index]]$latent_cols[[1]]
  value <- function(slot) {
    x <- fit[[slot]]
    if (is.null(x) || length(x) < latent_col) NA_real_ else as.numeric(x[latent_col])
  }
  data.frame(r_bm = value("r_cal_bm"), r_gnn = value("r_cal_gnn"),
             r_mean = value("r_cal_mean"))
}

.training_diagnostic <- function(fit, epochs_requested) {
  history <- fit$history
  if (is.null(history) || nrow(history) == 0L) {
    return(data.frame(
      epochs_requested = epochs_requested, epochs_run = epochs_requested,
      reached_epoch_cap = TRUE, recent_relative_improvement = NA_real_,
      suggests_more_epochs = FALSE
    ))
  }
  finite <- history[is.finite(history$val_loss), c("epoch", "val_loss"),
                    drop = FALSE]
  epochs_run <- max(history$epoch)
  reached_cap <- epochs_run >= epochs_requested
  recent_improvement <- NA_real_
  if (nrow(finite) >= 4L) {
    old <- finite$val_loss[nrow(finite) - 3L]
    new <- finite$val_loss[nrow(finite)]
    recent_improvement <- (old - new) / max(abs(old), .Machine$double.eps)
  }
  data.frame(
    epochs_requested = epochs_requested,
    epochs_run = epochs_run,
    reached_epoch_cap = reached_cap,
    recent_relative_improvement = recent_improvement,
    suggests_more_epochs = reached_cap && is.finite(recent_improvement) &&
      recent_improvement > 0.01
  )
}

.run_task <- function(task, dry_run = FALSE, controls_only = FALSE,
                      standard_smc_only = FALSE) {
  started <- Sys.time()
  dgp <- simulate_validation_dgp(
    task$dgp, task$regime, task$seed,
    missing_fraction = task$missing_fraction
  )
  oracle_fit <- fit_downstream(dgp$full, task$dgp)
  complete_case_fit <- fit_downstream(dgp$observed, task$dgp)
  fixed <- rbind(
    .fixed_effect_table(oracle_fit, "oracle", dgp$truth),
    .fixed_effect_table(complete_case_fit, "complete_case", dgp$truth)
  )
  variance <- rbind(
    .variance_components(list(oracle_fit), "oracle", dgp$truth_vc),
    .variance_components(list(complete_case_fit), "complete_case", dgp$truth_vc)
  )

  result <- list(
    meta = as.list(task),
    provenance = .provenance(),
    status = if (dry_run) "dry_run" else "running",
    error = NULL,
    fixed_effects = fixed,
    variance_components = variance,
    gate = data.frame(r_bm = NA_real_, r_gnn = NA_real_, r_mean = NA_real_),
    training = data.frame(
      epochs_requested = task$epochs, epochs_run = NA_integer_,
      reached_epoch_cap = NA, recent_relative_improvement = NA_real_,
      suggests_more_epochs = NA
    ),
    comparator_diagnostics = list(),
    downstream_failures = list(),
    missing_rows = dgp$missing_rows,
    started_at = started,
    completed_at = Sys.time(),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )
  if (dry_run) return(result)

  .load_package()
  standard_smc <- .draw_standard_smc_imputations(
    dgp, m = task$m, seed = task$seed + 90000L,
    smcfcs_numit = task$smcfcs_numit,
    jomo_nburn = task$jomo_nburn,
    jomo_nbetween = task$jomo_nbetween
  )
  standard_smc_analysis <- .analyse_method(
    standard_smc$datasets, task$dgp, "standard_smc",
    dgp$truth, dgp$truth_vc
  )
  oracle_conditional_analysis <- NULL
  if (!standard_smc_only) {
    oracle_conditional_datasets <- .draw_oracle_conditional_imputations(
      dgp, m = task$m, seed = task$seed + 80000L
    )
    oracle_conditional_analysis <- .analyse_method(
      oracle_conditional_datasets, task$dgp, "oracle_conditional",
      dgp$truth, dgp$truth_vc
    )
  }
  result$fixed_effects <- rbind(
    result$fixed_effects,
    if (is.null(oracle_conditional_analysis)) NULL else
      oracle_conditional_analysis$fixed_effects,
    standard_smc_analysis$fixed_effects
  )
  if (task$dgp == "lmer") {
    result$variance_components <- rbind(
      result$variance_components,
      if (is.null(oracle_conditional_analysis)) NULL else
        oracle_conditional_analysis$variance_components,
      standard_smc_analysis$variance_components
    )
  }
  result$comparator_diagnostics <- standard_smc$diagnostics
  result$downstream_failures <- list(
    standard_smc = list(
      fit_errors = standard_smc_analysis$fit_errors,
      pool_error = standard_smc_analysis$pool_error
    )
  )
  if (!is.null(oracle_conditional_analysis)) {
    result$downstream_failures$oracle_conditional <- list(
      fit_errors = oracle_conditional_analysis$fit_errors,
      pool_error = oracle_conditional_analysis$pool_error
    )
  }
  controls_ok <- c(
    standard_smc = nrow(standard_smc_analysis$fixed_effects) == length(dgp$truth)
  )
  if (!is.null(oracle_conditional_analysis)) {
    controls_ok <- c(
      controls_ok,
      oracle_conditional = nrow(oracle_conditional_analysis$fixed_effects) ==
        length(dgp$truth)
    )
  }
  if (controls_only) {
    result$status <- if (all(controls_ok)) "success" else if (any(controls_ok)) {
      "partial"
    } else {
      "analysis_error"
    }
    result$completed_at <- Sys.time()
    result$elapsed_seconds <- as.numeric(difftime(
      result$completed_at, started, units = "secs"
    ))
    return(result)
  }

  if (!requireNamespace("torch", quietly = TRUE) ||
      !isTRUE(torch::torch_is_installed())) {
    stop("A working torch R package and libtorch backend are required. ",
         "Use --dry-run to validate the manifest and DGPs without training.",
         call. = FALSE)
  }

  mi_mc <- pigauto::multi_impute(
    traits = dgp$traits,
    tree = dgp$tree,
    m = task$m,
    draws_method = "mc_dropout",
    species_col = dgp$species_col,
    log_transform = FALSE,
    missing_frac = 0.25,
    covariates = dgp$covariates,
    epochs = task$epochs,
    verbose = FALSE,
    seed = task$seed,
    eval_every = config$eval_every,
    patience = config$patience
  )

  mc_datasets <- lapply(mi_mc$datasets, complete_imputed_dataset,
                        dgp_object = dgp)
  conformal_datasets <- .conformal_datasets_from_fit(
    mi_mc, dgp$traits, dgp$species_col, task$m, task$seed
  )
  conformal_datasets <- lapply(conformal_datasets, complete_imputed_dataset,
                               dgp_object = dgp)
  pmm_datasets <- .pmm_datasets_from_fit(
    mi_mc, dgp$traits, dgp$species_col, task$m, task$seed
  )
  pmm_datasets <- lapply(pmm_datasets, complete_imputed_dataset,
                         dgp_object = dgp)

  conformal_analysis <- .analyse_method(
    conformal_datasets, task$dgp, "conformal", dgp$truth, dgp$truth_vc
  )
  mc_analysis <- .analyse_method(
    mc_datasets, task$dgp, "mc_dropout", dgp$truth, dgp$truth_vc
  )
  pmm_analysis <- .analyse_method(
    pmm_datasets, task$dgp, "pmm", dgp$truth, dgp$truth_vc
  )
  fixed <- rbind(
    result$fixed_effects,
    conformal_analysis$fixed_effects,
    mc_analysis$fixed_effects,
    pmm_analysis$fixed_effects
  )
  if (task$dgp == "lmer") {
    variance <- rbind(
      result$variance_components,
      conformal_analysis$variance_components,
      mc_analysis$variance_components,
      pmm_analysis$variance_components
    )
  }

  method_ok <- c(
    conformal = nrow(conformal_analysis$fixed_effects) == length(dgp$truth),
    mc_dropout = nrow(mc_analysis$fixed_effects) == length(dgp$truth),
    pmm = nrow(pmm_analysis$fixed_effects) == length(dgp$truth),
    controls_ok
  )
  result$status <- if (all(method_ok)) "success" else if (any(method_ok)) {
    "partial"
  } else {
    "analysis_error"
  }
  result$fixed_effects <- fixed
  result$variance_components <- variance
  result$gate <- .x_gate(mi_mc$fit)
  result$training <- .training_diagnostic(mi_mc$fit, task$epochs)
  result$downstream_failures <- c(result$downstream_failures, list(
    conformal = list(fit_errors = conformal_analysis$fit_errors,
                     pool_error = conformal_analysis$pool_error),
    mc_dropout = list(fit_errors = mc_analysis$fit_errors,
                      pool_error = mc_analysis$pool_error),
    pmm = list(fit_errors = pmm_analysis$fit_errors,
               pool_error = pmm_analysis$pool_error)
  ))
  result$completed_at <- Sys.time()
  result$elapsed_seconds <- as.numeric(difftime(result$completed_at, started,
                                                units = "secs"))
  result
}

task_id <- .task_from_args()
if (is.null(task_id) && !.has_flag("all-pending")) {
  stop("Specify --task=N (or PIGAUTO_MI_TASK_ID/SLURM_ARRAY_TASK_ID), ",
       "or use --all-pending for a sequential local run.", call. = FALSE)
}
if (!is.null(task_id) && task_id > nrow(manifest)) {
  stop("task ", task_id, " is outside manifest range 1..", nrow(manifest), ".",
       call. = FALSE)
}

selected <- if (!is.null(task_id)) task_id else manifest$task_id
dry_run <- .has_flag("dry-run") ||
  tolower(.env_value("PIGAUTO_MI_DRY_RUN", "false")) %in% c("true", "1", "yes")
force <- .has_flag("force")
controls_only <- .has_flag("controls-only") ||
  tolower(.env_value("PIGAUTO_MI_CONTROLS_ONLY", "false")) %in%
    c("true", "1", "yes")
standard_smc_only <- .has_flag("standard-smc-only") ||
  tolower(.env_value("PIGAUTO_MI_STANDARD_SMC_ONLY", "false")) %in%
    c("true", "1", "yes")
if (standard_smc_only) controls_only <- TRUE

for (id in selected) {
  path <- file.path(config$output_dir, sprintf("task-%06d.rds", id))
  if (file.exists(path) && !force) {
    cat("Skipping existing", path, "\n")
    next
  }
  task <- manifest[manifest$task_id == id, , drop = FALSE]
  cat(sprintf("Task %d: %s / %s / replicate %d%s\n", id, task$dgp,
              task$regime, task$replicate,
              if (dry_run) " [dry-run]" else ""))
  result <- tryCatch(
    .run_task(
      task, dry_run = dry_run, controls_only = controls_only,
      standard_smc_only = standard_smc_only
    ),
    error = function(e) {
      list(
        meta = as.list(task), status = "error",
        provenance = .provenance(),
        error = conditionMessage(e), fixed_effects = data.frame(),
        variance_components = data.frame(), gate = data.frame(),
        training = data.frame(),
        comparator_diagnostics = list(),
        downstream_failures = list(), missing_rows = integer(0),
        started_at = Sys.time(), completed_at = Sys.time(),
        elapsed_seconds = NA_real_
      )
    }
  )
  atomic_save_rds(result, path)
  cat("  ", result$status,
      if (!is.null(result$error)) paste0(": ", result$error) else "", "\n",
      sep = "")
}
