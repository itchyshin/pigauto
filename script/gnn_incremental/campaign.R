#!/usr/bin/env Rscript

# Resumable, provenance-bound runner for the GNN incremental-value study.
# Usage:
#   Rscript script/gnn_incremental/campaign.R init RUN_DIR N_REPS MASTER EPOCHS
#   Rscript script/gnn_incremental/campaign.R run RUN_DIR N_REPS MASTER EPOCHS
#   Rscript script/gnn_incremental/campaign.R task RUN_DIR TASK_INDEX

read_cli <- function(args) {
  if (length(args) < 1L || !args[1L] %in% c("init", "run", "task")) {
    stop("Use init, run, or task mode. See script header.", call. = FALSE)
  }
  mode <- args[1L]
  if (mode == "task") {
    if (length(args) != 3L) stop("task needs RUN_DIR and TASK_INDEX.")
    return(list(mode = mode, run_dir = args[2L], task_index = as.integer(args[3L])))
  }
  if (length(args) != 5L) stop(mode, " needs RUN_DIR N_REPS MASTER EPOCHS.")
  list(mode = mode, run_dir = args[2L], n_replicates = as.integer(args[3L]),
       master_seed = as.integer(args[4L]), epochs = as.integer(args[5L]))
}

script_dir <- "script/gnn_incremental"
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE)
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
} else if (requireNamespace("pigauto", quietly = TRUE)) {
  # Cluster workers may not have development tooling. In that case the job
  # installs this exact source bundle before invoking the runner.
  library("pigauto")
} else {
  stop("Install the source bundle, or install devtools/pkgload to load it.",
       call. = FALSE)
}
if (requireNamespace("torch", quietly = TRUE)) {
  try(torch::torch_set_num_threads(1L), silent = TRUE)
}
source(file.path(script_dir, "dgp.R"))
source(file.path(script_dir, "phylo_gls.R"))
source(file.path(script_dir, "ablation.R"))

hash_files <- function(paths) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The 'digest' package is required for SHA-256 source identity.",
         call. = FALSE)
  }
  paths <- sort(unique(paths[file.exists(paths)]))
  stats::setNames(vapply(paths, function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L)), paths)
}

source_identity <- function() {
  tracked <- list.files("R", pattern = "[.]R$", recursive = TRUE,
                        full.names = TRUE)
  runner <- list.files(script_dir, pattern = "[.](R|sh)$", full.names = TRUE)
  git_head <- tryCatch({
    out <- suppressWarnings(system2("git", c("rev-parse", "--verify", "HEAD"),
                                    stdout = TRUE, stderr = FALSE))
    if (length(out) == 1L && grepl("^[[:xdigit:]]{40}$", out)) out else NA_character_
  }, error = function(e) NA_character_)
  required_sources <- c("script/gnn_incremental/scenario_manifest.csv",
    "docs/dev-log/gnn-incremental/ADEMP.md",
    "tests/testthat/test-phylo-signal-training-mask.R")
  if (length(tracked) == 0L || length(runner) == 0L ||
      !all(file.exists(required_sources))) {
    stop("Required source files for campaign identity are missing.", call. = FALSE)
  }
  list(
    git_head = git_head,
    source_bundle = Sys.getenv("PIGAUTO_SOURCE_BUNDLE", unset = "working-tree"),
    source_sha256 = hash_files(c(tracked, runner,
      "tests/testthat/test-phylo-signal-training-mask.R",
      "script/gnn_incremental/scenario_manifest.csv",
      "docs/dev-log/gnn-incremental/ADEMP.md")),
    R = R.version.string,
    packages = vapply(c("pigauto", "torch", "ape", "phytools", "digest"),
      function(pkg) if (requireNamespace(pkg, quietly = TRUE)) {
        as.character(utils::packageVersion(pkg))
      } else NA_character_, character(1L)),
    library_paths = .libPaths(),
    session = utils::sessionInfo(),
    platform = Sys.info(),
    external_software = extSoftVersion(),
    thread_environment = Sys.getenv(c("OPENBLAS_NUM_THREADS",
      "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"))
  )
}

derive_seed <- function(master_seed, scenario_index, replicate, stream) {
  modulus <- 2147483646
  as.integer((as.double(master_seed) + scenario_index * 1000003 +
    replicate * 7919 + stream * 104729) %% modulus) + 1L
}

atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = ".pending-", tmpdir = dirname(path))
  saveRDS(object, tmp, compress = "xz")
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Atomic rename failed for ", path, call. = FALSE)
  }
  invisible(path)
}

configuration_fingerprint <- function(config) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The 'digest' package is required for result fingerprints.",
         call. = FALSE)
  }
  digest::digest(config, algo = "sha256", serialize = TRUE)
}

verify_config_fingerprint <- function(config) {
  if (is.null(config$fingerprint) || length(config$fingerprint) != 1L) {
    stop("Run config has no valid configuration fingerprint.", call. = FALSE)
  }
  supplied <- config$fingerprint
  unhashed <- config
  unhashed$fingerprint <- NULL
  if (!identical(supplied, configuration_fingerprint(unhashed))) {
    stop("Run configuration fingerprint does not match its settings.",
         call. = FALSE)
  }
  invisible(TRUE)
}

verify_task_manifest <- function(config, tasks) {
  n <- config$n_replicates
  expected <- do.call(rbind, lapply(seq_len(nrow(config$scenario_manifest)),
    function(i) data.frame(task_index = seq_len(n) + (i - 1L) * n,
      scenario_index = i,
      scenario_id = config$scenario_manifest$scenario_id[i],
      replicate = seq_len(n), stringsAsFactors = FALSE)))
  if (!identical(names(tasks), names(expected)) || nrow(tasks) != nrow(expected) ||
      !identical(as.integer(tasks$task_index), expected$task_index) ||
      !identical(as.integer(tasks$scenario_index), expected$scenario_index) ||
      !identical(as.character(tasks$scenario_id), expected$scenario_id) ||
      !identical(as.integer(tasks$replicate), expected$replicate)) {
    stop("Task manifest does not match the canonical scenario/replicate sequence.",
         call. = FALSE)
  }
  invisible(TRUE)
}

training_rows_from_splits <- function(data, splits, latent_col = 1L) {
  X <- data$X_scaled
  mask <- splits$mask
  if (!is.matrix(X) || !is.matrix(mask) || !identical(dim(X), dim(mask)) ||
      length(latent_col) != 1L || is.na(latent_col) ||
      latent_col < 1L || latent_col > ncol(X)) {
    stop("Training-row inputs must have matching X_scaled and split-mask matrices.",
         call. = FALSE)
  }
  # make_missing_splits() returns val/test indices and a mask (TRUE for
  # training-observed cells); it does not expose train_idx. Combine that mask
  # with the original observed values to exclude outer missing cells as well.
  which(is.finite(X[, latent_col]) & mask[, latent_col])
}

init_campaign <- function(run_dir, n_replicates, master_seed, epochs) {
  if (n_replicates < 1L || n_replicates > 500L || master_seed < 1L || epochs < 1L) {
    stop("N_REPS must be 1..500; MASTER and EPOCHS must be positive integers.", call. = FALSE)
  }
  if (master_seed != 20261019L) {
    stop("This frozen study uses master seed 20261019.", call. = FALSE)
  }
  if (!requireNamespace("phytools", quietly = TRUE)) {
    stop("The paired phylo_signal_gate = TRUE arm requires phytools.",
         call. = FALSE)
  }
  manifest <- utils::read.csv(file.path(script_dir, "scenario_manifest.csv"),
                               stringsAsFactors = FALSE)
  run_type <- if (all(manifest$n_replicates == n_replicates)) "production" else "pilot"
  tasks <- do.call(rbind, lapply(seq_len(nrow(manifest)), function(i) {
    data.frame(task_index = seq_len(n_replicates) + (i - 1L) * n_replicates,
      scenario_index = i, scenario_id = manifest$scenario_id[i],
      replicate = seq_len(n_replicates), stringsAsFactors = FALSE)
  }))
  config <- list(
    run_id = basename(normalizePath(run_dir, mustWork = FALSE)),
    run_dir = normalizePath(run_dir, mustWork = FALSE),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    n_replicates = n_replicates, run_type = run_type,
    master_seed = master_seed, epochs = epochs,
    scenario_manifest = manifest,
    gate_policies = c(FALSE, TRUE),
    fit_settings = list(
      gnn = TRUE, missing_frac = 0.25, safety_floor = FALSE,
      phylo_signal_threshold = 0.2, phylo_signal_method = "lambda",
      n_species = 500L, obs_per = 5L, n_covs = 10L, sigma_res = 0.5,
      k_eigen = 15L, predict_method = "per_column",
      gate_method = "cv_folds", gate_cv_folds = 5L,
      max_threads = 1L, hidden_dim = 64L, n_gnn_layers = 2L,
      gate_cap = 0.8, use_attention = TRUE,
      use_transformer_blocks = TRUE, n_heads = 4L, ffn_mult = 4L,
      use_trait_attention = FALSE, n_trait_heads = 2L,
      trait_embed_dim = 32L, dropout = 0.1, lr = 0.003,
      weight_decay = 1e-4, corruption_rate = 0.55,
      corruption_start = 0.20, corruption_ramp = 500L,
      refine_steps = 8L, lambda_shrink = 0.03, lambda_gate = 0.01,
      warmup_epochs = 200L, edge_dropout = 0.1, eval_every = 100L,
      patience = 10L, clip_norm = 1.0, min_val_cells = 20L,
      conformal_method = "split", conformal_split_val = FALSE,
      gate_splits_B = 31L, lambda_mode = "estimate",
      discrete_lambda = "estimate", joint_solver = "inhouse"
    ),
    source = source_identity(),
    session = utils::sessionInfo()
  )
  config$fingerprint <- configuration_fingerprint(config)
  dir.create(file.path(run_dir, "replicates"), recursive = TRUE,
             showWarnings = FALSE)
  config_path <- file.path(run_dir, "run_config.rds")
  tasks_path <- file.path(run_dir, "tasks.csv")
  if (file.exists(config_path) || file.exists(tasks_path)) {
    stop("Run directory already contains a configuration. Choose a new path.",
         call. = FALSE)
  }
  atomic_save_rds(config, config_path)
  utils::write.csv(tasks, tasks_path, row.names = FALSE)
  message("Initialized ", nrow(tasks), " datasets / ", 2L * nrow(tasks),
          " gate-specific fits at ", normalizePath(run_dir))
  invisible(config)
}

fit_prediction <- function(case, gate_policy, fit_seed, epochs, settings) {
  unattempted <- c(pigauto_full = "not_attempted",
    pigauto_neural_off = "not_attempted", mean = "not_attempted",
    lm_linear = "not_attempted", lm_poly_interaction = "not_attempted",
    phylo_gls = "not_attempted")
  traits <- case$df_observed[, c("species", "y"), drop = FALSE]
  covariates <- case$df_observed[, c("species", case$predictor_names),
                                  drop = FALSE]
  started <- proc.time()[["elapsed"]]
  result <- tryCatch(impute(
    traits = traits, tree = case$tree, species_col = "species",
    covariates = covariates, gnn = TRUE,
    phylo_signal_gate = gate_policy,
    phylo_signal_threshold = settings$phylo_signal_threshold,
    phylo_signal_method = settings$phylo_signal_method,
    safety_floor = settings$safety_floor,
    missing_frac = settings$missing_frac,
    epochs = epochs, k_eigen = settings$k_eigen,
    predict_method = settings$predict_method,
    gate_method = settings$gate_method,
    gate_cv_folds = settings$gate_cv_folds,
    hidden_dim = settings$hidden_dim,
    n_gnn_layers = settings$n_gnn_layers,
    gate_cap = settings$gate_cap,
    use_attention = settings$use_attention,
    use_transformer_blocks = settings$use_transformer_blocks,
    n_heads = settings$n_heads, ffn_mult = settings$ffn_mult,
    use_trait_attention = settings$use_trait_attention,
    n_trait_heads = settings$n_trait_heads,
    trait_embed_dim = settings$trait_embed_dim,
    dropout = settings$dropout, lr = settings$lr,
    weight_decay = settings$weight_decay,
    corruption_rate = settings$corruption_rate,
    corruption_start = settings$corruption_start,
    corruption_ramp = settings$corruption_ramp,
    refine_steps = settings$refine_steps,
    lambda_shrink = settings$lambda_shrink,
    lambda_gate = settings$lambda_gate,
    warmup_epochs = settings$warmup_epochs,
    edge_dropout = settings$edge_dropout,
    eval_every = settings$eval_every, patience = settings$patience,
    clip_norm = settings$clip_norm,
    min_val_cells = settings$min_val_cells,
    conformal_method = settings$conformal_method,
    conformal_split_val = settings$conformal_split_val,
    gate_splits_B = settings$gate_splits_B,
    lambda_mode = settings$lambda_mode,
    discrete_lambda = settings$discrete_lambda,
    joint_solver = settings$joint_solver,
    seed = fit_seed, verbose = FALSE
  ), error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(result, "error")) {
    method_status <- unattempted
    method_status["pigauto_full"] <- "failed_fit"
    return(list(status = "failed", failure_reason = conditionMessage(result),
                elapsed_seconds = elapsed, predictions = NULL,
                fit_diagnostics = list(fit_status = "failed_fit",
                  ablation_status = "not_attempted",
                  method_status = method_status)))
  }

  input_row_order <- result$data$input_row_order
  if (length(input_row_order) != nrow(case$df_observed) ||
      anyNA(input_row_order) || anyDuplicated(input_row_order)) {
    method_status <- unattempted
    method_status["pigauto_full"] <- "failed_row_alignment"
    return(list(status = "failed", failure_reason =
      "preprocess input-row mapping is incomplete or duplicated",
      elapsed_seconds = elapsed, predictions = NULL,
      fit_diagnostics = list(fit_status = "failed_row_alignment",
        ablation_status = "not_attempted",
        method_status = method_status)))
  }
  internal_for_input <- match(seq_len(nrow(case$df_observed)), input_row_order)
  pred_full <- result$prediction$imputed$y[internal_for_input]
  ablated <- tryCatch(ablate_gnn_in_fit(result$fit), error = function(e) e)
  ablation_status <- "ok"
  if (inherits(ablated, "error")) {
    ablation_status <- paste("failed_ablation_construction:", conditionMessage(ablated))
    pred_off <- rep(NA_real_, length(pred_full))
  } else {
    pred_off_result <- tryCatch(predict(ablated, return_se = FALSE),
                                error = function(e) e)
    if (inherits(pred_off_result, "error")) {
      ablation_status <- paste("failed_ablation_prediction:",
                               conditionMessage(pred_off_result))
      pred_off <- rep(NA_real_, length(pred_full))
    } else {
      pred_off <- pred_off_result$imputed$y[internal_for_input]
    }
  }
  if (length(pred_full) != nrow(case$df_observed) ||
      length(pred_off) != nrow(case$df_observed)) {
    method_status <- unattempted
    method_status["pigauto_full"] <- "failed_row_alignment"
    method_status["pigauto_neural_off"] <- if (ablation_status == "ok")
      "failed_row_alignment" else ablation_status
    return(list(status = "failed", failure_reason =
      "prediction row alignment failed", elapsed_seconds = elapsed,
      predictions = NULL, fit_diagnostics = list(
        fit_status = "ok", ablation_status = ablation_status,
        method_status = method_status)))
  }
  train_rows_internal <- training_rows_from_splits(result$data,
    result$splits, latent_col = 1L)
  train_rows_input <- input_row_order[train_rows_internal]
  train_rows_input <- train_rows_input[!is.na(train_rows_input)]
  train_y <- case$df_observed$y[train_rows_input]
  train_ok <- is.finite(train_y)
  train_rows_input <- train_rows_input[train_ok]
  train_y <- train_y[train_ok]
  scored <- which(!is.na(case$stratum))
  score_data <- case$df_complete[scored, , drop = FALSE]
  score_basis <- make_covariate_basis(score_data, case$predictor_names,
                                      nonlinear = TRUE)
  train_data <- case$df_complete[train_rows_input, , drop = FALSE]
  train_basis <- make_covariate_basis(train_data, case$predictor_names,
                                      nonlinear = TRUE)

  fit_lm <- function(x_train, x_new, method_name) {
    if (qr(x_train)$rank < ncol(x_train)) {
      return(list(prediction = rep(NA_real_, nrow(x_new)),
                  status = "failed: rank deficient"))
    }
    fit <- stats::lm.fit(x_train, train_y)
    list(prediction = drop(x_new %*% fit$coefficients), status = "ok")
  }
  linear_train <- train_basis[, seq_len(1L + length(case$predictor_names)),
                              drop = FALSE]
  linear_score <- score_basis[, colnames(linear_train), drop = FALSE]
  lm_linear <- fit_lm(linear_train, linear_score, "linear")
  lm_poly <- fit_lm(train_basis, score_basis, "polynomial_interactions")

  species_index <- match(case$df_complete$species, case$tree$tip.label)
  gls <- tryCatch(fit_phylo_random_intercept(
    y = train_y, X = train_basis,
    species_index = species_index[train_rows_input],
    V = case$phylo_correlation
  ), error = function(e) e)
  if (inherits(gls, "error")) {
    gls_pred <- rep(NA_real_, length(scored))
    gls_status <- paste0("failed: ", conditionMessage(gls))
    gls_converged <- NA
  } else {
    gls_pred <- gls$predict(score_basis, species_index[scored])
    gls_status <- if (gls$converged) "ok" else "nonconverged"
    gls_converged <- gls$converged
    if (!gls$converged) gls_pred[] <- NA_real_
  }

  fit <- result$fit
  hist_epoch <- if (nrow(fit$history)) max(fit$history$epoch) else 0L
  signal_values <- fit$phylo_signal_per_trait
  signal_status <- if (!isTRUE(gate_policy)) {
    "not_requested"
  } else if (!requireNamespace("phytools", quietly = TRUE)) {
    "dependency_missing"
  } else if (any(is.finite(signal_values))) {
    "estimated"
  } else {
    "unestimable"
  }
  paired_finite <- is.finite(pred_full[scored]) & is.finite(pred_off[scored])
  invalid_output <- !all(paired_finite)
  out <- list(
    status = if (invalid_output) "invalid_output" else "completed",
    failure_reason = if (invalid_output) sprintf(
      "non-finite paired predictions on %d/%d outer-scored rows",
      sum(!paired_finite), length(paired_finite)) else NA_character_,
    elapsed_seconds = elapsed,
    predictions = data.frame(
      row_id = case$df_complete$row_id[scored],
      species = case$df_complete$species[scored],
      stratum = case$stratum[scored], truth = case$df_complete$y[scored],
      pigauto_full = pred_full[scored], pigauto_neural_off = pred_off[scored],
      mean = rep(mean(train_y), length(scored)),
      lm_linear = lm_linear$prediction,
      lm_poly_interaction = lm_poly$prediction,
      phylo_gls = gls_pred, stringsAsFactors = FALSE
    ),
    fit_diagnostics = list(
      phylo_signal_per_trait = fit$phylo_signal_per_trait,
      phylo_gate_triggered = fit$phylo_gate_triggered,
      gate_signal_status = signal_status,
      r_cal_bm = fit$r_cal_bm, r_cal_gnn = fit$r_cal_gnn,
      r_cal_mean = fit$r_cal_mean,
      epochs_run = hist_epoch, epoch_limit = epochs,
      training_termination = if (hist_epoch < epochs) "early_stopped"
                             else "epoch_limit",
      gls_status = gls_status, gls_converged = gls_converged,
      linear_status = lm_linear$status, polynomial_status = lm_poly$status,
      fit_status = "ok", ablation_status = ablation_status,
      method_status = c(
        pigauto_full = if (all(is.finite(pred_full[scored]))) "ok" else "invalid_output",
        pigauto_neural_off = if (ablation_status != "ok") ablation_status else
          if (all(is.finite(pred_off[scored]))) "ok" else "invalid_output",
        mean = "ok", lm_linear = lm_linear$status,
        lm_poly_interaction = lm_poly$status, phylo_gls = gls_status
      )
    )
  )
  out$elapsed_seconds <- proc.time()[["elapsed"]] - started
  out
}

run_one_dataset <- function(config, task) {
  verify_config_fingerprint(config)
  if (!identical(source_identity(), config$source)) {
    stop("Current source/environment identity differs from initialized campaign; refusing to run or reuse results.",
         call. = FALSE)
  }
  run_dir <- config$run_dir
  result_path <- file.path(run_dir, "replicates",
    sprintf("%s__rep_%04d.rds", task$scenario_id, task$replicate))
  scenario <- config$scenario_manifest[task$scenario_index, , drop = FALSE]
  dgp_seed <- derive_seed(config$master_seed, task$scenario_index,
                          task$replicate, stream = 1L)
  fit_seed <- derive_seed(config$master_seed, task$scenario_index,
                          task$replicate, stream = 2L)
  expected_id <- list(
    fingerprint = config$fingerprint, task = task,
    dgp_seed = dgp_seed, fit_seed = fit_seed,
    f_type = scenario$f_type, variance_multiplier = scenario$variance_multiplier,
    beta_strength = scenario$beta_strength, epochs = config$epochs
  )
  fingerprint <- configuration_fingerprint(expected_id)
  if (file.exists(result_path)) {
    old <- readRDS(result_path)
    if (!identical(old$fingerprint, fingerprint)) {
      stop("Existing replicate has a stale or different fingerprint: ",
           result_path, call. = FALSE)
    }
    message("Reusing exact fingerprint: ", basename(result_path))
    return(invisible(old))
  }

  settings <- config$fit_settings
  started <- Sys.time()
  case <- tryCatch(make_gnn_case(
    n_species = settings$n_species, obs_per = settings$obs_per,
    n_covs = settings$n_covs,
    variance_multiplier = scenario$variance_multiplier,
    f_type = if (scenario$beta_strength == 0) "linear" else scenario$f_type,
    beta_strength = scenario$beta_strength, seed = dgp_seed,
    sigma_res = settings$sigma_res
  ), error = function(e) e)
  gate_results <- vector("list", length(config$gate_policies))
  names(gate_results) <- paste0("phylo_signal_gate_",
                                config$gate_policies)
  predictions <- NULL
  if (inherits(case, "error")) {
    gate_results[] <- list(list(status = "failed",
      failure_reason = paste("data generation:", conditionMessage(case)),
      elapsed_seconds = 0, predictions = NULL, fit_diagnostics = NULL))
  } else {
    for (j in seq_along(config$gate_policies)) {
      gate_started <- proc.time()[["elapsed"]]
      gate_results[[j]] <- tryCatch(fit_prediction(
        case, config$gate_policies[j], fit_seed, config$epochs, settings
      ), error = function(e) list(status = "failed",
        failure_reason = conditionMessage(e), elapsed_seconds = NA_real_,
        predictions = NULL, fit_diagnostics = NULL))
      gate_results[[j]]$elapsed_seconds <-
        proc.time()[["elapsed"]] - gate_started
      if (is.null(gate_results[[j]]$fit_diagnostics)) {
        gate_results[[j]]$fit_diagnostics <- list(
          fit_status = "failed_unhandled", ablation_status = "unknown",
          method_status = c(pigauto_full = "failed_unhandled",
            pigauto_neural_off = "failed_unhandled", mean = "not_attempted",
            lm_linear = "not_attempted", lm_poly_interaction = "not_attempted",
            phylo_gls = "not_attempted"))
      }
      if (is.null(predictions) && !is.null(gate_results[[j]]$predictions)) {
        predictions <- gate_results[[j]]$predictions[, "row_id", drop = FALSE]
      }
      if (!is.null(gate_results[[j]]$predictions) &&
          !identical(gate_results[[j]]$predictions$row_id,
                     predictions$row_id)) {
        stop("Gate-specific output rows are not aligned.", call. = FALSE)
      }
    }
  }
  record <- list(
    identity = expected_id, fingerprint = fingerprint,
    source = config$source,
    timestamps = list(started = format(started, tz = "UTC", usetz = TRUE),
                      finished = format(Sys.time(), tz = "UTC", usetz = TRUE)),
    elapsed_seconds_total = as.numeric(difftime(Sys.time(), started,
                                                 units = "secs")),
    case = if (inherits(case, "error")) NULL else list(
      meta = case$meta,
      fully_unobserved_species = case$fully_unobserved_species,
      n_scored = sum(!is.na(case$stratum)),
      stratum_counts = table(case$stratum, useNA = "no"),
      tree = case$tree
    ),
    gate_results = gate_results,
    predictions = NULL
  )
  atomic_save_rds(record, result_path)
  message(sprintf("%s rep %d: %s / %s", task$scenario_id, task$replicate,
    gate_results[[1L]]$status, gate_results[[2L]]$status))
  invisible(record)
}

audit_run <- function(run_dir) {
  config <- readRDS(file.path(run_dir, "run_config.rds"))
  tasks <- utils::read.csv(file.path(run_dir, "tasks.csv"),
                           stringsAsFactors = FALSE)
  results <- list.files(file.path(run_dir, "replicates"), pattern = "[.]rds$",
                        full.names = TRUE)
  found <- lapply(results, readRDS)
  expected_names <- sprintf("%s__rep_%04d.rds", tasks$scenario_id,
                            tasks$replicate)
  actual_names <- basename(results)
  if (anyDuplicated(actual_names)) stop("Duplicate result filenames found.")
  unexpected <- setdiff(actual_names, expected_names)
  missing <- setdiff(expected_names, actual_names)
  statuses <- unlist(lapply(found, function(x) {
    vapply(x$gate_results, `[[`, character(1L), "status")
  }), use.names = FALSE)
  counts <- table(factor(statuses,
    levels = c("completed", "invalid_output", "failed", "nonconverged", "infrastructure_skip")))
  data.frame(requested_datasets = nrow(tasks), completed_datasets = length(found),
    missing_datasets = length(missing), unexpected_results = length(unexpected),
    requested_fits = 2L * nrow(tasks), completed_fits = unname(counts["completed"]),
    failed_fits = unname(counts["failed"]),
    invalid_output_fits = unname(counts["invalid_output"]),
    nonconverged_fits = unname(counts["nonconverged"]),
    infrastructure_skips = unname(counts["infrastructure_skip"]),
    stringsAsFactors = FALSE
  )
}

if (sys.nframe() == 0L) {
  cli <- read_cli(commandArgs(trailingOnly = TRUE))
  if (cli$mode == "init") {
    init_campaign(cli$run_dir, cli$n_replicates, cli$master_seed, cli$epochs)
  } else if (cli$mode == "task") {
    config <- readRDS(file.path(cli$run_dir, "run_config.rds"))
    tasks <- utils::read.csv(file.path(cli$run_dir, "tasks.csv"),
                             stringsAsFactors = FALSE)
    verify_config_fingerprint(config)
    verify_task_manifest(config, tasks)
    if (is.na(cli$task_index) || cli$task_index < 1L ||
        cli$task_index > nrow(tasks)) stop("Invalid TASK_INDEX.")
    run_one_dataset(config, tasks[cli$task_index, , drop = FALSE])
  } else {
    init_campaign(cli$run_dir, cli$n_replicates, cli$master_seed, cli$epochs)
    config <- readRDS(file.path(cli$run_dir, "run_config.rds"))
    tasks <- utils::read.csv(file.path(cli$run_dir, "tasks.csv"),
                             stringsAsFactors = FALSE)
    verify_config_fingerprint(config)
    verify_task_manifest(config, tasks)
    for (i in seq_len(nrow(tasks))) {
      run_one_dataset(config, tasks[i, , drop = FALSE])
    }
    print(audit_run(cli$run_dir))
  }
}
