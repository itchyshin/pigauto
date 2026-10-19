# Adapters for the approved 2026-10-18 validation campaign.
# All routes are dispatched one manifest job at a time. Errors and unavailable
# routes are returned as records so callers can account for every attempt.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x

.vr_get <- function(job, name, default = NULL) {
  if (is.data.frame(job)) {
    if (!(name %in% names(job))) return(default)
    value <- job[[name]][[1L]]
  } else value <- job[[name]]
  absent <- is.null(value) || length(value) == 0L ||
    (is.atomic(value) && all(is.na(value)))
  if (absent) default else value
}

.vr_row <- function(job, method, estimand, status = "ok", truth = NA_real_,
                    estimate = NA_real_, se = NA_real_, lower = NA_real_,
                    upper = NA_real_, complete_reference = NA_real_,
                    failure_reason = NA_character_) {
  data.frame(
    job_id = as.character(.vr_get(job, "job_id", NA_character_)),
    scenario_id = as.character(.vr_get(job, "scenario_id", NA_character_)),
    rep = as.integer(.vr_get(job, "rep", .vr_get(job, "replicate", 1L))),
    seed = as.integer(.vr_get(job, "seed", NA_integer_)),
    slice = as.character(.vr_get(job, "slice", NA_character_)),
    adapter = as.character(.vr_get(job, "adapter", NA_character_)),
    method = as.character(method), estimand = as.character(estimand),
    status = as.character(status), truth = as.numeric(truth),
    estimate = as.numeric(estimate), se = as.numeric(se),
    lower = as.numeric(lower), upper = as.numeric(upper),
    complete_reference = as.numeric(complete_reference),
    failure_reason = as.character(failure_reason), stringsAsFactors = FALSE
  )
}

.vr_envelope <- function(job, status, results = NULL, error = NA_character_,
                         warnings = character(), metrics = list(), settings = list(),
                         diagnostics = list(), wall_sec = NA_real_, pilot = FALSE) {
  if (is.null(results)) results <- data.frame()
  list(
    job_id = as.character(.vr_get(job, "job_id", NA_character_)),
    scenario_id = as.character(.vr_get(job, "scenario_id", NA_character_)),
    rep = as.integer(.vr_get(job, "rep", .vr_get(job, "replicate", 1L))),
    seed = as.integer(.vr_get(job, "seed", NA_integer_)),
    slice = as.character(.vr_get(job, "slice", NA_character_)),
    adapter = as.character(.vr_get(job, "adapter", NA_character_)),
    status = status, error = error, warnings = warnings, metrics = metrics,
    source_identity = NULL, settings = settings, diagnostics = diagnostics,
    wall_sec = wall_sec, pilot = isTRUE(pilot), results = results
  )
}

.vr_load_package <- function(root) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required to load the frozen package root")
  pkgload::load_all(root, quiet = TRUE, helpers = FALSE)
  invisible(root)
}

.vr_prediction_metric_rows <- function(job, tab, run, truth = NULL) {
  if (!is.data.frame(tab) || !nrow(tab)) return(data.frame())
  rows <- list()
  for (i in seq_len(nrow(tab))) {
    arm <- as.character(tab$arm[[i]])
    metric <- as.character(tab$metric[[i]])
    metric_value <- as.numeric(tab$value[[i]])
    arm_failed <- isTRUE(run$failed[[arm]])
    arm_error <- run$errors[[arm]] %||% NA_character_
    trait <- if ("trait" %in% names(tab)) as.character(tab$trait[[i]]) else "trait"
    unsupported_ordinal_brier <- identical(metric, "brier") && !is.null(truth) &&
      trait %in% names(truth) && is.ordered(truth[[trait]])
    point_only_interval <- identical(arm, "floor") &&
      metric %in% c("width", "interval_score")
    metric_unavailable <- unsupported_ordinal_brier || point_only_interval
    score_status <- if (metric_unavailable) "unavailable" else
      if (arm_failed || !is.finite(metric_value)) "failed" else "ok"
    metric_reason <- if (unsupported_ordinal_brier) "Brier score is unsupported for ordinal traits" else
      if (point_only_interval) "point-only control no interval" else
      if (score_status == "ok") NA_character_ else "prediction metric was non-finite"
    rows[[length(rows) + 1L]] <- .vr_row(job, arm, paste(trait, metric, sep = "_"),
      status = score_status, estimate = if (metric_unavailable) NA_real_ else metric_value,
      failure_reason = if (!is.na(arm_error) && !metric_unavailable) arm_error else metric_reason)
    if (identical(metric, "zRMSE")) {
      coverage <- if ("coverage" %in% names(tab)) as.numeric(tab$coverage[[i]]) else NA_real_
      coverage_unavailable <- identical(arm, "floor")
      coverage_status <- if (coverage_unavailable) "unavailable" else
        if (arm_failed || !is.finite(coverage)) "failed" else "ok"
      coverage_reason <- if (coverage_unavailable) "point-only control no interval" else
        if (!is.na(arm_error)) arm_error else
        if (!is.finite(coverage)) "conformal coverage was absent or non-finite" else NA_character_
      rows[[length(rows) + 1L]] <- .vr_row(job, arm,
        paste(trait, "conformal_coverage", sep = "_"), status = coverage_status,
        estimate = if (coverage_unavailable) NA_real_ else coverage,
        failure_reason = coverage_reason)
    }
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

.vr_calibration_summary <- function(job, truth, mask, calibration, arms,
                                    arm_failed = list(), arm_errors = list(), n_bins = 10L) {
  bin_cols <- c("method", "trait", "bin", "bin_lower", "bin_upper", "n",
                "mean_confidence", "mean_correct", "absolute_gap", "status")
  empty_bins <- data.frame(method = character(), trait = character(), bin = integer(),
    bin_lower = numeric(), bin_upper = numeric(), n = integer(),
    mean_confidence = numeric(), mean_correct = numeric(), absolute_gap = numeric(),
    status = character(), stringsAsFactors = FALSE)
  traits <- names(truth)[vapply(truth, function(x) is.factor(x) || is.character(x), logical(1))]
  if (!length(traits)) return(list(results = data.frame(), bins = empty_bins))
  if (!is.matrix(mask) || is.null(rownames(mask)) || is.null(colnames(mask)))
    stop("calibration mask must be a named matrix")
  if (!is.null(calibration) && (!is.data.frame(calibration) ||
      !all(c("arm", "trait", "species", "confidence", "correct") %in% names(calibration))))
    stop("calibration evidence has an invalid schema")
  edges <- seq(0, 1, length.out = as.integer(n_bins) + 1L)
  bin_rows <- list(); result_rows <- list()
  for (arm in unique(as.character(arms))) for (trait in traits) {
    if (!(trait %in% colnames(mask))) stop("mask is missing discrete trait ", trait)
    expected_species <- rownames(mask)[which(mask[, trait])]
    raw <- if (is.data.frame(calibration)) calibration[
      calibration$arm == arm & calibration$trait == trait, , drop = FALSE] else data.frame()
    failed <- isTRUE(arm_failed[[arm]])
    err <- arm_errors[[arm]] %||% NA_character_
    if (!length(expected_species)) {
      status <- "unavailable"; reason <- "no masked cells requested for calibration"
    } else if (!nrow(raw)) {
      status <- if (failed) "failed" else "unavailable"
      reason <- if (failed && !is.na(err)) err else "probability calibration evidence unavailable"
    } else {
      confidence <- suppressWarnings(as.numeric(raw$confidence))
      correct <- suppressWarnings(as.numeric(raw$correct))
      species <- as.character(raw$species)
      in_scope <- species %in% expected_species
      valid <- in_scope & is.finite(confidence) & confidence >= 0 & confidence <= 1 &
        is.finite(correct) & correct %in% c(0, 1)
      complete_ids <- setequal(species[valid], expected_species) &&
        !anyDuplicated(species[valid]) && nrow(raw) == length(expected_species)
      if (!complete_ids || any(!valid) || failed) {
        status <- "failed"
        reason <- if (failed && !is.na(err)) err else sprintf(
          "calibration evidence covered %d/%d requested cells with valid unique confidence/correct pairs",
          length(unique(species[valid])), length(expected_species))
      } else {
        status <- "ok"; reason <- NA_character_
      }
    }
    if (nrow(raw)) {
      confidence <- suppressWarnings(as.numeric(raw$confidence))
      correct <- suppressWarnings(as.numeric(raw$correct))
      valid <- as.character(raw$species) %in% expected_species & is.finite(confidence) &
        confidence >= 0 & confidence <= 1 & is.finite(correct) & correct %in% c(0, 1)
    } else {
      confidence <- correct <- numeric(); valid <- logical()
    }
    bin_id <- findInterval(confidence[valid], edges, rightmost.closed = TRUE,
                           all.inside = TRUE)
    method_bins <- vector("list", as.integer(n_bins))
    for (b in seq_len(as.integer(n_bins))) {
      take <- bin_id == b
      n <- sum(take)
      mean_conf <- if (n) mean(confidence[valid][take]) else NA_real_
      mean_correct <- if (n) mean(correct[valid][take]) else NA_real_
      method_bins[[b]] <- data.frame(method = arm, trait = trait, bin = b,
        bin_lower = edges[[b]], bin_upper = edges[[b + 1L]], n = n,
        mean_confidence = mean_conf, mean_correct = mean_correct,
        absolute_gap = if (n) abs(mean_conf - mean_correct) else NA_real_,
        status = status, stringsAsFactors = FALSE)
    }
    bins <- do.call(rbind, method_bins)
    bin_rows[[length(bin_rows) + 1L]] <- bins
    ece <- if (status == "ok") sum(bins$n / length(expected_species) *
      ifelse(bins$n > 0L, bins$absolute_gap, 0)) else NA_real_
    result_rows[[length(result_rows) + 1L]] <- .vr_row(job, arm,
      paste(trait, "ece", sep = "_"), status = status, estimate = ece,
      failure_reason = reason)
  }
  list(results = do.call(rbind, result_rows),
       bins = if (length(bin_rows)) do.call(rbind, bin_rows) else empty_bins)
}

.vr_prediction_dgp_spec <- function(scenario, mechanism) {
  scenario <- as.character(scenario)
  miss <- tolower(as.character(mechanism))
  dgp <- if (grepl("^avonet", scenario)) "avonet" else
    if (grepl("^ou_mixed", scenario)) "ou_mixed" else
    if (grepl("^bm_mixed", scenario)) "bm_mixed" else "types_mixed"
  list(dgp = dgp, miss = miss, driver = identical(miss, "mar"))
}

.vr_prediction <- function(job, root, pilot) {
  n <- as.integer(.vr_get(job, "n", 100L)); seed <- as.integer(.vr_get(job, "seed", 1L))
  scenario <- as.character(.vr_get(job, "scenario", "types_mixed_BM_l1_r0_mcar0.3_n100"))
  miss_frac <- as.numeric(.vr_get(job, "frac", 0.30))
  lambda <- as.numeric(.vr_get(job, "lambda", 1))
  rho <- as.numeric(.vr_get(job, "rho", 0))
  evo <- as.character(.vr_get(job, "evo", "BM"))
  spec <- .vr_prediction_dgp_spec(scenario, .vr_get(job, "mechanism", "mcar"))
  dgp <- spec$dgp
  libenv <- new.env(parent = .GlobalEnv)
  sys.source(file.path(root, "script", "campaign_gnn_off_lib.R"), envir = libenv)
  cell <- libenv$make_cell(dgp, n, seed, miss_frac = miss_frac, miss = spec$miss,
                           lambda = lambda, rho = rho, evo = evo,
                           thresholds = "fixed", driver = spec$driver)
  # Reuse the approved baseline/floor configurations. GNN contrasts live in
  # the covariate workflow; the optional joint-BM comparator is explicit.
  arms <- c("gnn_off", "floor")
  if (isTRUE(.vr_get(job, "rphylopars", FALSE))) arms <- c(arms, "gnn_off_rphylopars")
  libenv$truth <- cell$truth
  libenv$mask <- cell$mask
  opts <- list(seed = seed, epochs = as.integer(.vr_get(job, "epochs", 2000L)),
               log_line = function(...) invisible(NULL))
  run <- libenv$run_arms(cell, arms, opts)
  tab <- run$results
  rows <- .vr_prediction_metric_rows(job, tab, run, truth = cell$truth)
  calibration <- .vr_calibration_summary(job, cell$truth, cell$mask,
    run$calib, arms, arm_failed = run$failed, arm_errors = run$errors)
  if (nrow(calibration$results)) rows <- rbind(rows, calibration$results)
  if (!nrow(rows)) rows <- .vr_row(job, "prediction", "masked_cell_score",
    status = "failed", failure_reason = "adapter returned no scored attempts")
  results <- rows; rownames(results) <- NULL
  failed <- any(results$status == "failed")
  list(results = results, status = if (failed) "failed" else "ok",
       error = if (length(run$errors)) paste(unlist(run$errors), collapse = " | ") else NA_character_,
       warnings = character(), metrics = list(raw_scores = tab,
         calibration_bins = calibration$bins),
       settings = list(n = n, seed = seed, dgp = dgp, missing_frac = miss_frac,
                       lambda = lambda, rho = rho, evo = evo, miss = spec$miss,
                       driver = spec$driver, arms = arms),
       diagnostics = list(walls = run$walls, errors = run$errors, paths = run$paths))
}

.vr_family <- function(job, root, composition = FALSE) {
  n <- as.integer(.vr_get(job, "n", 100L)); seed <- as.integer(.vr_get(job, "seed", 1L))
  scenario <- as.character(.vr_get(job, "scenario", "signal_0.6"))
  frac <- as.numeric(.vr_get(job, "frac", 0.25))
  set.seed(seed)
  tree <- ape::rtree(n)
  if (composition) {
    signal <- if (grepl("signal_", scenario)) as.numeric(sub("signal_", "", scenario)) else as.numeric(.vr_get(job, "rho", 0.6))
    K <- if (grepl("^K_", scenario)) as.integer(sub("K_", "", scenario)) else as.integer(.vr_get(job, "K", 5L))
    full <- pigauto:::simulate_multi_proportion_traits(tree, K = K, signal = signal, seed = seed)
    group <- setNames(list(names(full)), "comp")
  } else {
    zero_frac <- if (grepl("zf_", scenario)) as.numeric(sub("zf_", "", scenario)) else as.numeric(.vr_get(job, "rho", 0.5))
    mean_nz <- if (grepl("mean_nz_", scenario)) as.numeric(sub("mean_nz_", "", scenario)) else as.numeric(.vr_get(job, "evo", 20))
    full <- pigauto:::simulate_zi_count_traits(tree, n_traits = 3L, zero_frac = zero_frac,
                                               mean_nz = mean_nz, seed = seed)
    group <- NULL
  }
  trait_types <- if (!composition) setNames(rep("zi_count", ncol(full)), names(full)) else NULL
  pd <- pigauto::preprocess_traits(full, tree, trait_types = trait_types,
                                   multi_proportion_groups = group)
  splits <- pigauto::make_missing_splits(pd$X_scaled, missing_frac = frac,
    seed = seed, trait_map = pd$trait_map)
  graph <- pigauto::build_phylo_graph(tree, k_eigen = "auto")
  baseline <- pigauto::fit_baseline(pd, tree, splits = splits, graph = graph)
  mean_pred <- pd$X_scaled
  for (j in seq_len(ncol(mean_pred))) {
    train <- mean_pred[splits$mask[, j], j]
    center <- mean(train[is.finite(train)], na.rm = TRUE)
    if (!is.finite(center)) center <- 0
    mean_pred[!splits$mask[, j], j] <- center
  }
  fit <- pigauto::fit_pigauto(pd, tree, splits = splits, baseline = baseline,
    graph = graph, epochs = as.integer(.vr_get(job, "epochs", if (composition) 300L else 500L)),
    verbose = FALSE, seed = seed, gnn = FALSE)
  prediction <- stats::predict(fit, return_se = !composition, n_imputations = 1L)
  methods <- list(mean = mean_pred, baseline = baseline$mu, pigauto_gnn_off = prediction)
  rows <- list(); evals <- list()
  for (method in names(methods)) {
    ev <- pigauto::evaluate_imputation(methods[[method]], pd$X_scaled, splits,
                                        trait_map = pd$trait_map)
    if (is.null(ev) || !nrow(ev)) {
      rows[[length(rows) + 1L]] <- .vr_row(job, method, "test_metrics", "failed",
        failure_reason = "evaluation returned no held-out rows")
      next
    }
    ev <- ev[ev$split == "test", , drop = FALSE]
    evals[[method]] <- ev
    metric_names <- if (composition) c("aitchison", "rmse_clr", "simplex_mae") else
      c("brier", "zero_accuracy", "rmse")
    for (metric in metric_names) {
      value <- ev[[metric]]
      for (i in seq_len(nrow(ev))) {
        rows[[length(rows) + 1L]] <- .vr_row(job, method,
          paste(ev$trait[[i]], metric, sep = "_"),
          status = if (is.finite(value[[i]])) "ok" else "failed",
          estimate = value[[i]],
          failure_reason = if (is.finite(value[[i]])) NA_character_ else
            paste0("held-out ", metric, " was not finite"))
      }
    }
    if (!composition && method == "pigauto_gnn_off") {
      latent <- if (inherits(prediction, "pigauto_pred")) prediction$imputed_latent else prediction
      n_latent <- nrow(pd$X_scaled); idx <- splits$test_idx
      ri <- ((idx - 1L) %% n_latent) + 1L; cj <- ceiling(idx / n_latent)
      for (tm in pd$trait_map) {
        if (!identical(tm$type, "zi_count")) next
        ix <- ri[cj == tm$latent_cols[[1L]]]
        ix <- ix[is.finite(pd$X_scaled[ix, tm$latent_cols[[2L]]])]
        mag_estimand <- paste(tm$name, "positive_magnitude_latent_rmse", sep = "_")
        if (!length(ix)) {
          rows[[length(rows) + 1L]] <- .vr_row(job, method, mag_estimand, "failed",
            failure_reason = "no positive-count magnitudes in held-out cells")
        } else {
          err <- latent[ix, tm$latent_cols[[2L]]] - pd$X_scaled[ix, tm$latent_cols[[2L]]]
          finite_err <- is.finite(err)
          rows[[length(rows) + 1L]] <- .vr_row(job, method, mag_estimand,
            status = if (all(finite_err)) "ok" else "failed",
            estimate = if (all(finite_err)) sqrt(mean(err^2)) else NA_real_,
            failure_reason = if (all(finite_err)) NA_character_ else
              sprintf("only %d/%d positive magnitudes were finite", sum(finite_err), length(err)))
        }
      }
    }
  }
  results <- if (length(rows)) do.call(rbind, rows) else .vr_row(job, "family", "test_metrics", "failed",
    failure_reason = "no family metrics were produced")
  if (composition) metrics <- list(n_masked = length(splits$test_idx),
    test_metrics = evals, simulator = "pigauto::simulate_multi_proportion_traits") else
    metrics <- list(n_masked = length(splits$test_idx), test_metrics = evals,
      zero_frac = zero_frac, mean_nonzero = mean_nz,
      positive_magnitude_scale = "log1p z-scored latent scale",
      simulator = "pigauto::simulate_zi_count_traits", n_traits = 3L)
  list(results = results, status = if (any(results$status != "ok")) "failed" else "ok",
    error = if (any(results$status != "ok")) "one or more family metrics unavailable" else NA_character_,
    warnings = character(), metrics = metrics,
    settings = list(n = n, seed = seed, frac = frac, gnn = FALSE,
      composition = composition, simulator = if (composition) "multi_proportion" else "zi_count"),
    diagnostics = list(test_idx = splits$test_idx))
}

.vr_masked_score <- function(truth, prediction, continuous = TRUE) {
  if (continuous) {
    truth <- as.numeric(truth); prediction <- as.numeric(prediction)
    valid <- is.finite(truth) & is.finite(prediction)
    estimate <- if (all(valid) && length(valid)) sqrt(mean((truth - prediction)^2)) else NA_real_
  } else {
    valid <- !is.na(truth) & !is.na(prediction)
    estimate <- if (all(valid) && length(valid)) mean(as.character(truth) == as.character(prediction)) else NA_real_
  }
  list(status = if (is.finite(estimate)) "ok" else "failed", estimate = estimate,
       n_masked = length(valid), n_scored = sum(valid))
}

.vr_multiobs <- function(job, root, pilot) {
  n <- as.integer(.vr_get(job, "n", 200L)); seed <- as.integer(.vr_get(job, "seed", 1L))
  scenario <- as.character(.vr_get(job, "scenario", "ctmax_l0.5_b0.5_f0.5"))
  beta <- if (grepl("_b[^_]+", scenario)) as.numeric(sub(".*_b([^_]+)_f.*", "\\1", scenario)) else as.numeric(.vr_get(job, "rho", 0.5))
  lambda <- as.numeric(.vr_get(job, "lambda", 0.5))
  sp_frac <- as.numeric(.vr_get(job, "frac", 0.5)); within_frac <- 0.20
  set.seed(seed); tree <- ape::rtree(n); sp <- tree$tip.label
  phylo <- as.numeric(ape::rTraitCont(tree, model = "BM", sigma = sqrt(lambda), root.value = 0)[sp])
  names(phylo) <- sp
  nobs <- pmin(pmax(stats::rpois(n, 5L), 1L), 20L); names(nobs) <- sp
  species <- rep(sp, nobs); acclim <- stats::rnorm(length(species), 20, 5)
  y <- 30 + phylo[species] + beta * acclim + stats::rnorm(length(species), 0, 0.5)
  full <- data.frame(species = species, CTmax = y, acclim_temp = acclim)
  truth <- full$CTmax; set.seed(seed + 999L)
  missing_species <- sample(sp, floor(sp_frac * n))
  mask <- full$species %in% missing_species
  obs <- which(!mask); within <- sample(obs, floor(within_frac * length(obs)))
  mask[within] <- TRUE; full$CTmax[mask] <- NA_real_
  covs <- data.frame(acclim_temp = full$acclim_temp)
  covs$acclim_temp[!is.finite(covs$acclim_temp)] <- 20
  methods <- list(no_covariates = NULL, covariates = covs)
  rows <- list(); errors <- list()
  for (method in names(methods)) {
    fit <- tryCatch(pigauto::impute(full[, c("species", "CTmax")], tree,
      species_col = "species", covariates = methods[[method]], gnn = TRUE,
      epochs = as.integer(.vr_get(job, "epochs", 200L)), verbose = FALSE, seed = seed),
      error = identity)
    if (inherits(fit, "error")) {
      errors[[method]] <- conditionMessage(fit)
      rows[[length(rows) + 1L]] <- .vr_row(job, method, "masked_obs_rmse", "failed",
                                           failure_reason = errors[[method]])
    } else {
      pred <- fit$completed$CTmax
      score <- .vr_masked_score(truth[mask], pred[mask], continuous = TRUE)
      rows[[length(rows) + 1L]] <- .vr_row(job, method, "masked_obs_rmse", score$status,
        estimate = score$estimate,
        failure_reason = if (score$status == "ok") NA_character_ else
          sprintf("only %d/%d masked predictions were finite", score$n_scored, score$n_masked))
    }
  }
  results <- do.call(rbind, rows)
  results$status[!is.finite(results$estimate)] <- "failed"
  list(results = results, status = if (length(errors) || any(results$status != "ok")) "failed" else "ok",
       error = if (length(errors)) paste(names(errors), unlist(errors), collapse = " | ") else NA_character_,
       warnings = character(), metrics = list(n_obs = length(truth), masked = sum(mask)),
       settings = list(n_species = n, seed = seed, lambda = lambda, beta = beta,
                       species_missing_frac = sp_frac, within_obs_missing_frac = within_frac,
                       epochs = .vr_get(job, "epochs", 200L), generator = "bench_multi_obs.R::simulate_ctmax_data"),
       diagnostics = list(errors = errors))
}

.vr_multiobs_mixed_dgp <- function(n, seed, lambda, rho, miss_frac, within_frac = 0.20) {
  set.seed(seed); tree <- ape::rtree(n)
  if (abs(lambda - 1) >= 1e-6) {
    is_internal <- tree$edge[, 2L] > length(tree$tip.label)
    tree$edge.length[is_internal] <- tree$edge.length[is_internal] * lambda
  }
  sp <- tree$tip.label; L <- chol(ape::vcv(tree)); set.seed(seed)
  phy <- function() as.numeric(t(L) %*% stats::rnorm(n))
  x1s <- phy(); x2s <- rho * x1s + sqrt(1 - rho^2) * phy()
  ys <- rho * x1s + sqrt(1 - rho^2) * phy()
  zl <- matrix(0, nrow = n, ncol = 4L)
  for (k in seq_len(4L)) zl[, k] <- rho * x1s + sqrt(1 - rho^2) * phy()
  sp_df <- data.frame(x1 = x1s, x2 = x2s,
    y = factor(ifelse(ys > 0, "B", "A"), levels = c("A", "B")),
    z = factor(LETTERS[max.col(zl)], levels = LETTERS[1:4]), row.names = sp)
  set.seed(seed + 10000L); nobs <- pmin(pmax(stats::rpois(n, 5L), 1L), 15L)
  species <- rep(sp, nobs); no <- length(species)
  x1 <- rep(sp_df$x1, nobs) + stats::rnorm(no, 0, .3)
  x2 <- rep(sp_df$x2, nobs) + stats::rnorm(no, 0, .3)
  yy <- rep(as.character(sp_df$y), nobs); flip <- stats::runif(no) < .1
  yy[flip] <- ifelse(yy[flip] == "A", "B", "A")
  zz <- rep(as.character(sp_df$z), nobs); flip <- stats::runif(no) < .1
  if (any(flip)) zz[flip] <- vapply(zz[flip], function(z) sample(setdiff(LETTERS[1:4], z), 1L), character(1))
  truth <- data.frame(species = species, x1 = x1, x2 = x2,
    y = factor(yy, levels = c("A", "B")), z = factor(zz, levels = LETTERS[1:4]))
  full <- truth; set.seed(seed + 999L); miss_sp <- sample(sp, floor(miss_frac * n))
  mask <- matrix(FALSE, nrow(truth), 4L, dimnames = list(NULL, c("x1", "x2", "y", "z")))
  mask[truth$species %in% miss_sp, ] <- TRUE; obs <- which(!(truth$species %in% miss_sp))
  for (j in seq_len(ncol(mask))) mask[sample(obs, round(within_frac * length(obs))), j] <- TRUE
  for (j in seq_len(ncol(mask))) full[[colnames(mask)[j]]][mask[, j]] <- NA
  list(tree = tree, truth = truth, full = full, mask = mask,
       n_species = n, lambda = lambda, rho = rho, missing_frac = miss_frac,
       within_frac = within_frac)
}

.vr_multiobs_mixed <- function(job, root, pilot) {
  n <- as.integer(.vr_get(job, "n", 150L)); seed <- as.integer(.vr_get(job, "seed", 1L))
  scenario <- as.character(.vr_get(job, "scenario", "mod_phylo_correlated"))
  lambda <- as.numeric(.vr_get(job, "lambda", if (grepl("high", scenario)) 1 else if (grepl("low", scenario)) .2 else .5))
  rho <- as.numeric(.vr_get(job, "rho", if (grepl("high", scenario)) .7 else .6))
  miss_frac <- as.numeric(.vr_get(job, "frac", .3)); within_frac <- .20
  dgp <- .vr_multiobs_mixed_dgp(n, seed, lambda, rho, miss_frac, within_frac)
  tree <- dgp$tree; truth <- dgp$truth; full <- dgp$full; mask <- dgp$mask
  rows <- list(); errors <- list()
  for (aggregation in c("hard", "soft")) {
    fit <- tryCatch(pigauto::impute(full, tree, species_col = "species",
      multi_obs_aggregation = aggregation, gnn = FALSE, verbose = FALSE, seed = seed), error = identity)
    method <- paste0("pigauto_", aggregation)
    if (inherits(fit, "error")) {
      errors[[method]] <- conditionMessage(fit)
      rows[[length(rows) + 1L]] <- .vr_row(job, method, "x1_rmse", "failed", failure_reason = errors[[method]])
    } else for (v in c("x1", "x2", "y", "z")) {
      p <- fit$completed[[v]][mask[, v]]; t <- truth[[v]][mask[, v]]
      is_cont <- v %in% c("x1", "x2")
      valid <- if (is_cont) is.finite(as.numeric(p)) & is.finite(as.numeric(t)) else !is.na(p) & !is.na(t)
      val <- if (all(valid)) {
        if (is_cont) sqrt(mean((as.numeric(t) - as.numeric(p))^2)) else mean(as.character(t) == as.character(p))
      } else NA_real_
      rows[[length(rows) + 1L]] <- .vr_row(job, method,
        paste0(v, if (is_cont) "_rmse" else "_accuracy"),
        status = if (is.finite(val)) "ok" else "failed", estimate = val,
        failure_reason = if (is.finite(val)) NA_character_ else
          sprintf("only %d/%d masked predictions were finite", sum(valid), length(valid)))
    }
  }
  results <- do.call(rbind, rows)
  list(results = results, status = if (length(errors) || any(results$status != "ok")) "failed" else "ok",
    error = if (length(errors)) paste(unlist(errors), collapse = " | ") else NA_character_,
    warnings = character(), metrics = list(masked_observations = sum(mask)),
    settings = list(n_species = n, seed = seed, lambda = lambda, rho = rho,
      missing_frac = miss_frac, within_species_mask = within_frac, generator = "bench_multi_obs_mixed.R pure DGP"),
    diagnostics = list(errors = errors, masked_by_trait = colSums(mask)))
}

.vr_covariate <- function(job, root, pilot) {
  n <- as.integer(.vr_get(job, "n", 300L)); seed <- as.integer(.vr_get(job, "seed", 1L))
  scenario <- as.character(.vr_get(job, "scenario", "mod_phylo_strong_env"))
  frac <- as.numeric(.vr_get(job, "frac", .4)); lambda <- as.numeric(.vr_get(job, "lambda", .3))
  beta <- if (startsWith(scenario, "high_phylo_no_env")) 0 else if (startsWith(scenario, "high_phylo_mod_env")) .5 else if (startsWith(scenario, "mod_phylo_strong_env")) 1 else 1.5
  if (n == 300L) { e <- new.env(parent = baseenv()); utils::data("tree300", package = "pigauto", envir = e); tree <- e$tree300 }
  else tree <- ape::rtree(n)
  n <- length(tree$tip.label); sp <- tree$tip.label
  set.seed(seed); R <- stats::cov2cor(ape::vcv(tree)); bm <- t(chol(R)) %*% matrix(stats::rnorm(n * 3L), n, 3L)
  env <- matrix(stats::rnorm(n * 3L), n, 3L, dimnames = list(sp, paste0("env", 1:3)))
  rng_after_env <- .Random.seed
  set.seed(beta * 100 + 7); B <- matrix(stats::rnorm(9L), 3L, 3L)
  assign(".Random.seed", rng_after_env, envir = .GlobalEnv)
  eps <- matrix(stats::rnorm(n * 3L), n, 3L)
  linear_env <- beta * (env %*% t(B))
  nonlinear_env <- beta * (env %*% t(B)) + beta * 0.5 * (env[, 1L]^2 - 1)
  errors <- list(); rows <- list()
  variants <- list(linear = linear_env, nonlinear = nonlinear_env)
  for (variant in names(variants)) {
    y <- variants[[variant]] + sqrt(lambda) * bm + sqrt(1 - lambda) * eps
    colnames(y) <- paste0("trait", 1:3); rownames(y) <- sp
    full <- as.data.frame(y); covs <- as.data.frame(env)
    if (variant == "nonlinear") covs$env1_sq <- env[, 1L]^2 - 1
    set.seed(seed + 999L); masked <- full; mask <- matrix(FALSE, n, 3L, dimnames = dimnames(full))
    for (v in names(full)) { idx <- sample.int(n, round(n * frac)); mask[idx, v] <- TRUE; masked[idx, v] <- NA_real_ }
    methods <- list(pigauto_default = list(gnn = FALSE),
                    gnn_on = list(gnn = TRUE),
                    gnn_on_covariates = list(gnn = TRUE, covariates = covs))
    for (method in names(methods)) {
      opts <- methods[[method]]
      fit <- tryCatch(do.call(pigauto::impute, c(list(traits = masked, tree = tree,
        epochs = as.integer(.vr_get(job, "epochs", 300L)), verbose = FALSE, seed = seed), opts)), error = identity)
      estimand_prefix <- paste(variant, method, sep = "__")
      if (inherits(fit, "error")) {
        errors[[estimand_prefix]] <- conditionMessage(fit)
        rows[[length(rows) + 1L]] <- .vr_row(job, estimand_prefix, "trait1_rmse", "failed", failure_reason = errors[[estimand_prefix]])
      } else for (v in names(full)) {
        score <- .vr_masked_score(full[[v]][mask[, v]], fit$completed[[v]][mask[, v]], continuous = TRUE)
        rows[[length(rows) + 1L]] <- .vr_row(job, estimand_prefix, paste0(v, "_rmse"),
          status = score$status, estimate = score$estimate,
          failure_reason = if (score$status == "ok") NA_character_ else
            sprintf("only %d/%d masked predictions were finite", score$n_scored, score$n_masked))
      }
    }
  }
  results <- do.call(rbind, rows)
  list(results = results, status = if (length(errors) || any(results$status != "ok")) "failed" else "ok",
    error = if (length(errors)) paste(unlist(errors), collapse = " | ") else NA_character_,
    warnings = character(), metrics = list(masked_cells = sum(mask)),
    settings = list(n = n, seed = seed, scenario = scenario, lambda = lambda,
      beta = beta, frac = frac, generator = "bench_covariate_sim.R DGP plus env1-squared contrast",
      variants = c("linear", "nonlinear"), gnn_arms = c("gnn_on", "gnn_on_covariates")), diagnostics = list(errors = errors))
}

.vr_unavailable <- function(job, reason) {
  list(results = .vr_row(job, "unavailable", "attempt", "unavailable", failure_reason = reason),
       status = "unavailable", error = reason, warnings = character(), metrics = list(),
       settings = list(), diagnostics = list())
}

.validation_normalize_species <- function(x) {
  gsub("_", " ", gsub("[[:space:]]+", " ", trimws(as.character(x))))
}

.validation_posterior_interval_summary <- function(cell_interval, truth, expected_mask) {
  empty <- data.frame(trait = character(), n = integer(), n_expected = integer(),
    n_finite_scored = integer(), n_missing_or_nonfinite = integer(),
    coverage = numeric(), mean_width = numeric(), status = character(),
    failure_reason = character(), stringsAsFactors = FALSE)
  if (!is.list(expected_mask) || is.null(names(expected_mask)) ||
      any(!nzchar(names(expected_mask))))
    stop("expected_mask must be a named list of per-trait logical masks")
  if (!is.data.frame(cell_interval)) cell_interval <- data.frame()
  if (nrow(cell_interval)) {
  required <- c("trait", "row", "lower", "upper")
  if (!all(required %in% names(cell_interval)))
    stop("posterior cell_interval is missing required columns: ",
         paste(setdiff(required, names(cell_interval)), collapse = ", "))
  }
  interval_traits <- if (nrow(cell_interval)) unique(as.character(cell_interval$trait)) else character()
  traits <- unique(c(names(expected_mask), interval_traits))
  rows <- lapply(traits, function(trait) {
    ci <- cell_interval[cell_interval$trait == trait, , drop = FALSE]
    if (is.null(truth[[trait]])) stop("posterior interval trait absent from truth: ", trait)
    if (is.null(expected_mask[[trait]])) stop("posterior interval trait absent from expected mask: ", trait)
    expected_rows <- which(as.logical(expected_mask[[trait]]))
    expected_n <- length(expected_rows)
    idx <- as.integer(ci$row)
    requested <- !is.na(idx) & idx %in% expected_rows
    in_bounds <- requested & idx >= 1L & idx <= length(truth[[trait]])
    truth_i <- rep(NA_real_, nrow(ci))
    truth_i[in_bounds] <- as.numeric(truth[[trait]][idx[in_bounds]])
    valid <- in_bounds & is.finite(truth_i) & is.finite(ci$lower) & is.finite(ci$upper)
    n_finite <- length(unique(idx[valid]))
    n_missing <- max(0L, expected_n - n_finite)
    duplicate_rows <- anyDuplicated(idx[requested]) > 0L
    lost <- expected_n > 0L && (n_finite != expected_n || duplicate_rows)
    metric_status <- if (expected_n == 0L) "unavailable" else if (lost) "failed" else "ok"
    reason <- if (expected_n == 0L) "no masked cells requested for this trait" else
      if (lost) sprintf("only %d/%d requested masked cells had finite intervals%s",
        n_finite, expected_n, if (duplicate_rows) " (duplicate interval rows)" else "") else NA_character_
    data.frame(trait = trait, n = n_finite, n_expected = expected_n,
      n_finite_scored = n_finite, n_missing_or_nonfinite = n_missing,
      coverage = if (any(valid)) mean(truth_i[valid] >= ci$lower[valid] &
                                       truth_i[valid] <= ci$upper[valid]) else NA_real_,
      mean_width = if (any(valid)) mean(ci$upper[valid] - ci$lower[valid]) else NA_real_,
      status = metric_status, failure_reason = reason,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.validation_load_bien_cache <- function(traits_path, tree_path,
                                        subset_n = 300L,
                                        subset_seed = 20261018L) {
  if (!nzchar(traits_path) || !nzchar(tree_path) ||
      !file.exists(traits_path) || !file.exists(tree_path))
    return(list(available = FALSE,
      reason = "BIEN cache unavailable; set VALIDATION_BIEN_TRAITS and VALIDATION_BIEN_TREE to existing files"))
  trait_means <- tryCatch(readRDS(traits_path), error = identity)
  tree_pkg <- tryCatch(readRDS(tree_path), error = identity)
  if (inherits(trait_means, "error") || inherits(tree_pkg, "error"))
    return(list(available = FALSE, reason = "BIEN cache could not be read"))
  if (!is.list(trait_means) || is.null(names(trait_means)) ||
      !is.list(tree_pkg) || is.null(tree_pkg$scenario.3) ||
      !inherits(tree_pkg$scenario.3, "phylo"))
    return(list(available = FALSE, reason = "BIEN cache schema does not match bench_bien.R"))
  species_lists <- lapply(trait_means, function(d) {
    if (is.null(d) || !is.data.frame(d) || !("species" %in% names(d))) character(0)
    else .validation_normalize_species(d$species)
  })
  all_species <- Reduce(union, species_lists, init = character(0))
  wide <- data.frame(species = all_species, stringsAsFactors = FALSE)
  for (nm in names(trait_means)) {
    d <- trait_means[[nm]]
    if (is.null(d) || !is.data.frame(d)) { wide[[nm]] <- NA_real_; next }
    value_col <- intersect(c("mean_value", "trait_value", "value"), names(d))[1L]
    if (is.na(value_col) || !("species" %in% names(d))) { wide[[nm]] <- NA_real_; next }
    key <- .validation_normalize_species(d$species)
    wide[[nm]] <- suppressWarnings(as.numeric(d[[value_col]][match(wide$species, key)]))
  }
  finite_traits <- if (ncol(wide) > 1L)
    rowSums(is.finite(as.matrix(wide[, names(trait_means), drop = FALSE]))) else integer(nrow(wide))
  wide <- wide[finite_traits >= 2L, , drop = FALSE]
  tree <- tree_pkg$scenario.3
  tree$tip.label <- .validation_normalize_species(tree$tip.label)
  matched <- intersect(tree$tip.label, wide$species)
  if (length(matched) < subset_n)
    return(list(available = FALSE, reason = sprintf(
      "BIEN cache has only %d matched species with at least two traits; %d required",
      length(matched), subset_n), n_matched = length(matched)))
  set.seed(as.integer(subset_seed))
  selected <- sample(matched, as.integer(subset_n), replace = FALSE)
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, selected))
  wide <- wide[match(tree$tip.label, wide$species), , drop = FALSE]
  rownames(wide) <- wide$species
  wide$species <- NULL
  list(available = TRUE, data = wide, tree = tree, n_matched = length(matched),
       subset_n = nrow(wide), subset_seed = as.integer(subset_seed),
       trait_names = names(wide), traits_path = normalizePath(traits_path),
       tree_path = normalizePath(tree_path))
}

.validation_load_bundled_avonet <- function(include_trees = FALSE) {
  objects <- c("avonet300", "tree300", if (include_trees) "trees300")
  e <- new.env(parent = baseenv())
  utils::data(list = objects, package = "pigauto", envir = e)
  if (!exists("avonet300", envir = e, inherits = FALSE) ||
      !exists("tree300", envir = e, inherits = FALSE))
    return(list(available = FALSE, reason = "bundled AVONET300 data could not be loaded"))
  full <- e$avonet300
  if (!is.data.frame(full) || !("Species_Key" %in% names(full)) ||
      !inherits(e$tree300, "phylo"))
    return(list(available = FALSE, reason = "bundled AVONET300 data have an unexpected schema"))
  rownames(full) <- as.character(full$Species_Key)
  full$Species_Key <- NULL
  if (include_trees && (!exists("trees300", envir = e, inherits = FALSE) ||
      !is.list(e$trees300) || length(e$trees300) < 2L ||
      !all(vapply(e$trees300, inherits, logical(1), "phylo"))))
    return(list(available = FALSE, reason = "bundled AVONET tree posterior has an unexpected schema"))
  list(available = TRUE, full = full, tree = e$tree300,
       trees = if (include_trees) e$trees300 else NULL)
}

.validation_make_observed_mask <- function(full, traits, frac, seed) {
  if (!is.data.frame(full) || !all(traits %in% names(full)))
    stop("full data must contain all registered traits")
  observed <- as.matrix(full[, traits, drop = FALSE])
  storage.mode(observed) <- "numeric"
  mask <- matrix(FALSE, nrow(observed), ncol(observed),
                 dimnames = list(rownames(full), traits))
  set.seed(as.integer(seed))
  for (v in traits) {
    idx <- which(is.finite(observed[, v]))
    if (length(idx)) mask[sample(idx, ceiling(length(idx) * frac)), v] <- TRUE
  }
  masked <- full
  for (v in traits) masked[[v]][mask[, v]] <- NA_real_
  list(masked = masked, mask = mask, originally_observed = is.finite(observed))
}

.vr_real <- function(job, root, pilot) {
  dataset <- tolower(as.character(.vr_get(job, "scenario", "")))
  bien <- NULL
  if (dataset == "bien") {
    traits_path <- Sys.getenv("VALIDATION_BIEN_TRAITS", unset = "")
    tree_path <- Sys.getenv("VALIDATION_BIEN_TREE", unset = "")
    bien <- .validation_load_bien_cache(traits_path, tree_path, subset_n = 300L,
                                         subset_seed = 20261018L)
    if (!isTRUE(bien$available)) return(.vr_unavailable(job, bien$reason))
    full <- bien$data; tree <- bien$tree; traits <- names(full)
  } else if (dataset == "avonet") {
    avonet <- .validation_load_bundled_avonet()
    if (!isTRUE(avonet$available)) return(.vr_unavailable(job, avonet$reason))
    full <- avonet$full; tree <- avonet$tree
    traits <- intersect(c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length"), names(full))
  } else return(.vr_unavailable(job, paste0("no registered local source data for ", dataset)))
  if (length(traits) < 2L) return(.vr_unavailable(job, paste0(dataset, " data lack the registered traits")))
  seed <- as.integer(.vr_get(job, "seed", 1L)); frac <- as.numeric(.vr_get(job, "frac", .3))
  mask_out <- .validation_make_observed_mask(full, traits, frac, seed)
  masked <- mask_out$masked; mask <- mask_out$mask
  fit <- tryCatch(pigauto::impute(masked[, traits, drop = FALSE], tree,
    gnn = FALSE, verbose = FALSE, seed = seed), error = identity)
  dataset_label <- if (dataset == "bien") "BIEN300subset" else "AVONET300"
  if (inherits(fit, "error")) return(list(results = .vr_row(job, "pigauto_gnn_off", "Mass_rmse", "failed",
      failure_reason = conditionMessage(fit)), status = "failed", error = conditionMessage(fit),
      warnings = character(), metrics = list(), settings = list(dataset = dataset_label), diagnostics = list()))
  rows <- lapply(traits, function(v) {
    score <- .vr_masked_score(full[[v]][mask[, v]], fit$completed[[v]][mask[, v]], continuous = TRUE)
    val <- score$estimate
    .vr_row(job, "pigauto_gnn_off", paste0(v, "_rmse"),
      status = score$status, estimate = val,
      failure_reason = if (score$status == "ok") NA_character_ else
        sprintf("only %d/%d masked predictions were finite", score$n_scored, score$n_masked))
  })
  results <- do.call(rbind, rows)
  list(results = results, status = if (any(results$status != "ok")) "failed" else "ok",
    error = NA_character_, warnings = character(), metrics = list(masked_cells = sum(mask)),
    settings = list(dataset = if (dataset == "bien") "BIEN300subset" else "AVONET300",
      n = nrow(full), seed = seed, frac = frac, traits = traits, gnn = FALSE,
      subset_seed = if (dataset == "bien") bien$subset_seed else NA_integer_,
      n_matched = if (dataset == "bien") bien$n_matched else NA_integer_,
      cache_traits = if (dataset == "bien") bien$traits_path else NA_character_,
      cache_tree = if (dataset == "bien") bien$tree_path else NA_character_),
    diagnostics = list(masked_by_trait = colSums(mask),
      originally_observed_by_trait = colSums(mask_out$originally_observed)))
}

validation_analysis_dgp <- function(model = c("lm", "glm", "lmer"), n = 100L,
                                    seed = 1L, mechanism = "MCAR", miss_frac = 0.30) {
  model <- match.arg(model); n <- as.integer(n); seed <- as.integer(seed)
  set.seed(seed)
  x <- stats::rnorm(n); z <- stats::rnorm(n)
  group <- factor(rep(seq_len(ceiling(n / 5)), each = 5L)[seq_len(n)])
  eta <- 0.5 + 0.7 * x + 0.3 * z
  if (model == "glm") {
    y <- stats::rbinom(n, 1L, stats::plogis(eta))
  } else if (model == "lmer") {
    re <- stats::rnorm(nlevels(group), 0, 0.5)
    y <- eta + re[as.integer(group)] + stats::rnorm(n, 0, 1)
  } else y <- eta + stats::rnorm(n, 0, 1)
  dat <- data.frame(y = y, x = x, z = z, group = group)
  set.seed(seed + 4099L)
  if (identical(toupper(mechanism), "MCAR")) {
    hide <- sample(seq_len(n), ceiling(miss_frac * n))
  } else if (identical(toupper(mechanism), "MAR")) {
    driver <- as.numeric(scale(y)) + as.numeric(scale(z)) / 2
    f <- function(a) mean(stats::plogis(a + driver)) - miss_frac
    hide <- which(stats::runif(n) < stats::plogis(stats::uniroot(f, c(-30, 30))$root + driver))
    if (!length(hide)) hide <- which.max(driver)
  } else stop("analysis mechanism must be MCAR or MAR")
  full_data <- dat
  dat$x[hide] <- NA_real_
  list(data = dat, complete_data = full_data,
       missing = replace(rep(FALSE, n), hide, TRUE), truth = 0.7,
       model = model, mechanism = toupper(mechanism), seed = seed)
}

.vr_analysis <- function(job, pilot) {
  model <- tolower(as.character(.vr_get(job, "model", "lm")))
  if (model %in% c("analysislm", "analysis_lm")) model <- "lm"
  if (model %in% c("analysisglm", "analysis_glm")) model <- "glm"
  if (model %in% c("analysislmer", "analysis_lmer")) model <- "lmer"
  if (!(model %in% c("lm", "glm", "lmer"))) stop("unsupported analysis model: ", model, call. = FALSE)
  required <- switch(model, lm = character(), glm = "smcfcs", lmer = c("jomo", "lme4"))
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) return(.vr_unavailable(job,
    paste0("analysis model ", model, " requires unavailable package(s): ", paste(missing, collapse = ", "))))
  n <- as.integer(.vr_get(job, "n", 100L)); seed <- as.integer(.vr_get(job, "seed", 1L))
  mechanism <- as.character(.vr_get(job, "mechanism", "MCAR"))
  d <- validation_analysis_dgp(model, n, seed, mechanism, as.numeric(.vr_get(job, "frac", 0.30)))
  ref_fit <- if (model == "lm") stats::lm(y ~ x + z, data = d$complete_data) else if (model == "glm")
    stats::glm(y ~ x + z, data = d$complete_data, family = stats::binomial()) else
      lme4::lmer(y ~ x + z + (1 | group), data = d$complete_data, REML = FALSE)
  complete_reference <- if (model == "lmer") unname(lme4::fixef(ref_fit)[["x"]]) else unname(stats::coef(ref_fit)[["x"]])
  ref_tab <- summary(ref_fit)$coefficients
  ref_se <- unname(ref_tab["x", "Std. Error"])
  ref_df <- if (model == "lm") stats::df.residual(ref_fit) else NA_real_
  ref_crit <- if (model == "lm") stats::qt(.975, df = ref_df) else stats::qnorm(.975)
  m <- as.integer(.vr_get(job, "m", 50L))
  formula <- if (model == "lmer") y ~ x + z + (1 | group) else y ~ x + z
  mi <- pigauto::multi_impute_analysis(d$data, formula = formula, missing = "x",
       model = model, m = m, auxiliary = character(), seed = seed,
       control = .vr_get(job, "control", list()))
  fits <- pigauto::with_imputations(mi, function(data) {
    if (model == "lm") stats::lm(y ~ x + z, data = data)
    else if (model == "glm") stats::glm(y ~ x + z, data = data, family = stats::binomial())
    else lme4::lmer(y ~ x + z + (1 | group), data = data, REML = FALSE)
  })
  n_attempted <- length(fits)
  n_failed <- as.integer(attr(fits, "n_failed") %||%
    sum(vapply(fits, inherits, logical(1), "pigauto_mi_error")))
  n_ok <- n_attempted - n_failed
  ref_result <- .vr_row(job, "complete_data", "beta_x", truth = 0.7,
    estimate = complete_reference, se = ref_se,
    lower = complete_reference - ref_crit * ref_se,
    upper = complete_reference + ref_crit * ref_se,
    complete_reference = complete_reference)
  if (n_ok < 2L) {
    reason <- sprintf("only %d/%d downstream fits succeeded; at least two are required to pool", n_ok, n_attempted)
    failed_row <- .vr_row(job, paste0("analysis_", model), "beta_x", status = "failed",
      truth = 0.7, complete_reference = complete_reference, failure_reason = reason)
    return(list(results = rbind(ref_result, failed_row), status = "failed", error = reason,
      warnings = character(), metrics = list(n_attempted = n_attempted, n_ok = n_ok,
        n_failed = n_failed, m = m),
      settings = list(n = n, m = m, seed = seed, model = model, mechanism = mechanism,
        missing_frac = .vr_get(job, "frac", 0.30), pilot = isTRUE(pilot), smoke = FALSE),
      diagnostics = list(fit_failures = fits[vapply(fits, inherits, logical(1), "pigauto_mi_error")])))
  }
  pooled <- pigauto::pool_mi(fits)
  row <- pooled[pooled$term == "x", , drop = FALSE]
  if (nrow(row) != 1L) stop("pooled results did not contain exactly one x coefficient")
  est <- row$estimate[[1L]]; se <- row$std.error[[1L]]; df <- row$df[[1L]]
  crit <- stats::qt(0.975, df = df)
  mi_status <- if (n_failed == 0L && n_ok == m) "ok" else "failed"
  mi_failure <- if (mi_status == "ok") NA_character_ else
    sprintf("downstream fits succeeded for %d/%d imputations", n_ok, n_attempted)
  result <- .vr_row(job, paste0("analysis_", model), "beta_x", status = mi_status,
                    truth = 0.7,
                    estimate = est, se = se, lower = est - crit * se,
                    upper = est + crit * se, complete_reference = complete_reference,
                    failure_reason = mi_failure)
  list(results = rbind(ref_result, result),
       status = if (mi_status == "ok") "ok" else "failed",
       error = if (mi_status == "ok") NA_character_ else mi_failure,
       warnings = character(),
       metrics = list(n_attempted = n_attempted, n_ok = n_ok, n_failed = n_failed, m = m),
       settings = list(n = n, m = m, seed = seed,
       model = model, mechanism = mechanism, missing_frac = .vr_get(job, "frac", 0.30),
       pilot = isTRUE(pilot), smoke = isTRUE(pilot) && m < as.integer(.vr_get(job, "m", 50L)),
       complete_reference_df = ref_df,
       complete_reference_interval = if (model == "lm") "t" else "normal approximation"),
       diagnostics = list(pooled = row, mi_engine = mi$engine,
                          dgp = "y=0.5+0.7*x+0.3*z+noise; x,z iid N(0,1)"))
}

.vr_posterior <- function(job, root, pilot) {
  regime_id <- as.integer(.vr_get(job, "regime_id", .vr_get(job, "scenario", 25L)))
  rep_i <- as.integer(.vr_get(job, "rep", 1L))
  dgp_env <- new.env(parent = .GlobalEnv)
  sys.source(file.path(root, "script", "mi_gls", "regimes.R"), envir = dgp_env)
  dgp_expr <- parse(file.path(root, "script", "mi_gls", "dgp_v2.R"), keep.source = FALSE)
  # dgp_v2.R's first expression sources regimes.R into .GlobalEnv. It is
  # already loaded above in the private environment, so skip that expression.
  for (expr in as.list(dgp_expr)[-1L]) eval(expr, envir = dgp_env)
  simulate_regime_cell <- dgp_env$simulate_regime_cell
  cell <- simulate_regime_cell(regime_id, rep_i)
  ctl <- list(n_chains = as.integer(.vr_get(job, "n_chains", 4L)),
              keep_draws = as.integer(.vr_get(job, "keep_draws", 1000L)),
              param_uncertainty = if (cell$regime$missing[[1L]] == "both") "both" else "full",
              seed = cell$seed)
  for (nm in c("n_iter", "burnin", "thin")) {
    val <- .vr_get(job, nm, NULL); if (!is.null(val)) ctl[[nm]] <- as.integer(val)
  }
  m <- as.integer(.vr_get(job, "m", 50L))
  mi <- pigauto::multi_impute(cell$df, cell$tree, m = m, draws_method = "posterior",
       posterior_control = ctl, verbose = FALSE, seed = cell$seed)
  interval_metrics <- .validation_posterior_interval_summary(
    mi$posterior$cell_interval, cell$truth,
    expected_mask = list(x = cell$mask_x, y = cell$mask_y))
  full <- mi$datasets
  model <- tolower(as.character(.vr_get(job, "model", "gls")))
  fit_one <- function(dat) {
    dat$species <- rownames(dat)
    if (model == "gls") nlme::gls(y ~ x, correlation = ape::corBrownian(1, cell$tree, form = ~species), data = dat, method = "ML")
    else if (model == "lm") stats::lm(y ~ x, data = dat)
    else stop("posterior adapter supports model=gls or lm, got ", model)
  }
  ref <- fit_one(cell$truth)
  ref_est <- if (model == "gls") unname(summary(ref)$tTable["x", "Value"]) else unname(stats::coef(ref)[["x"]])
  ref_se <- if (model == "gls") unname(summary(ref)$tTable["x", "Std.Error"]) else unname(summary(ref)$coefficients["x", "Std. Error"])
  ref_df <- if (model == "gls") nrow(cell$truth) - length(stats::coef(ref)) else stats::df.residual(ref)
  ref_crit <- stats::qt(.975, df = ref_df)
  fits <- lapply(full, function(dat) tryCatch(fit_one(dat), error = identity))
  good <- !vapply(fits, inherits, logical(1), "error")
  pooled <- if (sum(good) >= 2L) pigauto::pool_mi(fits[good]) else stop("fewer than two downstream fits succeeded")
  row <- pooled[pooled$term == "x", , drop = FALSE]
  est <- row$estimate[[1L]]; se <- row$std.error[[1L]]; df <- row$df[[1L]]
  crit <- stats::qt(0.975, df)
  converged <- isTRUE(mi$posterior$converged)
  full_fit_ok <- all(good) && length(fits) == m
  full_status <- if (converged && full_fit_ok) "ok" else "failed"
  full_failure <- paste(c(
    if (!converged) "posterior sampler did not converge" else NULL,
    if (!full_fit_ok) sprintf("downstream fits succeeded for %d/%d imputations", sum(good), m) else NULL
  ), collapse = "; ")
  rows <- list(.vr_row(job, "complete_data", "beta_x", truth = cell$true_beta,
    estimate = ref_est, se = ref_se,
    lower = ref_est - ref_crit * ref_se,
    upper = ref_est + ref_crit * ref_se,
    complete_reference = ref_est),
    .vr_row(job, "posterior_full", "beta_x", status = full_status,
    truth = cell$true_beta,
    estimate = est, se = se, lower = est - crit * se, upper = est + crit * se,
    complete_reference = ref_est,
    failure_reason = if (full_status == "ok") NA_character_ else full_failure))
  if (!is.null(mi$posterior_improper) && length(mi$posterior_improper$datasets)) {
    fits_none <- lapply(mi$posterior_improper$datasets, function(dat) tryCatch(fit_one(dat), error = identity))
    ok_none <- !vapply(fits_none, inherits, logical(1), "error")
    none_fit_ok <- all(ok_none) && length(fits_none) == m
    none_status <- if (converged && none_fit_ok) "ok" else "failed"
    none_failure <- paste(c(
      if (!converged) "posterior sampler did not converge" else NULL,
      if (!none_fit_ok) sprintf("downstream fits succeeded for %d/%d imputations", sum(ok_none), m) else NULL
    ), collapse = "; ")
    if (sum(ok_none) >= 2L) {
      pnone <- pigauto::pool_mi(fits_none[ok_none]); rn <- pnone[pnone$term == "x", , drop = FALSE]
      sn <- rn$std.error[[1L]]; en <- rn$estimate[[1L]]; dn <- rn$df[[1L]]; cn <- stats::qt(0.975, dn)
      rows[[length(rows) + 1L]] <- .vr_row(job, "posterior_none", "beta_x", status = none_status,
        truth = cell$true_beta, estimate = en, se = sn, lower = en - cn * sn,
        upper = en + cn * sn, complete_reference = ref_est,
        failure_reason = if (none_status == "ok") NA_character_ else none_failure)
    } else rows[[length(rows) + 1L]] <- .vr_row(job, "posterior_none", "beta_x", "failed",
        truth = cell$true_beta, complete_reference = ref_est, failure_reason = "fewer than two downstream fits succeeded")
  }
  results <- do.call(rbind, rows)
  list(results = results, status = if (any(results$status != "ok")) "failed" else "ok",
       error = if (any(results$status != "ok")) paste(unique(results$failure_reason[results$status != "ok"]), collapse = "; ") else NA_character_,
       warnings = character(), metrics = list(n_attempted = length(fits), n_ok = sum(good), m = m,
         posterior_converged = converged, posterior_cell_interval = interval_metrics),
       settings = list(regime_id = regime_id, rep = rep_i, sampler_control = mi$posterior$control,
                       m = m, downstream = model, smoke = FALSE,
                       dgp_seed = cell$seed, sampler_seed = cell$seed,
                       interval_summary_scope = "replicate-level posterior cell intervals; descriptive",
                       complete_reference_df = ref_df,
                       complete_reference_interval = "t"),
       diagnostics = list(posterior = mi$posterior$diagnostics, complete_reference = ref_est))
}

.vr_trees <- function(job, root, pilot) {
  seed <- as.integer(.vr_get(job, "seed", 1L)); frac <- as.numeric(.vr_get(job, "frac", .3))
  avonet <- .validation_load_bundled_avonet(include_trees = TRUE)
  if (!isTRUE(avonet$available)) return(.vr_unavailable(job, avonet$reason))
  full <- avonet$full; tree <- avonet$tree; trees <- avonet$trees
  cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")
  missing <- full; mask <- matrix(FALSE, nrow(full), length(cols), dimnames = list(rownames(full), cols))
  set.seed(seed)
  for (v in cols) { idx <- sample.int(nrow(full), round(nrow(full) * frac)); mask[idx, v] <- TRUE; missing[idx, v] <- NA_real_ }
  T <- min(10L, length(trees)); m <- 5L
  mi <- pigauto::multi_impute_trees(missing, trees[seq_len(T)], m_per_tree = m,
    epochs = as.integer(.vr_get(job, "epochs", 300L)), missing_frac = .25,
    verbose = FALSE, seed = seed, eval_every = 50L, patience = 20L)
  rows <- list(); scored_by_trait <- setNames(integer(length(cols)), cols)
  for (v in cols) {
    vals <- unlist(lapply(mi$datasets, function(d) d[[v]][mask[, v]]))
    truth <- rep(full[[v]][mask[, v]], times = length(mi$datasets))
    valid <- is.finite(truth) & is.finite(vals)
    scored_by_trait[[v]] <- sum(valid)
    val <- if (all(valid) && length(valid)) sqrt(mean((truth - vals)^2)) else NA_real_
    rows[[length(rows) + 1L]] <- .vr_row(job, "multi_tree", paste0(v, "_rmse"),
      status = if (is.finite(val)) "ok" else "failed", estimate = val,
      failure_reason = if (is.finite(val)) NA_character_ else
        sprintf("only %d/%d masked predictions were finite", scored_by_trait[[v]], length(valid)))
  }
  results <- do.call(rbind, rows)
  list(results = results, status = if (any(results$status != "ok")) "failed" else "ok",
       error = NA_character_, warnings = character(), metrics = list(n_trees = T, m_per_tree = m),
       settings = list(n = nrow(full), seed = seed, missing_frac = frac, n_trees = T,
         m_per_tree = m, epochs = .vr_get(job, "epochs", 300L),
         data = "avonet300", trees = "trees300"), diagnostics = list(tree_index = mi$tree_index,
           masked_by_trait = colSums(mask), scored_by_trait_across_imputations = scored_by_trait))
}

validation_run_job <- function(job, root = ".", pilot = FALSE) {
  started <- proc.time()[["elapsed"]]
  adapter <- tolower(as.character(.vr_get(job, "adapter", "")))
  result <- tryCatch({
    .vr_load_package(root)
    if (adapter %in% c("prediction", "gnn_off", "prediction_reuse")) .vr_prediction(job, root, pilot)
    else if (adapter %in% c("familyzi", "family_zi", "zi_count", "zi")) .vr_family(job, root, FALSE)
    else if (adapter %in% c("composition", "multiproportion", "multi_proportion")) .vr_family(job, root, TRUE)
    else if (adapter %in% c("workflowmultiobs", "multiobs", "workflow_multiobs")) .vr_multiobs(job, root, pilot)
    else if (adapter == "multiobs_mixed") .vr_multiobs_mixed(job, root, pilot)
    else if (adapter == "covariate") .vr_covariate(job, root, pilot)
    else if (adapter %in% c("tree", "trees", "tree_uncertainty")) .vr_trees(job, root, pilot)
    else if (adapter %in% c("posterior", "mi_gls", "posterior_mi")) .vr_posterior(job, root, pilot)
    else if (adapter == "real") .vr_real(job, root, pilot)
    else if (adapter %in% c("analysislm", "analysis_lm", "analysisglm", "analysis_glm",
                            "analysislmer", "analysis_lmer", "analysis")) .vr_analysis(job, pilot)
    else structure(list(results = .vr_row(job, adapter %||% "unknown", "attempt", "unavailable",
         failure_reason = paste0("no adapter registered for '", adapter, "'")),
         status = "unavailable", error = paste0("no adapter registered for '", adapter, "'"),
         warnings = character(), metrics = list(), settings = list(), diagnostics = list()), class = "validation_adapter_result")
  }, error = function(e) structure(list(
      results = .vr_row(job, adapter, "attempt", "failed", failure_reason = conditionMessage(e)),
      status = "failed", error = conditionMessage(e), warnings = character(), metrics = list(),
      settings = list(), diagnostics = list()), class = "validation_adapter_result"))
  .vr_envelope(job, result$status, results = result$results, error = result$error,
    warnings = result$warnings %||% character(), metrics = result$metrics %||% list(),
    settings = result$settings %||% list(), diagnostics = result$diagnostics %||% list(),
    wall_sec = proc.time()[["elapsed"]] - started, pilot = pilot)
}
