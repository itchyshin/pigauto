#!/usr/bin/env Rscript

source(file.path("script", "validation_20261018", "runners.R"))

stopifnot(is.function(validation_run_job))
stopifnot(is.function(validation_analysis_dgp))

# The registered factorial MAR mechanism requires an always-observed driver.
# Exercise the borrowed campaign constructor and scorer directly; this does
# not fit pigauto or run a campaign arm.
mar_spec <- .vr_prediction_dgp_spec(
  "types_mixed_BM_l0.7_r0.5_mar0.3_n100", "MAR")
stopifnot(identical(mar_spec$dgp, "types_mixed"),
          identical(mar_spec$miss, "mar"), isTRUE(mar_spec$driver))
campaign_lib <- new.env(parent = .GlobalEnv)
sys.source("script/campaign_gnn_off_lib.R", envir = campaign_lib)
mar_cell <- campaign_lib$make_cell(mar_spec$dgp, 80L, 814L,
  miss_frac = 0.30, miss = mar_spec$miss, lambda = 0.7, rho = 0.5,
  evo = "BM", thresholds = "fixed", driver = mar_spec$driver)
stopifnot("d1" %in% names(mar_cell$truth), !any(mar_cell$mask[, "d1"]),
          identical(mar_cell$mechanism, "mar"), sum(mar_cell$mask) > 0L)
campaign_lib$truth <- mar_cell$truth
campaign_lib$mask <- mar_cell$mask
mar_score <- campaign_lib$score_arm("constructor_truth", mar_cell$truth)
mar_numeric <- mar_score[mar_score$metric == "zRMSE", , drop = FALSE]
stopifnot(nrow(mar_numeric) > 0L, all(is.finite(mar_numeric$value)),
          all(mar_numeric$value == 0))

# BIEN cache loading is deterministic, normalizes cache/tree species labels,
# requires at least two measured traits, and selects a bounded 300-tip subset.
bien_tmp <- tempfile("validation-bien-"); dir.create(bien_tmp)
bien_tree <- ape::rtree(305L)
bien_species <- gsub(" ", "_", bien_tree$tip.label)
bien_cache <- list(
  trait_a = data.frame(species = bien_species, mean_value = seq_along(bien_species)),
  trait_b = data.frame(species = bien_species, mean_value = seq_along(bien_species) / 2),
  trait_sparse = data.frame(species = bien_species[1:250], mean_value = 1:250))
traits_file <- file.path(bien_tmp, "traits.rds"); tree_file <- file.path(bien_tmp, "tree.rds")
saveRDS(bien_cache, traits_file); saveRDS(list(scenario.3 = bien_tree), tree_file)
bien_a <- .validation_load_bien_cache(traits_file, tree_file)
bien_b <- .validation_load_bien_cache(traits_file, tree_file)
stopifnot(isTRUE(bien_a$available), nrow(bien_a$data) == 300L,
          length(bien_a$tree$tip.label) == 300L,
          identical(bien_a$data, bien_b$data), identical(bien_a$tree, bien_b$tree),
          all(rownames(bien_a$data) == bien_a$tree$tip.label),
          identical(bien_a$subset_seed, 20261018L))

# Posterior cell-interval summaries retain a per-trait denominator and do not
# discard out-of-range or naturally missing truth rows as if they were scored.
interval_summary <- .validation_posterior_interval_summary(
  data.frame(trait = c("x", "x", "x", "y"), row = c(1L, 2L, 99L, 1L),
    lower = c(0, 0, 0, -1), upper = c(2, 2, 2, 1)),
  list(x = c(1, NA_real_), y = 2),
  expected_mask = list(x = c(TRUE, TRUE, FALSE), y = TRUE))
stopifnot(nrow(interval_summary) == 2L,
          interval_summary$n[interval_summary$trait == "x"] == 1L,
          interval_summary$n_expected[interval_summary$trait == "x"] == 2L,
          interval_summary$n_finite_scored[interval_summary$trait == "x"] == 1L,
          interval_summary$n_missing_or_nonfinite[interval_summary$trait == "x"] == 1L,
          interval_summary$coverage[interval_summary$trait == "x"] == 1,
          interval_summary$mean_width[interval_summary$trait == "x"] == 2,
          identical(interval_summary$status[interval_summary$trait == "x"], "failed"),
          interval_summary$n[interval_summary$trait == "y"] == 1L,
          interval_summary$coverage[interval_summary$trait == "y"] == 0,
          identical(interval_summary$status[interval_summary$trait == "y"], "ok"))

# The analysis generator must be deterministic for a fixed seed and preserve
# the known x coefficient, grouping structure, and requested mask mechanism.
a <- validation_analysis_dgp(model = "lm", n = 40L, seed = 101L,
                             mechanism = "MCAR", miss_frac = 0.30)
b <- validation_analysis_dgp(model = "lm", n = 40L, seed = 101L,
                             mechanism = "MCAR", miss_frac = 0.30)
stopifnot(identical(a$data, b$data), identical(a$missing, b$missing))
stopifnot(sum(a$missing) == ceiling(0.30 * 40), isTRUE(all.equal(a$truth, 0.7)))
stopifnot(all(is.finite(a$data$y)), all(is.finite(a$data$z)))
for (model in c("glm", "lmer")) {
  d <- validation_analysis_dgp(model = model, n = 40L, seed = 202L,
                               mechanism = "MAR", miss_frac = 0.30)
  stopifnot(sum(d$missing) > 0L, sum(d$missing) < 40L,
            all(is.finite(d$complete_data$x)), all(is.finite(d$data$y)),
            isTRUE(all.equal(d$truth, 0.7)))
  if (model == "lmer") stopifnot(nlevels(d$data$group) == 8L)
}

# This public API preflight catches formula/auxiliary overlap before doing any
# imputation work; y and z already occur in the analysis formula.
pkgload::load_all(".", quiet = TRUE, helpers = FALSE)
avonet_input <- .validation_load_bundled_avonet(include_trees = TRUE)
avonet_traits <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")
stopifnot(isTRUE(avonet_input$available), is.data.frame(avonet_input$full),
          all(avonet_traits %in% names(avonet_input$full)),
          inherits(avonet_input$tree, "phylo"), length(avonet_input$trees) >= 2L,
          setequal(rownames(avonet_input$full), avonet_input$tree$tip.label))
avonet_masked <- .validation_make_observed_mask(avonet_input$full, avonet_traits,
  frac = 0.30, seed = 815L)
for (v in avonet_traits) {
  original <- avonet_masked$originally_observed[, v]
  stopifnot(sum(avonet_masked$mask[, v]) == ceiling(sum(original) * 0.30),
            all(avonet_masked$mask[, v] <= original),
            all(is.na(avonet_masked$masked[[v]][avonet_masked$mask[, v]])))
}
preflight <- pigauto:::.analysis_mi_validate(
  data = a$data, formula = y ~ x + z, missing = "x", model = "lm",
  m = 5L, auxiliary = character(), seed = 303L, control = list()
)
stopifnot(identical(preflight$auxiliary, character()), preflight$m == 5L)
analysis_job <- list(job_id = "analysis-lm-n40-m5-test", scenario_id = "analysis_lm_mcar",
  rep = 1L, seed = 304L, slice = "analysis", adapter = "analysis", n = 40L,
  model = "lm", mechanism = "MCAR", frac = 0.30, m = 5L, claim_tier = "test")
analysis_run <- validation_run_job(analysis_job, root = ".", pilot = FALSE)
stopifnot(identical(analysis_run$status, "ok"), nrow(analysis_run$results) == 2L,
  all(analysis_run$results$status == "ok"),
  setequal(analysis_run$results$method, c("complete_data", "analysis_lm")),
  all(is.finite(analysis_run$results$complete_reference)))

# Named mixed multi-observation regimes must alter internal branch lengths,
# retain four mixed traits, and mask each trait independently.
mixed <- .vr_multiobs_mixed_dgp(30L, 811L, lambda = 0.2, rho = 0.6,
                                miss_frac = 0.30)
set.seed(811L); expected_tree <- ape::rtree(30L)
internal <- expected_tree$edge[, 2L] > length(expected_tree$tip.label)
expected_tree$edge.length[internal] <- expected_tree$edge.length[internal] * 0.2
stopifnot(identical(mixed$tree, expected_tree),
          identical(names(mixed$truth), c("species", "x1", "x2", "y", "z")),
          identical(colnames(mixed$mask), c("x1", "x2", "y", "z")),
          all(vapply(mixed$full[, c("x1", "x2", "y", "z")], function(x) anyNA(x), logical(1))),
          any(apply(mixed$mask, 1L, function(x) !all(x) && any(x))))

# Family generators use the approved three-trait ZI and multi-proportion
# designs; these checks are generator-only and do not fit models.
set.seed(812L); family_tree <- ape::rtree(40L)
zi <- pigauto:::simulate_zi_count_traits(family_tree, n_traits = 3L,
  zero_frac = 0.4, mean_nz = 20, seed = 812L)
comp <- pigauto:::simulate_multi_proportion_traits(family_tree, K = 5L,
  signal = 0.6, seed = 812L)
stopifnot(ncol(zi) == 3L, nrow(zi) == 40L,
          ncol(comp) == 5L, nrow(comp) == 40L,
          all(abs(rowSums(comp) - 1) < 1e-8))

# A zRMSE attempt must always retain its coverage denominator, even when the
# scoring helper did not return finite conformal coverage.
prediction_job <- list(job_id = "prediction-coverage-test", scenario_id = "prediction",
  rep = 1L, seed = 813L, slice = "prediction", adapter = "prediction")
prediction_tab <- data.frame(arm = "gnn_off", trait = "trait1", metric = "zRMSE",
  value = 1.2, coverage = NA_real_)
prediction_rows <- .vr_prediction_metric_rows(prediction_job, prediction_tab,
  list(failed = list(gnn_off = FALSE), errors = list()))
coverage_row <- prediction_rows[prediction_rows$estimand == "trait1_conformal_coverage", ]
stopifnot(nrow(coverage_row) == 1L, identical(coverage_row$status, "failed"),
  is.na(coverage_row$estimate), grepl("coverage was absent", coverage_row$failure_reason))

# Ten-bin calibration summaries use masked-cell pairs, retain every bin, and
# require complete valid evidence for the requested discrete mask.
cal_job <- list(job_id = "calibration-test", scenario_id = "prediction_cal",
  rep = 1L, seed = 816L, slice = "prediction", adapter = "prediction")
cal_n <- 200L
cal_truth <- data.frame(class = factor(rep(c("a", "b"), length.out = cal_n)),
                        row.names = paste0("sp", seq_len(cal_n)))
cal_mask <- matrix(TRUE, cal_n, 1L, dimnames = list(rownames(cal_truth), "class"))
centers <- (seq_len(10L) - 0.5) / 10
cal_pairs <- do.call(rbind, lapply(seq_len(10L), function(b) {
  correct_n <- as.integer(centers[[b]] * 20L)
  i <- (b - 1L) * 20L + seq_len(20L)
  data.frame(arm = "gnn_off", trait = "class", species = rownames(cal_truth)[i],
    confidence = centers[[b]], correct = c(rep(1, correct_n), rep(0, 20L - correct_n)))
}))
calibrated <- .vr_calibration_summary(cal_job, cal_truth, cal_mask, cal_pairs,
  arms = "gnn_off")
stopifnot(nrow(calibrated$bins) == 10L, all(calibrated$bins$n == 20L),
          calibrated$results$status == "ok",
          calibrated$results$estimate < 1e-12,
          all(abs(calibrated$bins$mean_confidence - calibrated$bins$mean_correct) < 1e-12))
miscal_pairs <- cal_pairs; miscal_pairs$correct <- 0
miscalibrated <- .vr_calibration_summary(cal_job, cal_truth, cal_mask, miscal_pairs,
  arms = "gnn_off")
stopifnot(abs(miscalibrated$results$estimate - 0.5) < 1e-12)
dropped <- .vr_calibration_summary(cal_job, cal_truth, cal_mask,
  cal_pairs[-1L, ], arms = "gnn_off")
stopifnot(identical(dropped$results$status, "failed"), is.na(dropped$results$estimate),
          grepl("covered 199/200", dropped$results$failure_reason, fixed = TRUE))
unsupported_cal <- .vr_calibration_summary(cal_job, cal_truth, cal_mask, NULL,
  arms = "gnn_off")
stopifnot(identical(unsupported_cal$results$status, "unavailable"),
          grepl("evidence unavailable", unsupported_cal$results$failure_reason, fixed = TRUE))

applicability_tab <- data.frame(
  arm = c("gnn_off", "floor", "floor", "floor", "gnn_off", "gnn_off"),
  trait = c("ord", "c1", "c1", "c1", "c1", "c1"),
  metric = c("brier", "width", "interval_score", "zRMSE", "zRMSE", "width"),
  value = c(NA_real_, NA_real_, NA_real_, 1.1, 1.1, NA_real_),
  coverage = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_))
applicability_rows <- .vr_prediction_metric_rows(prediction_job, applicability_tab,
  list(failed = as.list(setNames(rep(FALSE, 2L), c("gnn_off", "floor"))), errors = list()),
  truth = data.frame(ord = ordered(c("L1", "L2")), c1 = c(1, 2)))
status_of <- function(method, estimand) applicability_rows$status[
  applicability_rows$method == method & applicability_rows$estimand == estimand]
stopifnot(identical(status_of("gnn_off", "ord_brier"), "unavailable"),
  identical(status_of("floor", "c1_width"), "unavailable"),
  identical(status_of("floor", "c1_interval_score"), "unavailable"),
  identical(status_of("floor", "c1_conformal_coverage"), "unavailable"),
  identical(status_of("gnn_off", "c1_conformal_coverage"), "failed"),
  identical(status_of("gnn_off", "c1_width"), "failed"),
  !any(applicability_rows$status[applicability_rows$method == "floor"] == "failed"))

# Unsupported routes are retained as unavailable attempts, never converted
# into successful skips or silently dropped records.
job <- list(job_id = "adapter-test", scenario_id = "unknown", rep = 1L,
            seed = 1L, slice = "test", adapter = "no_such_adapter",
            n = 10L, scenario = "unknown", model = "unknown", reps = 1L,
            claim_tier = "test")
out <- validation_run_job(job, root = ".", pilot = TRUE)
stopifnot(identical(out$job_id, job$job_id), identical(out$status, "unavailable"))
stopifnot(is.data.frame(out$results), nrow(out$results) == 1L,
          identical(out$results$status, "unavailable"),
          nzchar(out$results$failure_reason))

# The real-data route must distinguish missing source data from a successful
# run. AVONET is handled separately by its registered local adapter.
job$adapter <- "real"; job$scenario <- "pantheria"
out <- validation_run_job(job, root = ".", pilot = TRUE)
stopifnot(identical(out$status, "unavailable"),
          grepl("no registered local source data", out$error, fixed = TRUE),
          identical(out$results$status, "unavailable"))

# BIEN is available only when both cache paths are explicitly supplied so the
# external input hashes can be bound into the frozen source identity.
old_bien_traits <- Sys.getenv("VALIDATION_BIEN_TRAITS", unset = "")
old_bien_tree <- Sys.getenv("VALIDATION_BIEN_TREE", unset = "")
Sys.unsetenv(c("VALIDATION_BIEN_TRAITS", "VALIDATION_BIEN_TREE"))
job$scenario <- "bien"
out <- validation_run_job(job, root = ".", pilot = TRUE)
if (nzchar(old_bien_traits)) Sys.setenv(VALIDATION_BIEN_TRAITS = old_bien_traits)
if (nzchar(old_bien_tree)) Sys.setenv(VALIDATION_BIEN_TREE = old_bien_tree)
stopifnot(identical(out$status, "unavailable"),
          grepl("set VALIDATION_BIEN_TRAITS and VALIDATION_BIEN_TREE", out$error, fixed = TRUE),
          identical(out$results$status, "unavailable"))

cat("RUNNERS_PASS\n")
