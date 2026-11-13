#!/usr/bin/env Rscript

# Flatten replicate-level posterior predictive interval metrics without
# treating the individual masked cells as independent simulation replicates.
validation_posterior_interval_tables <- function(records, jobs) {
  by_rep_cols <- c("scenario_id", "job_id", "rep", "seed", "status", "interval_status",
                   "trait", "n_expected", "n_finite_scored", "n_missing_or_nonfinite",
                   "n_cells", "coverage", "mean_width",
                   "scope")
  empty <- as.data.frame(setNames(replicate(length(by_rep_cols), logical(0),
                                           simplify = FALSE), by_rep_cols),
                        stringsAsFactors = FALSE)
  # Give the empty tables stable, useful column types.
  empty$scenario_id <- empty$job_id <- empty$status <- empty$interval_status <- empty$trait <- empty$scope <- character()
  empty$rep <- empty$seed <- empty$n_cells <- integer()
  empty$n_expected <- empty$n_finite_scored <- empty$n_missing_or_nonfinite <- integer()
  empty$coverage <- empty$mean_width <- numeric()
  summary_empty <- data.frame(scenario_id = character(), trait = character(),
    nrep_planned = integer(), nrep = integer(), nrep_missing = integer(),
    mean_coverage = numeric(), mcse_coverage = numeric(),
    mean_width = numeric(), mcse_mean_width = numeric(),
    mean_cells_per_replicate = numeric(),
    scope = character(), stringsAsFactors = FALSE)
  if (!is.data.frame(jobs) || !all(c("job_id", "scenario_id", "rep", "seed", "adapter") %in% names(jobs)))
    stop("jobs need job_id, scenario_id, rep, seed, and adapter")
  posterior_jobs <- jobs[tolower(as.character(jobs$adapter)) == "posterior", , drop = FALSE]
  if (!nrow(posterior_jobs))
    return(list(by_replicate = empty, summary = summary_empty))
  if (anyDuplicated(paste(posterior_jobs$scenario_id, posterior_jobs$rep, sep = "\034")))
    stop("posterior jobs contain duplicate scenario/replicate identities")
  if (!is.list(records)) stop("records must be a list of job envelopes")
  record_ids <- vapply(records, function(x) {
    if (is.list(x) && length(x$job_id) == 1L) as.character(x$job_id) else NA_character_
  }, character(1))
  rows <- list()
  for (i in seq_len(nrow(posterior_jobs))) {
    job <- posterior_jobs[i, , drop = FALSE]
    hit <- which(record_ids == as.character(job$job_id))
    if (length(hit) > 1L) stop("duplicate posterior envelopes for job ", job$job_id)
    env <- if (length(hit)) records[[hit]] else NULL
    metric <- if (is.list(env) && is.list(env$metrics)) env$metrics$posterior_cell_interval else NULL
    env_status <- if (is.list(env) && length(env$status) == 1L) as.character(env$status) else "missing"
    sampler_converged <- is.list(env) && is.list(env$metrics) &&
      isTRUE(env$metrics$posterior_converged)
    if (is.data.frame(metric) && nrow(metric)) {
      required <- c("trait", "n", "n_expected", "n_finite_scored",
                    "n_missing_or_nonfinite", "coverage", "mean_width", "status")
      if (!all(required %in% names(metric))) {
        legacy <- c("trait", "n", "coverage", "mean_width")
        if (!all(legacy %in% names(metric)))
          stop("posterior interval metrics missing columns for job ", job$job_id)
        # Older pilot metrics did not retain the requested-cell denominator.
        # Preserve their values but never promote them as valid interval evidence.
        for (nm in c("n_expected", "n_finite_scored", "n_missing_or_nonfinite"))
          metric[[nm]] <- NA_integer_
        metric$status <- "unavailable"
      }
      if (anyDuplicated(metric$trait)) stop("duplicate posterior interval trait for job ", job$job_id)
      for (j in seq_len(nrow(metric))) {
        rows[[length(rows) + 1L]] <- data.frame(
          scenario_id = as.character(job$scenario_id), job_id = as.character(job$job_id),
          rep = as.integer(job$rep), seed = as.integer(job$seed), status = env_status,
          interval_status = if (identical(as.character(metric$status[[j]]), "unavailable")) "unavailable" else
            if (sampler_converged && identical(as.character(metric$status[[j]]), "ok") &&
                is.finite(metric$coverage[[j]]) && is.finite(metric$mean_width[[j]])) "ok" else "failed",
          trait = as.character(metric$trait[[j]]),
          n_expected = as.integer(metric$n_expected[[j]]),
          n_finite_scored = as.integer(metric$n_finite_scored[[j]]),
          n_missing_or_nonfinite = as.integer(metric$n_missing_or_nonfinite[[j]]),
          n_cells = as.integer(metric$n[[j]]),
          coverage = as.numeric(metric$coverage[[j]]), mean_width = as.numeric(metric$mean_width[[j]]),
          scope = "descriptive; replicate-level interval metric", stringsAsFactors = FALSE)
      }
    } else {
      # Keep a visible row for attempts that failed or returned no interval
      # summary, so missing intervals cannot vanish from the replicate ledger.
      rows[[length(rows) + 1L]] <- data.frame(
        scenario_id = as.character(job$scenario_id), job_id = as.character(job$job_id),
        rep = as.integer(job$rep), seed = as.integer(job$seed), status = env_status,
        interval_status = if (is.null(env)) "missing" else "unavailable",
        trait = NA_character_, n_expected = NA_integer_, n_finite_scored = NA_integer_,
        n_missing_or_nonfinite = NA_integer_, n_cells = NA_integer_, coverage = NA_real_,
        mean_width = NA_real_, scope = "descriptive; replicate-level interval metric",
        stringsAsFactors = FALSE)
    }
  }
  by_rep <- if (length(rows)) do.call(rbind, rows) else empty
  rownames(by_rep) <- NULL
  traits_by_scenario <- split(by_rep$trait[!is.na(by_rep$trait)],
                              by_rep$scenario_id[!is.na(by_rep$trait)])
  summary_rows <- list()
  for (scenario in names(traits_by_scenario)) {
    for (trait in unique(traits_by_scenario[[scenario]])) {
      planned <- posterior_jobs[posterior_jobs$scenario_id == scenario, , drop = FALSE]
      sub <- by_rep[by_rep$scenario_id == scenario & by_rep$trait == trait, , drop = FALSE]
      valid <- sub$interval_status == "ok" & is.finite(sub$coverage) & is.finite(sub$mean_width)
      coverage <- sub$coverage[valid]
      width <- sub$mean_width[valid]
      cells <- sub$n_cells[valid]
      nrep <- length(coverage)
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        scenario_id = scenario, trait = trait, nrep_planned = nrow(planned),
        nrep = nrep, nrep_missing = nrow(planned) - nrep,
        mean_coverage = if (nrep) mean(coverage) else NA_real_,
        mcse_coverage = if (nrep >= 2L) stats::sd(coverage) / sqrt(nrep) else NA_real_,
        mean_width = if (nrep) mean(width) else NA_real_,
        mcse_mean_width = if (nrep >= 2L) stats::sd(width) / sqrt(nrep) else NA_real_,
        mean_cells_per_replicate = if (nrep) mean(cells) else NA_real_,
        scope = "descriptive; equal-weighted independent replicates; MCSE across replicates",
        stringsAsFactors = FALSE)
    }
  }
  summary <- if (length(summary_rows)) do.call(rbind, summary_rows) else summary_empty
  rownames(summary) <- NULL
  list(by_replicate = by_rep, summary = summary)
}

.validation_summarise_main <- function(args = commandArgs(TRUE)) {
  out <- if ("--out" %in% args) args[match("--out", args) + 1L] else
    file.path("script", "validation_20261018", "results",
              if ("--pilot" %in% args) "pilot" else "production")
  source("script/validation_20261018/evaluate.R")
  source("script/validation_20261018/freeze.R")
  if ("--pilot" %in% args) {
    jobs <- readRDS(file.path(out, "jobs.rds"))
    identity <- readRDS(file.path(out, "identity.rds"))
    paths <- file.path(out, paste0(jobs$job_id, ".rds"))
    records <- lapply(paths[file.exists(paths)], readRDS)
  } else {
    source("script/validation_20261018/manifest.R")
    jobs <- validation_jobs(validation_manifest("."))
    ip <- list.files(out, pattern = "^identity[.]rds$", recursive = TRUE,
                     full.names = TRUE)
    if (!length(ip)) stop("No production worker identities")
    identity <- readRDS(ip[1L])
    if (any(!vapply(ip, function(p) identical(readRDS(p)$hash, identity$hash), logical(1))))
      stop("Mixed source identities")
    paths <- list.files(out, pattern = "__[0-9]+[.]rds$", recursive = TRUE,
                        full.names = TRUE)
    records <- lapply(paths, readRDS)
  }
  if (any(!vapply(records, function(x)
      identical(x$runtime_identity, identity$runtime$hash), logical(1))))
    stop("Mixed or missing runtime identity")
  validated <- validate_records(records, jobs, identity$hash)
  summary <- summarise_records(validated)
  write.csv(validated$attempts, file.path(out, "attempts.csv"), row.names = FALSE)
  write.csv(validated$results, file.path(out, "result_rows.csv"), row.names = FALSE)
  posterior <- validation_posterior_interval_tables(records, jobs)
  paired <- paired_method_contrasts(validated)
  calibration <- lapply(records, function(e) {
    bins <- e$metrics$calibration_bins
    if (!is.data.frame(bins) || !nrow(bins)) return(NULL)
    cbind(job_id=e$job_id, scenario_id=e$scenario_id, rep=e$rep, seed=e$seed, bins)
  })
  calibration <- do.call(rbind, calibration)
  if(is.null(calibration))calibration <- data.frame(job_id=character())
  write.csv(calibration, file.path(out,"calibration_bins_by_replicate.csv"), row.names=FALSE)
  provenance <- list(raw_source_identity = identity$hash,
    raw_revision = identity$revision,
    evaluator_files = tools::md5sum(c("script/validation_20261018/evaluate.R",
                                     "script/validation_20261018/summarise.R")),
    evaluated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    evaluation_session = utils::sessionInfo())
  saveRDS(provenance, file.path(out, "evaluation_provenance.rds"))
  write.csv(paired, file.path(out, "paired_method_contrasts.csv"), row.names = FALSE)
  saveRDS(list(validated = validated, summary = summary, paired = paired,
               evaluation_provenance = provenance,
               posterior_intervals_by_replicate = posterior$by_replicate,
               posterior_intervals_summary = posterior$summary),
          file.path(out, "summary.rds"))
  for (nm in names(summary))
    if (is.data.frame(summary[[nm]]))
      write.csv(summary[[nm]], file.path(out, paste0("summary_", nm, ".csv")), row.names = FALSE)
  if (is.data.frame(summary)) write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
  write.csv(posterior$by_replicate,
            file.path(out, "posterior_intervals_by_replicate.csv"), row.names = FALSE)
  write.csv(posterior$summary,
            file.path(out, "posterior_intervals_summary.csv"), row.names = FALSE)
  cat("RESULT_FILES ", length(records), " EXPECTED ", nrow(jobs), "\n", sep = "")
  if (length(records) != nrow(jobs)) stop("Pilot incomplete: missing attempts remain")
  cat(if ("--pilot" %in% args) "PILOT_ACCOUNTED\n" else "PRODUCTION_ACCOUNTED\n")
}

if (sys.nframe() == 0L) .validation_summarise_main()
