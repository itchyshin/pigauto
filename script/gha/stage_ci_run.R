#!/usr/bin/env Rscript
# script/gha/stage_ci_run.R
#
# Stage a slim per-run summary into useful/ci_runs/<date>_<run-id>/
# for the draft-PR step at the end of the workflow.
#
# Usage (from CI):
#   Rscript script/gha/stage_ci_run.R "<github-run-id>"
#
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md §4.5.

DATASETS <- c("avonet","pantheria","amphibio","bien","globtherm","leptraits")

stage_ci_run <- function(results_root = file.path("script","gha","results"),
                          ci_runs_root = file.path("useful","ci_runs"),
                          run_id,
                          date_str = format(Sys.Date(), "%Y-%m-%d")) {
  stopifnot(is.character(run_id), nchar(run_id) > 0L)
  staged <- file.path(ci_runs_root, sprintf("%s_%s", date_str, run_id))
  dir.create(file.path(staged, "pigauto_per_dataset"),
              recursive = TRUE, showWarnings = FALSE)

  rep_src <- file.path(results_root, "cross_dataset", "report.md")
  if (file.exists(rep_src)) {
    file.copy(rep_src, file.path(staged, "report.md"), overwrite = TRUE)
  } else {
    writeLines("# Report unavailable (cross_dataset/report.md missing)",
                file.path(staged, "report.md"))
  }

  agg_timings <- list()
  for (d in DATASETS) {
    md_src <- file.path(results_root, "_artifacts",
                         paste0("bench-", d), "results.md")
    if (file.exists(md_src)) {
      file.copy(md_src,
                 file.path(staged, "pigauto_per_dataset", paste0(d, ".md")),
                 overwrite = TRUE)
    }
    t_src <- file.path(results_root, "_artifacts",
                        paste0("bench-", d), "timings.json")
    if (file.exists(t_src)) {
      agg_timings[[d]] <- jsonlite::read_json(t_src, simplifyVector = TRUE)
    }
  }
  jsonlite::write_json(agg_timings,
                       file.path(staged, "timings.json"),
                       auto_unbox = TRUE, pretty = TRUE)

  invisible(staged)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1L) stop("usage: stage_ci_run.R <run-id>")
  out <- stage_ci_run(run_id = args[[1L]])
  cat(sprintf("[stage] staged at %s\n", out))
}
