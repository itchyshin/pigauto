#!/usr/bin/env Rscript
# script/gha/make_headtohead_report.R
#
# Aggregate job's report generator. Reads pigauto's fresh per-dataset
# .rds + BACE pinned snapshot .rds, produces a combined head-to-head
# report.md + summary.rds under script/gha/results/cross_dataset/.
#
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md.

DATASETS <- c("avonet", "pantheria", "amphibio", "bien", "globtherm", "leptraits")
# Types whose primary metric is RMSE (lower is better). Ordinal is excluded
# because pigauto's CI evaluator treats ordered-factor truth values as
# factors and records accuracy (not RMSE), to match BACE's snapshot which
# also records accuracy for ordinal. Including ordinal here previously
# caused the h2h merge to read NA RMSE from both sides and produce NA
# winners for every ordinal trait in the report.
CONTINUOUS_TYPES <- c("continuous", "count", "proportion", "zi_count")

load_pigauto_results <- function() {
  base <- file.path("script", "gha", "results", "_artifacts")
  parts <- lapply(DATASETS, function(d) {
    f <- file.path(base, paste0("bench-", d), "results.rds")
    if (!file.exists(f)) {
      message(sprintf("[h2h] missing pigauto results for %s; skipping", d))
      return(NULL)
    }
    readRDS(f)
  })
  parts <- Filter(function(x) !is.null(x) && nrow(x) > 0L, parts)
  if (length(parts) == 0L) return(NULL)
  do.call(rbind, parts)
}

load_bace_snapshot <- function() {
  base <- file.path("useful", "bace_results_snapshot")
  parts <- lapply(DATASETS, function(d) {
    f <- file.path(base, paste0(d, ".rds"))
    if (!file.exists(f)) {
      message(sprintf("[h2h] missing BACE snapshot for %s; skipping", d))
      return(NULL)
    }
    readRDS(f)
  })
  parts <- Filter(function(x) !is.null(x) && nrow(x) > 0L, parts)
  if (length(parts) == 0L) return(NULL)
  do.call(rbind, parts)
}

# Decide winner per (dataset, trait, type) by primary metric.
.h2h_winner <- function(pig, bace, type) {
  if (is.na(pig) || is.na(bace)) return(NA_character_)
  if (type %in% CONTINUOUS_TYPES) {
    if (pig < bace) "pigauto" else if (bace < pig) "bace" else "tie"
  } else {
    if (pig > bace) "pigauto" else if (bace > pig) "bace" else "tie"
  }
}

build_h2h_report <- function(combined, out_dir) {
  stopifnot(is.data.frame(combined),
            all(c("dataset","trait","type","method","rmse","accuracy")
                %in% colnames(combined)))

  # coverage_95 and interval_width are pigauto-only columns; the BACE
  # snapshot has them as NA. Aggregate them too so the report can flag
  # pigauto's calibration credentials.
  for (col in c("coverage_95", "interval_width")) {
    if (!(col %in% colnames(combined))) combined[[col]] <- NA_real_
  }
  med <- stats::aggregate(
    cbind(rmse, mae, pearson_r, accuracy, brier,
          coverage_95, interval_width, time_sec) ~ dataset + trait + type + method,
    data = combined,
    FUN  = function(x) stats::median(x, na.rm = TRUE),
    na.action = stats::na.pass
  )

  pig <- med[med$method == "pigauto_ci",   , drop = FALSE]
  bac <- med[med$method == "BACE_snapshot", , drop = FALSE]
  m <- merge(pig, bac, by = c("dataset","trait","type"),
             suffixes = c("_pigauto","_bace"), all = TRUE)

  m$pigauto <- ifelse(m$type %in% CONTINUOUS_TYPES,
                      m$rmse_pigauto, m$accuracy_pigauto)
  m$bace    <- ifelse(m$type %in% CONTINUOUS_TYPES,
                      m$rmse_bace, m$accuracy_bace)
  # Label the metric and its preferred direction so the reader can
  # interpret the pigauto / bace columns at a glance. Continuous-family
  # types report RMSE (lower is better); discrete types report
  # accuracy (higher is better).
  m$metric <- ifelse(m$type %in% CONTINUOUS_TYPES, "rmse", "acc")
  m$better <- ifelse(m$type %in% CONTINUOUS_TYPES, "lower", "higher")
  m$winner <- vapply(seq_len(nrow(m)),
                      function(i) .h2h_winner(m$pigauto[i], m$bace[i], m$type[i]),
                      character(1))

  summary_tbl <- m[, c("dataset","trait","type","metric","better",
                       "pigauto","bace","winner")]

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(summary_tbl, file.path(out_dir, "summary.rds"))

  md <- c(
    "# pigauto vs BACE — cross-dataset head-to-head",
    "",
    sprintf("Run date: %s (UTC).", format(Sys.time(), tz = "UTC")),
    sprintf("Datasets: %s.", paste(DATASETS, collapse = ", ")),
    "",
    "## Wins / ties / losses per dataset",
    ""
  )
  if (nrow(summary_tbl) > 0L) {
    wtl <- as.data.frame(table(summary_tbl$dataset,
                                factor(summary_tbl$winner,
                                       levels = c("pigauto","tie","bace"))),
                         stringsAsFactors = FALSE)
    colnames(wtl) <- c("dataset", "winner", "n")
    md <- c(md, knitr::kable(wtl, format = "markdown"))
    md <- c(md, "", "## Per-trait detail (medians across imputations)", "",
            knitr::kable(summary_tbl, format = "markdown", digits = 3))

    # Brier comparison (categorical / binary / ordinal). Both methods
    # expose this; lower is better.
    brier_tbl <- m[!is.na(m$brier_pigauto) | !is.na(m$brier_bace),
                    c("dataset","trait","type","brier_pigauto","brier_bace"),
                    drop = FALSE]
    if (nrow(brier_tbl) > 0L) {
      brier_tbl$brier_winner <- vapply(seq_len(nrow(brier_tbl)),
        function(i) {
          p <- brier_tbl$brier_pigauto[i]; b <- brier_tbl$brier_bace[i]
          if (is.na(p) || is.na(b)) NA_character_
          else if (p < b) "pigauto" else if (b < p) "bace" else "tie"
        }, character(1))
      md <- c(md, "", "## Brier score (lower is better)", "",
              knitr::kable(brier_tbl, format = "markdown", digits = 3))
    }

    # Coverage / interval-width is pigauto-only. Show as pigauto's
    # calibration credentials; BACE has no comparable column.
    cov_tbl <- m[!is.na(m$coverage_95_pigauto),
                   c("dataset","trait","type","coverage_95_pigauto",
                     "interval_width_pigauto"),
                   drop = FALSE]
    if (nrow(cov_tbl) > 0L) {
      colnames(cov_tbl) <- c("dataset","trait","type",
                              "coverage_95","interval_width")
      md <- c(md, "",
              "## Pigauto 95% conformal coverage + interval width (BACE has none)",
              "",
              "Target coverage = 0.95. Coverage close to target with smaller",
              "interval width indicates well-calibrated, tight predictions.",
              "",
              knitr::kable(cov_tbl, format = "markdown", digits = 3))
    }
  } else {
    md <- c(md, "_No head-to-head rows: pigauto / BACE inputs both empty._")
  }
  writeLines(md, file.path(out_dir, "report.md"))

  invisible(list(summary = summary_tbl,
                 report_path = file.path(out_dir, "report.md")))
}

# When invoked as a script (not sourced), wire up I/O.
if (sys.nframe() == 0L) {
  source(file.path("script", "gha", "_ci_config.R"))
  out_dir <- file.path("script", "gha", "results", "cross_dataset")
  pig <- load_pigauto_results()
  bac <- load_bace_snapshot()
  combined <- rbind(pig, bac)
  if (is.null(combined) || nrow(combined) == 0L) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(c("# pigauto vs BACE — cross-dataset head-to-head",
                  "",
                  "No input rows from either pigauto or BACE; check matrix logs."),
               file.path(out_dir, "report.md"))
    quit(save = "no", status = 0L)
  }
  build_h2h_report(combined, out_dir)
  cat(sprintf("[h2h] wrote report to %s\n", file.path(out_dir, "report.md")))
}
