# Shared helpers for the mi-posterior real-data harness (G8,
# .unlazy/mi-posterior/GATES.md). Sourced by 00-03. No side effects at
# source time other than defining functions/constants.
#
# The Mondrian real-data harness (mask receipts, split/mondrian conformal
# receipts) lives on branch arc/mondrian-realdata, read-only, via
# `git show <ref>:<path>`. This worktree and that branch share one git
# object database (linked worktrees of the same repo), so `git show` works
# from here without a network fetch as long as origin/arc/mondrian-realdata
# has been fetched at least once (`git fetch origin arc/mondrian-realdata`).

mondrian_ref <- "origin/arc/mondrian-realdata"
mondrian_prefix <- "script/mondrian_confirmation/returned"

`%||%` <- function(x, y) if (is.null(x)) y else x

#' The 10 pre-registered (dataset, arm, seed) real-data cells (task spec).
#' PanTHERIA: mcar + structured x 3 seeds. AVONET: mcar x 3 seeds.
#' FishBase: structured x 1 seed. Matches exactly what exists under
#' script/mondrian_confirmation/returned/ on arc/mondrian-realdata.
planned_cells <- function() {
  seeds3 <- c(20260818L, 20260819L, 20260820L)
  rbind(
    data.frame(dataset = "pantheria", arm = "mcar", seed = seeds3),
    data.frame(dataset = "pantheria", arm = "structured", seed = seeds3),
    data.frame(dataset = "avonet", arm = "mcar", seed = seeds3),
    data.frame(dataset = "fishbase", arm = "structured", seed = 20260818L)
  )
}

cell_name <- function(dataset, arm, seed) sprintf("%s-%s-m%d", dataset, arm, as.integer(seed))

#' Extract one file from the mondrian-realdata branch via `git show` into
#' `out_file`. Fails loudly (does not silently write an empty/partial file)
#' if the ref or path does not exist.
git_show_to_file <- function(rel_path, out_file, ref = mondrian_ref) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  spec <- paste0(ref, ":", rel_path)
  status <- system2("git", c("show", spec), stdout = out_file, stderr = FALSE)
  if (!identical(status, 0L) || !file.exists(out_file) || file.size(out_file) == 0L) {
    unlink(out_file)
    stop("git show failed for ", spec,
         " (is arc/mondrian-realdata fetched? try `git fetch origin arc/mondrian-realdata`)",
         call. = FALSE)
  }
  out_file
}

#' Read an .rds straight from the mondrian-realdata branch into a temp file
#' (not persisted under inputs/ -- only mask receipts are persisted there;
#' split.rds/mondrian.rds are read on demand for the conformal comparison).
read_mondrian_rds <- function(dataset, arm, seed, leaf) {
  name <- cell_name(dataset, arm, seed)
  rel <- file.path(mondrian_prefix, name, leaf)
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  git_show_to_file(rel, tf)
  readRDS(tf)
}

#' Classify each column of a trait data.frame for the posterior method's
#' continuous-only restriction (design.md section 4: "any non-continuous
#' trait ... errors"). Returns a data.frame(trait, class, keep, reason).
classify_traits <- function(df) {
  cls <- vapply(df, function(x) paste(class(x), collapse = "/"), character(1))
  is_ordered <- vapply(df, is.ordered, logical(1))
  is_factor <- vapply(df, is.factor, logical(1)) & !is_ordered
  is_integer <- vapply(df, is.integer, logical(1))
  is_cont <- vapply(df, function(x) is.numeric(x) && !is.integer(x), logical(1))
  reason <- ifelse(is_cont, "kept: numeric, auto-detects as pigauto trait type 'continuous'",
             ifelse(is_ordered, "dropped: ordered factor -> pigauto type 'ordinal', not continuous",
             ifelse(is_factor, "dropped: (unordered) factor -> pigauto type 'binary'/'categorical', not continuous",
             ifelse(is_integer, "dropped: integer -> pigauto type 'count', not continuous",
                    "dropped: unrecognised column class, not continuous"))))
  data.frame(trait = names(df), class = cls, keep = is_cont, reason = reason,
             row.names = NULL, stringsAsFactors = FALSE)
}

#' Pre-registered downstream slope-check pairs applicable to one dataset.
pairs_for_dataset <- function(pairs, dataset) {
  Filter(function(p) identical(p$dataset, dataset), pairs)
}

# ---- Shared reporting / acceptance logic (used by 01 receipts, 02, 03) ----

#' Read the mi_posterior.rds receipt for each planned (dataset, arm, seed)
#' row, if present. Returns a named list (by cell_name) of either the
#' receipt (as written by 01_run.R) or list(status = "missing").
collect_receipts <- function(outdir, planned) {
  out <- list()
  for (i in seq_len(nrow(planned))) {
    nm <- cell_name(planned$dataset[[i]], planned$arm[[i]], planned$seed[[i]])
    f <- file.path(outdir, nm, "mi_posterior.rds")
    out[[nm]] <- if (file.exists(f)) readRDS(f) else list(status = "missing", name = nm)
  }
  out
}

#' One row of a `metrics` data.frame (trait, coverage/width columns under
#' a given name) for `trait`, or NA placeholders if unavailable.
metric_row <- function(conformal_obj, trait, cov_col, width_col) {
  if (is.null(conformal_obj) || !identical(conformal_obj$status, "ok") ||
        is.null(conformal_obj$metrics) || !(trait %in% conformal_obj$metrics$trait)) {
    return(list(coverage = NA_real_, width = NA_real_,
                status = if (is.null(conformal_obj)) "missing" else conformal_obj$status %||% "error"))
  }
  r <- conformal_obj$metrics[conformal_obj$metrics$trait == trait, , drop = FALSE][1L, ]
  list(coverage = r[[cov_col]], width = r[[width_col]], status = "ok")
}

#' Per-trait coverage table: model-based (posterior) vs split-conformal vs
#' Mondrian-conformal, for every kept trait of every "ok" receipt.
build_coverage_table <- function(receipts) {
  rows <- list()
  for (rc in receipts) {
    if (!identical(rc$status, "ok")) next
    mc <- rc$model_coverage
    for (i in seq_len(nrow(mc))) {
      tr <- mc$trait[[i]]
      sp <- metric_row(rc$split_conformal, tr, "coverage", "width")
      mo <- metric_row(rc$mondrian_conformal, tr, "coverage", "width")
      rows[[length(rows) + 1L]] <- data.frame(
        dataset = rc$dataset, arm = rc$arm, seed = rc$seed, trait = tr,
        n_masked = mc$n_masked[[i]],
        model_coverage = mc$coverage[[i]], model_width = mc$median_width[[i]],
        split_coverage = sp$coverage, split_width = sp$width, split_status = sp$status,
        mondrian_coverage = mo$coverage, mondrian_width = mo$width, mondrian_status = mo$status,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) {
    return(data.frame(dataset = character(), arm = character(), seed = integer(), trait = character(),
                       n_masked = integer(), model_coverage = numeric(), model_width = numeric(),
                       split_coverage = numeric(), split_width = numeric(), split_status = character(),
                       mondrian_coverage = numeric(), mondrian_width = numeric(), mondrian_status = character(),
                       stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

#' Downstream slope table: reference (phylolm lambda, originally-observed
#' rows) vs MI-pooled (gls corPagel + pool_mi, all rows), one row per pair
#' per successfully-run receipt that attempted it.
build_slope_table <- function(receipts) {
  rows <- list()
  for (rc in receipts) {
    if (!identical(rc$status, "ok") || is.null(rc$pairs)) next
    for (pr in rc$pairs) {
      ref <- pr$reference; mip <- pr$mi
      rows[[length(rows) + 1L]] <- data.frame(
        dataset = rc$dataset, arm = rc$arm, seed = rc$seed,
        response = pr$response, predictor = pr$predictor,
        ref_status = ref$status, ref_slope = ref$slope %||% NA_real_, ref_se = ref$se %||% NA_real_,
        ref_lambda = ref$lambda %||% NA_real_,
        mi_status = mip$status, mi_slope = mip$slope %||% NA_real_, mi_se = mip$se %||% NA_real_,
        mi_df = mip$df %||% NA_real_, mi_fmi = mip$fmi %||% NA_real_,
        rel_diff = pr$rel_diff, se_ratio = pr$se_ratio,
        within_5pct = if (is.finite(pr$rel_diff %||% NA_real_)) abs(pr$rel_diff) <= 0.05 else NA,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) {
    return(data.frame(dataset = character(), arm = character(), seed = integer(),
                       response = character(), predictor = character(),
                       ref_status = character(), ref_slope = numeric(), ref_se = numeric(),
                       ref_lambda = numeric(), mi_status = character(), mi_slope = numeric(),
                       mi_se = numeric(), mi_df = numeric(), mi_fmi = numeric(),
                       rel_diff = numeric(), se_ratio = numeric(), within_5pct = logical(),
                       stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

#' REALDATA_COMPLETE (G8) check: every planned cell has a receipt with
#' status "ok", and every pre-registered pair has an attempted result (a
#' `pairs` entry, `ok` or `error`, from at least one "ok" receipt for that
#' pair's dataset -- an entire dataset with no successful cell yields no
#' pair result and so is NOT complete). The 5% slope criterion (G8's
#' criterion) is reported, not required, for REALDATA_COMPLETE.
check_acceptance <- function(receipts, planned, pairs_list) {
  cell_status <- data.frame(
    name = vapply(seq_len(nrow(planned)), function(i)
      cell_name(planned$dataset[[i]], planned$arm[[i]], planned$seed[[i]]), character(1)),
    dataset = planned$dataset, arm = planned$arm, seed = planned$seed,
    status = vapply(seq_len(nrow(planned)), function(i) {
      nm <- cell_name(planned$dataset[[i]], planned$arm[[i]], planned$seed[[i]])
      receipts[[nm]]$status %||% "missing"
    }, character(1)),
    stringsAsFactors = FALSE
  )
  missing_cells <- cell_status[cell_status$status != "ok", , drop = FALSE]

  pair_status <- do.call(rbind, lapply(pairs_list, function(p) {
    found <- FALSE; ok_slope <- FALSE; within5 <- NA
    for (rc in receipts) {
      if (!identical(rc$status, "ok") || !identical(rc$dataset, p$dataset) || is.null(rc$pairs)) next
      for (pr in rc$pairs) {
        if (identical(pr$response, p$response) && identical(pr$predictor, p$predictor)) {
          found <- TRUE
          if (identical(pr$reference$status, "ok") && identical(pr$mi$status, "ok")) {
            ok_slope <- TRUE
            if (is.finite(pr$rel_diff %||% NA_real_)) within5 <- isTRUE(abs(pr$rel_diff) <= 0.05)
          }
        }
      }
    }
    data.frame(dataset = p$dataset, response = p$response, predictor = p$predictor,
               has_result = found, has_ok_slope = ok_slope, within_5pct = within5,
               stringsAsFactors = FALSE)
  }))

  list(
    complete = nrow(missing_cells) == 0L && all(pair_status$has_result),
    cell_status = cell_status, missing_cells = missing_cells, pair_status = pair_status
  )
}

#' Summarise mi$posterior$diagnostics down to the max split R-hat and min
#' bulk ESS over the Sigma_P, Sigma_E and lambda elements (design-review
#' addition, 2026-09-24: "record convergence for every real-data fit"),
#' plus the diagnostics object's own attr(,"converged"). Reads only fields
#' the frozen API (design.md section 4) already documents
#' (mi$posterior$diagnostics: parameter, rhat, ess_bulk; attr "converged")
#' -- the frozen API itself is unchanged.
summarize_convergence <- function(diagnostics) {
  if (is.null(diagnostics) || !is.data.frame(diagnostics) || !nrow(diagnostics)) {
    return(list(max_rhat = NA_real_, min_ess = NA_real_, converged = NA))
  }
  sel <- grepl("^(Sigma_P|Sigma_E|lambda)", diagnostics$parameter)
  if (!any(sel)) sel <- rep(TRUE, nrow(diagnostics))  # fallback: don't silently report NA
  list(
    max_rhat = suppressWarnings(max(diagnostics$rhat[sel], na.rm = TRUE)),
    min_ess = suppressWarnings(min(diagnostics$ess_bulk[sel], na.rm = TRUE)),
    converged = attr(diagnostics, "converged") %||% NA
  )
}

#' Per-cell convergence table (dataset, arm, seed, max_rhat, min_ess,
#' converged), one row per "ok" receipt. Consumed by 02_summarise.R
#' (carried through to convergence_table.csv / summary.md) and
#' 03_acceptance.R (NONCONVERGED flags, reported separately from the
#' coverage headline; never gates REALDATA_COMPLETE).
build_convergence_table <- function(receipts) {
  rows <- list()
  for (rc in receipts) {
    if (!identical(rc$status, "ok")) next
    cv <- rc$convergence %||% list(max_rhat = NA_real_, min_ess = NA_real_, converged = NA)
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = rc$dataset, arm = rc$arm, seed = rc$seed, name = rc$name,
      max_rhat = cv$max_rhat, min_ess = cv$min_ess, converged = cv$converged,
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) {
    return(data.frame(dataset = character(), arm = character(), seed = integer(), name = character(),
                       max_rhat = numeric(), min_ess = numeric(), converged = logical(),
                       stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

#' Build a small synthetic fixture (2 fake cells, 1 fake pair) for
#' --selftest in 02_summarise.R / 03_acceptance.R, entirely separate from
#' the real planned_cells()/mi_realdata_pairs so a selftest run can never
#' be mistaken for real-data acceptance.
build_selftest_fixture <- function(tmp_outdir) {
  planned <- data.frame(dataset = "synth", arm = "mcar", seed = c(1L, 2L))
  pairs_list <- list(list(dataset = "synth", response = "t2", predictor = "t1",
                           rationale = "synthetic selftest pair"))
  ok_seeds <- c(1L)  # seed 2 deliberately left without a receipt -> incomplete
  for (i in seq_len(nrow(planned))) {
    if (!(planned$seed[[i]] %in% ok_seeds)) next
    nm <- cell_name(planned$dataset[[i]], planned$arm[[i]], planned$seed[[i]])
    dir.create(file.path(tmp_outdir, nm), recursive = TRUE, showWarnings = FALSE)
    mc <- data.frame(trait = c("t1", "t2"), n_masked = c(10L, 10L), n_matched = c(10L, 10L),
                      coverage = c(0.94, 0.96), median_width = c(1.1, 0.9), stringsAsFactors = FALSE)
    sp_metrics <- data.frame(trait = c("t1", "t2"), n_masked = c(10L, 10L), n_interval = c(10L, 10L),
                              coverage = c(0.90, 0.93), width = c(1.3, 1.0), stringsAsFactors = FALSE)
    mo_metrics <- data.frame(trait = c("t1", "t2"), n_masked = c(10L, 10L), n_interval = c(10L, 10L),
                              coverage = c(0.95, 0.94), width = c(1.2, 0.95), stringsAsFactors = FALSE)
    ref <- list(status = "ok", n = 50L, slope = 2.0, se = 0.05, lambda = 0.8, y_log = FALSE, x_log = FALSE)
    mip <- list(status = "ok", m_used = 20L, m_total = 20L, slope = 2.02, se = 0.06, df = 15, fmi = 0.30,
                y_log = FALSE, x_log = FALSE)
    pr <- list(list(dataset = "synth", response = "t2", predictor = "t1",
                     rationale = "synthetic selftest pair", reference = ref, mi = mip,
                     rel_diff = (mip$slope - ref$slope) / ref$slope, se_ratio = mip$se / ref$se))
    receipt <- list(
      status = "ok", dataset = planned$dataset[[i]], arm = planned$arm[[i]], seed = planned$seed[[i]],
      name = nm, n_species = 50L, kept_traits = c("t1", "t2"),
      dropped_traits = data.frame(trait = "t3", class = "factor", keep = FALSE,
                                   reason = "dropped: synthetic factor trait", stringsAsFactors = FALSE),
      wall_time_s = 12.3, model_coverage = mc,
      split_conformal = list(status = "ok", error = NULL, metrics = sp_metrics),
      mondrian_conformal = list(status = "ok", error = NULL, metrics = mo_metrics),
      diagnostics = structure(
        data.frame(parameter = c("Sigma_P[1,1]", "Sigma_E[1,1]", "lambda[1]"),
                    rhat = c(1.09, 1.02, 1.01), ess_bulk = c(250, 900, 850),
                    stringsAsFactors = FALSE),
        converged = FALSE  # deliberately non-converged, to exercise 03's NONCONVERGED path
      ),
      pairs = pr
    )
    receipt$convergence <- summarize_convergence(receipt$diagnostics)
    saveRDS(receipt, file.path(tmp_outdir, nm, "mi_posterior.rds"))
  }
  list(planned = planned, pairs_list = pairs_list)
}
