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
#
# Offline runs (Totoro, DRAC compute nodes; 2026-09-24): machines without
# that git history read everything from inputs/<cell>/ instead.
# `00_fetch_masks.R --all --with-conformal` writes, per cell,
# mask_receipt.rds and conformal_metrics.rds (a small derived per-trait
# table for split and Mondrian conformal, see derive_conformal_metrics();
# never the raw split.rds/mondrian.rds objects). 01_run.R uses
# conformal_metrics.rds when present and, with MI_REALDATA_OFFLINE=1,
# refuses to fall back to git.

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

#' Read an .rds straight from the mondrian-realdata branch into a temp file.
#' The raw split.rds/mondrian.rds objects are never persisted under
#' inputs/; only their derived per-trait metrics are (see
#' derive_conformal_metrics()).
read_mondrian_rds <- function(dataset, arm, seed, leaf) {
  name <- cell_name(dataset, arm, seed)
  rel <- file.path(mondrian_prefix, name, leaf)
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  git_show_to_file(rel, tf)
  readRDS(tf)
}

#' `git rev-parse <spec>` or NA_character_ (provenance only, never fatal).
git_rev_parse <- function(spec, dir = NULL) {
  a <- c(if (!is.null(dir)) c("-C", dir), "rev-parse", spec)
  out <- tryCatch(suppressWarnings(system2("git", a, stdout = TRUE, stderr = FALSE)),
                  error = function(e) character(0))
  if (length(out) == 1L && nzchar(out) && is.null(attr(out, "status"))) out else NA_character_
}

#' Reduce one split.rds / mondrian.rds receipt from arc/mondrian-realdata
#' to a small per-trait table: trait, n_masked, n_interval, coverage,
#' mean_width, median_width. `coverage` and `mean_width` are the receipt's
#' own metrics$coverage and metrics$width (width there is a MEAN of hi - lo
#' over cells with a finite interval, 01_run_masked_confirmation.R);
#' median_width is the median of hi - lo over the same cells, from
#' receipt$cells. Fails (status "error") if the cells do not reproduce the
#' receipt's mean width, so the two summaries are provably like for like.
derive_conformal_metrics <- function(res) {
  if (inherits(res, "condition")) {
    return(list(status = "error", error = conditionMessage(res), metrics = NULL))
  }
  if (!identical(res$status, "ok") || is.null(res$metrics)) {
    return(list(status = res$status %||% "error", error = res$error %||% "no metrics in receipt",
                metrics = NULL))
  }
  m <- res$metrics; cells <- res$cells
  if (is.null(cells) || !nrow(cells)) {
    return(list(status = "error", metrics = NULL,
                error = "receipt has no per-cell intervals ($cells); median width not computable"))
  }
  w_by <- lapply(m$trait, function(tr) {
    sel <- cells$trait == tr
    cells$hi[sel] - cells$lo[sel]
  })
  mean_chk <- vapply(w_by, function(w) if (length(w)) mean(w) else NA_real_, numeric(1))
  n_chk <- vapply(w_by, length, integer(1))
  if (anyNA(mean_chk) || any(n_chk != m$n_interval) ||
        any(abs(mean_chk - m$width) > 1e-8 * pmax(1, abs(m$width)))) {
    return(list(status = "error", metrics = NULL,
                error = "receipt $cells do not reproduce $metrics (n_interval / mean width)"))
  }
  list(status = "ok", error = NULL, method = res$method %||% NA_character_,
       metrics = data.frame(trait = m$trait, n_masked = as.integer(m$n_masked),
                            n_interval = as.integer(m$n_interval), coverage = m$coverage,
                            mean_width = m$width,
                            median_width = vapply(w_by, stats::median, numeric(1)),
                            stringsAsFactors = FALSE))
}

#' Path of the offline conformal-metrics file for one cell.
conformal_metrics_path <- function(inputs_dir, name) {
  file.path(inputs_dir, name, "conformal_metrics.rds")
}

#' Split and Mondrian derived metrics for one cell, read from
#' arc/mondrian-realdata via git (needs the git history). Returns
#' list(split, mondrian, source). A git/read failure of either leaf is an
#' error (the caller decides whether to write anything); a receipt whose
#' own status is not "ok" is recorded as such, not hidden.
fetch_conformal_metrics <- function(dataset, arm, seed) {
  name <- cell_name(dataset, arm, seed)
  one <- function(leaf) derive_conformal_metrics(read_mondrian_rds(dataset, arm, seed, leaf))
  list(
    split = one("split.rds"), mondrian = one("mondrian.rds"),
    source = list(
      ref = mondrian_ref, commit = git_rev_parse(mondrian_ref),
      split_blob = git_rev_parse(paste0(mondrian_ref, ":", file.path(mondrian_prefix, name, "split.rds"))),
      mondrian_blob = git_rev_parse(paste0(mondrian_ref, ":", file.path(mondrian_prefix, name, "mondrian.rds"))),
      extracted_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    )
  )
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

#' One kept trait's conformal metrics (derived by derive_conformal_metrics())
#' or NA placeholders, with the reason, if unavailable.
metric_row <- function(conformal_obj, trait) {
  na <- list(n_masked = NA_integer_, n_interval = NA_integer_, coverage = NA_real_,
             mean_width = NA_real_, median_width = NA_real_)
  if (is.null(conformal_obj)) return(c(na, status = "missing"))
  if (!identical(conformal_obj$status, "ok")) return(c(na, status = conformal_obj$status %||% "error"))
  m <- conformal_obj$metrics
  if (is.null(m) || !all(c("mean_width", "median_width") %in% names(m))) {
    return(c(na, status = "no mean/median width columns"))
  }
  if (!(trait %in% m$trait)) return(c(na, status = "trait absent"))
  r <- m[m$trait == trait, , drop = FALSE][1L, ]
  list(n_masked = r$n_masked %||% NA_integer_, n_interval = r$n_interval, coverage = r$coverage,
       mean_width = r$mean_width, median_width = r$median_width, status = "ok")
}

#' Per-trait coverage table: model-based (posterior) vs split-conformal vs
#' Mondrian-conformal, for every kept trait of every "ok" receipt. Widths
#' are reported like for like: mean AND median of (upper - lower) over the
#' masked cells, on the original trait scale, for all three methods.
build_coverage_table <- function(receipts) {
  rows <- list()
  for (rc in receipts) {
    if (!identical(rc$status, "ok")) next
    mc <- rc$model_coverage
    if (is.null(mc)) next
    for (i in seq_len(nrow(mc))) {
      tr <- mc$trait[[i]]
      sp <- metric_row(rc$split_conformal, tr)
      mo <- metric_row(rc$mondrian_conformal, tr)
      rows[[length(rows) + 1L]] <- data.frame(
        dataset = rc$dataset, arm = rc$arm, seed = rc$seed, trait = tr,
        n_masked = mc$n_masked[[i]], n_matched = mc$n_matched[[i]],
        model_coverage = mc$coverage[[i]],
        split_coverage = sp$coverage, mondrian_coverage = mo$coverage,
        model_mean_width = mc$mean_width[[i]] %||% NA_real_,
        split_mean_width = sp$mean_width, mondrian_mean_width = mo$mean_width,
        model_median_width = mc$median_width[[i]] %||% NA_real_,
        split_median_width = sp$median_width, mondrian_median_width = mo$median_width,
        split_n_interval = sp$n_interval, mondrian_n_interval = mo$n_interval,
        split_status = sp$status, mondrian_status = mo$status,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) {
    return(data.frame(dataset = character(), arm = character(), seed = integer(), trait = character(),
                       n_masked = integer(), n_matched = integer(), model_coverage = numeric(),
                       split_coverage = numeric(), mondrian_coverage = numeric(),
                       model_mean_width = numeric(), split_mean_width = numeric(),
                       mondrian_mean_width = numeric(), model_median_width = numeric(),
                       split_median_width = numeric(), mondrian_median_width = numeric(),
                       split_n_interval = integer(), mondrian_n_interval = integer(),
                       split_status = character(), mondrian_status = character(),
                       stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

#' Paired contrast between the MI-pooled and the reference slope of one
#' pair in one cell: diff = mi - ref, rel_diff = diff / ref (NA when ref is
#' 0 or either leg failed), diff_ref_se = diff / ref_se (the difference in
#' units of the complete-rows slope's SE), se_ratio = mi_se / ref_se.
#' Single source of truth for 01_run.R, the slope table and the gate.
pair_contrast <- function(ref, mip) {
  out <- list(diff = NA_real_, rel_diff = NA_real_, diff_ref_se = NA_real_, se_ratio = NA_real_)
  if (!identical(ref$status, "ok") || !identical(mip$status, "ok")) return(out)
  if (!all(is.finite(c(ref$slope, mip$slope, ref$se, mip$se)))) return(out)
  out$diff <- mip$slope - ref$slope
  if (ref$slope != 0) out$rel_diff <- out$diff / ref$slope
  if (ref$se > 0) {
    out$diff_ref_se <- out$diff / ref$se
    out$se_ratio <- mip$se / ref$se
  }
  out
}

#' Downstream slope table: reference (phylolm lambda on the
#' originally-observed rows) vs MI-pooled (the same phylolm fit to each
#' completion on the same rows, Rubin-pooled), one row per pair per "ok"
#' receipt that attempted it. within_5pct is |rel_diff| <= 0.05 (reported,
#' not gated; see check_acceptance()). m_used / m_total / n_nonfinite show
#' how many completions entered the pool; check_acceptance() fails any
#' pair-cell that pooled fewer than all of them.
build_slope_table <- function(receipts) {
  rows <- list()
  for (rc in receipts) {
    if (!identical(rc$status, "ok") || is.null(rc$pairs)) next
    for (pr in rc$pairs) {
      ref <- pr$reference; mip <- pr$mi
      pc <- pair_contrast(ref, mip)
      rows[[length(rows) + 1L]] <- data.frame(
        dataset = rc$dataset, arm = rc$arm, seed = rc$seed,
        response = pr$response, predictor = pr$predictor,
        y_log = ref$y_log %||% mip$y_log %||% NA, x_log = ref$x_log %||% mip$x_log %||% NA,
        n = ref$n %||% mip$n %||% NA_integer_,
        ref_status = ref$status, ref_slope = ref$slope %||% NA_real_, ref_se = ref$se %||% NA_real_,
        ref_lambda = ref$lambda %||% NA_real_,
        mi_status = mip$status, m_used = mip$m_used %||% NA_integer_,
        m_total = mip$m_total %||% NA_integer_, n_nonfinite = mip$n_nonfinite %||% NA_integer_,
        mi_slope = mip$slope %||% NA_real_, mi_se = mip$se %||% NA_real_,
        mi_df = mip$df %||% NA_real_, mi_fmi = mip$fmi %||% NA_real_,
        diff = pc$diff, rel_diff = pc$rel_diff, diff_ref_se = pc$diff_ref_se, se_ratio = pc$se_ratio,
        within_5pct = if (is.finite(pc$rel_diff)) abs(pc$rel_diff) <= 0.05 else NA,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) {
    return(data.frame(dataset = character(), arm = character(), seed = integer(),
                       response = character(), predictor = character(),
                       y_log = logical(), x_log = logical(), n = integer(),
                       ref_status = character(), ref_slope = numeric(), ref_se = numeric(),
                       ref_lambda = numeric(), mi_status = character(), m_used = integer(),
                       m_total = integer(), n_nonfinite = integer(), mi_slope = numeric(),
                       mi_se = numeric(), mi_df = numeric(), mi_fmi = numeric(),
                       diff = numeric(), rel_diff = numeric(), diff_ref_se = numeric(),
                       se_ratio = numeric(), within_5pct = logical(),
                       stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

#' Receipt schema written by 01_run.R. Schema 3 (2026-09-24 repair round):
#' n_matched counts only masked cells with a finite model interval, and
#' coverage/width no longer drop NA. Receipts with an older schema (2:
#' n_matched counted rows, na.rm = TRUE; none/1: median-only model width,
#' raw conformal metrics, trait_map analysis scale) fail the gate and must
#' be re-run.
receipt_schema_current <- 3L

#' Code-provenance sources the gate accepts: MI_POST_SHA (set from
#' code_dir/SHA by 12_totoro_run.sh / 10_fir.sbatch) or the git HEAD of a
#' checkout whose R/ and DESCRIPTION had no uncommitted changes.
sha_sources_ok <- c("MI_POST_SHA", "git HEAD")

#' REALDATA_COMPLETE (G8) check, fail-closed (orchestrator decision R1,
#' 2026-09-24; tightened in the repair round the same day). Complete only
#' if ALL of:
#'   1. every planned cell has a receipt with status "ok" (a crashed or
#'      killed run leaves status "running", which fails);
#'   2. all ok receipts carry ONE non-NA code_sha whose code_sha_source is
#'      MI_POST_SHA or a clean git HEAD;
#'   3. every ok receipt has the current schema and a sampler record with
#'      no MI_POST_NITER / MI_POST_BURNIN override (pigauto defaults);
#'   4. every ok receipt records >= 1 kept trait, and for every trait in
#'      kept_traits or in the model coverage table: exactly one coverage
#'      row, n_matched == n_masked (n_matched counts masked cells with a
#'      finite interval), and, when n_masked > 0, a finite model coverage;
#'   5. every ok receipt read split AND Mondrian conformal results with
#'      status "ok" and, for every masked kept trait, a row with finite
#'      coverage and n_interval == n_masked == the model's n_masked;
#'   6. every pre-registered pair has at least one cell whose reference
#'      and MI legs both have status "ok", finite slopes and SEs, and an MI
#'      pool of ALL m completions (m_used == m_total), on the pre-registered
#'      analysis scale. A pair-cell with status "ok" that misses the
#'      finiteness or m_used == m_total condition is itself a failure (its
#'      number would otherwise be reported as if it were a full pool).
#' Each violation is named in `failures`. The 5% slope criterion is
#' REPORTED, NOT GATED (see 03_acceptance.R header and GATES.md G8): it is
#' aggregated over the counted pair-cells per pair (`pair_status`) and
#' overall (`slope_overall`), never overwritten by the last receipt read.
check_acceptance <- function(receipts, planned, pairs_list) {
  failures <- character(0)
  add_fail <- function(fmt, ...) failures <<- c(failures, sprintf(fmt, ...))
  fin <- function(v) length(v) == 1L && is.numeric(v) && is.finite(v)
  one_chr <- function(v) if (is.character(v) && length(v) == 1L) v else NA_character_

  nm_all <- vapply(seq_len(nrow(planned)), function(i)
    cell_name(planned$dataset[[i]], planned$arm[[i]], planned$seed[[i]]), character(1))
  cell_status <- data.frame(
    name = nm_all, dataset = planned$dataset, arm = planned$arm, seed = planned$seed,
    status = vapply(nm_all, function(nm) one_chr(receipts[[nm]]$status %||% "missing"), character(1)),
    stage = vapply(nm_all, function(nm) one_chr(receipts[[nm]]$stage %||% ""), character(1)),
    code_sha = vapply(nm_all, function(nm) one_chr(receipts[[nm]]$code_sha), character(1)),
    code_sha_source = vapply(nm_all, function(nm) one_chr(receipts[[nm]]$code_sha_source), character(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )
  cell_status$status[is.na(cell_status$status)] <- "unreadable"
  cell_status$stage[is.na(cell_status$stage)] <- ""
  missing_cells <- cell_status[cell_status$status != "ok", , drop = FALSE]
  for (i in seq_len(nrow(missing_cells))) {
    err <- receipts[[missing_cells$name[[i]]]]$error
    add_fail("cell %s: status '%s'%s%s", missing_cells$name[[i]], missing_cells$status[[i]],
             if (nzchar(missing_cells$stage[[i]])) paste0(" at stage '", missing_cells$stage[[i]], "'") else "",
             if (!is.null(err)) paste0(" (", substr(err, 1L, 160L), ")") else "")
  }

  ok_names <- cell_status$name[cell_status$status == "ok"]

  # ---- provenance: one clean, non-NA code SHA across all ok receipts ----
  ok_rows <- cell_status[cell_status$status == "ok", , drop = FALSE]
  for (i in seq_len(nrow(ok_rows))) {
    if (is.na(ok_rows$code_sha[[i]]) || !nzchar(ok_rows$code_sha[[i]])) {
      add_fail("provenance %s: code_sha is NA or empty (source '%s'); run with MI_POST_SHA set or from a clean git checkout",
               ok_rows$name[[i]], ok_rows$code_sha_source[[i]])
    } else if (!(ok_rows$code_sha_source[[i]] %in% sha_sources_ok)) {
      add_fail("provenance %s: code_sha_source '%s' is not MI_POST_SHA or a clean git HEAD",
               ok_rows$name[[i]], ok_rows$code_sha_source[[i]])
    }
  }
  shas <- unique(ok_rows$code_sha[!is.na(ok_rows$code_sha) & nzchar(ok_rows$code_sha)])
  if (length(shas) > 1L) {
    add_fail("provenance: ok receipts come from %d code SHAs (%s); rerun so every cell shares one SHA",
             length(shas), paste(vapply(shas, function(s) sprintf("%s: %s", substr(s, 1L, 12L),
               paste(ok_rows$name[ok_rows$code_sha %in% s], collapse = ", ")), character(1)), collapse = "; "))
  }

  for (nm in ok_names) {
    rc <- receipts[[nm]]
    if (!identical(rc$receipt_schema, receipt_schema_current)) {
      add_fail("receipt %s: schema %s, need %d (written before the 2026-09-24 repair); re-run 01_run.R",
               nm, format(rc$receipt_schema %||% 1L), receipt_schema_current)
      next
    }
    if (is.null(rc$sampler) || !("overrides" %in% names(rc$sampler))) {
      add_fail("sampler %s: no sampler record (receipt$sampler$overrides); cannot confirm pigauto defaults", nm)
    } else if (length(rc$sampler$overrides)) {
      ov <- rc$sampler$overrides
      add_fail("sampler %s: non-default posterior_control (%s); rerun without MI_POST_NITER/MI_POST_BURNIN",
               nm, paste(names(ov), unlist(ov), sep = " = ", collapse = ", "))
    }
    for (leg in c("split", "mondrian")) {
      co <- rc[[paste0(leg, "_conformal")]]
      if (!identical(co$status, "ok")) {
        add_fail("conformal %s: %s-conformal status '%s'%s", nm, leg, co$status %||% "missing",
                 if (!is.null(co$error)) paste0(" (", substr(co$error, 1L, 160L), ")") else "")
      }
    }
    mc <- rc$model_coverage
    if (is.null(mc) || !nrow(mc)) {
      add_fail("coverage %s: no model-coverage table", nm)
      next
    }
    if (!length(rc$kept_traits)) add_fail("coverage %s: no kept traits recorded in the receipt", nm)
    for (tr in union(rc$kept_traits, mc$trait)) {
      r <- mc[mc$trait == tr, , drop = FALSE]
      if (nrow(r) != 1L) {
        add_fail("coverage %s/%s: %d model-coverage rows, need 1", nm, tr, nrow(r))
        next
      }
      if (!(tr %in% rc$kept_traits)) add_fail("coverage %s/%s: coverage row for a trait not in kept_traits", nm, tr)
      if (!isTRUE(r$n_matched == r$n_masked)) {
        add_fail("coverage %s/%s: n_matched = %s != n_masked = %s", nm, tr,
                 format(r$n_matched), format(r$n_masked))
      }
      if (isTRUE(r$n_masked > 0L)) {
        if (!fin(r$coverage)) {
          add_fail("coverage %s/%s: model coverage is %s, need a finite value", nm, tr, format(r$coverage))
        }
        for (leg in c("split", "mondrian")) {
          co <- rc[[paste0(leg, "_conformal")]]
          if (!identical(co$status, "ok")) next  # already named above
          cr <- metric_row(co, tr)
          if (!identical(cr$status, "ok")) {
            add_fail("conformal %s/%s: no %s-conformal row (%s)", nm, tr, leg, cr$status)
            next
          }
          if (!fin(cr$coverage)) {
            add_fail("conformal %s/%s: %s-conformal coverage is %s, need a finite value",
                     nm, tr, leg, format(cr$coverage))
          }
          if (!(isTRUE(cr$n_interval == cr$n_masked) && isTRUE(cr$n_masked == r$n_masked))) {
            add_fail("conformal %s/%s: %s-conformal n_interval = %s, n_masked = %s; model n_masked = %s (need all equal)",
                     nm, tr, leg, format(cr$n_interval), format(cr$n_masked), format(r$n_masked))
          }
        }
      }
    }
    if (!any(mc$n_masked > 0L, na.rm = TRUE)) add_fail("coverage %s: no masked cells in any kept trait", nm)
  }

  pair_rows <- list(); w5_all <- logical(0)
  for (p in pairs_list) {
    lab <- sprintf("%s %s ~ %s", p$dataset, p$response, p$predictor)
    n_att <- 0L; n_ok <- 0L; n_deg <- 0L
    w5 <- logical(0); diffs <- numeric(0); dse <- numeric(0); rels <- numeric(0)
    detail <- character(0)
    for (nm in ok_names) {
      rc <- receipts[[nm]]
      if (!identical(rc$dataset, p$dataset) || is.null(rc$pairs)) next
      for (pr in rc$pairs) {
        if (!(identical(pr$response, p$response) && identical(pr$predictor, p$predictor))) next
        n_att <- n_att + 1L
        rs <- pr$reference$status %||% "missing"; ms <- pr$mi$status %||% "missing"
        for (leg in c("reference", "mi")) {
          lg <- pr[[leg]]
          if (identical(lg$status, "ok") &&
                !(identical(lg$y_log, p$y_log) && identical(lg$x_log, p$x_log))) {
            add_fail("pair %s: cell %s %s leg used y_log = %s, x_log = %s; pre-registered %s / %s",
                     lab, nm, leg, format(lg$y_log), format(lg$x_log), format(p$y_log), format(p$x_log))
          }
        }
        if (!(identical(rs, "ok") && identical(ms, "ok"))) {
          detail <- c(detail, sprintf("%s: ref=%s mi=%s", nm, rs, ms))
          next
        }
        est <- list(ref_slope = pr$reference$slope, ref_se = pr$reference$se,
                    mi_slope = pr$mi$slope, mi_se = pr$mi$se)
        nonfin <- names(est)[!vapply(est, fin, logical(1))]
        m_used <- pr$mi$m_used; m_total <- pr$mi$m_total
        bad <- character(0)
        if (length(nonfin)) bad <- c(bad, paste("non-finite", paste(nonfin, collapse = ", ")))
        if (!(fin(m_used) && fin(m_total) && m_total >= 1 && m_used == m_total)) {
          bad <- c(bad, sprintf("pooled %s of %s completions (%s non-finite on the analysis scale)",
                                format(m_used %||% NA), format(m_total %||% NA),
                                format(pr$mi$n_nonfinite %||% NA)))
        }
        if (length(bad)) {
          n_deg <- n_deg + 1L
          detail <- c(detail, sprintf("%s: ref=ok mi=ok but %s", nm, paste(bad, collapse = "; ")))
          add_fail("pair %s: cell %s has reference and MI status 'ok' but %s; not counted as ok",
                   lab, nm, paste(bad, collapse = "; "))
          next
        }
        n_ok <- n_ok + 1L
        detail <- c(detail, sprintf("%s: ref=ok mi=ok", nm))
        pc <- pair_contrast(pr$reference, pr$mi)
        if (is.finite(pc$diff)) diffs <- c(diffs, pc$diff)
        if (is.finite(pc$diff_ref_se)) dse <- c(dse, pc$diff_ref_se)
        if (is.finite(pc$rel_diff)) {
          rels <- c(rels, pc$rel_diff)
          w5 <- c(w5, abs(pc$rel_diff) <= 0.05)
        }
      }
    }
    if (n_ok == 0L) {
      add_fail("pair %s: no cell with both reference and MI slope status 'ok' (%s)", lab,
               if (length(detail)) paste(detail, collapse = "; ") else "never attempted by an ok cell")
    }
    w5_all <- c(w5_all, w5)
    mx <- function(v) if (length(v)) max(abs(v)) else NA_real_
    pair_rows[[length(pair_rows) + 1L]] <- data.frame(
      dataset = p$dataset, response = p$response, predictor = p$predictor,
      y_log = p$y_log %||% NA, x_log = p$x_log %||% NA,
      n_cells_attempted = n_att, n_cells_ok_slope = n_ok, n_cells_degraded = n_deg,
      has_ok_slope = n_ok > 0L,
      n_within_5pct = sum(w5), n_rel_finite = length(w5),
      all_within_5pct = if (length(w5)) all(w5) else NA,
      max_abs_rel_diff = mx(rels), max_abs_diff = mx(diffs), max_abs_diff_ref_se = mx(dse),
      stringsAsFactors = FALSE
    )
  }
  pair_status <- if (length(pair_rows)) do.call(rbind, pair_rows) else data.frame()

  list(
    complete = length(failures) == 0L,
    failures = failures,
    cell_status = cell_status, missing_cells = missing_cells, pair_status = pair_status,
    slope_overall = list(n_within_5pct = sum(w5_all), n_pair_cells = length(w5_all))
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

#' Synthetic fixture variants for --selftest in 02_summarise.R /
#' 03_acceptance.R, entirely separate from the real planned_cells() /
#' mi_realdata_pairs, so a selftest run can never be mistaken for
#' real-data acceptance. Each variant is c(expected complete, substrings
#' that must EACH start one of the gate's failures; NA for none). The
#' variants after `stale_schema` are the repair-round probes (2026-09-24):
#' each one printed REALDATA_COMPLETE before that repair.
selftest_variants <- function() {
  list(
    pass                 = list(complete = TRUE,  must_name = NA_character_),
    missing_cell         = list(complete = FALSE, must_name = "cell synth-mcar-m2: status 'missing'"),
    pair_error_only      = list(complete = FALSE, must_name = "pair synth t2 ~ t1: no cell with both reference and MI slope status 'ok'"),
    missing_conformal    = list(complete = FALSE, must_name = "conformal synth-mcar-m2: split-conformal status 'error'"),
    conformal_trait      = list(complete = FALSE, must_name = "conformal synth-mcar-m1/t2: no mondrian-conformal row"),
    n_matched_short      = list(complete = FALSE, must_name = "coverage synth-mcar-m1/t1: n_matched = 8 != n_masked = 10"),
    sampler_override     = list(complete = FALSE, must_name = "sampler synth-mcar-m1: non-default posterior_control"),
    wrong_scale          = list(complete = FALSE, must_name = "pair synth t2 ~ t1: cell synth-mcar-m2 reference leg used y_log = TRUE"),
    stale_schema         = list(complete = FALSE, must_name = "receipt synth-mcar-m1: schema 1, need 3"),
    # (a) model coverage NA although n_matched == n_masked
    model_coverage_na    = list(complete = FALSE, must_name = "coverage synth-mcar-m1/t2: model coverage is NA"),
    # (b) conformal coverage NA; conformal n_interval 0
    conformal_coverage_na = list(complete = FALSE, must_name = "conformal synth-mcar-m1/t1: split-conformal coverage is NA"),
    conformal_n_interval = list(complete = FALSE, must_name = "conformal synth-mcar-m2/t2: mondrian-conformal n_interval = 0, n_masked = 10"),
    # (c) both statuses ok but slope and SE NaN, in every cell
    pair_nan             = list(complete = FALSE, must_name = c(
      "pair synth t2 ~ t1: cell synth-mcar-m1 has reference and MI status 'ok' but non-finite mi_slope, mi_se",
      "pair synth t2 ~ t1: no cell with both reference and MI slope status 'ok'")),
    # (d) kept_traits NULL, with a short n_matched hidden behind it
    kept_traits_empty    = list(complete = FALSE, must_name = c(
      "coverage synth-mcar-m1: no kept traits recorded",
      "coverage synth-mcar-m1/t1: n_matched = 3 != n_masked = 10")),
    # (e) no sampler record at all
    sampler_missing      = list(complete = FALSE, must_name = "sampler synth-mcar-m2: no sampler record"),
    # (f) selective pooling: 2 of 20 completions pooled
    partial_pool         = list(complete = FALSE, must_name = "pair synth t2 ~ t1: cell synth-mcar-m1 has reference and MI status 'ok' but pooled 2 of 20 completions (18 non-finite"),
    # provenance: NA SHA, two SHAs, a dirty git HEAD, a crashed rerun
    sha_na               = list(complete = FALSE, must_name = "provenance synth-mcar-m1: code_sha is NA"),
    sha_mixed            = list(complete = FALSE, must_name = "provenance: ok receipts come from 2 code SHAs"),
    sha_dirty            = list(complete = FALSE, must_name = "provenance synth-mcar-m2: code_sha_source 'git HEAD, R/ or DESCRIPTION has uncommitted changes'"),
    running_placeholder  = list(complete = FALSE, must_name = "cell synth-mcar-m2: status 'running' at stage 'started'")
  )
}

#' Build one fixture variant (2 fake cells, 1 fake pair) under tmp_outdir.
#' "pass" is the positive control: every G8 condition holds, while cell 1
#' is deliberately non-converged and cell 2's slope is deliberately
#' outside 5%, so the fixture also shows that neither convergence nor the
#' 5% criterion gates REALDATA_COMPLETE. Every other variant breaks
#' exactly one condition and must fail, naming what failed.
build_selftest_fixture <- function(tmp_outdir, variant = "pass") {
  stopifnot(variant %in% names(selftest_variants()))
  planned <- data.frame(dataset = "synth", arm = "mcar", seed = c(1L, 2L))
  pairs_list <- list(list(dataset = "synth", response = "t2", predictor = "t1",
                           y_log = FALSE, x_log = FALSE, rationale = "synthetic selftest pair"))
  conf <- function(cov, w) {
    list(status = "ok", error = NULL, method = "synthetic",
         metrics = data.frame(trait = c("t1", "t2"), n_masked = c(10L, 10L), n_interval = c(10L, 10L),
                              coverage = cov, mean_width = w, median_width = 0.9 * w,
                              stringsAsFactors = FALSE))
  }
  mi_slopes <- c(2.02, 2.16)  # cell 1: +1% (within 5%); cell 2: +8% (outside 5%)
  for (i in seq_len(nrow(planned))) {
    if (variant == "missing_cell" && i == 2L) next
    nm <- cell_name(planned$dataset[[i]], planned$arm[[i]], planned$seed[[i]])
    dir.create(file.path(tmp_outdir, nm), recursive = TRUE, showWarnings = FALSE)
    if (variant == "running_placeholder" && i == 2L) {
      # what 01_run.R leaves behind when it is killed after its first write
      saveRDS(list(status = "running", stage = "started", dataset = "synth", arm = "mcar", seed = 2L,
                   name = nm, error = "01_run.R started but did not write a final receipt"),
              file.path(tmp_outdir, nm, "mi_posterior.rds"))
      next
    }
    mc <- data.frame(trait = c("t1", "t2"), n_masked = c(10L, 10L), n_matched = c(10L, 10L),
                      coverage = c(0.94, 0.96), mean_width = c(1.2, 1.0), median_width = c(1.1, 0.9),
                      stringsAsFactors = FALSE)
    if (variant == "n_matched_short" && i == 1L) mc$n_matched[[1L]] <- 8L
    if (variant == "kept_traits_empty" && i == 1L) mc$n_matched[[1L]] <- 3L
    if (variant == "model_coverage_na" && i == 1L) mc$coverage[[2L]] <- NA_real_
    sp <- conf(c(0.90, 0.93), c(1.3, 1.0))
    mo <- conf(c(0.95, 0.94), c(1.2, 0.95))
    if (variant == "missing_conformal" && i == 2L) {
      sp <- list(status = "error", error = "synthetic: split.rds unreadable", metrics = NULL)
    }
    if (variant == "conformal_trait" && i == 1L) mo$metrics <- mo$metrics[mo$metrics$trait != "t2", ]
    if (variant == "conformal_coverage_na" && i == 1L) sp$metrics$coverage[[1L]] <- NA_real_
    if (variant == "conformal_n_interval" && i == 2L) mo$metrics$n_interval[[2L]] <- 0L
    y_log <- if (variant == "wrong_scale" && i == 2L) TRUE else FALSE
    ref <- list(status = "ok", n = 50L, slope = 2.0, se = 0.05, lambda = 0.8, y_log = y_log, x_log = FALSE)
    mip <- list(status = "ok", m_used = 20L, m_total = 20L, n_nonfinite = 0L, n = 50L,
                slope = mi_slopes[[i]], se = 0.06, df = 15, fmi = 0.30, y_log = y_log, x_log = FALSE)
    if (variant == "pair_error_only") {
      mip <- list(status = "error", error = "synthetic: only 1/20 phylolm fits succeeded")
    }
    if (variant == "pair_nan") { mip$slope <- NaN; mip$se <- NaN }
    if (variant == "partial_pool" && i == 1L) { mip$m_used <- 2L; mip$n_nonfinite <- 18L }
    pr <- list(c(list(dataset = "synth", response = "t2", predictor = "t1",
                      rationale = "synthetic selftest pair", reference = ref, mi = mip),
                 pair_contrast(ref, mip)))
    receipt <- list(
      status = "ok", receipt_schema = receipt_schema_current,
      dataset = planned$dataset[[i]], arm = planned$arm[[i]], seed = planned$seed[[i]],
      name = nm, code_sha = "0123456789abcdef0123456789abcdef01234567", code_sha_source = "MI_POST_SHA",
      n_species = 50L, kept_traits = c("t1", "t2"),
      dropped_traits = data.frame(trait = "t3", class = "factor", keep = FALSE,
                                   reason = "dropped: synthetic factor trait", stringsAsFactors = FALSE),
      wall_time_s = 12.3, model_coverage = mc,
      split_conformal = sp, mondrian_conformal = mo, conformal_source = "selftest",
      sampler = list(overrides = if (variant == "sampler_override" && i == 1L) list(n_iter = 100L) else list()),
      diagnostics = structure(
        data.frame(parameter = c("Sigma_P[1,1]", "Sigma_E[1,1]", "lambda[1]"),
                    rhat = if (i == 1L) c(1.09, 1.02, 1.01) else c(1.01, 1.00, 1.00),
                    ess_bulk = if (i == 1L) c(250, 900, 850) else c(900, 1200, 1100),
                    stringsAsFactors = FALSE),
        converged = i != 1L  # cell 1 deliberately non-converged: exercises NONCONVERGED, must not gate
      ),
      pairs = pr
    )
    if (variant == "stale_schema" && i == 1L) receipt$receipt_schema <- NULL
    if (variant == "kept_traits_empty" && i == 1L) receipt$kept_traits <- NULL
    if (variant == "sampler_missing" && i == 2L) receipt$sampler <- NULL
    if (variant == "sha_na" && i == 1L) { receipt$code_sha <- NA_character_; receipt$code_sha_source <- "unknown" }
    if (variant == "sha_mixed" && i == 2L) receipt$code_sha <- "fedcba9876543210fedcba9876543210fedcba98"
    if (variant == "sha_dirty" && i == 2L) {
      receipt$code_sha_source <- "git HEAD, R/ or DESCRIPTION has uncommitted changes"
    }
    receipt$convergence <- summarize_convergence(receipt$diagnostics)
    saveRDS(receipt, file.path(tmp_outdir, nm, "mi_posterior.rds"))
  }
  list(planned = planned, pairs_list = pairs_list)
}
