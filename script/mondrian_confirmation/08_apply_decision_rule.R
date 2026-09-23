#!/usr/bin/env Rscript
# Usage:
#   Rscript 08_apply_decision_rule.R script/mondrian_confirmation/returned NEWS.md
#   Rscript 08_apply_decision_rule.R --selftest
#
# Implements the pre-registered decision rule verbatim from
# docs/dev-log/mondrian-realdata/00-preregistration.md "Decision rule for the
# default". Reads per-dataset summary .rds files (each the output of
# 02_summarise_masked_confirmation.R, i.e. list(stratum = data.frame(...))
# with columns method/trait/stratum/n_test/coverage/median_half_width/fallback,
# arm carried in the file name or an `arm` column) from `returned/`, applies
# the three conditions with Holm adjustment across datasets for condition 2's
# one-sided test, and prints RULE_VERDICT=KEEP_SPLIT or FLIP_MONDRIAN. It then
# checks whether NEWS.md already records that verdict in a line mentioning
# "mondrian" and prints NEWS_MATCHES=TRUE/FALSE.
#
# Traits whose mondrian fit fell back count as NO EVIDENCE, never a pass.

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Core: build a long paired table (one row per dataset x arm x trait x
# stratum) from a list of per-dataset "returned" summary objects, each
# carrying $stratum (per-method) and $paired (mondrian-split diffs), as
# produced by 02_summarise_masked_confirmation.R, plus an `arm` and
# `dataset` tag.
# ---------------------------------------------------------------------------
build_long_table <- function(returned) {
  rows <- lapply(names(returned), function(nm) {
    r <- returned[[nm]]
    p <- r$paired
    if (is.null(p) || !nrow(p)) return(NULL)
    s_mond <- r$stratum[r$stratum$method == "mondrian", , drop = FALSE]
    p$fallback <- s_mond$fallback[match(p$trait, s_mond$trait)]
    p$dataset <- r$dataset %||% nm
    p$arm <- r$arm %||% "unknown"
    p
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

one_sided_binom_z <- function(cov_mond, n_mond, cov_split, n_split, margin = -0.02) {
  # H0: p_mond - p_split <= margin  vs  H1: p_mond - p_split > margin.
  # Pooled-variance z on the observed gap minus the non-inferiority margin.
  se <- sqrt(cov_mond * (1 - cov_mond) / max(n_mond, 1) +
             cov_split * (1 - cov_split) / max(n_split, 1))
  if (!is.finite(se) || se <= 0) return(NA_real_)
  z <- ((cov_mond - cov_split) - margin) / se
  1 - stats::pnorm(z)
}

# ---------------------------------------------------------------------------
# apply_decision_rule(long): implements the three pre-registered conditions.
# `long` has one row per dataset/arm/trait/stratum with columns:
#   dataset, arm, trait, stratum, coverage_mondrian, coverage_split,
#   width_ratio, n_test_mondrian, n_test_split, fallback
# Returns list(verdict, detail).
# ---------------------------------------------------------------------------
apply_decision_rule <- function(long) {
  active <- long[!isTRUE_vec(long$fallback), , drop = FALSE]
  if (!nrow(active)) return(list(verdict = "KEEP_SPLIT", detail = "no activated trait"))

  # Condition 1: structured arm, far stratum.
  far_struct <- active[active$arm == "structured" & active$stratum == "far", , drop = FALSE]
  cond1 <- if (!nrow(far_struct)) FALSE else {
    gain <- stats::median(far_struct$coverage_gain)
    gain >= 0 && all(far_struct$coverage_mondrian >= 0.90)
  }

  # Condition 2: structured + MCAR arms, near stratum, non-inferiority,
  # Holm-adjusted across datasets on the pooled paired cells per dataset.
  near <- active[active$stratum == "near", , drop = FALSE]
  cond2 <- TRUE
  if (nrow(near)) {
    per_dataset <- split(near, near$dataset)
    pvals <- vapply(per_dataset, function(g) {
      n_mond <- sum(g$n_test_mondrian)
      n_split <- sum(g$n_test_split)
      cov_mond <- stats::weighted.mean(g$coverage_mondrian, g$n_test_mondrian)
      cov_split <- stats::weighted.mean(g$coverage_split, g$n_test_split)
      one_sided_binom_z(cov_mond, n_mond, cov_split, n_split, margin = -0.02)
    }, numeric(1))
    pvals <- pvals[is.finite(pvals)]
    if (length(pvals)) {
      adj <- stats::p.adjust(pvals, method = "holm")
      cond2 <- all(adj <= 0.05)
    }
  }

  # Condition 3: near stratum, paired half-width ratio.
  cond3 <- if (!nrow(near)) FALSE else stats::median(near$width_ratio) <= 1.10

  pass <- isTRUE(cond1) && isTRUE(cond2) && isTRUE(cond3)
  list(
    verdict = if (pass) "FLIP_MONDRIAN" else "KEEP_SPLIT",
    detail = list(cond1_far_coverage = cond1, cond2_near_noninferior = cond2,
                  cond3_near_width = cond3)
  )
}

isTRUE_vec <- function(x) vapply(x, isTRUE, logical(1))

news_matches <- function(news_path, verdict) {
  if (!file.exists(news_path)) return(FALSE)
  lines <- readLines(news_path, warn = FALSE)
  hit <- grepl("mondrian", lines, ignore.case = TRUE) & grepl(verdict, lines, fixed = TRUE)
  any(hit)
}

# ---------------------------------------------------------------------------
# Self-test: synthetic fixture exercising both verdicts, run in tempdir.
# ---------------------------------------------------------------------------
selftest <- function() {
  make_row <- function(dataset, arm, trait, stratum, cov_m, cov_s, wr, n = 100, fallback = FALSE) {
    data.frame(dataset = dataset, arm = arm, trait = trait, stratum = stratum,
               coverage_mondrian = cov_m, coverage_split = cov_s,
               coverage_gain = cov_m - cov_s, width_ratio = wr,
               n_test_mondrian = n, n_test_split = n, fallback = fallback,
               stringsAsFactors = FALSE)
  }

  keep_split <- rbind(
    make_row("d1", "structured", "t1", "far", 0.80, 0.85, 1.0),
    make_row("d1", "structured", "t1", "near", 0.93, 0.95, 1.05),
    make_row("d1", "mcar", "t1", "near", 0.93, 0.95, 1.05)
  )
  verdict_keep <- apply_decision_rule(keep_split)
  stopifnot(identical(verdict_keep$verdict, "KEEP_SPLIT"))

  flip <- rbind(
    make_row("d1", "structured", "t1", "far", 0.97, 0.85, 1.3, n = 2000),
    make_row("d1", "structured", "t1", "near", 0.96, 0.95, 1.02, n = 2000),
    make_row("d1", "mcar", "t1", "near", 0.96, 0.95, 1.02, n = 2000),
    make_row("d2", "structured", "t1", "far", 0.96, 0.85, 1.3, n = 2000),
    make_row("d2", "structured", "t1", "near", 0.96, 0.95, 1.02, n = 2000),
    make_row("d2", "mcar", "t1", "near", 0.96, 0.95, 1.02, n = 2000)
  )
  verdict_flip <- apply_decision_rule(flip)
  stopifnot(identical(verdict_flip$verdict, "FLIP_MONDRIAN"))

  tmp <- tempfile(fileext = ".md")
  writeLines(c("# NEWS", "* mondrian default: FLIP_MONDRIAN"), tmp)
  stopifnot(isTRUE(news_matches(tmp, "FLIP_MONDRIAN")))
  stopifnot(!isTRUE(news_matches(tmp, "KEEP_SPLIT")))
  unlink(tmp)

  cat("SELFTEST_OK\n")
  invisible(TRUE)
}

if (identical(Sys.getenv("PIGAUTO_08_SOURCE_ONLY"), "")) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 1L && identical(args[[1L]], "--selftest")) {
    selftest()
  } else if (length(args) == 2L) {
    returned_dir <- args[[1L]]
    news_path <- args[[2L]]
    files <- list.files(returned_dir, pattern = "\\.rds$", full.names = TRUE)
    if (!length(files)) stop("no .rds files found in ", returned_dir, call. = FALSE)
    returned <- stats::setNames(lapply(files, readRDS), tools::file_path_sans_ext(basename(files)))
    for (nm in names(returned)) {
      if (is.null(returned[[nm]]$dataset)) returned[[nm]]$dataset <- nm
    }
    long <- build_long_table(returned)
    if (is.null(long) || !nrow(long)) stop("no paired evidence found across returned receipts", call. = FALSE)
    res <- apply_decision_rule(long)
    cat(sprintf("RULE_VERDICT=%s\n", res$verdict))
    cat(sprintf("NEWS_MATCHES=%s\n", news_matches(news_path, res$verdict)))
  } else {
    stop("expected: returned_dir NEWS.md   OR   --selftest", call. = FALSE)
  }
}
