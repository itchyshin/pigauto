#!/usr/bin/env Rscript
# Usage: Rscript 12_build_results_doc.R [results_md] [results_csv]
# Env:   RESULTS_ROOT (default script/mondrian_confirmation/returned) -- one
#        subdirectory per confirmed cell, named <dataset>-<arm>-m<seed>, each
#        holding mask_receipt.rds, mondrian.rds, split.rds as written by
#        01_run_masked_confirmation.R.
#
# For every cell in the pre-registered registry (REGISTRY below) that has all
# three receipts, this drives 02_summarise_masked_confirmation.R (as an
# unmodified child Rscript process, the same way 09_gate_smoke.R drives
# 01/02) to get its per-trait/stratum table, pools across masks (seeds)
# within dataset x arm x trait x stratum for Table 1, and reduces further to
# dataset-level statistics with a between-mask SD for Table 2. Cells that are
# pre-registered but absent, or whose receipts are incomplete, are listed by
# name in results.md rather than silently skipped.
#
# docs/dev-log/mondrian-realdata/00-preregistration.md, Amendment 2
# (2026-09-23): decision-rule condition 1 (structured arm, far stratum) is
# evaluated only on traits with real missing fraction >= 5%. This script
# computes that fraction per trait from mask_receipt.rds$truth (is.na before
# masking) and carries it, plus a cond1_eligible flag, into results_table.csv
# so 08_apply_decision_rule.R applies the same filter without recomputing it.
#
# Run directly, this writes results_md (default
# docs/dev-log/mondrian-realdata/results.md) and results_csv (default
# results_table.csv next to it).
#
# Sourced by 07_check_results_doc.R (with PIGAUTO_12_SOURCE_ONLY set to a
# non-empty value first) to reuse build_results_table() / write_results_csv()
# for its ROWS_MATCH check, so the doc and the CSV can never silently drift
# from what the receipts say.

`%||%` <- function(a, b) if (is.null(a)) b else a

SCRIPT_DIR <- "script/mondrian_confirmation"
if (!file.exists(file.path(SCRIPT_DIR, "02_summarise_masked_confirmation.R"))) {
  stop("run this script from the pigauto repo root", call. = FALSE)
}

ALPHA <- 0.05
COND1_MIN_REAL_MISSING <- 0.05  # Amendment 2

# ---------------------------------------------------------------------------
# Registry: every (dataset, arm, seed) cell the pre-registration and its
# amendments call for (docs/dev-log/mondrian-realdata/00-preregistration.md):
#   - PanTHERIA: mcar and structured, seeds 20260818/19/20.
#   - AVONET: mcar only, seeds 20260818/19/20 (Amendment 1 -- too little real
#     missingness in the bundled subset for a structured propensity).
#   - FishBase: mcar and structured, seed 20260818 only ("1 for FishBase ...
#     if Shinichi approves the GPU campaign").
# ---------------------------------------------------------------------------
SEEDS_3 <- c(20260818L, 20260819L, 20260820L)
REGISTRY <- rbind(
  expand.grid(dataset = "pantheria", arm = c("mcar", "structured"), seed = SEEDS_3, stringsAsFactors = FALSE),
  expand.grid(dataset = "avonet", arm = "mcar", seed = SEEDS_3, stringsAsFactors = FALSE),
  expand.grid(dataset = "fishbase", arm = c("mcar", "structured"), seed = 20260818L, stringsAsFactors = FALSE)
)
REGISTRY <- REGISTRY[order(REGISTRY$dataset, REGISTRY$arm, REGISTRY$seed), ]
rownames(REGISTRY) <- NULL

cell_dir_name <- function(dataset, arm, seed) sprintf("%s-%s-m%d", dataset, arm, as.integer(seed))

# Same one-sided non-inferiority z-test as 08_apply_decision_rule.R's
# one_sided_binom_z(), kept in sync by hand (both are small and stable);
# used here only for Table 2's descriptive near-stratum p-value, never to
# gate anything -- 08 owns the verdict.
one_sided_binom_z <- function(cov_mond, n_mond, cov_split, n_split, margin = -0.02) {
  se <- sqrt(cov_mond * (1 - cov_mond) / max(n_mond, 1) +
             cov_split * (1 - cov_split) / max(n_split, 1))
  if (!is.finite(se) || se <= 0) return(NA_real_)
  z <- ((cov_mond - cov_split) - margin) / se
  1 - stats::pnorm(z)
}

# Same analytic MCSE as 02_summarise_masked_confirmation.R's per_stratum().
mcse_formula <- function(cov, n_test, n_s) {
  sqrt(cov * (1 - cov) / n_test + ALPHA * (1 - ALPHA) / (n_s + 2))
}

isTRUE_vec <- function(x) !is.na(x) & x

# ---------------------------------------------------------------------------
# summarise_one_cell(): drives 02_summarise_masked_confirmation.R against one
# receipt directory and returns its long per-trait/stratum table (mondrian
# and split merged into one row per trait x stratum), or a $status string
# explaining why the cell contributes nothing.
# ---------------------------------------------------------------------------
summarise_one_cell <- function(cell_dir) {
  receipts <- file.path(cell_dir, c("mask_receipt.rds", "mondrian.rds", "split.rds"))
  if (!dir.exists(cell_dir) || !all(file.exists(receipts))) {
    return(list(status = "not run (pre-registered/amended)"))
  }
  mask_receipt <- tryCatch(readRDS(receipts[[1L]]), error = function(e) NULL)
  if (is.null(mask_receipt) || is.null(mask_receipt$truth)) {
    return(list(status = "incomplete (unreadable mask_receipt.rds)"))
  }
  real_missing_frac <- vapply(mask_receipt$truth, function(x) mean(is.na(x)), numeric(1))

  summary_file <- tempfile(fileext = ".rds")
  on.exit(unlink(summary_file), add = TRUE)
  status <- system2("Rscript", shQuote(c(
    file.path(SCRIPT_DIR, "02_summarise_masked_confirmation.R"), cell_dir, summary_file
  )), stdout = FALSE, stderr = FALSE)
  if (!identical(status, 0L) || !file.exists(summary_file)) {
    return(list(status = "incomplete (02_summarise_masked_confirmation.R failed on this cell)"))
  }
  summ <- readRDS(summary_file)
  if (is.null(summ$stratum) || !nrow(summ$stratum)) {
    return(list(status = "incomplete (no stratum rows)"))
  }
  sp <- summ$stratum[summ$stratum$method == "split", , drop = FALSE]
  mo <- summ$stratum[summ$stratum$method == "mondrian", , drop = FALSE]
  m <- merge(mo, sp, by = c("trait", "stratum"), suffixes = c("_mondrian", "_split"))
  if (!nrow(m)) return(list(status = "incomplete (no matched split/mondrian rows)"))
  m$coverage_gain <- m$coverage_mondrian - m$coverage_split
  m$width_ratio <- m$median_half_width_mondrian / m$median_half_width_split
  m$fallback <- m$fallback_mondrian
  m$n_val <- m$n_val_mondrian
  m$n_near <- m$n_near_mondrian
  m$n_far <- m$n_far_mondrian
  m$real_missing_frac <- real_missing_frac[m$trait]
  list(status = "ok", long = m)
}

# ---------------------------------------------------------------------------
# build_raw_table(): walks REGISTRY, returns $raw (one row per dataset x arm
# x seed x trait x stratum, the finest grain) and $missing (character
# vector of "<cell>: <reason>" for every registered cell with no evidence).
# ---------------------------------------------------------------------------
build_raw_table <- function(results_root) {
  rows <- vector("list", nrow(REGISTRY))
  missing <- character(0)
  for (i in seq_len(nrow(REGISTRY))) {
    ds <- REGISTRY$dataset[i]; arm <- REGISTRY$arm[i]; seed <- REGISTRY$seed[i]
    cell_dir <- file.path(results_root, cell_dir_name(ds, arm, seed))
    res <- summarise_one_cell(cell_dir)
    if (identical(res$status, "ok")) {
      d <- res$long
      d$dataset <- ds; d$arm <- arm; d$seed <- seed
      rows[[i]] <- d
    } else {
      missing <- c(missing, sprintf("%s: %s", cell_dir_name(ds, arm, seed), res$status))
    }
  }
  list(raw = do.call(rbind, Filter(Negate(is.null), rows)), missing = missing)
}

TABLE1_COLS <- c("dataset", "arm", "trait", "stratum", "n_masks",
  "coverage_split", "coverage_mondrian", "coverage_gain", "mcse", "width_ratio",
  "winkler_split", "winkler_mondrian", "n_test_split", "n_test_mondrian",
  "n_val", "n_near", "n_far", "fallback", "real_missing_frac", "cond1_eligible")

empty_table1 <- function() {
  df <- as.data.frame(setNames(replicate(length(TABLE1_COLS), character(0), simplify = FALSE), TABLE1_COLS),
                       stringsAsFactors = FALSE)
  df
}

# ---------------------------------------------------------------------------
# pool_table1(): pools $raw across masks (seeds), one row per dataset x arm x
# trait x stratum. Coverage and Winkler pool exactly (n_test-weighted sums,
# since each per-mask value is itself a mean over that mask's cells); median
# half-width (and so width_ratio) pools as an n_test-weighted average of the
# per-mask medians, which is not a literal pooled median. mcse is the paired-
# difference MCSE, sqrt(mcse_mondrian^2 + mcse_split^2), each term the
# per_stratum() formula from 02 applied to the pooled coverage/n.
# ---------------------------------------------------------------------------
pool_table1 <- function(raw) {
  if (is.null(raw) || !nrow(raw)) return(empty_table1())
  key <- interaction(raw$dataset, raw$arm, raw$trait, raw$stratum, drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(raw)), key)
  out <- lapply(groups, function(idx) {
    g <- raw[idx, , drop = FALSE]
    n_test_mondrian <- sum(g$n_test_mondrian)
    n_test_split <- sum(g$n_test_split)
    cov_mondrian <- sum(g$coverage_mondrian * g$n_test_mondrian) / n_test_mondrian
    cov_split <- sum(g$coverage_split * g$n_test_split) / n_test_split
    hw_mondrian <- sum(g$median_half_width_mondrian * g$n_test_mondrian) / n_test_mondrian
    hw_split <- sum(g$median_half_width_split * g$n_test_split) / n_test_split
    wnk_mondrian <- sum(g$winkler_mondrian * g$n_test_mondrian) / n_test_mondrian
    wnk_split <- sum(g$winkler_split * g$n_test_split) / n_test_split
    n_val <- sum(g$n_val); n_near <- sum(g$n_near); n_far <- sum(g$n_far)
    n_s <- switch(g$stratum[1L], near = n_near, far = n_far, n_val)
    mcse_mondrian <- mcse_formula(cov_mondrian, n_test_mondrian, n_s)
    mcse_split <- mcse_formula(cov_split, n_test_split, n_s)
    real_missing_frac <- mean(g$real_missing_frac, na.rm = TRUE)
    data.frame(
      dataset = g$dataset[1L], arm = g$arm[1L], trait = g$trait[1L], stratum = g$stratum[1L],
      n_masks = length(unique(g$seed)),
      coverage_split = cov_split, coverage_mondrian = cov_mondrian,
      coverage_gain = cov_mondrian - cov_split,
      mcse = sqrt(mcse_mondrian^2 + mcse_split^2),
      width_ratio = hw_mondrian / hw_split,
      winkler_split = wnk_split, winkler_mondrian = wnk_mondrian,
      n_test_split = n_test_split, n_test_mondrian = n_test_mondrian,
      n_val = n_val, n_near = n_near, n_far = n_far,
      fallback = any(g$fallback),
      real_missing_frac = real_missing_frac,
      cond1_eligible = is.finite(real_missing_frac) && real_missing_frac >= COND1_MIN_REAL_MISSING,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$dataset, out$arm, out$trait, out$stratum), , drop = FALSE]
}

per_mask_stat <- function(sub, stat_fn) {
  if (is.null(sub) || !nrow(sub)) return(numeric(0))
  vals <- vapply(split(sub, sub$seed), stat_fn, numeric(1))
  vals[is.finite(vals)]
}
sd_or_na <- function(x) if (length(x) >= 2L) stats::sd(x) else NA_real_

# ---------------------------------------------------------------------------
# build_table2(): one row per dataset. median_far_gain_structured and
# min_mondrian_far_cov use only cond1_eligible traits (Amendment 2);
# ineligible traits present in the structured far stratum are named in
# ineligible_far_traits. Each dataset-level point estimate is paired with a
# between-mask SD computed from the same statistic evaluated on each mask
# (seed) alone (docs/dev-log/mondrian-realdata/00-preregistration.md:
# "The between-mask SD is reported beside it (2 df)").
# ---------------------------------------------------------------------------
build_table2 <- function(raw, table1) {
  datasets <- sort(unique(c(raw$dataset, table1$dataset)))
  rows <- lapply(datasets, function(ds) {
    t1_far <- table1[table1$dataset == ds & table1$arm == "structured" & table1$stratum == "far", , drop = FALSE]
    elig <- t1_far[isTRUE_vec(t1_far$cond1_eligible), , drop = FALSE]
    inelig_names <- t1_far$trait[!isTRUE_vec(t1_far$cond1_eligible)]

    far_gain <- if (nrow(elig)) stats::median(elig$coverage_gain) else NA_real_
    min_far_cov <- if (nrow(elig)) min(elig$coverage_mondrian) else NA_real_

    raw_far_elig <- raw[raw$dataset == ds & raw$arm == "structured" & raw$stratum == "far" &
                          raw$trait %in% elig$trait, , drop = FALSE]
    far_gain_by_mask <- per_mask_stat(raw_far_elig, function(g) stats::median(g$coverage_gain))
    min_cov_by_mask <- per_mask_stat(raw_far_elig, function(g) min(g$coverage_mondrian))

    t1_near <- table1[table1$dataset == ds & table1$stratum == "near", , drop = FALSE]
    n_mond <- sum(t1_near$n_test_mondrian); n_split <- sum(t1_near$n_test_split)
    cov_mond <- if (n_mond > 0) sum(t1_near$coverage_mondrian * t1_near$n_test_mondrian) / n_mond else NA_real_
    cov_split <- if (n_split > 0) sum(t1_near$coverage_split * t1_near$n_test_split) / n_split else NA_real_
    near_p <- if (isTRUE(n_mond > 0) && isTRUE(n_split > 0)) one_sided_binom_z(cov_mond, n_mond, cov_split, n_split) else NA_real_
    near_width <- if (nrow(t1_near)) stats::median(t1_near$width_ratio) else NA_real_

    raw_near <- raw[raw$dataset == ds & raw$stratum == "near", , drop = FALSE]
    near_p_by_mask <- per_mask_stat(raw_near, function(g) {
      nm <- sum(g$n_test_mondrian); ns <- sum(g$n_test_split)
      if (nm <= 0 || ns <= 0) return(NA_real_)
      cm <- sum(g$coverage_mondrian * g$n_test_mondrian) / nm
      cs <- sum(g$coverage_split * g$n_test_split) / ns
      one_sided_binom_z(cm, nm, cs, ns)
    })
    near_width_by_mask <- per_mask_stat(raw_near, function(g) stats::median(g$width_ratio))

    data.frame(
      dataset = ds,
      median_far_gain_structured = far_gain,
      far_gain_between_mask_sd = sd_or_na(far_gain_by_mask),
      far_gain_n_masks = length(far_gain_by_mask),
      min_mondrian_far_cov = min_far_cov,
      far_mincov_between_mask_sd = sd_or_na(min_cov_by_mask),
      near_noninferiority_p = near_p,
      near_p_between_mask_sd = sd_or_na(near_p_by_mask),
      near_width_ratio = near_width,
      near_width_between_mask_sd = sd_or_na(near_width_by_mask),
      ineligible_far_traits = paste(inelig_names, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

build_results_table <- function(results_root = Sys.getenv("RESULTS_ROOT", "script/mondrian_confirmation/returned")) {
  br <- build_raw_table(results_root)
  table1 <- pool_table1(br$raw)
  table2 <- build_table2(br$raw, table1)
  list(raw = br$raw, missing = br$missing, table1 = table1, table2 = table2, results_root = results_root)
}

write_results_csv <- function(table1, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(table1, path, row.names = FALSE)
}

fmt_num <- function(x, digits = 4) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
fmt_int <- function(x) ifelse(is.na(x), "NA", as.character(as.integer(x)))
fmt_bool <- function(x) ifelse(is.na(x), "NA", ifelse(x, "TRUE", "FALSE"))

md_table <- function(df, cols, headers) {
  lines <- c(paste0("| ", paste(headers, collapse = " | "), " |"),
             paste0("|", paste(rep("---", length(headers)), collapse = "|"), "|"))
  if (nrow(df)) {
    for (i in seq_len(nrow(df))) {
      vals <- vapply(cols, function(cn) as.character(df[[cn]][i]), character(1))
      lines <- c(lines, paste0("| ", paste(vals, collapse = " | "), " |"))
    }
  }
  lines
}

# ---------------------------------------------------------------------------
# render_results_md(): plain-text markdown, no verdict prose -- 08 owns the
# verdict. Every number that appears here also appears in results_csv.
# ---------------------------------------------------------------------------
render_results_md <- function(built, results_root, prereg_rel = "00-preregistration.md") {
  sha <- tryCatch({
    s <- suppressWarnings(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE))
    if (length(s) == 1L && nzchar(s)) s else NA_character_
  }, error = function(e) NA_character_)

  t1 <- built$table1
  t1p <- t1
  if (nrow(t1p)) {
    t1p$coverage_split <- fmt_num(t1p$coverage_split)
    t1p$coverage_mondrian <- fmt_num(t1p$coverage_mondrian)
    t1p$coverage_gain <- fmt_num(t1p$coverage_gain)
    t1p$mcse <- fmt_num(t1p$mcse)
    t1p$width_ratio <- fmt_num(t1p$width_ratio)
    t1p$winkler_split <- fmt_num(t1p$winkler_split)
    t1p$winkler_mondrian <- fmt_num(t1p$winkler_mondrian)
    t1p$n_test_split <- fmt_int(t1p$n_test_split)
    t1p$n_test_mondrian <- fmt_int(t1p$n_test_mondrian)
    t1p$n_val <- fmt_int(t1p$n_val)
    t1p$n_near <- fmt_int(t1p$n_near)
    t1p$n_far <- fmt_int(t1p$n_far)
    t1p$real_missing_frac <- fmt_num(t1p$real_missing_frac, 3)
    t1p$fallback <- fmt_bool(t1p$fallback)
    t1p$cond1_eligible <- fmt_bool(t1p$cond1_eligible)
  }

  t2 <- built$table2
  t2p <- t2
  if (nrow(t2p)) {
    for (cn in c("median_far_gain_structured", "far_gain_between_mask_sd", "min_mondrian_far_cov",
                 "far_mincov_between_mask_sd", "near_noninferiority_p", "near_p_between_mask_sd",
                 "near_width_ratio", "near_width_between_mask_sd")) {
      t2p[[cn]] <- fmt_num(t2p[[cn]])
    }
    t2p$far_gain_n_masks <- fmt_int(t2p$far_gain_n_masks)
    t2p$ineligible_far_traits[!nzchar(t2p$ineligible_far_traits)] <- "none"
  }

  lines <- c(
    "# Mondrian real-data confirmation: results",
    "",
    sprintf("Pre-registration: [%s](%s), Amendments 1 and 2.", prereg_rel, prereg_rel),
    sprintf("Source SHA: %s", if (is.na(sha)) "not recorded" else sha),
    sprintf("Generated: %s from RESULTS_ROOT=%s",
            format(Sys.time(), tz = "UTC", usetz = TRUE), results_root),
    "",
    "Method notes. Table 1 pools across masks (seeds) within each dataset x arm x",
    "trait x stratum. Coverage and Winkler score pool exactly, as n_test-weighted",
    "sums. Median half-width, and so width_ratio, pools as an n_test-weighted average",
    "of the per-mask medians, not a literal pooled median. mcse is the paired-",
    "difference MCSE, sqrt(mcse_mondrian^2 + mcse_split^2), each term the analytic",
    "formula from 02_summarise_masked_confirmation.R applied to the pooled",
    "coverage and n. Amendment 2: cond1_eligible marks traits with real_missing_frac",
    "at least 0.05; Table 2's structured-arm far-stratum aggregates use eligible",
    "traits only and name the traits excluded. Between-mask SD uses up to 3",
    "masks (2 df); NA when fewer than 2 masks contributed.",
    "",
    "Decision-rule script. 08_apply_decision_rule.R was rewritten at 5709a7f (after",
    "all results were read) to evaluate conditions 1 and 3 per dataset, as the",
    "pre-registration says, to fail closed on condition 2 with no evidence, and to",
    "gate on mask completeness. The verdict is KEEP_SPLIT under both the earlier",
    "pooled version and the current one: conditions 1 and 3 pass under both",
    "readings and condition 2 fails under both. Condition 1 also passes with and",
    "without Amendment 2 (FishBase median far gain 0.0164 either way; minimum",
    "Mondrian far coverage 0.923 over five traits without it, 0.963 over three",
    "with it; PanTHERIA 0.0110 and 0.922 unchanged).",
    "",
    "Condition-2 statistic. The one-sided non-inferiority test (one_sided_binom_z,",
    "fixed at 7af133a before any receipt) uses an unpaired pooled-binomial SE on",
    "n_test-weighted coverage per dataset and omits the calibration term of the",
    "pre-registered MCSE; both make it conservative for demonstrating",
    "non-inferiority. FishBase, pre-registered as descriptive, is nevertheless in",
    "the Holm family; the verdict is unchanged with it removed, because AVONET",
    "fails alone.",
    ""
  )

  if (length(built$missing)) {
    lines <- c(lines, "## Missing or incomplete cells", "")
    lines <- c(lines, paste0("- ", built$missing), "")
  }

  lines <- c(lines, "## Table 1: per dataset x arm x trait x stratum", "")
  headers1 <- c("dataset", "arm", "trait", "stratum", "n_masks", "split_cov", "mondrian_cov",
                "paired_diff", "mcse", "width_ratio", "winkler_split", "winkler_mondrian",
                "n_test_split", "n_test_mondrian", "n_val", "n_near", "n_far", "fallback",
                "real_missing_frac", "cond1_eligible")
  lines <- c(lines, md_table(t1p, TABLE1_COLS, headers1), "")

  lines <- c(lines, "## Table 2: dataset-level", "")
  headers2 <- c("dataset", "median_far_gain_structured", "far_gain_between_mask_sd", "far_gain_n_masks",
                "min_mondrian_far_cov", "far_mincov_between_mask_sd", "near_noninferiority_p",
                "near_p_between_mask_sd", "near_width_ratio", "near_width_between_mask_sd",
                "ineligible_far_traits")
  lines <- c(lines, md_table(t2p, headers2, headers2), "")

  lines
}

if (identical(Sys.getenv("PIGAUTO_12_SOURCE_ONLY"), "")) {
  args <- commandArgs(trailingOnly = TRUE)
  results_md <- if (length(args) >= 1L) args[[1L]] else "docs/dev-log/mondrian-realdata/results.md"
  results_csv <- if (length(args) >= 2L) args[[2L]] else file.path(dirname(results_md), "results_table.csv")
  results_root <- Sys.getenv("RESULTS_ROOT", "script/mondrian_confirmation/returned")

  built <- build_results_table(results_root)
  dir.create(dirname(results_md), recursive = TRUE, showWarnings = FALSE)
  writeLines(render_results_md(built, results_root), results_md)
  write_results_csv(built$table1, results_csv)
  cat(sprintf("Wrote %s (%d Table 1 rows, %d Table 2 rows, %d missing/incomplete cells)\n",
              results_md, nrow(built$table1), nrow(built$table2), length(built$missing)))
  cat(sprintf("Wrote %s\n", results_csv))
}
