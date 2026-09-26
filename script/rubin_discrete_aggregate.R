# script/rubin_discrete_aggregate.R
#
# Aggregates the discrete-trait re-run of the freq-vs-BACE Rubin study (rubin_cell.R --save_imp --discrete, on the
# continuous campaign's datasets) into the tables behind the "Discrete traits" section of
# script/rubin_study/study.qmd. Run from the worktree root:
#   Rscript script/rubin_discrete_aggregate.R [POOL] [TRUTH_CSV] [OUT_DIR] [CONT_POOL]
#     POOL       discrete pool, POOL/{bace,freq}/<host>/*.rds (default ~/pigauto_rubin_disc_pool)
#     TRUTH_CSV  Monte Carlo truth of the c1 ~ bin slope by n, lambda, rho (script/rubin_discrete_truth.R;
#                default script/rubin_study/data/discrete/truth_slope_c1_bin.csv)
#     OUT_DIR    default script/rubin_study/data/discrete
#     CONT_POOL  continuous pool, for the reproduction check (default ~/pigauto_rubin_pool)
#
# Unit of replication: one simulated dataset (the tag: cell + seed). BACE and frequentist files with the same tag
# hold the same dataset (checked here: identical mask and discrete traits, c1 within 1e-8), so arms pair by tag.
# The same fit found on two hosts is kept once per set x tag, as build_data.R does. Monte Carlo SEs are over
# datasets (sd / sqrt(count)), as in rubin_campaign_aggregate.R. Every expected fit x arm (x trait) gets a row with
# a status: "scored"; "one_class" (the trait, or bin for c1 ~ bin, has one realised class in the complete data, so
# there is nothing to score; score_discrete's and score_discrete_estimand's own rules); "undefined" (c1 ~ bin only:
# bin is constant in every imputed dataset, i.e. one observed class); "missing" (no row for another reason, e.g. an
# arm or scoring error; the fit's error messages are attached). Metric columns are NA unless scored.
# Arms: the fit's --arms plus the discrete-only comparison arms freqA_er and freqB_er (simulation v1's equal-rates
# castor models on the same continuous draws), which rubin_cell.R adds to frequentist fits.
# A malformed result file stops the run, listing every such file.
#
# Outputs (OUT_DIR):
#   fit_disc_cells.csv.gz      fit x arm x trait: status and the disc_cells scores, plus the reference arm
#                              "mode_floor" (the most frequent observed class predicted in every masked cell; ties
#                              fractional as in score_discrete; its other scores NA). A dataset in both sets carries
#                              its mode floor twice (once per fit); aggregates count it once.
#   fit_disc_estimands.csv.gz  fit x arm: status, the c1 ~ bin disc_estimands columns, truth (the population
#                              target, joined by n, lambda, rho), covered (lower <= truth <= upper),
#                              complete_covered (the same for the complete-data interval; NA in files written
#                              before 7a9d2ee); target_cond (the per-dataset target, rho x GLS slope of the bin
#                              liability on bin; NA before 4a4395a), covered_cond, complete_covered_cond
#   disc_fill.csv.gz           fit x arm x trait: missing cells filled by fill_degenerate (per dataset; summed over
#                              the M datasets in files written before 4a4395a), masked cells, classes observed
#   castor_diag.csv.gz         freq fit x arm x trait: rate model, castor's fitted rate, bootstrap median rate,
#                              resimulations, failed refits, the trait's error, wall time; one row with the error when
#                              castor did not run
#   agg_disc_cells.csv         n, lambda, rho, arm, trait: n_fits, mean and MC SE of each score, n_fits_trait_absent
#                              (one class in the complete data), n_fits_missing, n_na_draws (NA draws among scored
#                              cells, dropped cell by cell by score_discrete)
#   agg_disc_cells_l.csv       the same pooled over rho (the per-value tables of study.qmd)
#   agg_disc_down.csv          n, lambda, rho, arm: c1 ~ bin coverage of the population truth and of the per-dataset
#                              target, complete-data coverage of both on the same datasets, bias against both, RMSE,
#                              mean SE, empirical SD, mean width, mean FMI, estimate minus the complete-data estimate;
#                              arm "complete" is the complete-data analysis, once per dataset
#   agg_disc_paired.csv        per-dataset differences between arms (bace_chain - freqA, bace - freqA,
#                              bace_chain - bace, freqA - freqB, freqA - freqA_er, freqB - freqB_er) in accuracy and
#                              Brier score, by n, lambda, rho, trait
#   repro.csv                  fits in both pools (same set and tag): max and median |difference| of the continuous
#                              rows (pooled estimate and SE; per-value zRMSE and coverage) per set, host pair (host
#                              of the re-run, cont_host of the stored fit) and arm; "all" rows pool host pairs or
#                              arms; n_fits_identical counts fits with every difference < tol
#   sanity.txt                 counts and the cross-set checks

`%||%` <- function(x, y) if (is.null(x)) y else x

ID_COLS <- c("set", "host", "tag", "n", "lambda", "rho", "seed")
DISC_CELL_COLS <- c("n_cells", "accuracy", "brier", "ece", "set_coverage", "set_size", "mae_class", "frac_unanimous",
                    "n_na")
DISC_CELL_OPT <- "n_cells_all_na"                                      # added in 4a4395a
CELL_METRICS <- c("accuracy", "brier", "ece", "set_coverage", "set_size", "mae_class", "frac_unanimous")
DISC_EST_COLS <- c("estimate", "se", "lower", "upper", "df", "fmi", "m_ok", "complete_data")
DISC_EST_OPT <- c("complete_se", "complete_lower", "complete_upper",   # added in 7a9d2ee
                  "target_cond")                                        # added in 4a4395a
CASTOR_COLS <- c("rate_model", "degenerate", "q_hat", "q_star_med", "n_resim", "n_refit_fail", "n_reject", "n_fallback")
FREQ_DISC_ARMS <- c("freqA", "freqA_er", "freqB", "freqB_er")
CONTRASTS <- list(c("bace_chain", "freqA"), c("bace", "freqA"), c("bace_chain", "bace"), c("freqA", "freqB"),
                  c("freqA", "freqA_er"), c("freqB", "freqB_er"))
# "identical": freq A's refits reproduce across hosts to about 1e-6 (6e-7 measured, Mac against fir); a BACE chain
# on another build differs by 1e-2 or more
REPRO_TOL <- 1e-5

mc_se <- function(x) if (length(x) < 2L) NA_real_ else stats::sd(x) / sqrt(length(x))
mn <- function(x) if (length(x)) mean(x) else NA_real_
cell_key <- function(n, lambda, rho) paste(as.integer(n), as.character(round(as.numeric(lambda), 6)),
                                           as.character(round(as.numeric(rho), 6)))
tag_of_file <- function(f) sub("(_smoke)?(_dup[0-9]+)?\\.rds$", "", basename(f))
na_cols <- function(cols) as.data.frame(stats::setNames(rep(list(NA_real_), length(cols)), cols))
n_real <- function(f) nlevels(droplevels(f))                            # realised classes
#' The discrete arms of a fit: its --arms, the arms it scored, and for frequentist fits the equal-rates comparison
#' arms of each frequentist arm it ran (expected even when absent, so that a failure shows as "missing").
disc_arms_of <- function(x) {
  a <- unique(c(as.character(x$arms), names(x$disc_fill)))
  base <- intersect(c("freqA", "freqB"), a)
  if (length(base)) a <- unique(c(a, paste0(base, "_er")))    # paste0(character(0), "_er") would give "_er"
  c(intersect(FREQ_DISC_ARMS, a), setdiff(a, FREQ_DISC_ARMS))
}

#' Result files of a pool, one per set x tag: POOL/{bace,freq}/<host>/*.rds (a copy of the same fit on a second host,
#' or a _dupN copy, is dropped; the first in path order is kept, as in build_data.R).
list_pool <- function(pool) {
  f <- list.files(file.path(pool, c("bace", "freq")), "\\.rds$", recursive = TRUE, full.names = TRUE)
  d <- data.frame(file = f, set = basename(dirname(dirname(f))), host = basename(dirname(f)), stringsAsFactors = FALSE)
  bad <- d$file[!(d$set %in% c("bace", "freq"))]
  if (length(bad)) stop("files outside POOL/{bace,freq}/<host>/:\n", paste(bad, collapse = "\n"), call. = FALSE)
  key <- paste(d$set, tag_of_file(d$file))
  structure(d[!duplicated(key), , drop = FALSE], n_dup = sum(duplicated(key)))
}

#' NULL when x is a usable discrete result, otherwise the reason it is malformed.
check_fit <- function(x, file) {
  if (inherits(x, "error")) return(paste("unreadable:", conditionMessage(x)))
  if (!is.list(x)) return("not a list")
  need <- c("tag", "n", "lambda", "rho", "seed", "M", "arms", "truth", "mask")
  miss <- need[vapply(need, function(k) is.null(x[[k]]), logical(1))]
  if (length(miss)) return(paste("missing field(s):", paste(miss, collapse = ", ")))
  if (isTRUE(x$smoke)) return("a --smoke fit")
  if (!identical(tag_of_file(file), x$tag)) return(paste("file name does not match its tag", x$tag))
  if (!is.data.frame(x$truth)) return("truth is not a data.frame")
  if (!is.matrix(x$mask) || !is.logical(x$mask) || !identical(dim(x$mask), dim(x$truth)) ||
      !identical(dimnames(x$mask), dimnames(x$truth))) return("mask does not match truth")
  traits <- names(x$truth)[vapply(x$truth, is.factor, logical(1))]
  if (!length(traits)) return("no discrete trait in truth")
  if (is.null(x$disc_fill)) return("no discrete results (run without --discrete)")
  arms <- disc_arms_of(x)
  if (!is.list(x$disc_fill) || !all(names(x$disc_fill) %in% arms)) return("disc_fill names are not arms of the fit")
  dc <- x$disc_cells
  if (!is.null(dc)) {
    if (!is.data.frame(dc) || !all(c("arm", "trait", DISC_CELL_COLS) %in% names(dc))) return("disc_cells lacks columns")
    if (!all(dc$arm %in% arms) || !all(dc$trait %in% traits) || anyDuplicated(paste(dc$arm, dc$trait)))
      return("disc_cells rows do not match the fit's arms and traits")
  }
  de <- x$disc_estimands
  if (!is.null(de)) {
    if (!is.data.frame(de) || !all(c("arm", "estimand", DISC_EST_COLS) %in% names(de)))
      return("disc_estimands lacks columns")
    if (!all(de$arm %in% arms) || anyDuplicated(de$arm) || !all(de$estimand == "slope_c1_bin"))
      return("disc_estimands rows do not match the fit's arms")
  }
  NULL
}

#' Mode floor: predict the most frequent observed class in every masked cell; a tie of k classes scores 1/k when the
#' truth is among them (score_discrete's rule). Traits score_discrete skips (no masked cell, one class) are left out.
mode_floor <- function(truth, mask, traits) {
  rows <- lapply(traits, function(v) {
    idx <- which(mask[, v]); lev <- levels(truth[[v]])
    if (!length(idx) || n_real(truth[[v]]) < 2L) return(NULL)
    f <- tabulate(as.integer(truth[[v]][!mask[, v]]), length(lev))
    acc <- if (max(f) == 0) NA_real_ else {
      top <- which(f == max(f)); mean((as.integer(truth[[v]][idx]) %in% top) / length(top))
    }
    r <- cbind(data.frame(arm = "mode_floor", trait = v, stringsAsFactors = FALSE), na_cols(c(DISC_CELL_COLS, DISC_CELL_OPT)))
    r$n_cells <- length(idx); r$accuracy <- acc
    r
  })
  do.call(rbind, rows)
}

#' The per-fit rows of every output table, from one checked result list.
extract_fit <- function(x, set, host) {
  truth <- x$truth; mask <- x$mask; arms <- disc_arms_of(x)
  traits <- names(truth)[vapply(truth, is.factor, logical(1))]
  ids <- data.frame(set = set, host = host, tag = x$tag, n = as.integer(x$n), lambda = as.numeric(x$lambda),
                    rho = as.numeric(x$rho), seed = as.integer(x$seed), stringsAsFactors = FALSE)
  errs <- x$errors %||% list()
  err_for <- function(arm) {
    base <- sub("_er$", "", arm)                                         # an _er arm fails with its base arm
    k <- intersect(unique(c(base, arm, paste0(arm, c("_score", "_disc_score", "_disc_est", "_castor")),
                            if (startsWith(arm, "bace")) "bace_fit")), names(errs))
    if (!length(k)) return(NA_character_)
    paste(sprintf("%s: %s", k, vapply(errs[k], function(e) substr(paste(as.character(e), collapse = " "), 1, 200), "")),
          collapse = " | ")
  }
  n_masked <- colSums(mask[, traits, drop = FALSE])
  one_class <- vapply(traits, function(v) n_real(truth[[v]]) < 2L, logical(1))   # score_discrete's rule
  n_obs_cls <- vapply(traits, function(v) length(unique(truth[[v]][!mask[, v]])), integer(1))
  # expected (arm, trait) rows with a status; scored rows come from `scored` (arm, trait + DISC_CELL_COLS)
  status_rows <- function(arm_set, scored) {
    if (!is.null(scored)) for (k in DISC_CELL_OPT) if (!(k %in% names(scored))) scored[[k]] <- NA_real_
    do.call(rbind, lapply(arm_set, function(arm) do.call(rbind, lapply(traits, function(v) {
      r <- if (!is.null(scored)) scored[scored$arm == arm & scored$trait == v, c("arm", "trait", DISC_CELL_COLS,
                                                                              DISC_CELL_OPT), drop = FALSE]
      if (!is.null(r) && nrow(r) == 1L) return(cbind(ids, r, status = "scored", errors = NA_character_))
      st <- if (one_class[[v]] && n_masked[[v]] > 0) "one_class" else "missing"
      er <- if (st == "missing") (if (n_masked[[v]] == 0) "no masked cell" else err_for(arm)) else NA_character_
      cbind(ids, data.frame(arm = arm, trait = v, stringsAsFactors = FALSE), na_cols(c(DISC_CELL_COLS, DISC_CELL_OPT)),
            status = st, errors = er)
    }))))
  }
  cells <- rbind(status_rows(arms, x$disc_cells), status_rows("mode_floor", mode_floor(truth, mask, traits)))

  bin_one_class <- "bin" %in% names(truth) && nlevels(droplevels(truth$bin)) < 2L   # score_discrete_estimand's rule
  est_status <- x$disc_est_status %||% character(0)                   # written from 4a4395a
  de <- x$disc_estimands
  if (!is.null(de)) for (k in DISC_EST_OPT) if (!(k %in% names(de))) de[[k]] <- NA_real_
  est <- do.call(rbind, lapply(arms, function(arm) {
    r <- if (!is.null(de)) de[de$arm == arm, c("arm", "estimand", DISC_EST_COLS, DISC_EST_OPT), drop = FALSE]
    if (!is.null(r) && nrow(r) == 1L) return(cbind(ids, r, status = "scored", errors = NA_character_))
    st <- if (bin_one_class) "one_class" else if (identical(unname(est_status[arm]), "undefined")) "undefined" else "missing"
    cbind(ids, data.frame(arm = arm, estimand = "slope_c1_bin", stringsAsFactors = FALSE),
          na_cols(c(DISC_EST_COLS, DISC_EST_OPT)), status = st, errors = if (st == "missing") err_for(arm) else NA_character_)
  }))

  fill <- do.call(rbind, lapply(arms, function(arm) {
    f <- x$disc_fill[[arm]]
    data.frame(ids, arm = arm, trait = traits, n_filled = if (is.null(f)) NA_integer_ else as.integer(f[traits]),
               n_masked = as.integer(n_masked), n_classes_observed = n_obs_cls, row.names = NULL,
               stringsAsFactors = FALSE)
  }))

  castor <- do.call(rbind, lapply(intersect(FREQ_DISC_ARMS, arms), function(arm) {
    d <- x$diag[[paste0(arm, "_castor")]]
    w <- unname(x$walls[paste0(arm, "_castor")]); if (!length(w)) w <- NA_real_
    if (is.data.frame(d) && nrow(d)) {
      for (k in c(CASTOR_COLS, "error")) if (!(k %in% names(d))) d[[k]] <- NA
      return(data.frame(ids, arm = arm, trait = as.character(d$trait), d[, CASTOR_COLS], wall_s = w,
                        error = as.character(d$error), row.names = NULL, stringsAsFactors = FALSE))
    }
    e <- err_for(arm)
    cbind(ids, data.frame(arm = arm, trait = NA_character_, stringsAsFactors = FALSE), na_cols(CASTOR_COLS),
          data.frame(wall_s = w, error = if (is.na(e)) "castor did not run" else e, stringsAsFactors = FALSE))
  }))

  keep_df <- function(d, cols) if (is.data.frame(d) && all(cols %in% names(d))) d[, cols, drop = FALSE] else NULL
  list(set = set, host = host, tag = x$tag, cells = cells, est = est, fill = fill, castor = castor,
       cont = list(est = keep_df(x$estimands, c("arm", "estimand", "estimate", "se")),
                   cells = keep_df(x$cells, c("arm", "trait", "zRMSE", "coverage"))),
       fp = list(rn = rownames(truth), mask = mask,
                 disc = vapply(traits, function(v) as.character(truth[[v]]), character(nrow(truth))),
                 c1 = truth$c1),
       has_imp = !is.null(x$imputations),
       has_ci = if (is.null(x$disc_estimands)) NA else all(DISC_EST_OPT %in% names(x$disc_estimands)))
}

#' Monte Carlo truth table (script/rubin_discrete_truth.R): one row per n, lambda, rho.
read_truth <- function(truth_csv) {
  if (!file.exists(truth_csv)) stop("truth table not found: ", truth_csv, call. = FALSE)
  tt <- utils::read.csv(truth_csv, stringsAsFactors = FALSE)
  need <- c("n", "lambda", "rho", "truth", "mc_se")
  if (!all(need %in% names(tt))) stop("truth table lacks columns: ", paste(setdiff(need, names(tt)), collapse = ", "),
                                      call. = FALSE)
  tt$key <- cell_key(tt$n, tt$lambda, tt$rho)
  if (anyDuplicated(tt$key)) stop("truth table has duplicated n, lambda, rho rows", call. = FALSE)
  tt
}

#' Per-value aggregates by `by`; the mode floor is counted once per dataset.
summ_cells <- function(FC, by) {
  mf <- FC$arm == "mode_floor"
  d <- rbind(FC[!mf, ], FC[mf, ][!duplicated(paste(FC$tag[mf], FC$trait[mf])), ])
  out <- do.call(rbind, lapply(split(d, d[by], drop = TRUE), function(s) {
    r <- s[1, by, drop = FALSE]; ok <- s[s$status == "scored", ]
    r$n_fits <- nrow(ok); r$n_fits_trait_absent <- sum(s$status == "one_class")
    r$n_fits_missing <- sum(s$status == "missing")
    for (v in CELL_METRICS) { r[[v]] <- mn(ok[[v]]); r[[paste0(v, "_se")]] <- mc_se(ok[[v]]) }
    r$n_na_draws <- if (all(is.na(ok$n_na))) NA_real_ else sum(ok$n_na)   # NA for the mode floor (no draws)
    r
  }))
  rownames(out) <- NULL
  out[do.call(order, out[by]), ]
}

#' Downstream (c1 ~ bin) aggregates by n, lambda, rho, arm, with the complete-data analysis as arm "complete"
#' (once per dataset; excluded when bin has one class, missing when no arm of the dataset has a row).
summ_down <- function(FE) {
  by <- c("n", "lambda", "rho", "arm")
  cp <- do.call(rbind, lapply(split(FE, FE$tag), function(s) {
    sc <- s[s$status == "scored", ]
    st <- if (nrow(sc)) "scored" else if (any(s$status == "one_class")) "one_class" else
      if (any(s$status == "undefined")) "undefined" else "missing"
    r <- if (nrow(sc)) sc[1, ] else s[1, ]
    data.frame(r[, c(ID_COLS, "truth", "truth_mc_se")], arm = "complete", status = st, estimate = r$complete_data,
               se = r$complete_se, lower = r$complete_lower, upper = r$complete_upper, fmi = 0, m_ok = NA_real_,
               complete_data = r$complete_data, covered = r$complete_covered, complete_covered = r$complete_covered,
               target_cond = r$target_cond, covered_cond = r$complete_covered_cond,
               complete_covered_cond = r$complete_covered_cond, stringsAsFactors = FALSE)
  }))
  cols <- names(cp)
  d <- rbind(FE[, cols], cp)
  out <- do.call(rbind, lapply(split(d, d[by], drop = TRUE), function(s) {
    ok <- s[s$status == "scored", ]
    err <- ok$estimate - ok$truth; dcd <- ok$estimate - ok$complete_data; errc <- ok$estimate - ok$target_cond
    data.frame(s[1, by, drop = FALSE], n_fits = nrow(ok), n_excluded = sum(s$status == "one_class"),
               n_undefined = sum(s$status == "undefined"),
               n_missing = sum(s$status == "missing"), truth = s$truth[1], truth_mc_se = s$truth_mc_se[1],
               coverage = mn(ok$covered), coverage_se = mc_se(as.numeric(ok$covered)),
               complete_coverage = mn(ok$complete_covered), complete_coverage_se = mc_se(as.numeric(ok$complete_covered)),
               coverage_cond = mn(ok$covered_cond), coverage_cond_se = mc_se(as.numeric(ok$covered_cond)),
               complete_coverage_cond = mn(ok$complete_covered_cond),
               bias_cond = mn(errc), bias_cond_se = mc_se(errc),
               bias = mn(err), bias_se = mc_se(err), rmse = sqrt(mn(err^2)), se_mean = mn(ok$se),
               emp_sd = if (nrow(ok) > 1L) stats::sd(ok$estimate) else NA_real_, width = mn(ok$upper - ok$lower),
               fmi_mean = mn(ok$fmi), m_ok_min = if (nrow(ok)) min(ok$m_ok) else NA_real_,
               diff_cd_mean = mn(dcd), diff_cd_sd = if (nrow(ok) > 1L) stats::sd(dcd) else NA_real_,
               diff_cd_se = mc_se(dcd), stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out[do.call(order, out[by]), ]
}

#' Paired per-dataset differences between arms in accuracy and Brier score, by n, lambda, rho, trait.
pair_disc <- function(FC) {
  sc <- FC[FC$status == "scored" & FC$arm != "mode_floor", ]
  out <- do.call(rbind, lapply(CONTRASTS, function(p) {
    m <- merge(sc[sc$arm == p[1], c("tag", "n", "lambda", "rho", "trait", "accuracy", "brier")],
               sc[sc$arm == p[2], c("tag", "trait", "accuracy", "brier")], by = c("tag", "trait"))
    if (!nrow(m)) return(NULL)
    do.call(rbind, lapply(c("accuracy", "brier"), function(v) {
      m$d <- m[[paste0(v, ".x")]] - m[[paste0(v, ".y")]]
      do.call(rbind, lapply(split(m, m[c("n", "lambda", "rho", "trait")], drop = TRUE), function(g) data.frame(
        contrast = paste(p[1], "-", p[2]), measure = v, n = g$n[1], lambda = g$lambda[1], rho = g$rho[1],
        trait = g$trait[1], diff = mean(g$d), se = mc_se(g$d), datasets = nrow(g), stringsAsFactors = FALSE)))
    }))
  }))
  if (is.null(out)) return(data.frame(contrast = character(0), measure = character(0), n = integer(0),
                                      lambda = numeric(0), rho = numeric(0), trait = character(0), diff = numeric(0),
                                      se = numeric(0), datasets = integer(0)))
  rownames(out) <- NULL
  out[order(match(out$contrast, vapply(CONTRASTS, paste, "", collapse = " - ")), out$measure, out$n, out$lambda,
            out$rho, out$trait), ]
}

#' Reproduction of the stored continuous rows: each fit of the discrete pool against the fit of record with the same
#' set and tag in the continuous pool (first in path order, as build_data.R keeps).
repro_check <- function(fits, cont_pool) {
  empty <- data.frame(set = character(0), host = character(0), cont_host = character(0), arm = character(0),
                      n_fits = integer(0),
                      n_fits_identical = integer(0), est_rows = integer(0), max_abs_estimate = numeric(0),
                      median_abs_estimate = numeric(0), max_abs_se = numeric(0), median_abs_se = numeric(0),
                      cell_rows = integer(0), max_abs_zRMSE = numeric(0), median_abs_zRMSE = numeric(0),
                      max_abs_coverage = numeric(0), median_abs_coverage = numeric(0), tol = numeric(0))
  if (is.null(cont_pool) || !dir.exists(cont_pool)) return(list(table = empty, n_matched = 0L, n_unreadable = 0L))
  cf <- list.files(file.path(cont_pool, c("bace", "freq")), "\\.rds$", recursive = TRUE, full.names = TRUE)
  ckey <- paste(basename(dirname(dirname(cf))), tag_of_file(cf))
  cf <- cf[!duplicated(ckey)]; ckey <- ckey[!duplicated(ckey)]
  rows <- list(); n_unreadable <- 0L
  for (r in fits) {
    j <- match(paste(r$set, r$tag), ckey); if (is.na(j)) next
    y <- tryCatch(readRDS(cf[j]), error = function(e) NULL)
    if (is.null(y)) { n_unreadable <- n_unreadable + 1L; next }
    base <- data.frame(set = r$set, host = r$host, cont_host = basename(dirname(cf[j])), tag = r$tag,
                       stringsAsFactors = FALSE)
    if (!is.null(r$cont$est) && is.data.frame(y$estimands)) {
      m <- merge(r$cont$est, y$estimands[, c("arm", "estimand", "estimate", "se")], by = c("arm", "estimand"))
      if (nrow(m)) rows[[length(rows) + 1L]] <- data.frame(base, arm = m$arm, kind = "estimand",
        d1 = abs(m$estimate.x - m$estimate.y), d2 = abs(m$se.x - m$se.y), stringsAsFactors = FALSE)
    }
    if (!is.null(r$cont$cells) && is.data.frame(y$cells)) {
      m <- merge(r$cont$cells, y$cells[, c("arm", "trait", "zRMSE", "coverage")], by = c("arm", "trait"))
      if (nrow(m)) rows[[length(rows) + 1L]] <- data.frame(base, arm = m$arm, kind = "cell",
        d1 = abs(m$zRMSE.x - m$zRMSE.y), d2 = abs(m$coverage.x - m$coverage.y), stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(list(table = empty, n_matched = 0L, n_unreadable = n_unreadable))
  D <- do.call(rbind, rows)
  mx <- function(x) if (length(x)) max(x) else NA_real_       # NA propagates: a non-finite stored value shows
  md <- function(x) if (length(x)) stats::median(x) else NA_real_
  summ <- function(s, set, host, cont_host, arm) {
    e <- s[s$kind == "estimand", ]; ce <- s[s$kind == "cell", ]
    per_fit <- tapply(pmax(s$d1, s$d2), s$tag, max)             # a fit with an NA difference is not identical
    data.frame(set = set, host = host, cont_host = cont_host, arm = arm, n_fits = length(per_fit),
               n_fits_identical = sum(per_fit < REPRO_TOL, na.rm = TRUE),
               est_rows = nrow(e), max_abs_estimate = mx(e$d1), median_abs_estimate = md(e$d1),
               max_abs_se = mx(e$d2), median_abs_se = md(e$d2), cell_rows = nrow(ce), max_abs_zRMSE = mx(ce$d1),
               median_abs_zRMSE = md(ce$d1), max_abs_coverage = mx(ce$d2), median_abs_coverage = md(ce$d2),
               tol = REPRO_TOL, stringsAsFactors = FALSE)
  }
  # an exact reproduction is expected only where the same cluster (same R and BACE builds) ran both fits, so the
  # rows break down by the pair (host of the re-run, host of the stored fit); host = cont_host = "all" pools them
  out <- list()
  for (st in sort(unique(D$set))) {
    Ds <- D[D$set == st, ]
    pairs <- unique(Ds[order(Ds$host, Ds$cont_host), c("host", "cont_host")])
    for (k in c(0L, seq_len(nrow(pairs)))) {
      Dh <- if (k == 0L) Ds else Ds[Ds$host == pairs$host[k] & Ds$cont_host == pairs$cont_host[k], ]
      h <- if (k == 0L) "all" else pairs$host[k]; ch <- if (k == 0L) "all" else pairs$cont_host[k]
      for (a in c("all", sort(unique(Dh$arm))))
        out[[length(out) + 1L]] <- summ(if (a == "all") Dh else Dh[Dh$arm == a, ], st, h, ch, a)
    }
  }
  list(table = do.call(rbind, out), n_matched = length(unique(paste(D$set, D$tag))), n_unreadable = n_unreadable)
}

#' Read the pool, check every file, write every table. Returns the tables invisibly.
aggregate_disc <- function(pool, truth_csv, out, cont_pool = NULL) {
  fl <- list_pool(pool)
  if (!nrow(fl)) stop("no result files under ", file.path(pool, "{bace,freq}", "<host>"), call. = FALSE)
  tt <- read_truth(truth_csv)
  fits <- vector("list", nrow(fl)); bad <- character(0)
  seen <- new.env(); mism <- character(0); both <- 0L; max_dc1 <- 0
  for (i in seq_len(nrow(fl))) {
    x <- tryCatch(readRDS(fl$file[i]), error = function(e) e)
    msg <- check_fit(x, fl$file[i])
    if (!is.null(msg)) { bad <- c(bad, sprintf("%s: %s", fl$file[i], msg)); next }
    r <- extract_fit(x, fl$set[i], fl$host[i]); rm(x)
    o <- seen[[r$tag]]
    if (is.null(o)) assign(r$tag, r$fp, envir = seen) else {
      both <- both + 1L
      dc1 <- if (!is.null(o$c1) && !is.null(r$fp$c1) && length(o$c1) == length(r$fp$c1)) max(abs(o$c1 - r$fp$c1)) else Inf
      max_dc1 <- max(max_dc1, dc1)
      if (!identical(o$rn, r$fp$rn) || !identical(o$mask, r$fp$mask) || !identical(o$disc, r$fp$disc) || dc1 > 1e-8)
        mism <- c(mism, r$tag)
    }
    r$fp <- NULL; fits[[i]] <- r
  }
  if (length(bad)) stop(sprintf("%d malformed result file(s):\n%s", length(bad), paste(bad, collapse = "\n")),
                        call. = FALSE)
  if (length(mism)) stop(sprintf("%d dataset(s) differ between the bace and freq files (mask, discrete traits or c1):\n%s",
                                 length(mism), paste(mism, collapse = "\n")), call. = FALSE)
  pick <- function(k) do.call(rbind, lapply(fits, `[[`, k))
  FC <- pick("cells"); FE <- pick("est"); FILL <- pick("fill"); CAS <- pick("castor")
  arm_dup <- function(d, cols) { s <- d[d$status == "scored" & d$arm != "mode_floor", ]
    unique(s$tag[duplicated(do.call(paste, s[cols]))]) }
  dup <- c(arm_dup(FC, c("tag", "arm", "trait")), arm_dup(FE, c("tag", "arm")))
  if (length(dup)) stop("an arm is scored in two fits of the same dataset: ", paste(unique(dup), collapse = ", "),
                        call. = FALSE)

  m <- match(cell_key(FE$n, FE$lambda, FE$rho), tt$key)
  FE$truth <- tt$truth[m]; FE$truth_mc_se <- tt$mc_se[m]
  FE$covered <- FE$lower <= FE$truth & FE$truth <= FE$upper
  FE$complete_covered <- FE$complete_lower <= FE$truth & FE$truth <= FE$complete_upper
  FE$covered_cond <- FE$lower <= FE$target_cond & FE$target_cond <= FE$upper
  FE$complete_covered_cond <- FE$complete_lower <= FE$target_cond & FE$target_cond <= FE$complete_upper
  ord_fit <- function(d, extra) d[do.call(order, d[c("n", "lambda", "rho", "seed", "set", extra)]), ]
  FC <- ord_fit(FC[, c(ID_COLS, "arm", "trait", "status", DISC_CELL_COLS, DISC_CELL_OPT, "errors")], c("arm", "trait"))
  FE <- ord_fit(FE[, c(ID_COLS, "arm", "status", "estimand", DISC_EST_COLS, DISC_EST_OPT, "truth", "truth_mc_se",
                       "covered", "complete_covered", "covered_cond", "complete_covered_cond", "errors")], "arm")
  FILL <- ord_fit(FILL, c("arm", "trait"))
  if (!is.null(CAS)) CAS <- ord_fit(CAS, c("arm", "trait"))

  agg_cells <- summ_cells(FC, c("n", "lambda", "rho", "arm", "trait"))
  agg_cells_l <- summ_cells(FC, c("n", "lambda", "arm", "trait"))
  agg_down <- summ_down(FE)
  agg_paired <- pair_disc(FC)
  rp <- repro_check(fits, cont_pool)

  # construction checks: BACE + residual perturbs Gaussian traits only, so its discrete scores are BACE's; the
  # mode floor and the complete-data c1 ~ bin estimate are properties of the dataset, so equal across its fits
  sc <- FC[FC$status == "scored", ]
  rb <- merge(sc[sc$arm == "bace", c("tag", "trait", CELL_METRICS)], sc[sc$arm == "bace_resid", c("tag", "trait", CELL_METRICS)],
              by = c("tag", "trait"))
  d_resid <- if (nrow(rb)) max(abs(as.matrix(rb[paste0(CELL_METRICS, ".x")]) - as.matrix(rb[paste0(CELL_METRICS, ".y")])),
                              na.rm = TRUE) else NA
  spread <- function(v, g) if (length(v)) max(tapply(v, g, function(a) diff(range(a)))) else NA
  mf <- sc[sc$arm == "mode_floor", ]; ce <- FE[FE$status == "scored", ]
  na_d <- sc[sc$arm != "mode_floor" & sc$n_na > 0, ]
  na_d <- if (nrow(na_d)) tapply(na_d$n_na, paste(na_d$arm, na_d$trait), sum) else integer(0)
  hc <- vapply(fits, function(r) r$has_ci, NA)
  n_set <- table(factor(vapply(fits, `[[`, "", "set"), levels = c("bace", "freq")))
  host_tab <- table(paste(fl$set, fl$host))
  no_truth <- FE$status == "scored" & is.na(FE$truth)
  sanity <- c(
    sprintf("files %d after dropping %d duplicate(s) (same set and tag on another host); bace %d, freq %d",
            nrow(fl), attr(fl, "n_dup"), n_set[["bace"]], n_set[["freq"]]),
    sprintf("by set and host: %s", paste(names(host_tab), host_tab, sep = " = ", collapse = "; ")),
    sprintf("datasets in both sets %d: mask and discrete traits identical, max |c1 difference| %.2e", both,
            if (both) max_dc1 else NA),
    sprintf("same dataset across fits: max range of the mode floor %.2e, of the complete-data c1 ~ bin estimate %.2e",
            spread(mf$accuracy, paste(mf$tag, mf$trait)), spread(ce$complete_data, ce$tag)),
    sprintf("BACE + residual vs BACE as shipped, discrete per-value scores: max |difference| %s",
            format(signif(d_resid, 3))),
    sprintf("per-value rows (fit x arm x trait, mode floor included): scored %d, one class in the complete data %d, missing %d",
            sum(FC$status == "scored"), sum(FC$status == "one_class"), sum(FC$status == "missing")),
    sprintf("c1 ~ bin rows (fit x arm): scored %d, bin one class in the complete data %d, undefined (bin constant in every imputed dataset) %d, missing %d",
            sum(FE$status == "scored"), sum(FE$status == "one_class"), sum(FE$status == "undefined"),
            sum(FE$status == "missing")),
    sprintf("scored c1 ~ bin rows without the per-dataset target (files written before 4a4395a): %d",
            sum(FE$status == "scored" & is.na(FE$target_cond))),
    sprintf("NA draws among scored cells (dropped cell by cell by score_discrete): %s",
            if (length(na_d)) paste(names(na_d), na_d, collapse = "; ") else "none"),
    sprintf("scored c1 ~ bin rows without a truth row: %d%s", sum(no_truth),
            if (any(no_truth)) paste0(" (n lambda rho: ", paste(unique(cell_key(FE$n, FE$lambda, FE$rho)[no_truth]),
                                                                collapse = "; "), ")") else ""),
    sprintf("files without the complete-data interval (written before 7a9d2ee): %d; files without saved imputations: %d",
            sum(hc %in% FALSE), sum(!vapply(fits, `[[`, NA, "has_imp"))),
    sprintf("reproduction check: %d fit(s) matched in the continuous pool %s (%d unreadable)", rp$n_matched,
            if (is.null(cont_pool)) "(none given)" else cont_pool, rp$n_unreadable))

  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(FC, gzfile(file.path(out, "fit_disc_cells.csv.gz")), row.names = FALSE)
  utils::write.csv(FE, gzfile(file.path(out, "fit_disc_estimands.csv.gz")), row.names = FALSE)
  utils::write.csv(FILL, gzfile(file.path(out, "disc_fill.csv.gz")), row.names = FALSE)
  utils::write.csv(CAS %||% data.frame(set = character(0), host = character(0), tag = character(0), arm = character(0),
                                       trait = character(0)), gzfile(file.path(out, "castor_diag.csv.gz")), row.names = FALSE)
  utils::write.csv(agg_cells, file.path(out, "agg_disc_cells.csv"), row.names = FALSE)
  utils::write.csv(agg_cells_l, file.path(out, "agg_disc_cells_l.csv"), row.names = FALSE)
  utils::write.csv(agg_down, file.path(out, "agg_disc_down.csv"), row.names = FALSE)
  utils::write.csv(agg_paired, file.path(out, "agg_disc_paired.csv"), row.names = FALSE)
  utils::write.csv(rp$table, file.path(out, "repro.csv"), row.names = FALSE)
  writeLines(sanity, file.path(out, "sanity.txt"))
  cat(sanity, sep = "\n")
  invisible(list(fit_cells = FC, fit_estimands = FE, fill = FILL, castor = CAS, agg_cells = agg_cells,
                 agg_cells_l = agg_cells_l, agg_down = agg_down, agg_paired = agg_paired, repro = rp$table,
                 sanity = sanity))
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  home <- Sys.getenv("HOME")
  aggregate_disc(pool = if (length(a) >= 1) a[1] else file.path(home, "pigauto_rubin_disc_pool"),
                 truth_csv = if (length(a) >= 2) a[2] else
                   file.path("script", "rubin_study", "data", "discrete", "truth_slope_c1_bin.csv"),
                 out = if (length(a) >= 3) a[3] else file.path("script", "rubin_study", "data", "discrete"),
                 cont_pool = if (length(a) >= 4) a[4] else file.path(home, "pigauto_rubin_pool"))
}
