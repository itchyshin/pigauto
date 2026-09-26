# Discrete aggregation (script/rubin_discrete_aggregate.R) on hand-made result lists: the mode floor and its ties, the
# truth join and coverage, paired differences between arms of the same dataset, the status of traits and estimands
# that cannot be scored, duplicates across hosts, the reproduction check, and loud failures.
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
suppressMessages(source(file.path(root, "script", "rubin_discrete_aggregate.R")))

sp <- paste0("t", 1:10)
# rows 8-10 are masked in every discrete trait; observed rows 1-7 give bin "yes" 4:3, cat3 A 3 = B 3 > C 1 (a tie),
# ord L1 2 = L2 2 = L3 2 > L4 1 (a three-way tie)
mk_truth <- function(bin = c(rep("yes", 4), rep("no", 3), "yes", "no", "yes")) data.frame(
  row.names = sp, c1 = seq(-1, 1, length.out = 10), bin = factor(bin),
  ord = factor(c("L1", "L1", "L2", "L2", "L3", "L3", "L4", "L1", "L4", "L3"), levels = paste0("L", 1:4), ordered = TRUE),
  cat3 = factor(c("A", "A", "A", "B", "B", "B", "C", "A", "C", "B"), levels = c("A", "B", "C")))
mk_mask <- function(truth) {
  m <- matrix(FALSE, nrow(truth), ncol(truth), dimnames = dimnames(truth)); m[8:10, c("bin", "ord", "cat3")] <- TRUE; m
}
cells_df <- function(arms, acc, brier, traits = c("bin", "ord", "cat3")) {
  g <- expand.grid(trait = traits, arm = arms, stringsAsFactors = FALSE)
  data.frame(arm = g$arm, trait = g$trait, n_cells = 3L, accuracy = acc, brier = brier, ece = 0, set_coverage = 1,
             set_size = 1, mae_class = ifelse(g$trait == "ord", 0.5, NA), frac_unanimous = 1, n_na = 0L,
             stringsAsFactors = FALSE)
}
est_df <- function(arms, estimate, lower, upper, cd, clo, cup) data.frame(
  arm = arms, estimand = "slope_c1_bin", estimate = estimate, se = 0.1, lower = lower, upper = upper, df = 20, fmi = 0.3,
  m_ok = 20L, complete_data = cd, complete_se = 0.05, complete_lower = clo, complete_upper = cup, stringsAsFactors = FALSE)
mk_fit <- function(seed, arms, disc_cells, disc_estimands, truth = mk_truth(), lambda = 0.3, rho = 0.5, errors = list(),
                   estimands = NULL) {
  tag <- sprintf("rubin_types_mixed_BM_l%s_r%s_mcar0.3_n10_M20_s%d", format(lambda), format(rho), seed)
  list(tag = tag, n = 10L, seed = as.integer(seed), M = 20L, arms = arms, smoke = FALSE, lambda = lambda, rho = rho,
       truth = truth, mask = mk_mask(truth), errors = errors, diag = list(), walls = c(bace_fit = 1),
       cells = NULL, estimands = estimands, disc_cells = disc_cells, disc_estimands = disc_estimands,
       disc_fill = stats::setNames(rep(list(c(bin = 0L, ord = 0L, cat3 = 0L)), length(arms)), arms))
}
put <- function(pool, set, host, x, name = paste0(x$tag, ".rds")) {
  d <- file.path(pool, set, host); dir.create(d, recursive = TRUE, showWarnings = FALSE); saveRDS(x, file.path(d, name))
}
truth_csv <- function(dir) {
  f <- file.path(dir, "truth.csv")
  utils::write.csv(data.frame(n = 10, lambda = 0.3, rho = 0.5, R = 100, R_used = 100, truth = 0.5, mc_se = 0.01,
                              sd_complete = 0.1, liability_ref = 0.8, wall_s = 1), f, row.names = FALSE)
  f
}
# the shared synthetic pool: two datasets (seeds 1, 2) in both sets, a third (seed 3, lambda = 1, no truth row) in freq
build_pool <- function(dir) {
  pool <- file.path(dir, "pool")
  fq <- c("freqA", "freqB"); bc <- c("bace", "bace_chain", "bace_resid")
  # seed 1: freqA covers the truth at its lower bound, freqB does not; the complete-data interval misses
  f1 <- mk_fit(1, fq, cells_df(fq, acc = c(0.6, 0.6, 0.6, 0.5, 0.5, 0.5), brier = c(0.35, 0.35, 0.35, 0.4, 0.4, 0.4)),
               est_df(fq, c(0.6, 0.7), c(0.5, 0.51), c(0.7, 0.9), 0.42, 0.40, 0.45))
  b1 <- mk_fit(1, bc, cells_df(bc, acc = c(0.7, 0.7, 0.7, 0.8, 0.8, 0.8, 0.7, 0.7, 0.7), brier = c(rep(0.2, 3), rep(0.3, 3), rep(0.2, 3))),
               est_df(bc, c(0.55, 0.52, 0.58), c(0.3, 0.3, 0.3), c(0.8, 0.8, 0.8), 0.42, 0.40, 0.45))
  # seed 2: freqA misses, freqB covers; the complete-data interval covers
  f2 <- mk_fit(2, fq, cells_df(fq, acc = c(0.7, 0.7, 0.7, 0.4, 0.4, 0.4), brier = c(0.30, 0.30, 0.30, 0.5, 0.5, 0.5)),
               est_df(fq, c(0.3, 0.45), c(0.2, 0.3), c(0.4, 0.6), 0.45, 0.30, 0.60))
  b2 <- mk_fit(2, bc, cells_df(bc, acc = c(0.9, 0.9, 0.9, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9), brier = c(rep(0.1, 3), rep(0.4, 3), rep(0.1, 3))),
               est_df(bc, c(0.5, 0.5, 0.5), c(0.2, 0.2, 0.2), c(0.9, 0.9, 0.9), 0.45, 0.30, 0.60))
  # seed 3: bin has one class in the complete data (skipped everywhere); freqB's discrete scoring failed
  t3 <- mk_truth(bin = rep("no", 10))
  f3 <- mk_fit(3, fq, cells_df("freqA", acc = c(0.9, 0.9), brier = c(0.1, 0.1), traits = c("ord", "cat3")), NULL,
               truth = t3, lambda = 1, errors = list(freqB_disc_score = "boom"))
  put(pool, "freq", "fir", f1); put(pool, "bace", "nibi", b1); put(pool, "freq", "fir", f2); put(pool, "bace", "rorqual", b2)
  put(pool, "freq", "totoro", f3)
  # the same freq fit of seed 1 again, from a second host and as a _dup copy: both dropped
  put(pool, "freq", "totoro", f1); put(pool, "freq", "totoro", f1, name = paste0(f1$tag, "_dup1.rds"))
  pool
}

testthat::test_that("mode floor: the most frequent observed class in every masked cell, ties fractional", {
  tr <- mk_truth()
  mf <- mode_floor(tr, mk_mask(tr), c("bin", "ord", "cat3"))
  acc <- stats::setNames(mf$accuracy, mf$trait)
  testthat::expect_equal(acc[["bin"]], 2 / 3)                        # "yes" predicted; truth yes, no, yes
  testthat::expect_equal(acc[["cat3"]], (1 / 2 + 0 + 1 / 2) / 3)     # A-B tie; truth A, C, B
  testthat::expect_equal(acc[["ord"]], (1 / 3 + 0 + 1 / 3) / 3)      # L1-L2-L3 tie; truth L1, L4, L3
  testthat::expect_true(all(is.na(mf$brier)) && all(mf$n_cells == 3L))
  t1 <- mk_truth(bin = rep("no", 10))                               # one class in the complete data: nothing to score
  testthat::expect_false("bin" %in% mode_floor(t1, mk_mask(t1), c("bin", "ord", "cat3"))$trait)
})

testthat::test_that("end to end: duplicates dropped, statuses, truth join, coverage, the mode floor once per dataset", {
  dir <- withr::local_tempdir()
  pool <- build_pool(dir)
  testthat::expect_equal(attr(list_pool(pool), "n_dup"), 2L)
  utils::capture.output(aggregate_disc(pool, truth_csv(dir), file.path(dir, "out"), cont_pool = NULL))
  out <- file.path(dir, "out")
  for (f in c("fit_disc_cells.csv.gz", "fit_disc_estimands.csv.gz", "disc_fill.csv", "castor_diag.csv", "agg_disc_cells.csv",
              "agg_disc_cells_l.csv", "agg_disc_down.csv", "agg_disc_paired.csv", "repro.csv", "sanity.txt"))
    testthat::expect_true(file.exists(file.path(out, f)), info = f)
  fc <- utils::read.csv(file.path(out, "fit_disc_cells.csv.gz"), stringsAsFactors = FALSE)
  testthat::expect_equal(anyDuplicated(fc[c("set", "tag", "arm", "trait")]), 0L)
  s3 <- fc[fc$seed == 3, ]
  testthat::expect_setequal(s3$status[s3$trait == "bin"], "one_class")          # freqA, freqB and the mode floor
  testthat::expect_equal(s3$status[s3$arm == "freqB" & s3$trait == "ord"], "missing")
  testthat::expect_match(s3$errors[s3$arm == "freqB" & s3$trait == "ord"], "freqB_disc_score: boom")

  fe <- utils::read.csv(file.path(out, "fit_disc_estimands.csv.gz"), stringsAsFactors = FALSE)
  k <- function(seed, arm, col) fe[[col]][fe$seed == seed & fe$arm == arm]
  testthat::expect_equal(k(1, "freqA", "truth"), 0.5)
  testthat::expect_true(k(1, "freqA", "covered"))                                # lower == truth counts as covered
  testthat::expect_false(k(1, "freqB", "covered"))
  testthat::expect_false(k(1, "freqA", "complete_covered")); testthat::expect_true(k(2, "freqA", "complete_covered"))
  testthat::expect_setequal(fe$status[fe$seed == 3], "one_class")
  testthat::expect_true(all(is.na(fe$truth[fe$seed == 3])))                      # lambda = 1: no truth row

  ad <- utils::read.csv(file.path(out, "agg_disc_down.csv"), stringsAsFactors = FALSE)
  g <- function(arm, col, lambda = 0.3) ad[[col]][ad$arm == arm & ad$lambda == lambda]
  testthat::expect_equal(g("freqA", "coverage"), 0.5); testthat::expect_equal(g("freqA", "coverage_se"), 0.5)
  testthat::expect_equal(g("freqB", "coverage"), 0.5)
  testthat::expect_equal(g("freqA", "bias"), mean(c(0.6, 0.3)) - 0.5)
  testthat::expect_equal(g("freqA", "diff_cd_mean"), mean(c(0.6 - 0.42, 0.3 - 0.45)))
  testthat::expect_equal(g("freqA", "complete_coverage"), 0.5)
  testthat::expect_equal(g("complete", "n_fits"), 2L)                            # once per dataset, not per set
  testthat::expect_equal(g("complete", "coverage"), 0.5); testthat::expect_equal(g("complete", "bias"), mean(c(0.42, 0.45)) - 0.5)
  testthat::expect_equal(g("freqA", "n_excluded", lambda = 1), 1L); testthat::expect_equal(g("freqA", "n_fits", lambda = 1), 0L)

  ac <- utils::read.csv(file.path(out, "agg_disc_cells.csv"), stringsAsFactors = FALSE)
  h <- function(arm, trait, col, lambda = 0.3) ac[[col]][ac$arm == arm & ac$trait == trait & ac$lambda == lambda]
  testthat::expect_equal(h("mode_floor", "bin", "n_fits"), 2L)                   # seeds 1, 2 are in both sets
  testthat::expect_equal(h("mode_floor", "bin", "accuracy"), 2 / 3)
  testthat::expect_equal(h("freqA", "bin", "accuracy"), 0.65); testthat::expect_equal(h("freqA", "bin", "accuracy_se"), 0.05)
  testthat::expect_equal(h("freqA", "bin", "n_fits_trait_absent", lambda = 1), 1L)
  testthat::expect_equal(h("mode_floor", "bin", "n_fits_trait_absent", lambda = 1), 1L)
  testthat::expect_equal(h("freqB", "ord", "n_fits_missing", lambda = 1), 1L)
  testthat::expect_true(is.na(h("freqA", "bin", "mae_class"))); testthat::expect_equal(h("freqA", "ord", "mae_class"), 0.5)
  acl <- utils::read.csv(file.path(out, "agg_disc_cells_l.csv"), stringsAsFactors = FALSE)
  testthat::expect_false("rho" %in% names(acl))
})

testthat::test_that("paired differences are taken per dataset between arms of the two sets", {
  dir <- withr::local_tempdir()
  pool <- build_pool(dir)
  utils::capture.output(res <- aggregate_disc(pool, truth_csv(dir), file.path(dir, "out")))
  p <- res$agg_paired
  q <- function(contrast, measure, trait = "bin", lambda = 0.3)
    p[p$contrast == contrast & p$measure == measure & p$trait == trait & p$lambda == lambda, ]
  a <- q("bace_chain - freqA", "accuracy")                  # 0.8 - 0.6 and 0.5 - 0.7
  testthat::expect_equal(a$diff, 0); testthat::expect_equal(a$se, 0.2); testthat::expect_equal(a$datasets, 2L)
  b <- q("bace_chain - freqA", "brier")                     # 0.3 - 0.35 and 0.4 - 0.3
  testthat::expect_equal(b$diff, 0.025); testthat::expect_equal(b$se, stats::sd(c(-0.05, 0.1)) / sqrt(2))
  testthat::expect_equal(q("bace - freqA", "accuracy")$diff, mean(c(0.7 - 0.6, 0.9 - 0.7)))
  testthat::expect_equal(q("bace_chain - bace", "accuracy")$diff, mean(c(0.8 - 0.7, 0.5 - 0.9)))
  testthat::expect_equal(q("freqA - freqB", "accuracy")$diff, mean(c(0.6 - 0.5, 0.7 - 0.4)))
  testthat::expect_equal(nrow(q("freqA - freqB", "accuracy", trait = "bin", lambda = 1)), 0L)   # bin absent at seed 3
  testthat::expect_equal(nrow(q("freqA - freqB", "accuracy", trait = "ord", lambda = 1)), 0L)   # freqB missing there
})

testthat::test_that("reproduction check compares continuous rows of the same set and tag", {
  dir <- withr::local_tempdir()
  pool <- build_pool(dir)
  est <- function(delta) data.frame(arm = c("complete", "bace"), estimand = "slope", estimate = c(0.5, 0.5 + delta),
                                    se = 0.1, stringsAsFactors = FALSE)
  b1 <- readRDS(file.path(pool, "bace", "nibi", "rubin_types_mixed_BM_l0.3_r0.5_mcar0.3_n10_M20_s1.rds"))
  b1$estimands <- est(0.02); saveRDS(b1, file.path(pool, "bace", "nibi", paste0(b1$tag, ".rds")))
  cont <- file.path(dir, "cont"); b1$estimands <- est(0); put(cont, "bace", "nibi", b1)
  utils::capture.output(res <- aggregate_disc(pool, truth_csv(dir), file.path(dir, "out"), cont_pool = cont))
  rp <- res$repro
  all_row <- rp[rp$set == "bace" & rp$host == "all" & rp$arm == "all", ]
  testthat::expect_equal(all_row$n_fits, 1L); testthat::expect_equal(all_row$n_fits_identical, 0L)
  testthat::expect_equal(all_row$max_abs_estimate, 0.02)
  testthat::expect_equal(rp$n_fits_identical[rp$host == "all" & rp$arm == "complete"], 1L)
  testthat::expect_equal(nrow(rp[rp$host == "nibi" & rp$cont_host == "nibi" & rp$arm == "all", ]), 1L)  # host pair row
})

testthat::test_that("malformed files and datasets that differ between sets stop the run, naming the files", {
  dir <- withr::local_tempdir()
  pool <- build_pool(dir)
  bad <- mk_fit(9, c("freqA", "freqB"), NULL, NULL); bad$disc_fill <- NULL
  put(pool, "freq", "fir", bad)
  testthat::expect_error(utils::capture.output(aggregate_disc(pool, truth_csv(dir), file.path(dir, "out"))),
                         "malformed.*_s9\\.rds: no discrete results")
  dir2 <- withr::local_tempdir()
  pool2 <- build_pool(dir2)
  f <- file.path(pool2, "freq", "fir", "rubin_types_mixed_BM_l0.3_r0.5_mcar0.3_n10_M20_s2.rds")
  x <- readRDS(f); x$mask[1, "bin"] <- TRUE; saveRDS(x, f)                     # the freq copy of seed 2 differs
  testthat::expect_error(utils::capture.output(aggregate_disc(pool2, truth_csv(dir2), file.path(dir2, "out"))),
                         "differ between the bace and freq files.*_s2")
})

testthat::test_that("disc_arms_of adds the equal-rates arms to frequentist fits only (no phantom '_er')", {
  testthat::expect_identical(disc_arms_of(list(arms = c("bace", "bace_chain", "bace_resid"))), c("bace", "bace_chain", "bace_resid"))
  testthat::expect_identical(disc_arms_of(list(arms = c("freqA", "freqB"))), c("freqA", "freqA_er", "freqB", "freqB_er"))
})
