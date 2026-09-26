# script/rubin_cell.R
#
# One cell of the freq-vs-BACE multiple-imputation study (docs/dev-log/arc/2026-09-24-rubin-freq-bace-plan.md).
# Builds the cell with the v1 make_cell() (script/campaign_gnn_off_lib.R), produces M completed datasets per
# arm, and scores everything under Rubin's rules:
#   per cell   : mean of the M imputations, interval mean +/- t_{M-1} sqrt((1 + 1/M) B)
#   slope      : PGLS slope of c2 on c1 (corPagel), pooled W + (1 + 1/M) B, Barnard-Rubin df
#   correlation: GLS-whitened phylogenetic correlation of c1 and c2, pooled on Fisher's z
# In this DGP the population slope and the population phylogenetic correlation both equal rho, so coverage is
# scored against rho. Each cell also carries a "complete" row (the same estimands on the complete data, same df
# convention), because the complete-data coverage itself falls short of 0.95 at small n (Meng review N3).
#
# Arms: freqA (parametric bootstrap then joint conditional draw), freqB (joint conditional draws at one fixed
# Rphylopars-lambda fit), bace (as shipped: the installed build draws one posterior iteration plus a residual per
# final run, but every final run starts from the same converged dataset; Meng review B2), bace_resid (BACE plus
# a second, post hoc residual draw: a negative control, not a fix), bace_chain (BACE's own final step run M times,
# each run starting from the previous draw instead of the shared converged dataset; paired with bace on one fit).
# Each arm reseeds with seed + a fixed offset, so an arm's draws do not depend on which arms ran before it.
#
# Usage:
#   Rscript script/rubin_cell.R --n 100 --seed 1 --lambda 0.7 --rho 0.5 --out results/ \
#           [--arms freqA,freqB,bace,bace_chain,bace_resid] [--M 20] [--miss mcar --frac 0.30] \
#           [--thresholds fixed] [--no_driver] [--save_imp] [--discrete] \
#           [--bace_nitt 50000 --bace_burnin 10000 --bace_thin 25 --bace_runs 10] [--smoke]
#
# Output: rubin_types_mixed_BM_l<lambda>_r<rho>_<miss><frac>_n<n>_M<M>_s<seed>[_smoke].rds
# Resumable: a cell whose rds already exists is skipped.

suppressPackageStartupMessages({ library(ape) })
RNGkind("L'Ecuyer-CMRG")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1L], "--")) return(TRUE)
  args[i + 1L]
}
n        <- as.integer(get_arg("--n", 100L))
seed     <- as.integer(get_arg("--seed", 1L))
out      <- get_arg("--out", "results")
arms     <- strsplit(get_arg("--arms", "freqA,freqB,bace,bace_chain,bace_resid"), ",")[[1]]
M        <- as.integer(get_arg("--M", 20L))
lambda   <- as.numeric(get_arg("--lambda", 0.7))
rho      <- as.numeric(get_arg("--rho", 0.5))
miss     <- get_arg("--miss", "mcar")
frac     <- as.numeric(get_arg("--frac", 0.30))
# v1's core stage ran every cell with --driver --thresholds fixed (script/campaign_sim_totoro.sh); match it so
# the MI results sit beside v1's. --no_driver / --thresholds sample reproduce make_cell()'s defaults instead.
thresholds <- get_arg("--thresholds", "fixed")
driver     <- !isTRUE(get_arg("--no_driver", FALSE))
smoke    <- isTRUE(get_arg("--smoke", FALSE))
# --save_imp keeps each arm's M completed datasets in the rds (res$imputations), so discrete traits can be
# scored later; off by default, which leaves the continuous campaign's output unchanged.
save_imp <- isTRUE(get_arg("--save_imp", FALSE))
# --discrete adds the discrete traits (script/rubin_discrete.R): castor Mk with flexible rates fills the freq A
# (bootstrapped rates) and freq B (fixed rates) datasets after their continuous draws, under their own seeds (505,
# 606), so the continuous draws and scores do not change. Two comparison arms, freqA_er and freqB_er (seeds 707,
# 808), fill the same continuous draws with simulation v1's equal-rates models and are scored on discrete traits
# only. Every arm is scored on bin, ord, cat3 and on c1 ~ bin.
discrete <- isTRUE(get_arg("--discrete", FALSE))
bace_nitt   <- as.integer(get_arg("--bace_nitt", 50000L))
bace_burnin <- as.integer(get_arg("--bace_burnin", 10000L))
bace_thin   <- as.integer(get_arg("--bace_thin", 25L))
bace_runs   <- as.integer(get_arg("--bace_runs", 10L))
if (smoke) { bace_nitt <- 6000L; bace_burnin <- 1000L; bace_thin <- 5L; bace_runs <- 3L }

dir.create(out, showWarnings = FALSE, recursive = TRUE)
tag <- sprintf("rubin_types_mixed_BM_l%s_r%s_%s%s_n%d_M%d_s%d", format(lambda), format(rho), miss,
               format(frac), n, M, seed)
out_path <- file.path(out, paste0(tag, if (smoke) "_smoke" else "", ".rds"))
log_line <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))
if (file.exists(out_path)) { log_line("cell %s: rds exists, skipping (resume)", tag); quit(save = "no", status = 0) }

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
source(file.path(here, "campaign_gnn_off_lib.R"))
source(file.path(here, "rubin_lib.R"))
source(file.path(here, "rubin_freq.R"))
source(file.path(here, "rubin_bace.R"))
if (discrete) source(file.path(here, "rubin_discrete.R"))

git_hash <- tryCatch({
  h <- system2("git", c("-C", here, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)
  if (length(h) == 1 && grepl("^[0-9a-f]{40}$", h)) h else NA_character_   # a repo without commits prints "HEAD"
}, error = function(e) NA_character_)
# cluster copies are rsynced, not git checkouts: fall back to the commit recorded at sync time (script/DISC_COMMIT)
if (is.na(git_hash) && file.exists(file.path(here, "DISC_COMMIT"))) git_hash <- readLines(file.path(here, "DISC_COMMIT"))[1]

cell <- make_cell("types_mixed", n, seed, miss_frac = frac, miss = miss, lambda = lambda, rho = rho,
                  thresholds = thresholds, driver = driver)
truth <- cell$truth; tree <- cell$tree; mask <- cell$mask; df_miss <- cell$df_miss
log_line("cell %s: n=%d masked=%d realised_frac=%.3f", tag, nrow(truth), sum(mask), cell$realised_frac)

# ---- scoring -----------------------------------------------------------------------------------
num_traits <- names(truth)[vapply(truth, is.numeric, logical(1))]

score_cells <- function(arm, sets) {
  rows <- list()
  for (v in num_traits) {
    idx <- which(mask[, v]); if (!length(idx)) next
    draws <- t(vapply(sets, function(s) as.numeric(s[[v]][idx]), numeric(length(idx))))
    if (length(idx) == 1L) draws <- matrix(draws, ncol = 1L)
    if (all(is.na(draws))) next                 # trait not imputed by this arm (e.g. outside the freq block)
    n_na <- sum(is.na(draws))                   # a partial failure is recorded, never hidden (Meng N14)
    if (n_na) draws <- draws[stats::complete.cases(draws), , drop = FALSE]
    sc <- score_scale(v, draws); draws <- sc$x          # prp scored on the logit scale (Meng N10)
    ci <- rubin_cell_intervals(draws)
    tr <- score_scale(v, truth[[v]][idx])$x
    sdt <- stats::sd(score_scale(v, truth[[v]][!mask[, v]])$x)
    lo <- ci$lower; hi <- ci$upper; a <- 0.05
    rows[[length(rows) + 1L]] <- data.frame(
      arm = arm, trait = v, n_cells = length(idx),
      zRMSE = sqrt(mean(((tr - ci$mean) / sdt)^2)),
      coverage = mean(tr >= lo & tr <= hi),
      width = mean((hi - lo) / sdt),
      interval_score = mean(((hi - lo) + (2 / a) * (lo - tr) * (tr < lo) + (2 / a) * (tr - hi) * (tr > hi)) / sdt),
      frac_B0 = mean(ci$B == 0), n_na = n_na, m_used = nrow(draws), n_oob = sc$n_oob,
      scale = if (v == "prp") "logit" else "raw")
  }
  do.call(rbind, rows)
}

# The fast eigenbasis estimator (est_pgls_slope_fast; equal to nlme::gls REML with corPagel to ~1e-5, gate
# test-lib-fast.R) replaces gls here: every dataset in a cell shares the tree, so one eigendecomposition
# serves them all, and a bounded-lambda fit at n = 1000 drops from about a minute to milliseconds.
eig <- pagel_eigen(tree, rownames(truth))
ref_slope <- est_pgls_slope_fast(truth, tree, eig = eig)
ref_cor   <- est_phylo_cor(truth, tree, lambda = ref_slope$lambda_hat)

score_estimands <- function(arm, sets) {
  sl <- lapply(sets, est_pgls_slope_fast, tree = tree, eig = eig)
  ok <- vapply(sl, function(s) is.finite(s$estimate) && is.finite(s$variance), logical(1))
  # the correlation reuses each dataset's slope-fit lambda instead of refitting (Meng N4)
  co <- Map(function(s, fit) est_phylo_cor(s, tree, lambda = fit$lambda_hat), sets, sl)
  okc <- vapply(co, function(s) is.finite(s$z), logical(1))
  if (sum(ok) < 2L || sum(okc) < 2L) stop(sprintf("%s: fewer than 2 analysable imputations", arm))
  ps <- rubin_pool(vapply(sl[ok], `[[`, numeric(1), "estimate"), vapply(sl[ok], `[[`, numeric(1), "variance"),
                   df_com = n - 2)
  pc <- pool_cor(vapply(co[okc], `[[`, numeric(1), "z"), n)
  rbind(
    data.frame(arm = arm, estimand = "slope", estimate = ps$estimate, se = ps$se, lower = ps$lower,
               upper = ps$upper, df = ps$df, fmi = ps$fmi, m_ok = sum(ok), truth = rho,
               complete_data = ref_slope$estimate, covered = ps$lower <= rho & rho <= ps$upper),
    # se, df and fmi of the correlation are on Fisher's z scale; estimate and CI are back-transformed
    data.frame(arm = arm, estimand = "cor", estimate = pc$r, se = pc$z_pool$se, lower = pc$lower,
               upper = pc$upper, df = pc$z_pool$df, fmi = pc$z_pool$fmi, m_ok = sum(okc), truth = rho,
               complete_data = ref_cor$r, covered = pc$lower <= rho & rho <= pc$upper))
}

# ---- arms ----------------------------------------------------------------------------------------
cells_tab <- list(); est_tab <- list(); walls <- list(); errors <- list(); diag <- list(); imps <- list()
disc_cells_tab <- list(); disc_est_tab <- list(); disc_fill <- list(); disc_est_status <- list()
L5 <- if (discrete) cell$L[rownames(truth), 5] else NULL   # bin's liability, for the per-dataset c1 ~ bin target
run_arm <- function(arm, expr) {
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(expr, error = function(e) e)
  walls[[arm]] <<- proc.time()[["elapsed"]] - t0
  if (inherits(res, "error")) { errors[[arm]] <<- conditionMessage(res); log_line("arm %s ERROR: %s", arm, conditionMessage(res)); return(NULL) }
  log_line("arm %s done in %.1fs", arm, walls[[arm]])
  res
}

# Scoring is guarded per arm (Meng B1): a failed draw or estimand records an error and the cell still writes its
# rds. NULL datasets (freq A draws whose refits failed twice) are dropped and counted.
score_disc <- function(arm, sets) {
  f <- fill_degenerate(sets, df_miss, disc_traits_of(truth)); disc_fill[[arm]] <<- f$n_filled
  dc <- tryCatch(score_discrete(arm, f$sets, truth, mask), error = function(e) e)
  if (inherits(dc, "error")) {
    errors[[paste0(arm, "_disc_score")]] <<- conditionMessage(dc)
    log_line("arm %s discrete scoring ERROR: %s", arm, conditionMessage(dc))
  } else disc_cells_tab[[arm]] <<- dc
  de <- tryCatch(score_discrete_estimand(arm, f$sets, truth, tree, eig, n, L5 = L5, lambda = lambda, rho = rho),
                 error = function(e) e)
  if (inherits(de, "error")) {
    errors[[paste0(arm, "_disc_est")]] <<- conditionMessage(de); disc_est_status[[arm]] <<- "error"
    log_line("arm %s c1 ~ bin ERROR: %s", arm, conditionMessage(de))
  } else if (is.null(de)) {
    disc_est_status[[arm]] <<- "undefined"      # bin one class in the complete data, or constant in every dataset
  } else { disc_est_tab[[arm]] <<- de; disc_est_status[[arm]] <<- "scored" }
}
score_arm_sets <- function(arm, sets, continuous = TRUE) {
  sets <- Filter(Negate(is.null), sets)
  if (save_imp) imps[[arm]] <<- sets
  if (continuous) {
    res <- tryCatch(list(cells = score_cells(arm, sets), est = score_estimands(arm, sets)),
                    error = function(e) e)
    if (inherits(res, "error")) {
      errors[[paste0(arm, "_score")]] <<- conditionMessage(res)
      log_line("arm %s scoring ERROR: %s", arm, conditionMessage(res))
    } else { cells_tab[[arm]] <<- res$cells; est_tab[[arm]] <<- res$est }
  }
  if (discrete) score_disc(arm, sets)
}
castor_fill <- function(arm, sets, offset, proper, models = CASTOR_FLEX) {
  if (!discrete) return(sets)
  arm_seed(offset)
  cf <- run_arm(paste0(arm, "_castor"), mi_castor(cell, length(sets), proper = proper, base_sets = sets, models = models))
  if (is.null(cf)) return(sets)
  diag[[paste0(arm, "_castor")]] <<- cf$diag
  cf$datasets
}
arm_seed <- function(offset) set.seed(seed + offset)

est_tab$complete <- tryCatch(score_estimands("complete", list(truth, truth)), error = function(e) NULL)
if (!is.null(est_tab$complete)) {
  # two copies of the truth give B = 0, so the pooled interval is the complete-data interval with the same
  # df convention as the MI arms; FMI is 0 by construction
  est_tab$complete$m_ok <- 1L
  est_tab$complete$fmi <- 0   # mice's small-sample fmi formula gives 2/(df+3) at B = 0, not 0
}

if ("freqA" %in% arms) {
  arm_seed(101L); a <- run_arm("freqA", mi_freq_A(cell, M))
  if (!is.null(a)) {
    a_sets <- Filter(Negate(is.null), a$datasets)
    score_arm_sets("freqA", castor_fill("freqA", a_sets, 505L, proper = TRUE))
    if (discrete) score_arm_sets("freqA_er", castor_fill("freqA_er", a_sets, 707L, proper = TRUE, models = CASTOR_V1),
                                 continuous = FALSE)
    diag$freqA <- list(n_fail = a$n_fail, n_degenerate = a$n_degenerate, m_used = sum(!vapply(a$datasets, is.null, logical(1))),
                       lambda_star = vapply(Filter(Negate(is.null), a$pars_star), function(p) p$lambda, numeric(1)))
  }
}
if ("freqB" %in% arms) {
  arm_seed(202L); b <- run_arm("freqB", mi_freq_B(cell, M))
  if (!is.null(b)) {
    b_sets <- Filter(Negate(is.null), b$datasets)
    score_arm_sets("freqB", castor_fill("freqB", b_sets, 606L, proper = FALSE))
    if (discrete) score_arm_sets("freqB_er", castor_fill("freqB_er", b_sets, 808L, proper = FALSE, models = CASTOR_V1),
                                 continuous = FALSE)
  }
}
if (any(c("bace", "bace_resid", "bace_chain") %in% arms)) {
  arm_seed(303L)
  fb <- run_arm("bace_fit", fit_bace_mi(cell, M = M, nitt = bace_nitt, burnin = bace_burnin, thin = bace_thin,
                                        runs = bace_runs))
  if (!is.null(fb)) {
    diag$bace <- fb$diag
    if ("bace" %in% arms) score_arm_sets("bace", mi_bace_shipped(fb$outb, df_miss))
    if ("bace_resid" %in% arms) score_arm_sets("bace_resid", mi_bace_resid(fb$outb, df_miss, seed = seed))
    if ("bace_chain" %in% arms) {
      arm_seed(404L); ch <- run_arm("bace_chain", mi_bace_chain(fb, df_miss, M = M))
      if (!is.null(ch)) score_arm_sets("bace_chain", ch$datasets)
    }
  }
}

cells_df <- do.call(rbind, cells_tab); est_df <- do.call(rbind, est_tab)
res <- list(tag = tag, n = n, seed = seed, M = M, arms = arms, smoke = smoke, lambda = lambda, rho = rho,
            miss = miss, miss_frac = frac, realised_frac = cell$realised_frac,
            thresholds = thresholds, driver = driver,
            bace = c(nitt = bace_nitt, burnin = bace_burnin, thin = bace_thin, runs = bace_runs),
            cells = cells_df, estimands = est_df,
            reference = list(slope = ref_slope, cor = ref_cor),
            walls = unlist(walls), errors = errors, diag = diag, truth = truth, mask = mask,
            imputations = if (save_imp) imps else NULL,
            disc_cells = if (discrete) do.call(rbind, disc_cells_tab) else NULL,
            disc_estimands = if (discrete) do.call(rbind, disc_est_tab) else NULL,
            disc_fill = if (discrete) disc_fill else NULL,
            disc_est_status = if (discrete) unlist(disc_est_status) else NULL,
            RNGkind = RNGkind(), sessionInfo = utils::sessionInfo(), git_hash = git_hash,
            host = Sys.info()[["nodename"]], time = Sys.time())
saveRDS(res, out_path)
print(cells_df, digits = 3); print(est_df, digits = 3); print(round(unlist(walls), 1))
if (length(errors)) { cat("ERRORS:\n"); print(errors) }
