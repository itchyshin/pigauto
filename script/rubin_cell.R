# script/rubin_cell.R
#
# One cell of the freq-vs-BACE multiple-imputation study (docs/dev-log/arc/2026-09-24-rubin-freq-bace-plan.md).
# Builds the cell with the v1 make_cell() (script/campaign_gnn_off_lib.R), produces M completed datasets per
# arm, and scores everything under Rubin's rules:
#   per cell   : mean of the M imputations, interval mean +/- t_{M-1} sqrt((1 + 1/M) B)
#   slope      : PGLS slope of c2 on c1 (corPagel), pooled W + (1 + 1/M) B, Barnard-Rubin df
#   correlation: GLS-whitened phylogenetic correlation of c1 and c2, pooled on Fisher's z
# In this DGP the population slope and the population phylogenetic correlation both equal rho, so coverage is
# scored against rho; bias is also reported against the complete-data estimate ("ref_full").
#
# Arms: freqA (parametric bootstrap then joint conditional draw), freqB (joint conditional draws at one fixed
# Rphylopars-lambda fit), bace (as shipped: continuous imputations are posterior means), bace_resid (BACE plus a
# post hoc residual draw from each final fit's units posterior).
#
# Usage:
#   Rscript script/rubin_cell.R --n 100 --seed 1 --lambda 0.7 --rho 0.5 --out results/ \
#           [--arms freqA,freqB,bace,bace_resid] [--M 20] [--miss mcar --frac 0.30] \
#           [--thresholds fixed] [--no_driver] \
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
arms     <- strsplit(get_arg("--arms", "freqA,freqB,bace,bace_resid"), ",")[[1]]
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

git_hash <- tryCatch({
  h <- system2("git", c("-C", here, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)
  if (length(h) == 1 && nzchar(h)) h else NA_character_
}, error = function(e) NA_character_)

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
    if (anyNA(draws)) next                      # trait not imputed by this arm (e.g. outside the freq block)
    ci <- rubin_cell_intervals(draws)
    tr <- truth[[v]][idx]
    sdt <- stats::sd(truth[[v]][!mask[, v]])
    lo <- ci$lower; hi <- ci$upper; a <- 0.05
    rows[[length(rows) + 1L]] <- data.frame(
      arm = arm, trait = v, n_cells = length(idx),
      zRMSE = sqrt(mean(((tr - ci$mean) / sdt)^2)),
      coverage = mean(tr >= lo & tr <= hi),
      width = mean((hi - lo) / sdt),
      interval_score = mean(((hi - lo) + (2 / a) * (lo - tr) * (tr < lo) + (2 / a) * (tr - hi) * (tr > hi)) / sdt),
      frac_B0 = mean(ci$B == 0))
  }
  do.call(rbind, rows)
}

ref_slope <- est_pgls_slope(truth, tree)
ref_cor   <- est_phylo_cor(truth, tree)

score_estimands <- function(arm, sets) {
  sl <- lapply(sets, est_pgls_slope, tree = tree)
  ok <- vapply(sl, function(s) is.finite(s$estimate) && is.finite(s$variance), logical(1))
  ps <- rubin_pool(vapply(sl[ok], `[[`, numeric(1), "estimate"), vapply(sl[ok], `[[`, numeric(1), "variance"),
                   df_com = n - 2)
  co <- lapply(sets, est_phylo_cor, tree = tree)
  okc <- vapply(co, function(s) is.finite(s$z), logical(1))
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
cells_tab <- list(); est_tab <- list(); walls <- list(); errors <- list(); diag <- list()
run_arm <- function(arm, expr) {
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(expr, error = function(e) e)
  walls[[arm]] <<- proc.time()[["elapsed"]] - t0
  if (inherits(res, "error")) { errors[[arm]] <<- conditionMessage(res); log_line("arm %s ERROR: %s", arm, conditionMessage(res)); return(NULL) }
  log_line("arm %s done in %.1fs", arm, walls[[arm]])
  res
}

if ("freqA" %in% arms) {
  a <- run_arm("freqA", mi_freq_A(cell, M))
  if (!is.null(a)) {
    cells_tab$freqA <- score_cells("freqA", a$datasets); est_tab$freqA <- score_estimands("freqA", a$datasets)
    diag$freqA <- list(n_fail = a$n_fail, lambda_star = vapply(a$pars_star, function(p) p$lambda, numeric(1)))
  }
}
if ("freqB" %in% arms) {
  b <- run_arm("freqB", mi_freq_B(cell, M))
  if (!is.null(b)) { cells_tab$freqB <- score_cells("freqB", b$datasets); est_tab$freqB <- score_estimands("freqB", b$datasets) }
}
if (any(c("bace", "bace_resid") %in% arms)) {
  fb <- run_arm("bace_fit", fit_bace_mi(cell, M = M, nitt = bace_nitt, burnin = bace_burnin, thin = bace_thin,
                                        runs = bace_runs))
  if (!is.null(fb)) {
    diag$bace <- fb$diag
    if ("bace" %in% arms) {
      s <- mi_bace_shipped(fb$outb, df_miss)
      cells_tab$bace <- score_cells("bace", s); est_tab$bace <- score_estimands("bace", s)
    }
    if ("bace_resid" %in% arms) {
      s <- mi_bace_resid(fb$outb, df_miss, seed = seed)
      cells_tab$bace_resid <- score_cells("bace_resid", s); est_tab$bace_resid <- score_estimands("bace_resid", s)
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
            RNGkind = RNGkind(), sessionInfo = utils::sessionInfo(), git_hash = git_hash,
            host = Sys.info()[["nodename"]], time = Sys.time())
saveRDS(res, out_path)
print(cells_df, digits = 3); print(est_df, digits = 3); print(round(unlist(walls), 1))
if (length(errors)) { cat("ERRORS:\n"); print(errors) }
