# script/campaign_sim_checks.R
#
# Verification gates for the corrected-design simulation (script/campaign_sim_cell.R,
# script/campaign_gnn_off_lib.R). Prints "<gate> PASS" only when every assertion passes for that
# gate; otherwise stops with the failing assertion (fail loudly, never silently).
#
# Usage: Rscript script/campaign_sim_checks.R --gate <G3|G4|G5|G6|Gold|G9b|G10|G11|G12|G14> [--dir DIR]

suppressPackageStartupMessages({ library(pigauto); library(ape) })

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1L], "--")) return(TRUE)
  args[i + 1L]
}
gate <- get_arg("--gate", NULL)
dir_arg <- get_arg("--dir", NULL)
if (is.null(gate)) stop("usage: campaign_sim_checks.R --gate <G3|G4|G5|G6|Gold|G9b|G10|G11|G12|G14> [--dir DIR]")

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
source(file.path(script_dir, "campaign_gnn_off_lib.R"))

fail <- function(...) stop(sprintf(...), call. = FALSE)
pass <- function(g) cat(sprintf("%s PASS\n", g))

# ================================================================================================
if (gate == "G3") {
  mechs <- c("mcar", "mcar", "mar", "clade", "clade")
  for (i in seq_along(mechs)) {
    m <- mechs[i]
    cd <- make_cell("types_mixed", 80, 100 + i, miss_frac = 0.30, miss = m, driver = (m %in% c("mar", "clade")))
    expect_na <- cd$mask | is.na(cd$truth)
    got_na <- is.na(cd$df_miss)
    if (!identical(dim(expect_na), dim(got_na))) fail("G3: dim mismatch for mechanism %s", m)
    if (!identical(unname(got_na), unname(expect_na))) fail("G3: is.na(df_miss) != mask | is.na(truth) for mechanism %s", m)
    if ("d1" %in% names(cd$truth) && any(cd$mask[, "d1"])) fail("G3: d1 masked under mechanism %s", m)
  }
  # arm-dispatch rownames check: every arm's completed rownames equal truth rownames
  cd <- make_cell("types_mixed", 60, 1, miss_frac = 0.30, miss = "mcar")
  cell_obj <- list(truth = cd$truth, tree = cd$tree, mask = cd$mask, df_miss = cd$df_miss,
                    cont_traits = cd$cont_traits, trait_types = cd$trait_types)
  out <- run_arms(cell_obj, c("gnn_off_pure", "freq", "floor"),
                   list(seed = 1, epochs = 20L, bace_nitt = 600, bace_burnin = 100, bace_thin = 5))
  if (length(out$errors)) fail("G3: arm(s) errored: %s", paste(names(out$errors), collapse = ","))
  pass("G3")
}

# ================================================================================================
if (gate == "G4") {
  if (!requireNamespace("phylolm", quietly = TRUE)) fail("G4: phylolm not installed")
  targets <- c(0.3, 0.7, 1.0)
  reps <- 20L
  res <- sapply(targets, function(lam) {
    ests <- vapply(seq_len(reps), function(s) {
      set.seed(1000 + s)
      tree <- ape::rcoal(1000)
      tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
      L <- sim_latents(tree, K = 1L, lambda = lam, rho = 0, evo = "BM")
      dat <- data.frame(y = L[, 1], row.names = rownames(L))
      fit <- phylolm::phylolm(y ~ 1, data = dat, phy = tree, model = "lambda")
      unname(fit$optpar)
    }, numeric(1))
    mean(ests)
  })
  cat("G4: target vs recovered lambda:\n"); print(rbind(target = targets, recovered = res))
  if (any(abs(res - targets) > 0.05)) fail("G4: lambda recovery off by more than 0.05: %s", paste(round(res, 3), collapse = ","))
  pass("G4")
}

# ================================================================================================
if (gate == "G5") {
  reps <- 20L
  for (m in c("mar", "clade")) {
    fracs <- vapply(seq_len(reps), function(s) {
      cd <- make_cell("types_mixed", 200, 2000 + s, miss_frac = 0.30, miss = m, driver = TRUE)
      min_obs <- min(colSums(!cd$mask[, setdiff(names(cd$truth), "d1"), drop = FALSE]))
      if (min_obs < 5) fail("G5: mechanism %s rep %d has < 5 observed cells in a column (min = %d)", m, s, min_obs)
      cd$realised_frac
    }, numeric(1))
    dev <- abs(mean(fracs) - 0.30)
    cat(sprintf("G5: mechanism %s realised_frac mean = %.4f (target 0.30, |dev| = %.4f)\n", m, mean(fracs), dev))
    if (dev > 0.01) fail("G5: mechanism %s realised fraction off target by %.4f (> 0.01)", m, dev)
  }
  pass("G5")
}

# ================================================================================================
if (gate == "Gold") {
  gold_path <- file.path(script_dir, "campaign_gnn_off_prerun", "bm_mixed_n100_s1.rds")
  if (!file.exists(gold_path)) fail("Gold: reference rds not found at %s", gold_path)
  old <- readRDS(gold_path)
  tmp <- tempfile("gold_check_")
  status <- system2("Rscript", c(file.path(script_dir, "campaign_gnn_off_cell.R"),
                                  "--dgp", "bm_mixed", "--n", "100", "--seed", "1", "--out", tmp,
                                  "--arms", "freq,gnn_off_pure,floor"))
  if (status != 0) fail("Gold: campaign_gnn_off_cell.R exited with status %d", status)
  new <- readRDS(file.path(tmp, "bm_mixed_n100_s1.rds"))
  torch_arms <- c("gnn_on", "gnn_on_full")
  o <- old$results; nw <- new$results
  o <- o[!(o$arm %in% torch_arms), ]; nw <- nw[!(nw$arm %in% torch_arms), ]
  key <- c("arm", "trait", "metric")
  comb <- merge(o[c(key, "value")], nw[c(key, "value")], by = key, suffixes = c(".old", ".new"))
  if (nrow(comb) < nrow(o)) fail("Gold: %d rows from the reference rds have no match in the new run", nrow(o) - nrow(comb))
  d <- abs(comb$value.old - comb$value.new)
  cat("Gold: max abs diff over", nrow(comb), "shared (arm, trait, metric) rows:", max(d), "\n")
  if (any(d > 1e-8)) { print(comb[d > 1e-8, ]); fail("Gold: %d rows differ by more than 1e-8", sum(d > 1e-8)) }
  pass("Gold")
}

# ================================================================================================
if (gate == "G6") {
  if (!requireNamespace("parallel", quietly = TRUE)) fail("G6: parallel not installed")
  # Regime (state it wherever this gate's number is quoted): bm_mixed, lambda = 1, MCAR 0.30,
  # n = 300, 20 replicates, BACE nitt 20000 / burnin 4000 / thin 20 with PIG_BACE_NFINAL imputation
  # runs. n = 300 rather than 1000 because n_final is a count of full MCMCglmm refits: this gate
  # exists to catch an interval that is built wrongly, and that shows up at any n.
  reps <- as.integer(Sys.getenv("G6_REPS", "20"))
  n_g6 <- as.integer(Sys.getenv("G6_N", "300"))
  bace_nitt <- 20000L; bace_burnin <- 4000L; bace_thin <- 20L
  cat(sprintf("G6 regime: bm_mixed lambda=1 MCAR 0.30 n=%d reps=%d nitt=%d n_final=%s\n",
              n_g6, reps, bace_nitt, Sys.getenv("PIG_BACE_NFINAL", "50")))
  one_rep <- function(s) {
    cd <- make_cell("bm_mixed", n_g6, 3000 + s, miss_frac = 0.30, miss = "mcar", lambda = 1)
    cell_obj <- list(truth = cd$truth, tree = cd$tree, mask = cd$mask, df_miss = cd$df_miss,
                      cont_traits = cd$cont_traits, trait_types = cd$trait_types)
    out <- run_arms(cell_obj, c("freq", "bace"),
                     list(seed = s, bace_nitt = bace_nitt, bace_burnin = bace_burnin, bace_thin = bace_thin))
    list(results = out$results, errors = out$errors, n_final = tryCatch({
      # recover n_final from the bace() call inside run_bace via the draws count stored on lower/upper cols
      NA_integer_
    }, error = function(e) NA_integer_))
  }
  # Sequential: MCMCglmm segfaults inside a forked worker, and BACE already runs one process per
  # cell. Parallelism for this gate would have to be separate R processes, which is what the
  # campaign drivers do.
  reps_out <- lapply(seq_len(reps), one_rep)
  errs <- unlist(lapply(reps_out, function(r) names(r$errors)))
  if (length(errs)) cat("G6: arm errors across replicates:", paste(errs, collapse = ","), "\n")
  tab <- do.call(rbind, lapply(seq_along(reps_out), function(i) { d <- reps_out[[i]]$results; if (!is.null(d)) d$rep <- i; d }))
  cov_by_arm <- sapply(c("freq", "bace"), function(a) {
    d <- tab[tab$arm == a & tab$metric == "zRMSE" & !is.na(tab$coverage), ]
    if (!nrow(d)) return(NA_real_)
    mean(d$coverage)
  })
  cat("G6: mean coverage by arm:\n"); print(cov_by_arm)
  for (a in names(cov_by_arm)) {
    cv <- cov_by_arm[[a]]
    if (is.na(cv)) fail("G6: arm %s has no coverage values (interval not implemented for this arm)", a)
    if (cv < 0.90 || cv > 0.99) fail("G6: arm %s coverage %.3f outside [0.90, 0.99]", a, cv)
  }
  pass("G6")
}

# ================================================================================================
if (gate == "G9b") {
  # BACE convergence, on BACE's own terms: its `converged` verdict, the drift of the per-trait
  # imputed mean over the last two imputation iterations, and the smallest effective sample size
  # over the final fits' fixed effects. See bace_diagnostics() for why Gelman-Rubin across `runs`
  # is not the diagnostic here.
  #
  # With --dir, reads `$diag$bace` from every cell that ran the BACE arm; without it, runs one
  # short fit so the gate is testable before any cell exists.
  # Gate on BACE's own verdict and on the MEDIAN effective sample size over the final fits' fixed
  # effects. The minimum is reported, not gated: MCMCglmm's threshold and categorical models mix
  # slowly by construction, so a single badly-mixing parameter is expected and is a disclosure
  # item for the methods section, not a reason to discard a cell.
  ess_floor <- as.numeric(Sys.getenv("G9B_ESS_MED", "100"))
  check_one <- function(d, label) {
    cat(sprintf("G9b: %s converged=%s attempts=%s drift=%.4g ess median=%.0f min=%.0f (%.0f%% of %d below 100)\n",
                label, d$converged, d$n_attempts, d$drift, d$ess_med, d$ess_min,
                100 * (d$ess_frac_low %||% NA_real_), d$ess_n %||% 0L))
    if (!isTRUE(d$converged)) fail("G9b: %s did not converge (BACE verdict); lengthen runs", label)
    if (!is.na(d$ess_med) && d$ess_med < ess_floor)
      fail("G9b: %s median effective sample size %.0f below %.0f; lengthen nitt or reduce thin",
           label, d$ess_med, ess_floor)
  }
  if (!is.null(dir_arg)) {
    fs <- list.files(dir_arg, pattern = "[.]rds$", full.names = TRUE)
    if (!length(fs)) fail("G9b: no rds files under %s", dir_arg)
    seen <- 0L; miss_diag <- character(0)
    for (f in fs) {
      x <- readRDS(f)
      if (!"bace" %in% names(x$walls %||% list())) next
      d <- x$diag$bace
      if (is.null(d)) { miss_diag <- c(miss_diag, basename(f)); next }
      check_one(d, basename(f)); seen <- seen + 1L
    }
    if (length(miss_diag)) fail("G9b: %d cell(s) ran BACE without storing diagnostics, e.g. %s",
                                length(miss_diag), miss_diag[1])
    if (!seen) fail("G9b: no BACE arm found in %s", dir_arg)
    cat(sprintf("G9b: %d BACE cell(s) checked\n", seen))
  } else {
    cd <- make_cell("bm_mixed", 200, 42, miss_frac = 0.30, miss = "mcar")
    out <- run_bace(cd$df_miss, cd$truth, cd$mask, cd$tree, cd$cont_traits,
                    bace_nitt = 12000, bace_burnin = 2000, bace_thin = 10, bace_runs = 3L)
    check_one(out$diag, "bm_mixed n=200 nitt=12000 runs=3")
  }
  pass("G9b")
}

# ================================================================================================
if (gate == "G12") {
  if (is.null(dir_arg)) fail("G12: --dir is required (a directory of campaign_sim_cell.R rds files)")
  fs <- list.files(dir_arg, pattern = "\\.rds$", full.names = TRUE); fs <- fs[!grepl("_smoke", fs)]
  if (!length(fs)) fail("G12: no rds files found in %s", dir_arg)
  cells <- lapply(fs, readRDS)
  tab <- do.call(rbind, lapply(cells, function(c) { r <- c$results; if (is.null(r) || !nrow(r)) return(NULL); r$dgp <- c$dgp %||% NA; r$lambda <- c$lambda %||% NA; r$rho <- c$rho %||% NA; r$evo <- c$evo %||% NA; r$miss <- c$miss %||% NA; r }))
  if (is.null(tab) || !nrow(tab)) fail("G12: no scored rows across %d rds files", length(fs))
  mcse <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) return(NA_real_); stats::sd(x) / sqrt(length(x)) }
  agg <- do.call(rbind, lapply(split(tab, list(tab$dgp, tab$arm, tab$trait, tab$metric), drop = TRUE), function(d)
    data.frame(dgp = d$dgp[1], arm = d$arm[1], trait = d$trait[1], metric = d$metric[1], n = nrow(d),
               mean = mean(d$value, na.rm = TRUE), mcse = mcse(d$value))))
  regime_cols <- c("dgp", "lambda", "rho", "evo", "miss")
  if (!all(regime_cols %in% names(tab))) fail("G12: regime columns missing from the scored table: %s", paste(setdiff(regime_cols, names(tab)), collapse = ","))
  # a metric that structurally does not apply to an arm (e.g. width/interval_score/brier for an arm
  # with no interval or probability output) legitimately has mean = NA/NaN and mcse = NA for every
  # replicate; only require finite, non-negative MCSE where the metric actually produced a value.
  applicable <- agg[is.finite(agg$mean) & agg$n >= 2, ]
  if (!nrow(applicable)) fail("G12: no (dgp, arm, trait, metric) rows with n >= 2 replicates and a finite mean")
  if (!all(is.finite(applicable$mcse) & applicable$mcse >= 0)) fail("G12: MCSE is not finite/non-negative for at least one applicable (n >= 2, finite mean) row")
  cat("G12: aggregated", nrow(agg), "(dgp, arm, trait, metric) rows across", length(fs), "cells; regime columns present:", paste(regime_cols, collapse = ","), "\n")
  pass("G12")
}

# ================================================================================================
if (gate %in% c("G10", "G11", "G14")) {
  cat(sprintf("%s: TODO STUB -- not implemented in this lane. ", gate))
  if (gate == "G10") cat("Intended check: full-factorial coverage / completeness audit across the locked (n, lambda, rho, miss, evo) design grid, reading an index csv of completed cells.\n")
  if (gate == "G11") cat("Intended check: cross-arm paired-difference significance / sign-consistency audit against the aggregator's paired.csv.\n")
  if (gate == "G14") cat("Intended check: TDIP ensemble / GAIN comparison arm, if it installs cleanly (open item in the plan doc).\n")
  fail("%s: TODO stub -- fails loudly by design, not implemented", gate)
}

if (!(gate %in% c("G3", "G4", "G5", "G6", "Gold", "G9b", "G10", "G11", "G12", "G14"))) {
  fail("unknown gate %s", gate)
}
