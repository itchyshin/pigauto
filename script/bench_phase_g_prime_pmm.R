#!/usr/bin/env Rscript
#
# script/bench_phase_g_prime_pmm.R
#
# Phase G' acceptance bench: PMM vs clamp vs none on AVONET 4
# log-transformed continuous traits + a heavy-tail simulation.
#
# Pre-registered question:
#   Does match_observed = "pmm" (a) close the AVONET Mass tail
#   blow-up at least as well as clamp_outliers; (b) preserve
#   accuracy on traits where Phase G's clamp had no effect; (c)
#   work cleanly on a synthetic heavy-tail dataset (extra-AVONET
#   evidence)?
#
# Pre-registered acceptance criteria for default-flip in v0.9.2:
#   1. PMM Mass RMSE on seed 2030 <= clamp RMSE on seed 2030.
#   2. PMM does not regress more than 5 % on any (seed, trait, dataset)
#      cell vs the no-clamp baseline.
#   3. On the synthetic heavy-tail dataset, PMM imputed values are
#      ALWAYS in the observed range (programmatic check).
#
# Bench design:
#   AVONET n=1500, 3 seeds (2030/2031/2032), N_IMP=20.
#   Configurations: {none, clamp@5, pmm@K=5} -- 3 fits per seed,
#   ~9 min each = ~80 min total.
#   Plus synthetic heavy-tail simulation (n=200, 3 reps, 1 trait,
#   ~3 min total).
#
# Output:
#   script/bench_phase_g_prime_pmm.{rds,md}

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH", unset = "")
  if (nzchar(pkg_path) && dir.exists(pkg_path) &&
      file.exists(file.path(pkg_path, "DESCRIPTION"))) {
    devtools::load_all(pkg_path, quiet = TRUE)
  } else {
    library(pigauto)
  }
})

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "bench_phase_g_prime_pmm.rds")
out_md  <- file.path(here, "script", "bench_phase_g_prime_pmm.md")

SEEDS     <- c(2030L, 2031L, 2032L)
MISS_FRAC <- 0.30
N_SUB     <- 1500L
N_IMP     <- 20L
CONT_TRAITS <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - t0),
      ..., "\n", sep = "")
  flush.console()
}

# ---- AVONET part -----------------------------------------------------------

run_avonet_one <- function(seed, config) {
  log_line(sprintf("AVONET seed %d config=%s ...", seed, config))
  e <- new.env(parent = emptyenv())
  utils::data("avonet_full", package = "pigauto", envir = e)
  utils::data("tree_full",   package = "pigauto", envir = e)
  df   <- e$avonet_full
  tree <- e$tree_full
  rownames(df) <- df$Species_Key
  df$Species_Key <- NULL

  set.seed(seed)
  keep <- sample(rownames(df), N_SUB)
  df   <- df[keep, , drop = FALSE]
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
  df   <- df[tree$tip.label, , drop = FALSE]

  set.seed(seed)
  df_truth  <- df
  mask_test <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df),
                      dimnames = list(rownames(df), names(df)))
  for (v in names(df)) {
    obs_idx <- which(!is.na(df[[v]]))
    to_hide <- sample(obs_idx, ceiling(length(obs_idx) * MISS_FRAC))
    mask_test[to_hide, v] <- TRUE
  }
  df_miss <- df
  for (v in names(df)) df_miss[[v]][mask_test[, v]] <- NA

  args <- switch(config,
    none  = list(clamp_outliers = FALSE, match_observed = "none"),
    clamp = list(clamp_outliers = TRUE,  match_observed = "none",
                 clamp_factor = 5),
    pmm   = list(clamp_outliers = FALSE, match_observed = "pmm",
                 pmm_K = 5L)
  )
  res <- do.call(pigauto::impute, c(
    list(df_miss, tree,
         epochs = 200L, n_imputations = N_IMP,
         verbose = FALSE, seed = seed),
    args
  ))

  out <- vector("list", length(CONT_TRAITS))
  for (i in seq_along(CONT_TRAITS)) {
    tr <- CONT_TRAITS[i]
    idx <- which(mask_test[, tr])
    truth <- as.numeric(df_truth[[tr]][idx])
    pred  <- as.numeric(res$completed[[tr]][idx])
    rmse  <- sqrt(mean((truth - pred)^2, na.rm = TRUE))
    obs_max <- max(as.numeric(df_truth[[tr]]), na.rm = TRUE)
    obs_min <- min(as.numeric(df_truth[[tr]]), na.rm = TRUE)
    in_range <- mean(pred >= obs_min & pred <= obs_max, na.rm = TRUE)
    out[[i]] <- data.frame(
      dataset = "AVONET",
      seed    = seed,
      config  = config,
      trait   = tr,
      rmse    = rmse,
      pct_in_observed_range = 100 * in_range,
      max_pred = max(pred, na.rm = TRUE),
      obs_max  = obs_max,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

# ---- Synthetic heavy-tail simulation ---------------------------------------

run_synthetic_one <- function(rep_id, config) {
  log_line(sprintf("Synthetic rep %d config=%s ...", rep_id, config))
  set.seed(rep_id * 100L + 1L)
  n <- 200L
  tree <- ape::rtree(n)
  # Heavy-tail: mostly mass ~ 1000g, but 5 % of species are 10x-50x
  # larger (an isolated clade analogue).  Phylogenetically structured
  # via tip index modulo.
  log_mass <- stats::rnorm(n, 7, 0.5)        # ~1000g
  outlier_idx <- sample.int(n, ceiling(0.05 * n))
  log_mass[outlier_idx] <- log_mass[outlier_idx] + stats::rnorm(length(outlier_idx),
                                                                 mean = 2,
                                                                 sd = 0.5)
  df <- data.frame(mass = exp(log_mass), row.names = tree$tip.label)
  df_truth <- df
  obs_max <- max(df$mass)
  obs_min <- min(df$mass)
  set.seed(rep_id * 100L + 2L)
  mask_idx <- sample.int(n, ceiling(MISS_FRAC * n))
  df$mass[mask_idx] <- NA

  args <- switch(config,
    none  = list(clamp_outliers = FALSE, match_observed = "none"),
    clamp = list(clamp_outliers = TRUE,  match_observed = "none",
                 clamp_factor = 5),
    pmm   = list(clamp_outliers = FALSE, match_observed = "pmm",
                 pmm_K = 5L)
  )
  res <- do.call(pigauto::impute, c(
    list(df, tree,
         epochs = 100L, n_imputations = 10L,
         verbose = FALSE, seed = rep_id * 100L + 3L),
    args
  ))
  truth <- df_truth$mass[mask_idx]
  pred  <- as.numeric(res$completed$mass[mask_idx])
  rmse  <- sqrt(mean((truth - pred)^2, na.rm = TRUE))
  in_range <- mean(pred >= obs_min & pred <= obs_max, na.rm = TRUE)
  data.frame(
    dataset = "synthetic_heavy_tail",
    seed    = rep_id,
    config  = config,
    trait   = "mass",
    rmse    = rmse,
    pct_in_observed_range = 100 * in_range,
    max_pred = max(pred, na.rm = TRUE),
    obs_max  = obs_max,
    stringsAsFactors = FALSE
  )
}

# ---- Sweep -----------------------------------------------------------------

CONFIGS <- c("none", "clamp", "pmm")

results <- list()

# AVONET
for (s in SEEDS) {
  for (cfg in CONFIGS) {
    results[[length(results) + 1L]] <- run_avonet_one(s, cfg)
  }
}

# Synthetic heavy-tail
for (rp in 1:3) {
  for (cfg in CONFIGS) {
    results[[length(results) + 1L]] <- run_synthetic_one(rp, cfg)
  }
}

all_rows <- do.call(rbind, results)
saveRDS(all_rows, out_rds)

# ---- Acceptance check ------------------------------------------------------

# Criterion 1: PMM Mass RMSE on seed 2030 <= clamp RMSE on seed 2030
mass_2030 <- subset(all_rows, dataset == "AVONET" & trait == "Mass" &
                              seed == 2030L)
crit_1 <- mass_2030$rmse[mass_2030$config == "pmm"] <=
          mass_2030$rmse[mass_2030$config == "clamp"]

# Criterion 2: PMM doesn't regress > 5 % vs none on any (dataset, seed, trait)
none_rows <- subset(all_rows, config == "none")
pmm_rows  <- subset(all_rows, config == "pmm")
key <- paste(none_rows$dataset, none_rows$seed, none_rows$trait, sep = "/")
pkey <- paste(pmm_rows$dataset, pmm_rows$seed, pmm_rows$trait, sep = "/")
ord <- match(key, pkey)
pct_change <- 100 * (pmm_rows$rmse[ord] - none_rows$rmse) / none_rows$rmse
crit_2_violators <- which(pct_change > 5)
crit_2 <- length(crit_2_violators) == 0L

# Criterion 3: PMM imputed values 100 % in observed range
pmm_in_range <- subset(all_rows, config == "pmm")$pct_in_observed_range
crit_3 <- all(pmm_in_range > 99.5)

verdict <- if (crit_1 && crit_2 && crit_3) {
  "**PHASE G' PASSES all 3 criteria** -- recommend default flip to match_observed = 'pmm' for at-risk types in v0.9.2"
} else {
  parts <- c(
    if (!crit_1) sprintf("FAIL 1: PMM seed-2030 Mass RMSE %.0f > clamp %.0f",
                         mass_2030$rmse[mass_2030$config == "pmm"],
                         mass_2030$rmse[mass_2030$config == "clamp"]),
    if (!crit_2) sprintf("FAIL 2: PMM regresses > 5 %% on %d cells",
                         length(crit_2_violators)),
    if (!crit_3) sprintf("FAIL 3: PMM in-range only %.1f %%-%.1f %%",
                         min(pmm_in_range), max(pmm_in_range))
  )
  paste("**PHASE G' PARTIAL** --", paste(parts, collapse = "; "))
}

# ---- Markdown report -------------------------------------------------------

md <- c(
  "# Phase G' acceptance bench: PMM vs clamp vs default",
  "",
  sprintf("Run: %s", format(Sys.time())),
  sprintf("AVONET n=%d, miss_frac=%.2f, seeds=%s, N_IMP=%d.  Synthetic heavy-tail: 3 reps.",
          N_SUB, MISS_FRAC,
          paste(SEEDS, collapse = ","), N_IMP),
  "",
  "## Pre-registered acceptance criteria",
  "",
  "1. PMM Mass RMSE on seed 2030 <= clamp RMSE on seed 2030 (PMM at least as good as clamp).",
  "2. PMM doesn't regress > 5 % vs no-clamp on any (dataset, seed, trait) cell (no harm).",
  "3. PMM imputed values 100 % in observed range (by design; programmatic check).",
  "",
  "## Per-cell results",
  "",
  "```",
  capture.output(print(all_rows[order(all_rows$dataset, all_rows$trait,
                                       all_rows$seed, all_rows$config),
                                 c("dataset", "trait", "seed", "config",
                                   "rmse", "pct_in_observed_range",
                                   "max_pred", "obs_max")],
                       row.names = FALSE)),
  "```",
  "",
  sprintf("Crit 1 (PMM seed-2030 Mass <= clamp): %s",
          if (crit_1) "PASS" else "FAIL"),
  sprintf("Crit 2 (PMM no-regress vs none, threshold 5%%): %s",
          if (crit_2) "PASS"
          else sprintf("FAIL on %d cells", length(crit_2_violators))),
  sprintf("Crit 3 (PMM 100%% in observed range): %s",
          if (crit_3)
            sprintf("PASS (min %.1f %%, max %.1f %%)",
                    min(pmm_in_range), max(pmm_in_range))
          else
            sprintf("FAIL (min %.1f %%, max %.1f %%)",
                    min(pmm_in_range), max(pmm_in_range))),
  "",
  "## Verdict",
  "",
  verdict,
  ""
)
writeLines(md, out_md)
log_line("DONE -- ", out_rds)
log_line("DONE -- ", out_md)
