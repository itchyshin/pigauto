#!/usr/bin/env Rscript
#
# script/bench_amphibio_discrete.R
#
# Real-data discrete verification of the strict val-floor v3 fix.
#
# AVONET's discrete columns (Trophic.Level, Primary.Lifestyle, Migration)
# are heavily phylo-conserved within bird families, so the calibration
# pushes their gates to corner values and the strict-discrete check
# never fires.  AmphiBIO's habitat indicators (Fos / Ter / Aqu / Arb)
# are more variable, so they exercise the discrete strict-floor path on
# real data.
#
# Trait subset (mixed continuous + binary):
#   Body_size_mm    numeric continuous (auto-log)
#   Body_mass_g     numeric continuous (auto-log)
#   Aqu             binary (aquatic 0/1)
#   Arb             binary (arboreal 0/1)
#   Fos             binary (fossorial 0/1)
#   Ter             binary (terrestrial 0/1)
#
# Default N_SUBSET = 1500.
#
# KNOWN BLOCKER (2026-04-29): the original bench_amphibio.R header
# noted that AmphiBIO binary columns trigger a Rphylopars internal
# "Not compatible with requested type: [type=character; target=double]"
# error inside fit_joint_threshold_baseline -> Rphylopars::phylopars().
# This bug reproduces at n=200 with ONE binary column (Aqu) on the
# 2026-04-29 codebase, so it is NOT a scale-dependent issue as the
# original comment suggested. The Rphylopars source path is
# `Rphylopars::phylopars` -> `estim_pars` -> `tp` (C++ core).
# Triage notes:
#   * AmphiBIO's tax-formula tree has very short branch lengths from
#     compute.brlen(method = "Grafen") — likely the C++ optimizer
#     hits a degenerate state when the threshold-joint path tries to
#     fit a binary liability column on this tree.
#   * Continuous columns alone (the production AmphiBIO bench) do
#     not trigger the bug, so the issue is specific to the threshold-
#     joint baseline's call signature for binary on this dataset.
# This script is preserved for the day the bug is resolved; today it
# errors out inside Rphylopars and cannot run end-to-end.
#
# Compares to mean / mode baseline to confirm the strict val-floor
# keeps pigauto >= baseline at the cell level for the binary columns.

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

SEED      <- 2026L
MISS_FRAC <- 0.30
N_IMP <- {
  v <- Sys.getenv("PIGAUTO_N_IMPUTATIONS", unset = "")
  if (nzchar(v)) as.integer(v) else 5L
}
N_SUB <- {
  v <- Sys.getenv("PIGAUTO_N_SUBSET", unset = "")
  if (nzchar(v)) as.integer(v) else 1500L
}

here <- "/Users/z3437171/Dropbox/Github Local/pigauto"
cache_dir <- file.path(here, "script", "data-cache")
csv_path  <- file.path(cache_dir, "AmphiBIO_v1.csv")

out_rds <- file.path(here, "script", "bench_amphibio_discrete.rds")
out_md  <- file.path(here, "script", "bench_amphibio_discrete.md")

script_start <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start),
      ..., "\n", sep = "")
  flush.console()
}

stopifnot(file.exists(csv_path))
amphi <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
log_line(sprintf("Loaded AmphiBIO: %d rows x %d cols", nrow(amphi), ncol(amphi)))

# Trait selection
cont_cols <- c("Body_size_mm", "Body_mass_g")
bin_cols  <- c("Aqu", "Arb", "Fos", "Ter")
pick_traits <- intersect(c(cont_cols, bin_cols), names(amphi))
log_line("Selected traits: ", paste(pick_traits, collapse = ", "))

tax_cols <- intersect(c("Order", "Family", "Genus", "Species"),
                       names(amphi))
stopifnot(all(c("Order", "Family", "Genus", "Species") %in% tax_cols))

df <- amphi[, c(tax_cols, pick_traits), drop = FALSE]
df <- df[nzchar(df$Species) & !is.na(df$Species), , drop = FALSE]
df <- df[!duplicated(df$Species), , drop = FALSE]

for (v in cont_cols) if (v %in% names(df))
  df[[v]] <- suppressWarnings(as.numeric(df[[v]]))

for (v in bin_cols) {
  if (v %in% names(df)) {
    raw <- df[[v]]
    raw_chr <- as.character(raw)
    raw_chr[raw_chr %in% c("1", "TRUE", "T", "yes", "Yes")] <- "yes"
    raw_chr[raw_chr %in% c("0", "FALSE", "F", "no", "No")]  <- "no"
    raw_chr[!raw_chr %in% c("yes", "no")] <- NA
    df[[v]] <- factor(raw_chr, levels = c("no", "yes"))
  }
}

trait_cols <- setdiff(names(df), tax_cols)
has_any <- rowSums(!is.na(df[, trait_cols, drop = FALSE])) > 0L
df <- df[has_any, , drop = FALSE]
tax_ok <- rowSums(is.na(df[, tax_cols, drop = FALSE])) == 0 &
           nzchar(df$Order) & nzchar(df$Family) & nzchar(df$Genus)
df <- df[tax_ok, , drop = FALSE]

if (!is.na(N_SUB) && N_SUB < nrow(df)) {
  set.seed(SEED)
  df <- df[sample.int(nrow(df), N_SUB), , drop = FALSE]
  log_line(sprintf("N_SUBSET = %d -> %d species", N_SUB, nrow(df)))
}
log_line(sprintf("Final: %d species x %d traits (%d cont + %d bin)",
                  nrow(df), length(trait_cols), length(cont_cols),
                  length(bin_cols)))

# Build taxonomic tree (Order/Family/Genus/Species).  ape::as.phylo.formula
# requires the response variable Species to be a factor.
log_line("Building taxonomic tree (Grafen branch lengths) ...")
df$Order   <- factor(df$Order)
df$Family  <- factor(df$Family)
df$Genus   <- factor(df$Genus)
df$Species <- factor(df$Species)
tax_formula <- ~Order/Family/Genus/Species
tree <- ape::as.phylo.formula(tax_formula, data = df)
tree <- ape::compute.brlen(tree, method = "Grafen")
tree$tip.label <- df$Species[match(tree$tip.label, df$Species)]
df  <- df[match(tree$tip.label, df$Species), , drop = FALSE]
rownames(df) <- df$Species
df_traits <- df[, trait_cols, drop = FALSE]
log_line(sprintf("Tree: %d tips, %d nodes", length(tree$tip.label),
                  tree$Nnode))

# -------------------------------------------------------------------------
# Hold out 30% MCAR test cells per trait
# -------------------------------------------------------------------------
set.seed(SEED)
df_truth <- df_traits
mask_test <- matrix(FALSE, nrow = nrow(df_traits),
                     ncol = ncol(df_traits),
                     dimnames = list(rownames(df_traits),
                                     colnames(df_traits)))
for (v in colnames(df_traits)) {
  obs_idx <- which(!is.na(df_traits[[v]]))
  to_hide <- sample(obs_idx, ceiling(length(obs_idx) * MISS_FRAC))
  mask_test[to_hide, v] <- TRUE
}
df_miss <- df_traits
for (v in colnames(df_miss)) df_miss[[v]][mask_test[, v]] <- NA

log_line("Held-out test cells per trait:")
print(colSums(mask_test))

# -------------------------------------------------------------------------
# Method 1: mean / mode baseline
# -------------------------------------------------------------------------
log_line("[STAGE] mean_mode baseline starting")
t_mb <- proc.time()[["elapsed"]]
df_mb <- df_miss
for (v in cont_cols) {
  if (!v %in% names(df_mb)) next
  df_mb[[v]][is.na(df_mb[[v]])] <- mean(df_traits[[v]], na.rm = TRUE)
}
for (v in bin_cols) {
  if (!v %in% names(df_mb)) next
  modal <- names(which.max(table(df_traits[[v]])))
  df_mb[[v]][is.na(df_mb[[v]])] <- factor(modal, levels = levels(df_traits[[v]]))
}
mb_wall <- proc.time()[["elapsed"]] - t_mb
log_line(sprintf("[STAGE] mean_mode done in %.1fs", mb_wall))

# -------------------------------------------------------------------------
# Method 2: pigauto with v3 strict val-floor
# -------------------------------------------------------------------------
log_line("[STAGE] pigauto starting (n_imputations = ", N_IMP, ")")
t_pg <- proc.time()[["elapsed"]]
res <- pigauto::impute(df_miss, tree,
                        epochs = 500L, n_imputations = N_IMP,
                        verbose = FALSE, seed = SEED)
pg_wall <- proc.time()[["elapsed"]] - t_pg
log_line(sprintf("[STAGE] pigauto done in %.1fs", pg_wall))
df_pg <- res$completed

# -------------------------------------------------------------------------
# Score: per-trait RMSE (continuous) and accuracy + Brier (binary)
# -------------------------------------------------------------------------
results <- data.frame(method = character(0), trait = character(0),
                      metric = character(0), value = numeric(0),
                      n_cells = integer(0), wall_s = numeric(0))

for (v in cont_cols) {
  if (!v %in% colnames(mask_test)) next
  truth <- df_truth[[v]][mask_test[, v]]
  ok <- is.finite(truth)
  if (sum(ok) == 0L) next
  for (m in c("mean_baseline", "pigauto")) {
    pred <- if (m == "mean_baseline") df_mb[[v]] else df_pg[[v]]
    pv   <- pred[mask_test[, v]][ok]
    rmse <- sqrt(mean((pv - truth[ok])^2, na.rm = TRUE))
    rho  <- suppressWarnings(stats::cor(pv, truth[ok], use = "complete.obs"))
    wall <- if (m == "mean_baseline") mb_wall else pg_wall
    results <- rbind(results,
      data.frame(method = m, trait = v, metric = "rmse", value = rmse,
                  n_cells = sum(ok), wall_s = wall),
      data.frame(method = m, trait = v, metric = "pearson_r",
                  value = ifelse(is.finite(rho), rho, NA),
                  n_cells = sum(ok), wall_s = wall))
  }
}

for (v in bin_cols) {
  if (!v %in% colnames(mask_test)) next
  truth <- df_truth[[v]][mask_test[, v]]
  ok <- !is.na(truth)
  if (sum(ok) == 0L) next
  for (m in c("mean_baseline", "pigauto")) {
    pred <- if (m == "mean_baseline") df_mb[[v]] else df_pg[[v]]
    pv   <- pred[mask_test[, v]][ok]
    acc  <- mean(as.character(pv) == as.character(truth[ok]))
    p_yes <- mean(as.character(pv) == "yes")
    truth_yes <- as.numeric(as.character(truth[ok]) == "yes")
    pred_yes  <- as.numeric(as.character(pv) == "yes")
    brier <- mean((pred_yes - truth_yes)^2)
    wall <- if (m == "mean_baseline") mb_wall else pg_wall
    results <- rbind(results,
      data.frame(method = m, trait = v, metric = "accuracy",
                  value = acc, n_cells = sum(ok), wall_s = wall),
      data.frame(method = m, trait = v, metric = "brier",
                  value = brier, n_cells = sum(ok), wall_s = wall))
  }
}

log_line("=== Per-trait results ===")
print(results)

# -------------------------------------------------------------------------
# Save
# -------------------------------------------------------------------------
saveRDS(list(results = results,
              seed = SEED, miss_frac = MISS_FRAC,
              n_species = nrow(df), n_subset = N_SUB,
              n_imputations = N_IMP), out_rds)
log_line(sprintf("Wrote %s", out_rds))

# Summary md
md <- c(
  "# AmphiBIO discrete-trait verification (post v3 strict val-floor)",
  "",
  sprintf("n = %d species x %d traits (cont: %s; binary: %s).",
           nrow(df), length(trait_cols),
           paste(cont_cols, collapse = ", "),
           paste(bin_cols, collapse = ", ")),
  sprintf("Seed = %d, miss_frac = %.2f, n_imputations = %d.",
           SEED, MISS_FRAC, N_IMP),
  "",
  "Purpose: real-data verification of the strict-discrete val-floor",
  "(commit 1ac34b1) on binary trait columns where the calibrated gate",
  "is more likely to be non-corner than AVONET's bird categoricals.",
  "",
  "## Per-trait metrics",
  "",
  "```",
  capture.output(print(results)),
  "```"
)
writeLines(md, out_md)
log_line(sprintf("Wrote %s", out_md))
log_line("=== DONE ===")
