#!/usr/bin/env Rscript
#
# dev/bace_timing_probe.R
#
# D-139 pre-run timing probe: ONE BACE::bace() call on a small subset of
# the bundled AVONET 300 data, with reduced MCMC settings, to estimate
# the wall-time of the full head-to-head benches before committing to a
# full run. Single-threaded, single R process.
#
# Not part of the package build (dev/ is .Rbuildignore'd). Scratch only.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(".", quiet = TRUE)
})

SEED      <- 2026L
MISS_FRAC <- 0.30
N_SUB     <- 100L
TRAITS    <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length")

cat(sprintf("[probe] pigauto loaded from source (worktree). BACE installed: %s\n",
            requireNamespace("BACE", quietly = TRUE)))
if (!requireNamespace("BACE", quietly = TRUE)) {
  stop("BACE not installed -- cannot run timing probe.", call. = FALSE)
}

# -------------------------------------------------------------------------
# 1. Load AVONET 300, subset to N_SUB species x 2-3 continuous traits
# -------------------------------------------------------------------------

e <- new.env(parent = emptyenv())
utils::data("avonet300", package = "pigauto", envir = e)
utils::data("tree300",   package = "pigauto", envir = e)
df   <- e$avonet300
tree <- e$tree300
rownames(df) <- df$Species_Key
df$Species_Key <- NULL

df <- df[, TRAITS, drop = FALSE]
stopifnot(all(rownames(df) == tree$tip.label))

set.seed(SEED)
keep <- sample(rownames(df), N_SUB)
df   <- df[keep, , drop = FALSE]
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
df   <- df[tree$tip.label, , drop = FALSE]
cat(sprintf("[probe] Subset: %d species x %d continuous traits (%s)\n",
            nrow(df), ncol(df), paste(TRAITS, collapse = ", ")))

# -------------------------------------------------------------------------
# 2. MCAR 30% mask
# -------------------------------------------------------------------------

set.seed(SEED)
mask_test <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df),
                     dimnames = list(rownames(df), names(df)))
for (v in names(df)) {
  obs_idx <- which(!is.na(df[[v]]))
  to_hide <- sample(obs_idx, ceiling(length(obs_idx) * MISS_FRAC))
  mask_test[to_hide, v] <- TRUE
}
df_miss <- df
for (v in names(df)) df_miss[[v]][mask_test[, v]] <- NA
cat("[probe] Held-out test cells per trait:\n")
print(colSums(mask_test))

# -------------------------------------------------------------------------
# 3. Patch zero-length edges (BACE/MCMCglmm requirement)
# -------------------------------------------------------------------------

tree_b <- tree
if (any(tree_b$edge.length == 0, na.rm = TRUE)) {
  n_zero <- sum(tree_b$edge.length == 0, na.rm = TRUE)
  cat(sprintf("[probe] Tree has %d zero-length edges -- patching to 1e-8\n", n_zero))
  tree_b$edge.length[tree_b$edge.length == 0] <- 1e-8
}

# -------------------------------------------------------------------------
# 4. Build fixformula (chained-equations, one formula per trait)
# -------------------------------------------------------------------------

df_b <- df_miss
df_b$Species <- rownames(df_miss)
all_traits <- setdiff(names(df_b), "Species")
fixformula <- lapply(all_traits, function(v) {
  others <- setdiff(all_traits, v)
  paste0(v, " ~ ", paste(others, collapse = " + "))
})
cat("[probe] fixformula:\n")
print(fixformula)

# -------------------------------------------------------------------------
# 5. Time ONE bace() call, reduced MCMC
# -------------------------------------------------------------------------

cat("[probe] Calling BACE::bace() with nitt=3000, burnin=500, thin=5, runs=1, n_final=5 ...\n")
t0 <- proc.time()[["elapsed"]]
bace_error <- NULL
out <- tryCatch({
  BACE::bace(
    fixformula     = fixformula,
    ran_phylo_form = "~ 1 |Species",
    phylo          = tree_b,
    data           = df_b,
    nitt           = 3000L,
    burnin         = 500L,
    thin           = 5L,
    runs           = 1L,
    n_final        = 5L,
    verbose        = FALSE,
    skip_conv      = TRUE
  )
}, error = function(e) {
  bace_error <<- conditionMessage(e)
  message("[probe] BACE run failed: ", bace_error)
  NULL
})
wall <- proc.time()[["elapsed"]] - t0

cat(sprintf("[probe] Wall time: %.1f s\n", wall))

if (is.null(out)) {
  cat(sprintf("[probe] RESULT: bace() errored: %s\n", bace_error))
} else {
  cat("[probe] RESULT: bace() returned successfully.\n")
  cat("[probe] Top-level names: ", paste(names(out), collapse = ", "), "\n")
  imputed_sets <- if ("imputed_datasets" %in% names(out)) out$imputed_datasets
    else if ("imputed_data" %in% names(out)) out$imputed_data
    else if ("data" %in% names(out)) list(out$data)
    else NULL
  if (is.null(imputed_sets) || !length(imputed_sets)) {
    cat("[probe] Could not locate imputed datasets in output shape.\n")
  } else {
    cat(sprintf("[probe] imputed_sets: M = %d draws, each %d x %d\n",
                length(imputed_sets), nrow(imputed_sets[[1]]), ncol(imputed_sets[[1]])))
    cat(sprintf("[probe] Got %d rows (expected %d species).\n",
                nrow(imputed_sets[[1]]), N_SUB))
  }
}

# -------------------------------------------------------------------------
# 6. Extrapolate full bench_avonet_bace.R wall-time
# -------------------------------------------------------------------------
#
# bench_avonet_bace.R (bundled AVONET 300 mode) uses:
#   n_species = 300 (vs probe's 100)
#   n_traits  = 7 (4 continuous + 2 categorical + 1 ordinal; vs probe's 3)
#   nitt = 2000, burnin = 500, thin = 5, runs = 2, n_final = 2
#     (vs probe's nitt = 3000, runs = 1, n_final = 5)
#
# MCMCglmm cost scales roughly linearly in nitt and in runs (independent
# chains), and each fixformula in the chained-equations list contributes
# roughly one MCMCglmm call per trait, so cost also scales ~linearly in
# n_traits. Per-iteration cost scales with n_species (through the phylo
# random-effect covariance structure), sub-linearly for small n but
# treated here as linear for a conservative (upper-bound) estimate.

n_iter_probe <- 3000L; runs_probe <- 1L
n_iter_full  <- 2000L; runs_full  <- 2L

scale_traits  <- 7 / length(TRAITS)
scale_species <- 300 / N_SUB
scale_iter    <- (n_iter_full * runs_full) / (n_iter_probe * runs_probe)

extrap <- wall * scale_traits * scale_species * scale_iter

cat("\n[probe] === Extrapolation to full bench_avonet_bace.R (bundled AVONET 300) ===\n")
cat(sprintf("[probe] probe wall time      : %.1f s (%d species, %d traits, nitt*runs = %d)\n",
            wall, N_SUB, length(TRAITS), n_iter_probe * runs_probe))
cat(sprintf("[probe] full-bench settings  : %d species, 7 traits, nitt*runs = %d\n",
            300L, n_iter_full * runs_full))
cat(sprintf("[probe] scale factors        : traits x%.2f, species x%.2f, iter*runs x%.2f\n",
            scale_traits, scale_species, scale_iter))
cat(sprintf("[probe] EXTRAPOLATED full BACE wall time: %.1f s (%.1f min)\n",
            extrap, extrap / 60))
cat("[probe] Note: this is the bace_default stage only; pigauto_default +\n")
cat("[probe]       mean_baseline stages add on top (bench_avonet_bace.R full\n")
cat("[probe]       script header estimates ~5-10 min total without BACE).\n")

cat("\n[probe] === DONE ===\n")
