#!/usr/bin/env Rscript
# script/regress.R
# Full-canary regression check for safety_floor = TRUE (v0.9.1.9002+).
#
# Runs pigauto on 5 taxa at n = 1000 each with safety_floor = TRUE vs
# FALSE, enforces the hard regression gates from the safety-floor design
# spec (specs/2026-04-23-safety-floor-mean-gate-design.md §5.2):
#
#   - Continuous RMSE:     RMSE_on  <= 1.02 * RMSE_off
#   - Discrete accuracy:   acc_on   >= acc_off - 0.01
#   - Plants safety:       RMSE_on  <= 1.02 * mean_baseline RMSE
#   - Conformal coverage:  coverage_on in [0.90, 1.00]
#
# Writes script/regress_result.json with per-bench per-trait pass/fail.
# Exits non-zero if any bench fails.
#
# Invocation:
#   Rscript script/regress.R
#
# Wall: ~20 min on Apple Silicon MPS with all caches present.
# ~5 min (AVONET birds only) without any caches.

options(warn = 1, stringsAsFactors = FALSE)

# Silence R-level warnings. Armadillo's "solve(): system is singular"
# stream from Rphylopars EM iterations comes from C++ stderr and can
# overflow shell pipe buffers — redirect stderr to /dev/null at the
# shell level to suppress:
#
#   Rscript script/regress.R > regress.log 2>/dev/null
#
# The approx-solution path is still mathematically valid; the warnings
# are purely informational.
options(warn = -1L)

suppressPackageStartupMessages({
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH", unset = "")
  if (nzchar(pkg_path) && dir.exists(pkg_path) &&
      file.exists(file.path(pkg_path, "DESCRIPTION"))) {
    devtools::load_all(pkg_path, quiet = TRUE)
  } else {
    # Fallback: if pigauto isn't installed system-wide, try to load from cwd.
    tryCatch(library(pigauto),
              error = function(e) devtools::load_all(".", quiet = TRUE))
  }
  library(jsonlite)
  library(ape)
})

set.seed(2026L)

N_PER_TAXON  <- 1000L
SEED         <- 2026L
MISS_FRAC    <- 0.30
N_IMPS       <- 10L
EPOCHS       <- 300L
# Tolerances are configurable via PIGAUTO_REGRESS_MODE.
#
# "smoke" (default): calibrated for the n=1000 subset used below. At this
# size ~30-100 test cells per continuous trait, the val-to-test RMSE
# sampling noise is ~15-30% on heavy-tailed traits (plant seed_mass, bird
# Tarsus). Enforcing the full design-spec 2% tolerance here would flag
# noise as regression.
#
# "strict": design-spec production thresholds, intended for the manual
# pre-release canary at n >= 5000 species where the sampling noise drops
# to ~2-5%. See spec §5.2.
mode <- Sys.getenv("PIGAUTO_REGRESS_MODE", unset = "smoke")
if (identical(mode, "strict")) {
  TOL_RMSE   <- 1.02
  TOL_ACC    <- 0.01
  TOL_COV_LO <- 0.90
  TOL_COV_HI <- 1.00
} else {
  TOL_RMSE   <- 1.30
  TOL_ACC    <- 0.03
  TOL_COV_LO <- 0.85
  TOL_COV_HI <- 1.00
}
cat_line_mode <- function() cat(sprintf("[mode=%s] TOL_RMSE=%.2f TOL_ACC=%.2f TOL_COV=[%.2f,%.2f]\n",
                                          mode, TOL_RMSE, TOL_ACC, TOL_COV_LO, TOL_COV_HI))
cat_line_mode()

cat_line <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

taxa_results <- list()

# ---------------------------------------------------------------------------
# bench_one(): runs pigauto twice (safety_floor off vs on) and scores
# ---------------------------------------------------------------------------
bench_one <- function(taxon, df_truth, tree, cont_cols, disc_cols) {
  cat_line(sprintf("== %s: n = %d, p_cont = %d, p_disc = %d ==",
                    taxon, nrow(df_truth), length(cont_cols), length(disc_cols)))
  if (is.null(rownames(df_truth)) ||
      all(rownames(df_truth) == as.character(seq_len(nrow(df_truth))))) {
    stop("rownames(df_truth) must be species labels matching tree tip labels")
  }

  df <- df_truth
  all_cols <- c(cont_cols, disc_cols)
  mask <- matrix(FALSE, nrow = nrow(df), ncol = length(all_cols),
                 dimnames = list(rownames(df), all_cols))
  for (v in all_cols) {
    ok  <- which(!is.na(df_truth[[v]]))
    n_m <- round(MISS_FRAC * length(ok))
    if (n_m < 3L) next
    idx <- sample(ok, size = n_m)
    mask[idx, v] <- TRUE
    if (v %in% cont_cols) df[[v]][idx] <- NA_real_ else df[[v]][idx] <- NA
  }

  # Fit safety_floor = FALSE
  t0 <- proc.time()[["elapsed"]]
  res_off <- tryCatch(
    pigauto::impute(df, tree, safety_floor = FALSE,
                    epochs = EPOCHS, n_imputations = N_IMPS,
                    verbose = FALSE, seed = SEED),
    error = function(e) { cat_line("ERROR (off): ", conditionMessage(e)); NULL }
  )
  wall_off <- proc.time()[["elapsed"]] - t0

  # Fit safety_floor = TRUE
  t0 <- proc.time()[["elapsed"]]
  res_on <- tryCatch(
    pigauto::impute(df, tree, safety_floor = TRUE,
                    epochs = EPOCHS, n_imputations = N_IMPS,
                    verbose = FALSE, seed = SEED),
    error = function(e) { cat_line("ERROR (on): ", conditionMessage(e)); NULL }
  )
  wall_on <- proc.time()[["elapsed"]] - t0

  if (is.null(res_off) || is.null(res_on)) {
    return(list(taxon = taxon, n = nrow(df_truth),
                error = "impute() failed for at least one condition",
                pass = FALSE))
  }

  out <- list(taxon = taxon, n = nrow(df_truth),
              wall_off_s = wall_off, wall_on_s = wall_on,
              traits = list(), pass = TRUE)

  for (v in cont_cols) {
    if (!any(mask[, v])) next
    truth_v <- df_truth[[v]][mask[, v]]
    ok <- is.finite(truth_v)
    if (sum(ok) < 3L) next

    pred_off_v <- res_off$completed[[v]][mask[, v]][ok]
    pred_on_v  <- res_on$completed[[v]][mask[, v]][ok]
    tv_ok      <- truth_v[ok]

    rmse_off  <- sqrt(mean((pred_off_v - tv_ok)^2))
    rmse_on   <- sqrt(mean((pred_on_v  - tv_ok)^2))
    mean_pred <- mean(df[[v]], na.rm = TRUE)
    rmse_mean <- sqrt(mean((mean_pred - tv_ok)^2))

    # Conformal coverage for safety-on
    cov_on <- NA_real_
    cli <- res_on$prediction$conformal_lower
    cui <- res_on$prediction$conformal_upper
    if (!is.null(cli) && !is.null(cui) && v %in% colnames(cli)) {
      lo   <- cli[mask[, v], v][ok]
      hi   <- cui[mask[, v], v][ok]
      covv <- is.finite(lo) & is.finite(hi)
      if (any(covv)) {
        cov_on <- mean(tv_ok[covv] >= lo[covv] & tv_ok[covv] <= hi[covv])
      }
    }

    pass_rmse <- isTRUE(rmse_on <= rmse_off * TOL_RMSE)
    pass_mean <- isTRUE(rmse_on <= rmse_mean * TOL_RMSE)
    pass_cov  <- is.na(cov_on) || (cov_on >= TOL_COV_LO & cov_on <= TOL_COV_HI)
    pass_all  <- pass_rmse && pass_mean && pass_cov

    out$traits[[v]] <- list(
      type = "continuous",
      rmse_off = rmse_off, rmse_on = rmse_on, rmse_mean = rmse_mean,
      coverage_on = cov_on,
      pass_rmse = pass_rmse, pass_mean = pass_mean, pass_cov = pass_cov,
      pass = pass_all)
    if (!pass_all) out$pass <- FALSE
  }

  for (v in disc_cols) {
    if (!any(mask[, v])) next
    truth_v <- as.character(df_truth[[v]][mask[, v]])
    ok       <- !is.na(truth_v)
    if (sum(ok) < 3L) next
    acc_off <- mean(as.character(res_off$completed[[v]][mask[, v]])[ok] == truth_v[ok],
                    na.rm = TRUE)
    acc_on  <- mean(as.character(res_on$completed[[v]][mask[, v]])[ok]  == truth_v[ok],
                    na.rm = TRUE)
    pass_acc <- isTRUE(acc_on >= acc_off - TOL_ACC)
    out$traits[[v]] <- list(
      type = "discrete",
      acc_off = acc_off, acc_on = acc_on,
      pass = pass_acc)
    if (!pass_acc) out$pass <- FALSE
  }

  out
}

# ---------------------------------------------------------------------------
# Birds: AVONET full (bundled; always exercised)
# ---------------------------------------------------------------------------

cat_line("Loading avonet_full + tree_full from pigauto ...")
e <- new.env(parent = emptyenv())
utils::data("avonet_full", package = "pigauto", envir = e)
utils::data("tree_full",   package = "pigauto", envir = e)
av_all   <- e$avonet_full
av_tree0 <- e$tree_full

# Species_Key is the row identifier and matches tree$tip.label
rownames(av_all) <- av_all$Species_Key
av_all$Species_Key <- NULL

set.seed(SEED)
av_sp <- sample(rownames(av_all), min(N_PER_TAXON, nrow(av_all)))
av_df   <- av_all[av_sp, , drop = FALSE]
av_tree <- ape::keep.tip(av_tree0, av_sp)
av_df   <- av_df[av_tree$tip.label, , drop = FALSE]

taxa_results$birds <- bench_one(
  "birds", av_df, av_tree,
  cont_cols = c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length"),
  disc_cols = c("Trophic.Level", "Primary.Lifestyle", "Migration"))

# ---------------------------------------------------------------------------
# Mammals: PanTHERIA (cache: script/data-cache/pantheria_matched.rds)
# Expected schema: list(df = data.frame [rownames = species],
#                       tree = phylo,
#                       cont_cols = character,
#                       disc_cols = character)
# ---------------------------------------------------------------------------
cache_ph <- "script/data-cache/pantheria_matched.rds"
if (file.exists(cache_ph)) {
  ph  <- readRDS(cache_ph)
  set.seed(SEED)
  idx <- sample(seq_len(nrow(ph$df)), min(N_PER_TAXON, nrow(ph$df)))
  ph_sp   <- rownames(ph$df)[idx]
  ph_df   <- ph$df[idx, , drop = FALSE]
  ph_tree <- ape::keep.tip(ph$tree, ph_sp)
  ph_df   <- ph_df[ph_tree$tip.label, , drop = FALSE]
  taxa_results$mammals <- bench_one(
    "mammals", ph_df, ph_tree,
    cont_cols = ph$cont_cols, disc_cols = ph$disc_cols)
} else {
  cat_line("SKIP mammals (no cache at ", cache_ph, ")")
  taxa_results$mammals <- list(taxon = "mammals", skipped = TRUE,
                                reason = "no cache", pass = TRUE)
}

# ---------------------------------------------------------------------------
# Fish: FishBase + fishtree (cache: script/data-cache/fishbase_matched.rds)
# Expected schema: same as pantheria — list(df, tree, cont_cols, disc_cols)
# ---------------------------------------------------------------------------
cache_fs <- "script/data-cache/fishbase_matched.rds"
if (file.exists(cache_fs)) {
  fs  <- readRDS(cache_fs)
  set.seed(SEED)
  idx <- sample(seq_len(nrow(fs$df)), min(N_PER_TAXON, nrow(fs$df)))
  fs_sp   <- rownames(fs$df)[idx]
  fs_df   <- fs$df[idx, , drop = FALSE]
  fs_tree <- ape::keep.tip(fs$tree, fs_sp)
  fs_df   <- fs_df[fs_tree$tip.label, , drop = FALSE]
  taxa_results$fish <- bench_one(
    "fish", fs_df, fs_tree,
    cont_cols = fs$cont_cols, disc_cols = fs$disc_cols)
} else {
  cat_line("SKIP fish (no cache at ", cache_fs, ")")
  taxa_results$fish <- list(taxon = "fish", skipped = TRUE,
                              reason = "no cache", pass = TRUE)
}

# ---------------------------------------------------------------------------
# Amphibians: AmphiBIO (cache: script/data-cache/amphibio_matched.rds)
# Expected schema: same as pantheria — list(df, tree, cont_cols, disc_cols)
# Note: bench_amphibio.R builds from raw CSV; this cache is for regress.R.
#       Run script/bench_amphibio.R with PIGAUTO_SAVE_CACHE=1 to build it,
#       or create it manually from the bench output.
# ---------------------------------------------------------------------------
cache_am <- "script/data-cache/amphibio_matched.rds"
if (file.exists(cache_am)) {
  am  <- readRDS(cache_am)
  set.seed(SEED)
  idx <- sample(seq_len(nrow(am$df)), min(N_PER_TAXON, nrow(am$df)))
  am_sp   <- rownames(am$df)[idx]
  am_df   <- am$df[idx, , drop = FALSE]
  am_tree <- ape::keep.tip(am$tree, am_sp)
  am_df   <- am_df[am_tree$tip.label, , drop = FALSE]
  am_disc <- if (is.null(am$disc_cols)) character(0) else am$disc_cols
  taxa_results$amphibians <- bench_one(
    "amphibians", am_df, am_tree,
    cont_cols = am$cont_cols, disc_cols = am_disc)
} else {
  cat_line("SKIP amphibians (no cache at ", cache_am, ")")
  taxa_results$amphibians <- list(taxon = "amphibians", skipped = TRUE,
                                    reason = "no cache", pass = TRUE)
}

# ---------------------------------------------------------------------------
# Plants: BIEN + V.PhyloMaker2 (cache: bien_trait_means.rds + bien_tree.rds)
#
# Cache schema (as saved by bench_bien.R):
#   bien_trait_means.rds : named list of per-trait data.frames,
#                          each with columns $species and $mean_value.
#                          Names: height_m, leaf_area, sla, seed_mass,
#                                 wood_density.
#   bien_tree.rds        : list(scenario.3 = phylo, species.list = df)
#                          tip.label uses underscores ("Genus_species").
# ---------------------------------------------------------------------------
cache_bt <- "script/data-cache/bien_trait_means.rds"
cache_br <- "script/data-cache/bien_tree.rds"

if (file.exists(cache_bt) && file.exists(cache_br)) {
  cat_line("Loading BIEN trait means cache ...")
  trait_means <- readRDS(cache_bt)

  # Pivot from per-trait list -> species x trait wide data.frame
  species_lists <- lapply(trait_means, function(d) if (is.null(d)) character(0) else d$species)
  all_species <- Reduce(union, species_lists)
  bwide <- data.frame(species = all_species, stringsAsFactors = FALSE)
  for (nm in names(trait_means)) {
    d <- trait_means[[nm]]
    if (is.null(d)) { bwide[[nm]] <- NA_real_; next }
    val_col <- intersect(c("mean_value", "trait_value", "value"), names(d))[1]
    if (is.na(val_col)) { bwide[[nm]] <- NA_real_; next }
    m <- match(bwide$species, d$species)
    bwide[[nm]] <- suppressWarnings(as.numeric(d[[val_col]][m]))
  }
  trait_cols_bien <- names(trait_means)
  n_obs_per_sp <- rowSums(!is.na(bwide[, trait_cols_bien, drop = FALSE]))
  bwide <- bwide[n_obs_per_sp >= 1L, , drop = FALSE]
  rownames(bwide) <- bwide$species
  bwide$species   <- NULL

  cat_line("Loading BIEN tree cache ...")
  tree_pkg <- readRDS(cache_br)
  btree    <- tree_pkg$scenario.3
  # V.PhyloMaker2 stores tips as "Genus_species" -> convert to "Genus species"
  btree$tip.label <- gsub("_", " ", btree$tip.label)

  matched <- intersect(rownames(bwide), btree$tip.label)
  if (length(matched) < 50L) {
    cat_line(sprintf("SKIP plants (only %d matched species after cache load)", length(matched)))
    taxa_results$plants <- list(taxon = "plants", skipped = TRUE,
                                  reason = sprintf("only %d matched species", length(matched)),
                                  pass = TRUE)
  } else {
    set.seed(SEED)
    sp <- sample(matched, min(N_PER_TAXON, length(matched)))
    b_df   <- bwide[sp, , drop = FALSE]
    b_tree <- ape::keep.tip(btree, sp)
    b_df   <- b_df[b_tree$tip.label, , drop = FALSE]
    cont_bien <- colnames(b_df)[vapply(b_df, is.numeric, logical(1))]
    taxa_results$plants <- bench_one(
      "plants", b_df, b_tree,
      cont_cols = cont_bien,
      disc_cols = character(0))
  }
} else {
  cat_line("SKIP plants (missing cache at ", cache_bt, " or ", cache_br, ")")
  taxa_results$plants <- list(taxon = "plants", skipped = TRUE,
                                reason = "no cache", pass = TRUE)
}

# ---------------------------------------------------------------------------
# Summary + JSON output
# ---------------------------------------------------------------------------

any_fail <- any(!vapply(taxa_results, function(r) isTRUE(r$pass), logical(1)))
taxa_results$ALL_BENCHES_PASS <- !any_fail
taxa_results$generated_at     <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

json_path <- "script/regress_result.json"
jsonlite::write_json(taxa_results, json_path, pretty = TRUE, auto_unbox = TRUE)

cat_line(if (any_fail) "FAIL" else "PASS",
          " --- see ", json_path)
quit(status = if (any_fail) 1L else 0L)
