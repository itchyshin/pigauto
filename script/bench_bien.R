#!/usr/bin/env Rscript
#
# script/bench_bien.R
#
# Plant trait imputation -- the kingdom jump for the pigauto
# vertebrate-breadth series (birds + mammals + fish + amphibians +
# now plants).
#
# Data
#   BIEN (Botanical Information and Ecology Network) v4 via the
#   `BIEN` R package on CRAN. ~400k species globally, curated traits.
#   We pull species-level means via BIEN_trait_mean for a small set
#   of well-populated continuous traits.
#
# Tree
#   V.PhyloMaker2 (Jin & Qian 2022 J. Plant Ecology) from GitHub
#   `jinyizju/V.PhyloMaker2`.  Builds a megaphylogeny-backed
#   phylogeny for any seed-plant species list using the Smith &
#   Brown 2018 backbone (~74k taxa).
#
# Traits (continuous-only v1, mirroring amphibian bench rationale --
# threshold-joint baseline character/double issue with binaries
# parked for v0.9.2):
#   maximum whole plant height                    (m)
#   leaf area                                      (mm^2)
#   leaf area per leaf dry mass                    (SLA, mm^2 / mg)
#   seed mass                                      (mg)
#   stem wood density                              (g / cm^3)
#
# Run
#   PIGAUTO_PKG_PATH=$(pwd) Rscript script/bench_bien.R
# Optional env vars:
#   PIGAUTO_BIEN_N_SPECIES=NNNN   limit species count after intersect
#   PIGAUTO_N_IMPUTATIONS=K       defaults to 1 (predict-leak workaround)
#
# Output:
#   script/bench_bien.rds
#   script/bench_bien.md

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  library(BIEN)
  library(V.PhyloMaker2)
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
N_LIMIT <- {
  v <- Sys.getenv("PIGAUTO_BIEN_N_SPECIES", unset = "")
  if (nzchar(v)) as.integer(v) else NA_integer_
}
N_IMP <- {
  v <- Sys.getenv("PIGAUTO_N_IMPUTATIONS", unset = "")
  if (nzchar(v)) as.integer(v) else 1L
}

here <- "/Users/z3437171/Dropbox/Github Local/pigauto"
cache_dir <- file.path(here, "script", "data-cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
cache_traits <- file.path(cache_dir, "bien_trait_means.rds")
cache_tree   <- file.path(cache_dir, "bien_tree.rds")

out_rds <- file.path(here, "script", "bench_bien.rds")
out_md  <- file.path(here, "script", "bench_bien.md")

script_start <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start),
      ..., "\n", sep = "")
  flush.console()
}

# -------------------------------------------------------------------------
# 1. Pull species-level trait means from BIEN
# -------------------------------------------------------------------------

trait_names <- c(
  "maximum whole plant height",
  "leaf area",
  "leaf area per leaf dry mass",
  "seed mass",
  "stem wood density"
)
trait_short <- c(height_m = "maximum whole plant height",
                  leaf_area = "leaf area",
                  sla = "leaf area per leaf dry mass",
                  seed_mass = "seed mass",
                  wood_density = "stem wood density")

if (file.exists(cache_traits)) {
  log_line("Loading cached BIEN trait means from ", cache_traits)
  trait_means <- readRDS(cache_traits)
} else {
  log_line("Querying BIEN_trait_trait() for 5 traits (slow on first run) ...")
  trait_dfs <- lapply(trait_names, function(tn) {
    log_line(sprintf("  fetching '%s' ...", tn))
    obs <- tryCatch(BIEN::BIEN_trait_trait(trait = tn),
                     error = function(e) {
                       log_line("    error: ", conditionMessage(e))
                       NULL
                     })
    if (is.null(obs) || !nrow(obs)) return(NULL)
    log_line(sprintf("    got %d obs", nrow(obs)))
    # Aggregate observation-level rows to species means.
    # BIEN obs cols typically include: scrubbed_species_binomial,
    # trait_value, trait_name, ...
    sp_col <- intersect(c("scrubbed_species_binomial", "species",
                            "verbatim_species"), names(obs))[1]
    val_col <- intersect(c("trait_value", "value"), names(obs))[1]
    if (is.na(sp_col) || is.na(val_col)) {
      log_line(sprintf("    cols missing; have: %s",
                        paste(names(obs), collapse = ", ")))
      return(NULL)
    }
    val_num <- suppressWarnings(as.numeric(obs[[val_col]]))
    keep <- is.finite(val_num) & nzchar(obs[[sp_col]])
    if (!any(keep)) return(NULL)
    means <- tapply(val_num[keep], obs[[sp_col]][keep], mean,
                     na.rm = TRUE)
    data.frame(species = names(means),
                mean_value = as.numeric(means),
                stringsAsFactors = FALSE)
  })
  names(trait_dfs) <- names(trait_short)
  saveRDS(trait_dfs, cache_traits)
  trait_means <- trait_dfs
}

# -------------------------------------------------------------------------
# 2. Pivot to species x trait wide table
# -------------------------------------------------------------------------

species_lists <- lapply(trait_means, function(d) {
  if (is.null(d)) return(character(0))
  d$species
})
all_species <- Reduce(union, species_lists)
log_line(sprintf("Union of all trait species: %d", length(all_species)))

# Build wide matrix
wide <- data.frame(species = all_species)
for (nm in names(trait_means)) {
  d <- trait_means[[nm]]
  if (is.null(d)) {
    wide[[nm]] <- NA_real_
    next
  }
  # BIEN_trait_mean returns columns: species, mean_value, sd_value, ...
  val_col <- intersect(c("mean_value", "trait_value", "value"), names(d))[1]
  if (is.na(val_col)) {
    log_line(sprintf("  %s: no value column found (cols: %s)",
                      nm, paste(names(d), collapse = ", ")))
    wide[[nm]] <- NA_real_
    next
  }
  m <- match(wide$species, d$species)
  wide[[nm]] <- suppressWarnings(as.numeric(d[[val_col]][m]))
}

# Require at least one trait observed per species
trait_cols <- names(trait_short)
n_obs_per_sp <- rowSums(!is.na(wide[, trait_cols, drop = FALSE]))
wide <- wide[n_obs_per_sp >= 1L, , drop = FALSE]
log_line(sprintf("After filtering (>=1 trait): %d species", nrow(wide)))

# Optionally limit
if (!is.na(N_LIMIT) && N_LIMIT < nrow(wide)) {
  set.seed(SEED)
  wide <- wide[sample.int(nrow(wide), N_LIMIT), , drop = FALSE]
  log_line(sprintf("PIGAUTO_BIEN_N_SPECIES = %d: random subset",
                    N_LIMIT))
}

# -------------------------------------------------------------------------
# 3. Build phylogeny via V.PhyloMaker2
# -------------------------------------------------------------------------

if (file.exists(cache_tree)) {
  log_line("Loading cached V.PhyloMaker2 tree ...")
  tree_pkg <- readRDS(cache_tree)
} else {
  log_line("Building V.PhyloMaker2 tree (slow on first run) ...")
  # V.PhyloMaker2 needs a species_list df with columns: species, genus, family
  sp_clean <- gsub("\\s+", " ", trimws(wide$species))
  parts <- strsplit(sp_clean, " ")
  genus <- vapply(parts, function(x) x[1], character(1))
  spec_df <- data.frame(species = sp_clean, genus = genus,
                         family = NA_character_)  # phylo.maker fills family
  tree_pkg <- V.PhyloMaker2::phylo.maker(sp.list = spec_df,
                                           tree = V.PhyloMaker2::GBOTB.extended.TPL,
                                           nodes = V.PhyloMaker2::nodes.info.1.TPL,
                                           scenarios = "S3")
  saveRDS(tree_pkg, cache_tree)
}

tree <- tree_pkg$scenario.3
# V.PhyloMaker2 returns species labels with underscores -- convert
tree$tip.label <- gsub("_", " ", tree$tip.label)
log_line(sprintf("Tree: %d tips after V.PhyloMaker2", length(tree$tip.label)))

# Match wide table to tree
keep <- intersect(wide$species, tree$tip.label)
wide <- wide[wide$species %in% keep, , drop = FALSE]
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
wide <- wide[match(tree$tip.label, wide$species), , drop = FALSE]
rownames(wide) <- wide$species
df <- wide[, trait_cols, drop = FALSE]
log_line(sprintf("Matched: %d species x %d traits", nrow(df), ncol(df)))

# -------------------------------------------------------------------------
# 4. MCAR mask + eval helpers (same pattern as fish/amphibian)
# -------------------------------------------------------------------------

set.seed(SEED)
df_truth <- df
mask_test <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df),
                     dimnames = list(rownames(df), names(df)))
for (v in names(df)) {
  obs_idx <- which(!is.na(df[[v]]))
  if (!length(obs_idx)) next
  to_hide <- sample(obs_idx, ceiling(length(obs_idx) * MISS_FRAC))
  mask_test[to_hide, v] <- TRUE
}
df_miss <- df
for (v in names(df)) df_miss[[v]][mask_test[, v]] <- NA

log_line("Held-out test cells per trait:")
print(colSums(mask_test))

safe_cor <- function(x, y) {
  idx <- which(is.finite(x) & is.finite(y))
  if (length(idx) < 2L) return(NA_real_)
  if (stats::sd(x[idx]) == 0 || stats::sd(y[idx]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[idx], y[idx]))
}

eval_completed <- function(completed, truth, mask, method, wall_s,
                            res_obj = NULL) {
  if (is.null(completed)) return(NULL)
  lo      <- if (!is.null(res_obj)) res_obj$prediction$conformal_lower else NULL
  hi      <- if (!is.null(res_obj)) res_obj$prediction$conformal_upper else NULL
  mi_list <- if (!is.null(res_obj)) res_obj$prediction$imputed_datasets else NULL
  rows <- list()
  for (v in colnames(mask)) {
    idx <- which(mask[, v])
    if (!length(idx)) next
    t_v <- truth[[v]][idx]; c_v <- completed[[v]][idx]
    rmse <- sqrt(mean((as.numeric(t_v) - as.numeric(c_v))^2, na.rm = TRUE))
    pear <- safe_cor(as.numeric(t_v), as.numeric(c_v))
    rows[[length(rows) + 1L]] <- data.frame(
      method = method, trait = v, metric = "rmse",
      value = rmse, n_cells = length(idx), wall_s = wall_s)
    rows[[length(rows) + 1L]] <- data.frame(
      method = method, trait = v, metric = "pearson_r",
      value = pear, n_cells = length(idx), wall_s = wall_s)
    t_num <- as.numeric(t_v)
    if (!is.null(lo) && !is.null(hi) && v %in% colnames(lo)) {
      lo_v <- lo[idx, v]; hi_v <- hi[idx, v]
      valid <- is.finite(lo_v) & is.finite(hi_v) & is.finite(t_num)
      if (any(valid)) {
        hits <- t_num[valid] >= lo_v[valid] & t_num[valid] <= hi_v[valid]
        rows[[length(rows) + 1L]] <- data.frame(
          method = method, trait = v, metric = "coverage95_conformal",
          value = mean(hits), n_cells = sum(valid), wall_s = wall_s)
      }
    }
    if (!is.null(mi_list) && length(mi_list) > 1L &&
        v %in% names(mi_list[[1]])) {
      draws_mat <- vapply(mi_list, function(d) as.numeric(d[idx, v]),
                           numeric(length(idx)))
      if (!is.matrix(draws_mat))
        draws_mat <- matrix(draws_mat, ncol = length(mi_list))
      if (nrow(draws_mat) == length(idx) && ncol(draws_mat) > 1L) {
        q_lo <- apply(draws_mat, 1L, stats::quantile, probs = 0.025,
                       na.rm = TRUE)
        q_hi <- apply(draws_mat, 1L, stats::quantile, probs = 0.975,
                       na.rm = TRUE)
        valid <- is.finite(q_lo) & is.finite(q_hi) & is.finite(t_num)
        if (any(valid)) {
          hits <- t_num[valid] >= q_lo[valid] & t_num[valid] <= q_hi[valid]
          rows[[length(rows) + 1L]] <- data.frame(
            method = method, trait = v, metric = "coverage95_mcdropout",
            value = mean(hits), n_cells = sum(valid), wall_s = wall_s)
        }
      }
    }
  }
  if (!length(rows)) NULL else do.call(rbind, rows)
}

run_mean <- function() {
  out <- df_miss
  for (v in names(out)) {
    out[[v]][mask_test[, v]] <- mean(out[[v]], na.rm = TRUE)
  }
  list(completed = out, res = NULL)
}

run_pigauto <- function() {
  res <- pigauto::impute(df_miss, tree,
                           log_transform = TRUE,
                           missing_frac  = 0.20,
                           n_imputations = N_IMP,
                           epochs        = 500L,
                           verbose       = TRUE,
                           seed          = SEED)
  list(completed = res$completed, res = res)
}

# -------------------------------------------------------------------------
# 5. Run
# -------------------------------------------------------------------------

timed <- function(tag, fn) {
  log_line("[STAGE] ", tag, " starting")
  t0 <- proc.time()[["elapsed"]]
  val <- fn()
  wall <- proc.time()[["elapsed"]] - t0
  log_line("[STAGE] ", tag, " done in ", round(wall, 1), " s")
  list(val = val, wall = wall)
}

r_mean <- timed("mean_baseline",   run_mean)
r_pig  <- timed("pigauto_default", run_pigauto)

ev_mean <- eval_completed(r_mean$val$completed, df_truth, mask_test,
                            "mean_baseline",   r_mean$wall, NULL)
ev_pig  <- eval_completed(r_pig$val$completed,  df_truth, mask_test,
                            "pigauto_default", r_pig$wall, r_pig$val$res)

all_rows <- do.call(rbind, Filter(Negate(is.null),
                                     list(ev_mean, ev_pig)))

saveRDS(list(
  results   = all_rows,
  seed      = SEED,
  miss_frac = MISS_FRAC,
  n_species = nrow(df),
  n_imputations = N_IMP
), out_rds)

md <- c(
  "# BIEN x V.PhyloMaker2 -- pigauto on plants (kingdom jump)",
  "",
  sprintf("n = %d plant species x %d traits.", nrow(df), ncol(df)),
  sprintf("Seed = %d, miss_frac = %.2f, n_imputations = %d.",
          SEED, MISS_FRAC, N_IMP),
  "",
  "## Per-trait metrics",
  "",
  "```",
  capture.output(print(all_rows, row.names = FALSE, max = 500)),
  "```"
)
writeLines(md, out_md)

log_line("=== DONE ===")
log_line("  rds: ", out_rds)
log_line("  md : ", out_md)
