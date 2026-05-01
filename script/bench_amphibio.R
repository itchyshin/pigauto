#!/usr/bin/env Rscript
#
# script/bench_amphibio.R
#
# Fourth real-data benchmark in the vertebrate-breadth series:
# amphibians. Completes the tetrapod story (birds via AVONET,
# mammals via PanTHERIA, fish via FishBase, amphibians via
# AmphiBIO).
#
# Data
#   AmphiBIO (Oliveira et al. 2017, Sci. Data 4: 170123)
#   DOI: 10.1038/sdata.2017.123
#   CSV: https://ndownloader.figshare.com/files/8828866
#   ~6,800 species with body size, activity, diet, habitat,
#   reproductive mode, etc.  Columns are a mix of continuous and
#   categorical.
#
# Phylogeny
#   Taxonomic tree built from AmphiBIO's Order / Family / Genus /
#   Species columns via ape::as.phylo.formula(), with Grafen
#   rank-based branch lengths via ape::compute.brlen(method = "Grafen").
#   This mirrors the PanTHERIA bench (feature/bench-pantheria) --
#   we hit an Rcpp version conflict with rotl so we fall back to
#   taxonomic tree construction.  For pigauto's phylogenetic-graph
#   baseline this is adequate (the graph is built from cophenetic
#   distances, not specific speciation times).
#
# Traits (6 mixed types, modeled after AmphiBIO's curated columns):
#   Body_size_mm    numeric (SVL or total length, mm) continuous (auto-log)
#   Body_mass_g     numeric                            continuous (auto-log)
#   Diu             factor (Diurnal T/F)               binary
#   Noc             factor (Nocturnal T/F)             binary
#   Diet_breadth    integer (1-5)                      ordinal
#   Habitat         factor (Aquatic/Terrestrial/...)   categorical
#
# Config: 30% MCAR held-out, seed 2026, n_imputations = 20,
# epochs = 500. Median-pool MI fix (dc8cffa) active.
#
# Output:
#   script/bench_amphibio.rds  tidy results + timings
#   script/bench_amphibio.md   human-readable summary

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
N_IMP     <- {
  v <- Sys.getenv("PIGAUTO_N_IMPUTATIONS", unset = "")
  if (nzchar(v)) as.integer(v) else 20L
}
N_SUB <- {
  v <- Sys.getenv("PIGAUTO_N_SUBSET", unset = "")
  if (nzchar(v)) as.integer(v) else NA_integer_
}

here <- "/Users/z3437171/Dropbox/Github Local/pigauto"
cache_dir <- file.path(here, "script", "data-cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
zip_path <- file.path(cache_dir, "amphibio_v1.zip")
csv_path <- file.path(cache_dir, "AmphiBIO_v1.csv")

out_rds <- file.path(here, "script", "bench_amphibio.rds")
out_md  <- file.path(here, "script", "bench_amphibio.md")

script_start <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start),
      ..., "\n", sep = "")
  flush.console()
}

# -------------------------------------------------------------------------
# 1. Load AmphiBIO CSV (download on first run, cached thereafter)
# -------------------------------------------------------------------------

if (!file.exists(csv_path)) {
  if (!file.exists(zip_path)) {
    log_line("Downloading AmphiBIO_v1.zip from figshare ...")
    # figshare article 4644424, file id 8828578 (AmphiBIO_v1.zip, 1.4 MB).
    # See figshare API: https://api.figshare.com/v2/articles/4644424
    download.file("https://ndownloader.figshare.com/files/8828578",
                  zip_path, mode = "wb", quiet = TRUE)
  }
  log_line("Extracting AmphiBIO_v1.csv from zip ...")
  utils::unzip(zip_path, files = "AmphiBIO_v1.csv", exdir = cache_dir)
  if (!file.exists(csv_path)) {
    # Try listing zip contents in case filename differs
    zc <- utils::unzip(zip_path, list = TRUE)
    stop("AmphiBIO_v1.csv not found in zip. Contents:\n",
         paste(zc$Name, collapse = "\n"), call. = FALSE)
  }
}

amphi <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
log_line(sprintf("Loaded AmphiBIO: %d rows x %d cols", nrow(amphi), ncol(amphi)))
log_line("First 20 col names:")
print(head(names(amphi), 20))

# -------------------------------------------------------------------------
# 2. Pick a subset of traits (mixed types)
# -------------------------------------------------------------------------

# AmphiBIO canonical columns (Oliveira et al. 2017 Table 1):
#   Order, Family, Genus, Species, Diu, Noc, Fos, Ter, Aqu, Arb, Ver, Mar,
#   Body_size_mm, Body_mass_g, Reproductive_output_y, Offspring_size_mm,
#   Breeding_interval, Habitat, Diet, Diet_categ, Observed_max_age,
#   Age_at_maturity_min_y, Age_at_maturity_max_y, Egg_size_mm, ...
#
# Pick a tractable mixed-type subset:

pick_traits <- c("Body_size_mm", "Body_mass_g")
# Continuous-only for v1. Binary traits (Diu/Noc/Fos/Ter/Aqu/Arb)
# and diet items currently trigger a character/double type mismatch
# inside pigauto's threshold-joint baseline at 5000+ species with
# the AmphiBIO data -- likely a factor encoding edge case. Tracked
# for investigation; continuous-only is a clean first pass.

# If any picked trait is missing from the CSV, drop it with a warning.
pick_traits <- intersect(pick_traits, names(amphi))
log_line("Selected traits: ", paste(pick_traits, collapse = ", "))

tax_cols <- intersect(c("Order", "Family", "Genus", "Species"),
                       names(amphi))
stopifnot(all(c("Order", "Family", "Genus", "Species") %in% tax_cols))

# -------------------------------------------------------------------------
# 3. Clean + coerce types
# -------------------------------------------------------------------------

df <- amphi[, c(tax_cols, pick_traits), drop = FALSE]

# Drop rows with no species name
df <- df[nzchar(df$Species) & !is.na(df$Species), , drop = FALSE]

# Collapse intraspecific duplicates (AmphiBIO has some): keep first.
df <- df[!duplicated(df$Species), , drop = FALSE]

# Continuous columns -> numeric
for (v in c("Body_size_mm", "Body_mass_g")) {
  if (v %in% names(df)) df[[v]] <- suppressWarnings(as.numeric(df[[v]]))
}

# Binary columns (Diu, Noc) -> factor with 2 levels
for (v in c("Diu", "Noc", "Fos", "Ter", "Aqu", "Arb")) {
  if (v %in% names(df)) {
    # AmphiBIO uses 0/1 integer or T/F strings -- normalise to factor c("no","yes")
    raw <- df[[v]]
    raw_chr <- as.character(raw)
    raw_chr[raw_chr %in% c("1", "TRUE", "T", "yes", "Yes")] <- "yes"
    raw_chr[raw_chr %in% c("0", "FALSE", "F", "no", "No")]  <- "no"
    raw_chr[!raw_chr %in% c("yes", "no")] <- NA
    df[[v]] <- factor(raw_chr, levels = c("no", "yes"))
  }
}

# Ordinal: Diet_breadth (integer, typically 1-5)
if ("Diet_breadth" %in% names(df)) {
  raw <- suppressWarnings(as.integer(df$Diet_breadth))
  raw[!is.finite(raw) | raw < 1] <- NA
  df$Diet_breadth <- ordered(raw)
}

# Categorical: Habitat (string categories)
if ("Habitat" %in% names(df)) {
  h <- as.character(df$Habitat)
  h[!nzchar(h) | h %in% c("NA", "na", "-", "?")] <- NA
  # Collapse rare levels (< 50 species) to NA
  tab <- sort(table(h), decreasing = TRUE)
  keep <- names(tab)[tab >= 50L]
  h[!h %in% keep] <- NA
  df$Habitat <- factor(h, levels = keep)
}

# Keep only species with >= 1 observed trait (excluding taxonomy cols)
trait_cols <- setdiff(names(df), tax_cols)
has_any <- rowSums(!is.na(df[, trait_cols, drop = FALSE])) > 0L
df <- df[has_any, , drop = FALSE]

# Drop rows with incomplete taxonomy (needed for tree)
tax_ok <- rowSums(is.na(df[, tax_cols, drop = FALSE])) == 0 &
           nzchar(df$Order) & nzchar(df$Family) & nzchar(df$Genus)
df <- df[tax_ok, , drop = FALSE]

# Subset sample if requested
if (!is.na(N_SUB) && N_SUB < nrow(df)) {
  set.seed(SEED)
  df <- df[sample.int(nrow(df), N_SUB), , drop = FALSE]
  log_line(sprintf("N_SUBSET = %d: random subset -> %d species",
                    N_SUB, nrow(df)))
}

rownames(df) <- df$Species
log_line(sprintf("After cleaning: %d species x %d traits", nrow(df),
                  length(trait_cols)))

# -------------------------------------------------------------------------
# 4. Build taxonomic tree via as.phylo.formula
# -------------------------------------------------------------------------

log_line("Building taxonomic tree from Order/Family/Genus/Species ...")
tax_df <- df[, tax_cols, drop = FALSE]
tax_df[] <- lapply(tax_df, factor)
tree <- as.phylo(~Order/Family/Genus/Species, data = tax_df,
                  collapse = FALSE)
# Rphylopars (BM baseline) requires a fully dichotomous, rooted tree.
# as.phylo.formula() produces polytomies and singleton (single-child)
# internal nodes.  Pipeline:
#   1. collapse.singles() -- remove single-child internal nodes
#   2. root() if not already rooted
#   3. multi2di(random=TRUE) -- randomly resolve polytomies to binary
#   4. compute.brlen("Grafen") -- rank-based branch lengths
set.seed(SEED)
tree <- ape::collapse.singles(tree)
if (!ape::is.rooted(tree)) {
  tree <- ape::root.phylo(tree, outgroup = 1L, resolve.root = TRUE)
}
tree <- ape::multi2di(tree, random = TRUE)
tree <- compute.brlen(tree, method = "Grafen")
# Some internal edges may still be zero after Grafen if the tree has
# degenerate structure; bump them to tiny positive.
tree$edge.length[tree$edge.length <= 0] <- 1e-8
log_line(sprintf("Tree after resolve: rooted=%s, dichotomous=%s, n_edges=%d, any_zero_edge=%s",
                  ape::is.rooted(tree), ape::is.binary(tree),
                  nrow(tree$edge),
                  any(tree$edge.length <= 0)))
stopifnot(all(rownames(df) %in% tree$tip.label))
tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(df)))
df   <- df[match(tree$tip.label, rownames(df)), , drop = FALSE]
log_line(sprintf("Tree: %d tips, max depth = %.3f",
                  length(tree$tip.label), max(node.depth.edgelength(tree))))

# Keep only trait columns for pigauto
df <- df[, trait_cols, drop = FALSE]

# -------------------------------------------------------------------------
# 5. MCAR 30% held-out mask
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

# -------------------------------------------------------------------------
# 6. Eval helper
# -------------------------------------------------------------------------

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
    t_v <- truth[[v]][idx]
    c_v <- completed[[v]][idx]
    if (is.factor(t_v) || is.ordered(t_v) || is.character(t_v)) {
      acc <- mean(as.character(c_v) == as.character(t_v), na.rm = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(
        method = method, trait = v, metric = "accuracy",
        value = acc, n_cells = length(idx), wall_s = wall_s)
    } else {
      rmse <- sqrt(mean((as.numeric(t_v) - as.numeric(c_v))^2, na.rm = TRUE))
      pear <- safe_cor(as.numeric(t_v), as.numeric(c_v))
      rows[[length(rows) + 1L]] <- data.frame(
        method = method, trait = v, metric = "rmse",
        value = rmse, n_cells = length(idx), wall_s = wall_s)
      rows[[length(rows) + 1L]] <- data.frame(
        method = method, trait = v, metric = "pearson_r",
        value = pear, n_cells = length(idx), wall_s = wall_s)
      t_num <- as.numeric(t_v)
      # Conformal coverage (from pigauto's $conformal_lower / upper,
      # present on all pigauto fits regardless of n_imputations).
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
      # MC-dropout coverage (needs n_imputations > 1)
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
  }
  if (!length(rows)) NULL else do.call(rbind, rows)
}

# -------------------------------------------------------------------------
# 7. Methods
# -------------------------------------------------------------------------

run_mean <- function() {
  out <- df_miss
  for (v in names(out)) {
    if (is.factor(out[[v]])) {
      mode <- names(sort(table(out[[v]], useNA = "no"), decreasing = TRUE))[1]
      if (!length(mode)) next
      out[[v]][mask_test[, v]] <- factor(mode, levels = levels(out[[v]]),
                                           ordered = is.ordered(out[[v]]))
    } else {
      out[[v]][mask_test[, v]] <- mean(out[[v]], na.rm = TRUE)
    }
  }
  list(completed = out)
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
# 8. Run
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
  n_subset  = N_SUB,
  n_imputations = N_IMP
), out_rds)

md <- c(
  "# AmphiBIO x pigauto",
  "",
  sprintf("n = %d amphibian species x %d traits (AmphiBIO + taxonomic tree).",
          nrow(df), ncol(df)),
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
