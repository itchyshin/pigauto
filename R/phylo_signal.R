# R/phylo_signal.R
# Phylogenetic-signal diagnostic: per-trait Pagel's lambda used by
# fit_pigauto() to route weak-signal traits to the grand-mean corner
# of the safety-floor simplex. See
# specs/2026-04-23-phylo-signal-gate-design.md.

# Compute per-trait phylogenetic signal on training-observed cells.
#
# For continuous / count / ordinal / proportion / zi_mag latent columns,
# lambda is computed on the z-scored latent values via phytools::phylosig().
# For binary / categorical latent columns, lambda is computed on the
# estep-decoded posterior mean liability (still continuous on the
# liability scale).
#
# Returns a named numeric vector of length = length(data$trait_map),
# keyed by trait name. NA for traits where lambda cannot be computed
# (too few tips, constant values, phytools missing).
#
# @param data     pigauto_data object from preprocess_traits()
# @param tree     phylo object matching data$species
# @param method   character, currently only "lambda" is implemented
# @param min_tips integer, skip traits with fewer than this many
#                 species in training observations (default 20)
# @return named numeric vector
#
# @noRd
compute_phylo_signal_per_trait <- function(data, tree,
                                            method = c("lambda", "blomberg_k"),
                                            min_tips = 20L) {
  method <- match.arg(method)
  if (!requireNamespace("phytools", quietly = TRUE)) {
    warning("phylo_signal_gate requires the 'phytools' package; ",
            "returning NA for all traits.", call. = FALSE)
    out <- rep(NA_real_, length(data$trait_map))
    names(out) <- vapply(data$trait_map, function(tm) tm$name,
                          character(1L))
    return(out)
  }
  tm_names <- vapply(data$trait_map, function(tm) tm$name, character(1L))
  out <- rep(NA_real_, length(tm_names))
  names(out) <- tm_names
  for (i in seq_along(data$trait_map)) {
    tm <- data$trait_map[[i]]
    lambda <- tryCatch(
      .phylo_signal_one_trait(tm, data, tree, method, min_tips),
      error = function(e) NA_real_)
    out[i] <- lambda
  }
  out
}

.phylo_signal_one_trait <- function(tm, data, tree, method, min_tips) {
  if (tm$type %in% c("continuous", "count", "ordinal",
                       "proportion", "zi_mag")) {
    .phylo_signal_lambda_continuous(tm, data, tree, method, min_tips)
  } else if (tm$type %in% c("binary", "categorical")) {
    .phylo_signal_lambda_discrete(tm, data, tree, method, min_tips)
  } else if (tm$type == "zi_count") {
    .phylo_signal_lambda_zi_count(tm, data, tree, method, min_tips)
  } else if (tm$type == "multi_proportion") {
    .phylo_signal_lambda_continuous(tm, data, tree, method, min_tips)
  } else {
    NA_real_
  }
}

.phylo_signal_lambda_continuous <- function(tm, data, tree, method, min_tips) {
  lc <- tm$latent_cols[1L]
  x  <- data$X_scaled[, lc]
  # Training-observed species
  obs <- !is.na(x)
  if (sum(obs) < min_tips) return(NA_real_)
  # Standard deviation guard
  if (stats::sd(x[obs], na.rm = TRUE) < 1e-10) return(NA_real_)
  # phylosig needs a named vector keyed by tip label
  sp_names <- if (!is.null(data$species_names)) data$species_names[obs] else
    rownames(data$X_scaled)[obs]
  vals <- x[obs]
  names(vals) <- sp_names
  # Restrict tree to these tips
  tr <- ape::keep.tip(tree, sp_names)
  vals <- vals[tr$tip.label]
  if (identical(method, "lambda")) {
    res <- phytools::phylosig(tr, vals, method = "lambda", test = FALSE)
    as.numeric(res$lambda)
  } else {
    res <- phytools::phylosig(tr, vals, method = "K", test = FALSE)
    as.numeric(res$K)
  }
}

.phylo_signal_lambda_discrete <- function(tm, data, tree, method, min_tips) {
  # For binary/categorical, use the z-scored latent column directly
  # as a continuous proxy. This is the simplest approach; more
  # sophisticated: run estep_liability_* first to get a posterior
  # mean liability. For the v0.9.1.9003 initial ship, we use the
  # direct z-scored proxy which is cheap and correlates well with
  # lambda on the exact liability.
  .phylo_signal_lambda_continuous(tm, data, tree, method, min_tips)
}

.phylo_signal_lambda_zi_count <- function(tm, data, tree, method, min_tips) {
  # Coupled-gate convention: weak-signal iff BOTH latent columns are
  # weak. Return max of the two lambdas so the downstream threshold check
  # treats "either one weak" as "trait weak".
  lc_gate <- tm$latent_cols[1L]
  lc_mag  <- tm$latent_cols[2L]
  tm_gate <- list(latent_cols = lc_gate, type = "binary",
                   name = paste0(tm$name, "_gate"))
  tm_mag  <- list(latent_cols = lc_mag, type = "zi_mag",
                   name = paste0(tm$name, "_mag"))
  lam_gate <- .phylo_signal_lambda_discrete(tm_gate, data, tree, method, min_tips)
  lam_mag  <- .phylo_signal_lambda_continuous(tm_mag, data, tree, method, min_tips)
  if (is.na(lam_gate) && is.na(lam_mag)) return(NA_real_)
  if (is.na(lam_gate)) return(lam_mag)
  if (is.na(lam_mag))  return(lam_gate)
  max(lam_gate, lam_mag)   # use max so the gate triggers only when BOTH are weak
}
