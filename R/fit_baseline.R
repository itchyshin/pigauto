#' Fit a phylogenetic BM baseline using Rphylopars
#'
#' Fits a Brownian Motion (or OU) model with \pkg{Rphylopars} and returns
#' imputed means and standard errors for every species.
#'
#' @details
#' When \code{splits} is supplied the val and test cells are masked to
#' \code{NA} before fitting, so the baseline is evaluated under the same
#' conditions as \code{\link{fit_pigauto}}.
#'
#' For mixed-type data, Rphylopars is applied only to continuous, count,
#' and ordinal traits (which are continuous in latent space).  Binary
#' and categorical traits use phylogenetic label propagation: each
#' species receives a personalised baseline computed as the
#' phylo-similarity-weighted average of observed values, using a
#' Gaussian kernel on cophenetic distances.
#'
#' @param data object of class \code{"pigauto_data"}.
#' @param tree object of class \code{"phylo"}.
#' @param splits list (output of \code{\link{make_missing_splits}}) or
#'   \code{NULL}.
#' @param model character. Evolutionary model: \code{"BM"} (default) or
#'   \code{"OU"}.
#' @return A list with:
#'   \describe{
#'     \item{mu}{Numeric matrix (n_species x p_latent), baseline means in
#'       latent scale.}
#'     \item{se}{Numeric matrix (n_species x p_latent), standard errors.}
#'   }
#' @examples
#' \dontrun{
#' data(avonet300, tree300, package = "pigauto")
#' traits <- avonet300; rownames(traits) <- traits$Species_Key
#' traits$Species_Key <- NULL
#' pd     <- preprocess_traits(traits, tree300)
#' splits <- make_missing_splits(pd$X_scaled, trait_map = pd$trait_map)
#' bl     <- fit_baseline(pd, tree300, splits)
#' }
#' @importFrom stats complete.cases
#' @export
fit_baseline <- function(data, tree, splits = NULL, model = "BM") {
  if (!inherits(data, "pigauto_data")) {
    stop("'data' must be a pigauto_data object (output of preprocess_traits).")
  }
  if (!inherits(tree, "phylo")) stop("'tree' must be a phylo object.")

  X   <- data$X_scaled
  p   <- ncol(X)

  multi_obs <- isTRUE(data$multi_obs)
  if (multi_obs) {
    # Multi-obs: X is n_obs x p, species_names is n_species
    n_obs     <- data$n_obs
    n_species <- data$n_species
    spp       <- data$species_names          # unique species (n_species)
    obs_spp   <- data$obs_species            # species per obs (n_obs)
    obs_to_sp <- data$obs_to_species         # integer mapping (n_obs)
  } else {
    n_obs     <- nrow(X)
    n_species <- n_obs
    spp       <- data$species_names
    obs_spp   <- spp
    obs_to_sp <- NULL
  }

  # ---- Phylogenetic similarity for discrete-trait label propagation ------
  D_phylo  <- ape::cophenetic.phylo(tree)
  # Reorder to match species order
  D_phylo  <- D_phylo[spp, spp]
  sigma_lp <- stats::median(D_phylo) * 0.5
  sim_phylo <- exp(-(D_phylo^2) / (2 * sigma_lp^2))
  diag(sim_phylo) <- 0  # exclude self for label propagation

  # Mask val + test cells before fitting
  if (!is.null(splits)) {
    X[splits$val_idx]  <- NA
    X[splits$test_idx] <- NA
  }

  trait_map <- data$trait_map

  # Output at species level (n_species x p)
  mu <- matrix(0, nrow = n_species, ncol = p)
  se <- matrix(0, nrow = n_species, ncol = p)
  dimnames(mu) <- list(spp, data$latent_names)
  dimnames(se) <- list(spp, data$latent_names)

  # ---- Identify BM-eligible columns (continuous in latent space) -----------
  bm_cols <- integer(0)
  for (tm in trait_map) {
    if (tm$type %in% c("continuous", "count", "ordinal")) {
      bm_cols <- c(bm_cols, tm$latent_cols)
    }
  }

  # ---- Fit Rphylopars on BM-eligible columns --------------------------------
  if (length(bm_cols) > 0) {
    X_bm <- X[, bm_cols, drop = FALSE]
    df_ph <- as.data.frame(X_bm)
    # Use observation-level species names (allows duplicates for multi-obs)
    df_ph <- cbind(species = obs_spp, df_ph)

    fit_rphylo <- function(tr) {
      suppressWarnings(
        Rphylopars::phylopars(
          trait_data  = df_ph,
          tree        = tr,
          model       = model,
          pheno_error = FALSE
        )
      )
    }

    fit <- tryCatch(
      fit_rphylo(tree),
      error = function(e) {
        message("Rphylopars failed (", conditionMessage(e),
                "). Retrying with branch-length jitter.")
        tree_j <- tree
        tree_j$edge.length <- tree_j$edge.length + 1e-6
        tryCatch(
          fit_rphylo(tree_j),
          error = function(e2) {
            stop("Rphylopars failed after jitter retry: ",
                 conditionMessage(e2))
          }
        )
      }
    )

    bm_names <- colnames(X_bm)
    mu_bm  <- as.matrix(fit$anc_recon[spp, bm_names, drop = FALSE])
    var_bm <- as.matrix(fit$anc_var[spp,  bm_names, drop = FALSE])
    se_bm  <- sqrt(pmax(var_bm, 0))

    mu[, bm_cols] <- mu_bm
    se[, bm_cols] <- se_bm
  }

  # ---- Binary baseline: phylogenetic label propagation -------------------
  for (tm in trait_map) {
    if (tm$type != "binary") next
    lc   <- tm$latent_cols

    # Get species-level observations
    if (multi_obs) {
      sp_vals <- tapply(X[, lc], obs_spp, function(v) {
        v <- v[!is.na(v)]
        if (length(v) == 0) NA_real_ else mean(v)
      })
      vals_species <- rep(NA_real_, n_species)
      names(vals_species) <- spp
      vals_species[names(sp_vals)] <- as.numeric(sp_vals)
    } else {
      vals_species <- X[, lc]
      names(vals_species) <- spp
    }

    observed <- !is.na(vals_species)
    if (sum(observed) == 0) {
      mu[, lc] <- logit(0.5)
      se[, lc] <- 0
      next
    }

    # Phylo-weighted probability for each species
    sim_obs <- sim_phylo[, observed, drop = FALSE]
    row_weights <- rowSums(sim_obs)
    row_weights[row_weights < 1e-10] <- 1e-10
    probs <- as.numeric(sim_obs %*% vals_species[observed]) / row_weights
    probs <- pmin(pmax(probs, 0.01), 0.99)  # clip for stability

    mu[, lc] <- logit(probs)
    se[, lc] <- 0
  }

  # ---- Categorical baseline: phylogenetic label propagation ---------------
  for (tm in trait_map) {
    if (tm$type != "categorical") next
    K    <- tm$n_latent
    lc   <- tm$latent_cols
    oh   <- X[, lc, drop = FALSE]  # n_obs x K one-hot (with NAs)

    # Get species-level one-hot observations
    if (multi_obs) {
      # Average one-hot within species (handles multiple obs)
      oh_species <- matrix(NA_real_, n_species, K)
      for (s in seq_len(n_species)) {
        rows <- which(obs_to_sp == s)
        obs_rows <- rows[complete.cases(oh[rows, , drop = FALSE])]
        if (length(obs_rows) > 0) {
          oh_species[s, ] <- colMeans(oh[obs_rows, , drop = FALSE])
        }
      }
    } else {
      oh_species <- oh
    }

    # Which species have observed values
    observed <- complete.cases(oh_species)
    if (sum(observed) == 0) {
      freqs <- rep(1 / K, K)
      log_freqs <- log(freqs)
      for (k in seq_len(K)) mu[, lc[k]] <- log_freqs[k]
      se[, lc] <- 0
      next
    }

    # Phylo-weighted category probabilities per species
    sim_obs <- sim_phylo[, observed, drop = FALSE]
    row_weights <- rowSums(sim_obs)
    row_weights[row_weights < 1e-10] <- 1e-10

    # Weighted category probs: (n_species x K)
    weighted_probs <- (sim_obs %*% oh_species[observed, , drop = FALSE]) /
      row_weights
    # Add small floor and renormalise
    weighted_probs <- pmax(weighted_probs, 1e-6)
    weighted_probs <- weighted_probs / rowSums(weighted_probs)

    for (k in seq_len(K)) {
      mu[, lc[k]] <- log(weighted_probs[, k])
    }
    se[, lc] <- 0
  }

  list(mu = mu, se = se)
}
