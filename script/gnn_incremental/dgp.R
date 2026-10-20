# Frozen repeated-observation DGP for the GNN incremental-value study.
# The historical `lambda` variable is a Brownian variance multiplier here,
# not Pagel's lambda.

make_gnn_case <- function(n_species = 500L, obs_per = 5L, n_covs = 10L,
                          variance_multiplier = 0.3,
                          f_type = c("nonlinear", "interactive", "linear"),
                          beta_strength = 1,
                          seed,
                          sp_miss_frac = 0.5,
                          within_miss = 0.2,
                          sigma_res = 0.5) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required to generate the phylogeny.", call. = FALSE)
  }
  f_type <- match.arg(f_type)
  scalar_ok <- function(x) is.numeric(x) && length(x) == 1L &&
    is.finite(x)
  if (!scalar_ok(n_species) || n_species < 20L ||
      n_species != as.integer(n_species)) {
    stop("n_species must be an integer of at least 20.", call. = FALSE)
  }
  if (!scalar_ok(obs_per) || obs_per <= 0 ||
      !scalar_ok(n_covs) || n_covs < 5L || n_covs != as.integer(n_covs)) {
    stop("obs_per must be positive and n_covs must be an integer >= 5.",
         call. = FALSE)
  }
  if (!scalar_ok(variance_multiplier) || variance_multiplier < 0 ||
      !scalar_ok(beta_strength) || !scalar_ok(sigma_res) || sigma_res <= 0) {
    stop("Variance, signal, or residual settings are invalid.", call. = FALSE)
  }
  if (!scalar_ok(sp_miss_frac) || sp_miss_frac <= 0 || sp_miss_frac >= 1 ||
      !scalar_ok(within_miss) || within_miss < 0 || within_miss >= 1) {
    stop("Missingness fractions are outside their allowed ranges.",
         call. = FALSE)
  }
  if (!scalar_ok(seed) || seed < 1 || seed > .Machine$integer.max ||
      seed != as.integer(seed)) {
    stop("seed must be a positive integer.", call. = FALSE)
  }

  set.seed(as.integer(seed))
  tree <- ape::rtree(as.integer(n_species))
  species <- tree$tip.label
  C_tree <- ape::vcv(tree)[species, species, drop = FALSE]
  R_tree <- stats::cov2cor(C_tree)
  phylo_raw <- ape::rTraitCont(
    tree, model = "BM", sigma = sqrt(variance_multiplier), root.value = 0
  )
  phylo_raw <- as.numeric(phylo_raw[species])
  phylo <- phylo_raw / sqrt(diag(C_tree))
  names(phylo) <- species

  n_per_species <- pmin(pmax(stats::rpois(n_species, obs_per), 1L), 15L)
  names(n_per_species) <- species
  species_id <- rep(species, n_per_species)
  species_index <- match(species_id, species)
  n_obs <- length(species_id)
  Z <- matrix(stats::rnorm(n_obs * n_covs), nrow = n_obs, ncol = n_covs)
  colnames(Z) <- paste0("z", seq_len(n_covs))

  f <- switch(f_type,
    linear = {
      coefficients <- c(0.5, 0.3, 0.2, 0.1, -0.1, 0.15, -0.2, 0.05,
                        0.1, -0.15)[seq_len(n_covs)]
      drop(Z %*% coefficients)
    },
    nonlinear = {
      sin(Z[, 1L]) * exp(0.3 * Z[, 2L]) + 0.5 * tanh(Z[, 3L]) +
        0.3 * sin(2 * Z[, 4L]) - 0.2 * cos(Z[, 5L])
    },
    interactive = {
      Z[, 1L] * Z[, 2L] + 0.5 * Z[, 3L]^2 +
        0.3 * Z[, 1L] * Z[, 4L] - 0.2 * Z[, 2L] * Z[, 5L]
    }
  )
  f_sd <- stats::sd(f)
  if (!is.finite(f_sd) || f_sd <= 0) {
    stop("Covariate contribution has zero or invalid variation.",
         call. = FALSE)
  }
  f <- as.numeric((f - mean(f)) / f_sd)
  y_mean <- unname(phylo[species_id]) + beta_strength * f
  y <- y_mean + stats::rnorm(n_obs, mean = 0, sd = sigma_res)

  whole_missing_species <- sample(
    species, size = floor(n_species * sp_miss_frac), replace = FALSE
  )
  full_mask <- species_id %in% whole_missing_species
  partial_species <- setdiff(species, whole_missing_species)
  within_pool <- which(!full_mask)
  within_capacity <- sum(pmax(n_per_species[partial_species] - 1L, 0L))
  n_within <- min(floor(within_miss * sum(!full_mask)), within_capacity)
  within_mask <- rep(FALSE, n_obs)
  within_by_species <- integer(n_species)
  if (n_within > 0L) {
    for (idx in sample(within_pool, length(within_pool), replace = FALSE)) {
      sp <- species_index[idx]
      if (within_by_species[sp] < n_per_species[sp] - 1L) {
        within_mask[idx] <- TRUE
        within_by_species[sp] <- within_by_species[sp] + 1L
      }
      if (sum(within_mask) == n_within) break
    }
  }
  mask <- full_mask | within_mask
  stratum <- ifelse(full_mask, "fully_unobserved_species",
                    ifelse(within_mask, "partly_observed_species", NA_character_))

  df_complete <- data.frame(
    row_id = seq_len(n_obs), species = species_id, y = y,
    Z, check.names = FALSE, stringsAsFactors = FALSE
  )
  df_observed <- df_complete
  df_observed$y[mask] <- NA_real_

  list(
    tree = tree,
    df_complete = df_complete,
    df_observed = df_observed,
    species = species_id,
    whole_missing_mask = full_mask,
    within_missing_mask = within_mask,
    phylo_effect = stats::setNames(phylo, species),
    phylo_effect_raw = stats::setNames(phylo_raw, species),
    phylo_correlation = R_tree,
    mask = mask,
    stratum = stratum,
    fully_unobserved_species = whole_missing_species,
    predictor_names = colnames(Z),
    truth_mean = y_mean,
    meta = list(
      n_species = as.integer(n_species), obs_per = obs_per,
      n_covs = as.integer(n_covs),
      variance_multiplier = variance_multiplier, f_type = f_type,
      beta_strength = beta_strength, sigma_res = sigma_res,
      seed = as.integer(seed), sp_miss_frac = sp_miss_frac,
      within_miss = within_miss
    )
  )
}

make_covariate_basis <- function(data, predictor_names, nonlinear = TRUE) {
  if (!all(predictor_names %in% names(data))) {
    stop("Some requested predictors are missing from data.", call. = FALSE)
  }
  Z <- data[, predictor_names, drop = FALSE]
  if (!all(vapply(Z, is.numeric, logical(1)))) {
    stop("All covariates must be numeric.", call. = FALSE)
  }
  X <- cbind(`(Intercept)` = 1, as.matrix(Z))
  if (isTRUE(nonlinear)) {
    sq <- as.matrix(Z)^2
    colnames(sq) <- paste0(names(Z), "^2")
    X <- cbind(X, sq)
    if (ncol(Z) >= 2L) {
      pairs <- utils::combn(seq_len(ncol(Z)), 2L)
      inter <- vapply(seq_len(ncol(pairs)), function(j) {
        Z[[pairs[1L, j]]] * Z[[pairs[2L, j]]]
      }, numeric(nrow(Z)))
      if (is.null(dim(inter))) inter <- matrix(inter, ncol = 1L)
      colnames(inter) <- paste0(names(Z)[pairs[1L, ]], ":",
                                names(Z)[pairs[2L, ]])
      X <- cbind(X, inter)
    }
  }
  storage.mode(X) <- "double"
  X
}
