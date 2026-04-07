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
#' traits receive the population-level logit proportion as baseline;
#' categorical traits receive log marginal frequencies.
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
  n   <- nrow(X)
  p   <- ncol(X)
  spp <- data$species_names

  # Mask val + test cells before fitting
  if (!is.null(splits)) {
    X[splits$val_idx]  <- NA
    X[splits$test_idx] <- NA
  }

  trait_map <- data$trait_map
  mu <- matrix(0, nrow = n, ncol = p)
  se <- matrix(0, nrow = n, ncol = p)
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
    df_ph <- cbind(species = spp, df_ph)

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

  # ---- Binary baseline: logit of observed proportion ------------------------
  for (tm in trait_map) {
    if (tm$type != "binary") next
    lc   <- tm$latent_cols
    vals <- X[, lc]
    obs  <- vals[!is.na(vals)]
    prop <- if (length(obs) > 0) mean(obs) else 0.5
    mu[, lc] <- logit(prop)
    se[, lc] <- 0
  }

  # ---- Categorical baseline: log marginal frequencies -----------------------
  for (tm in trait_map) {
    if (tm$type != "categorical") next
    K    <- tm$n_latent
    lc   <- tm$latent_cols
    oh   <- X[, lc, drop = FALSE]  # n x K one-hot (with NAs)

    # Compute marginal frequencies from observed species
    obs_rows <- which(complete.cases(oh))
    if (length(obs_rows) > 0) {
      freqs <- colMeans(oh[obs_rows, , drop = FALSE])
      freqs <- pmax(freqs, 1e-6)  # avoid log(0)
      freqs <- freqs / sum(freqs)
    } else {
      freqs <- rep(1 / K, K)
    }
    log_freqs <- log(freqs)

    for (k in seq_len(K)) {
      mu[, lc[k]] <- log_freqs[k]
    }
    se[, lc] <- 0
  }

  list(mu = mu, se = se)
}
