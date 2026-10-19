#!/usr/bin/env Rscript

# Deterministic acceptance checks for the simulation generator and phylogenetic
# comparator. This file is written before its implementation as the TDD oracle.

source_if_present <- function(path) {
  if (file.exists(path)) sys.source(path, envir = .GlobalEnv)
}
source_if_present("script/gnn_incremental/dgp.R")
source_if_present("script/gnn_incremental/phylo_gls.R")

required <- c("make_gnn_case", "make_covariate_basis",
              "predict_phylo_gls_conditional")
missing <- required[!vapply(required, exists, logical(1), where = .GlobalEnv, mode = "function", inherits = FALSE)]
if (length(missing)) {
  stop("Expected implementation missing: ", paste(missing, collapse = ", "),
       call. = FALSE)
}

case1 <- make_gnn_case(
  n_species = 40L, obs_per = 4L, n_covs = 10L,
  variance_multiplier = 0.3, f_type = "nonlinear",
  beta_strength = 1, seed = 10301L
)
case2 <- make_gnn_case(
  n_species = 40L, obs_per = 4L, n_covs = 10L,
  variance_multiplier = 0.3, f_type = "nonlinear",
  beta_strength = 1, seed = 10301L
)
small_basis_data <- case1$df_complete[seq_len(12L), , drop = FALSE]
covariate_basis <- make_covariate_basis(
  small_basis_data, case1$predictor_names, nonlinear = TRUE
)
expected_basis_columns <- 1L + length(case1$predictor_names) * 2L +
  choose(length(case1$predictor_names), 2L)
stopifnot(all(dim(covariate_basis) ==
                c(nrow(small_basis_data), expected_basis_columns)))
stopifnot(all(c("(Intercept)", "z1", "z1^2", "z1:z2") %in%
                colnames(covariate_basis)))
stopifnot(isTRUE(all.equal(unname(covariate_basis[, "z1^2"]),
                           unname(small_basis_data$z1^2))))
stopifnot(isTRUE(all.equal(unname(covariate_basis[, "z1:z2"]),
                           unname(small_basis_data$z1 * small_basis_data$z2))))
stopifnot(identical(case1$tree$tip.label, case2$tree$tip.label))
stopifnot(isTRUE(all.equal(case1$df_complete, case2$df_complete)))
stopifnot(identical(case1$mask, case2$mask))
stopifnot(length(case1$mask) == nrow(case1$df_complete))
stopifnot(all(case1$df_complete$row_id == seq_len(nrow(case1$df_complete))))
stopifnot(all(is.finite(as.matrix(case1$df_complete[paste0("z", 1:10)]))))
stopifnot(all(is.na(case1$df_observed$y) == case1$mask))
stopifnot(all(case1$df_observed$y[!case1$mask] ==
              case1$df_complete$y[!case1$mask]))
stopifnot(!any(case1$whole_missing_mask & case1$within_missing_mask))
stopifnot(isTRUE(all.equal(
  unname(case1$phylo_effect),
  unname(case1$phylo_effect_raw / sqrt(diag(
    ape::vcv(case1$tree)[case1$tree$tip.label,
                         case1$tree$tip.label, drop = FALSE]
  ))), tolerance = 1e-12
)))
stopifnot(isTRUE(all.equal(
  unname(case1$phylo_correlation),
  unname(stats::cov2cor(ape::vcv(case1$tree)[case1$tree$tip.label,
                                            case1$tree$tip.label,
                                            drop = FALSE])),
  tolerance = 1e-12
)))
stopifnot(all(c("fully_unobserved_species", "partly_observed_species") %in%
                unique(case1$stratum[case1$mask])))
partial <- !case1$species %in% case1$fully_unobserved_species
observed_per_partial <- tapply(!case1$mask[partial], case1$species[partial], sum)
stopifnot(all(observed_per_partial >= 1L))
stopifnot(all(tapply(case1$mask[!partial], case1$species[!partial], all)))
for (seed_i in 10302:10341) {
  case_i <- make_gnn_case(
    n_species = 40L, obs_per = 4L, n_covs = 10L,
    variance_multiplier = 0.3, f_type = "nonlinear",
    beta_strength = 1, seed = seed_i
  )
  partial_i <- !case_i$species %in% case_i$fully_unobserved_species
  keep_i <- tapply(!case_i$mask[partial_i], case_i$species[partial_i], sum)
  stopifnot(all(keep_i >= 1L))
  stopifnot(!any(case_i$whole_missing_mask & case_i$within_missing_mask))
  stopifnot(all(case_i$stratum[case_i$within_missing_mask] ==
                  "partly_observed_species"))
}

# Independent dense covariance oracle with repeated observations from a species.
V <- matrix(c(1.4, 0.5, 0.2,
              0.5, 1.1, 0.3,
              0.2, 0.3, 0.9), nrow = 3L, byrow = TRUE)
sp_obs <- c(1L, 1L, 2L, 2L)
sp_new <- c(1L, 2L, 3L, 3L)
X_obs <- cbind(1, c(-1.0, 0.2, 0.7, 1.1))
X_new <- cbind(1, c(-0.4, 0.6, -0.2, 1.3))
y_obs <- c(0.1, 0.7, 1.4, 1.2)
tau2 <- 0.65
sigma2 <- 0.30

pred <- predict_phylo_gls_conditional(
  y_obs = y_obs, X_obs = X_obs, X_new = X_new, V = V,
  species_obs = sp_obs, species_new = sp_new,
  tau2 = tau2, sigma2 = sigma2
)
C_oo <- tau2 * V[sp_obs, sp_obs, drop = FALSE] +
  diag(sigma2, length(sp_obs))
C_no <- tau2 * V[sp_new, sp_obs, drop = FALSE]
C_nn <- tau2 * V[sp_new, sp_new, drop = FALSE] +
  diag(sigma2, length(sp_new))
C_oo_inv <- solve(C_oo)
beta_oracle <- solve(crossprod(X_obs, C_oo_inv %*% X_obs),
                     crossprod(X_obs, C_oo_inv %*% y_obs))
resid_oracle <- y_obs - drop(X_obs %*% beta_oracle)
pred_oracle <- drop(X_new %*% beta_oracle +
                      C_no %*% C_oo_inv %*% resid_oracle)
stopifnot(isTRUE(all.equal(as.numeric(pred), pred_oracle, tolerance = 1e-10)))

# The optimized observation-level fit must agree with the dense conditional
# prediction at its fitted covariance parameters.
small_fit <- fit_phylo_random_intercept(
  y = y_obs, X = X_obs, species_index = sp_obs, V = V, maxit = 50L
)
pred_efficient <- small_fit$predict(X_new, sp_new)
pred_dense <- predict_phylo_gls_conditional(
  y_obs = y_obs, X_obs = X_obs, X_new = X_new, V = V,
  species_obs = sp_obs, species_new = sp_new,
  tau2 = small_fit$tau2, sigma2 = small_fit$sigma2
)
stopifnot(isTRUE(all.equal(as.numeric(pred_efficient),
                            as.numeric(pred_dense), tolerance = 1e-8)))

# A tree-derived, named correlation fixture exercises species-to-tip alignment
# in addition to the arbitrary positive-definite matrix above.
tree_fixture <- ape::read.tree(text = "((sp1:1,sp2:1):1,sp3:2);")
V_tree <- stats::cov2cor(ape::vcv(tree_fixture))
stopifnot(identical(rownames(V_tree), tree_fixture$tip.label))
stopifnot(identical(colnames(V_tree), tree_fixture$tip.label))
pred_tree <- predict_phylo_gls_conditional(
  y_obs = y_obs, X_obs = X_obs, X_new = X_new, V = V_tree,
  species_obs = sp_obs, species_new = sp_new,
  tau2 = tau2, sigma2 = sigma2
)
C_tree_oo <- tau2 * V_tree[sp_obs, sp_obs, drop = FALSE] +
  diag(sigma2, length(sp_obs))
C_tree_no <- tau2 * V_tree[sp_new, sp_obs, drop = FALSE]
C_tree_inv <- solve(C_tree_oo)
beta_tree <- solve(crossprod(X_obs, C_tree_inv %*% X_obs),
                   crossprod(X_obs, C_tree_inv %*% y_obs))
pred_tree_oracle <- drop(X_new %*% beta_tree + C_tree_no %*%
  C_tree_inv %*% (y_obs - drop(X_obs %*% beta_tree)))
stopifnot(isTRUE(all.equal(as.numeric(pred_tree), unname(pred_tree_oracle),
                            tolerance = 1e-10)))

# Negative control: dropping the phylogenetic cross-covariance must be detected.
pred_broken <- drop(X_new %*% beta_oracle)
stopifnot(max(abs(as.numeric(pred_broken) - pred_oracle)) > 1e-3)
stopifnot(inherits(try(predict_phylo_gls_conditional(
  y_obs = y_obs, X_obs = X_obs, X_new = X_new[-1L, , drop = FALSE],
  V = V, species_obs = sp_obs, species_new = sp_new,
  tau2 = tau2, sigma2 = sigma2
), silent = TRUE), "try-error"))
cat("HARNESS_VALIDATION_PASS\n")
