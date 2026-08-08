# D-Blk1: covariates must follow species / observation identity, not row position.
# Traits are reordered (and padded) to tree tip order; matching nrow alone used
# to silently pair the wrong environment with each species.

recover_cov <- function(pd, col = "temp") {
  as.numeric(pd$covariates[, col] * pd$cov_sds[[col]] + pd$cov_means[[col]])
}

test_that("shuffled trait rows pair named covariates by species, not position", {
  n <- 8L
  tree <- ape::rtree(n)
  sp <- tree$tip.label
  temp_by_sp <- setNames(as.numeric(seq_len(n) - 1L), sp)

  # Same shuffled input order for traits and covs (documented nrow contract).
  # After tip-align, pairing must still be by species — not leftover position.
  traits <- data.frame(row.names = rev(sp), y = rnorm(n))
  covs <- data.frame(row.names = rev(sp), temp = unname(temp_by_sp[rev(sp)]))

  pd <- preprocess_traits(traits, tree, covariates = covs, log_transform = FALSE)

  expect_equal(pd$species_names, sp)
  expect_equal(recover_cov(pd), unname(temp_by_sp[sp]))
  # Unpermuted covs would remain in the reversed input order.
  expect_false(isTRUE(all.equal(recover_cov(pd), unname(temp_by_sp[rev(sp)]))))
})

test_that("shuffled covariate rows pair to tree-ordered traits by species", {
  n <- 8L
  tree <- ape::rtree(n)
  sp <- tree$tip.label
  temp_by_sp <- setNames(as.numeric(seq_len(n) - 1L), sp)

  traits <- data.frame(row.names = sp, y = rnorm(n))
  covs <- data.frame(row.names = rev(sp), temp = unname(temp_by_sp[rev(sp)]))

  pd <- preprocess_traits(traits, tree, covariates = covs, log_transform = FALSE)

  expect_equal(pd$species_names, sp)
  expect_equal(recover_cov(pd), unname(temp_by_sp[sp]))
})

test_that("nrow-matching covariates with wrong species names error", {
  n <- 6L
  tree <- ape::rtree(n)
  sp <- tree$tip.label
  traits <- data.frame(row.names = sp, y = rnorm(n))
  covs <- data.frame(row.names = paste0("not_", sp), temp = seq_len(n))

  expect_error(
    preprocess_traits(traits, tree, covariates = covs, log_transform = FALSE),
    regexp = "match|species|rownames"
  )
})

test_that("multi-obs shuffled traits pair covs by observation species, not position", {
  n <- 5L
  k <- 2L
  tree <- ape::rtree(n)
  sp <- tree$tip.label
  # Distinct temp per observation; within-species order is stable across inputs.
  df <- data.frame(
    species = rep(sp, each = k),
    y = rnorm(n * k),
    temp = rep(as.numeric(seq_len(n) - 1L), each = k) + c(0, 0.25)
  )
  # Reverse species *blocks* on both tables (within-species order unchanged).
  block_ord <- unlist(lapply(rev(seq_len(n)), function(i) {
    ((i - 1L) * k + 1L):(i * k)
  }), use.names = FALSE)
  traits <- df[block_ord, c("species", "y")]
  covs <- df[block_ord, c("species", "temp")]

  pd <- preprocess_traits(
    traits, tree, species_col = "species",
    covariates = covs, log_transform = FALSE
  )

  expect_equal(pd$obs_species, df$species)
  expect_equal(recover_cov(pd), df$temp)
  expect_false(isTRUE(all.equal(recover_cov(pd), df$temp[block_ord])))
})

test_that("multi-obs nrow match but unmatched species names error", {
  n <- 4L
  k <- 2L
  tree <- ape::rtree(n)
  sp <- tree$tip.label
  traits <- data.frame(
    species = rep(sp, each = k),
    y = rnorm(n * k)
  )
  covs <- data.frame(
    species = paste0("not_", rep(sp, each = k)),
    temp = seq_len(n * k)
  )

  expect_error(
    preprocess_traits(
      traits, tree, species_col = "species",
      covariates = covs, log_transform = FALSE
    ),
    regexp = "match|species"
  )
})

test_that("predict cov tensor follows species identity after trait shuffle", {
  skip_if_not_installed("torch")

  n <- 4L
  tree <- ape::rtree(n)
  sp <- tree$tip.label
  temp_by_sp <- setNames(as.numeric(seq_len(n) - 1L), sp)

  traits <- data.frame(row.names = rev(sp), y = rnorm(n))
  covs <- data.frame(row.names = rev(sp), temp = unname(temp_by_sp[rev(sp)]))
  pd <- preprocess_traits(traits, tree, covariates = covs, log_transform = FALSE)

  temp_aligned <- recover_cov(pd)
  expect_equal(temp_aligned, unname(temp_by_sp[sp]))

  model <- ResidualPhyloDAE(
    input_dim = 1L, hidden_dim = 4L, coord_dim = 1L, cov_dim = 3L,
    per_column_rs = TRUE, n_gnn_layers = 1L, gate_cap = 0.8,
    use_attention = FALSE, n_user_cov = 1L, dropout = 0,
    use_transformer_blocks = FALSE, n_heads = 1L, ffn_mult = 1L
  )
  torch::with_no_grad({
    model$cov_linear$weight$copy_(
      torch::torch_tensor(matrix(2, nrow = 1L, ncol = 1L))
    )
    model$cov_linear$bias$copy_(torch::torch_tensor(1))
  })

  fit <- structure(
    list(
      model_state = lapply(model$state_dict(), function(t) t$detach()$cpu()),
      model_config = list(
        input_dim = 1L, hidden_dim = 4L, k_eigen = 1L, cov_dim = 3L,
        per_column_rs = TRUE, n_gnn_layers = 1L, gate_cap = 0.8,
        use_attention = FALSE, n_user_cov = 1L, dropout = 0,
        refine_steps = 1L, use_transformer_blocks = FALSE,
        n_heads = 1L, ffn_mult = 1L
      ),
      graph = list(coords = matrix(0, n, 1), adj = diag(n), D_sq = NULL),
      baseline = list(
        mu = matrix(0, n, 1, dimnames = list(sp, "y")),
        se = matrix(0, n, 1, dimnames = list(sp, "y"))
      ),
      species_names = pd$species_names,
      obs_species = sp,
      obs_to_species = NULL,
      X_scaled = matrix(NA_real_, n, 1, dimnames = list(sp, "y")),
      n_species = n,
      n_obs = n,
      multi_obs = FALSE,
      trait_names = "y",
      latent_names = "y",
      trait_map = list(list(
        name = "y", type = "continuous", latent_cols = 1L,
        mean = 0, sd = 1, log_transform = FALSE
      )),
      calibrated_gates = c(y = 0),
      r_cal_bm = c(y = 0),
      r_cal_gnn = c(y = 0),
      r_cal_mean = c(y = 0),
      mean_baseline_per_col = c(y = 0),
      conformal_scores = NULL,
      covariates = matrix(temp_aligned, ncol = 1, dimnames = list(sp, "temp")),
      cov_means = 0,
      cov_sds = 1,
      cov_names = "temp"
    ),
    class = "pigauto_fit"
  )

  pred <- predict(fit, return_se = FALSE)
  expected <- 1 + 2 * unname(temp_by_sp[sp])
  expect_equal(as.numeric(pred$imputed_latent[, "y"]), expected, tolerance = 1e-6)
  # Positional pairing with reversed traits would flip the cov_linear surface.
  expect_false(isTRUE(all.equal(
    as.numeric(pred$imputed_latent[, "y"]),
    1 + 2 * unname(temp_by_sp[rev(sp)]),
    tolerance = 1e-6
  )))
})
