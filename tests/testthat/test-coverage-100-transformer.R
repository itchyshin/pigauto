test_that("tiny model gates blend real CPU predictions and backpropagate", {
  skip_if_no_libtorch()
  torch::torch_manual_seed(42)
  x <- torch::torch_tensor(matrix(seq_len(8) / 10, 4, 2))
  coords <- torch::torch_zeros(c(4L, 2L))
  covs <- torch::torch_zeros(c(4L, 4L))
  adj <- torch::torch_eye(4L)
  baseline <- torch::torch_ones(c(4L, 2L))
  for (transformer in c(FALSE, TRUE)) {
    for (attention in if (transformer) FALSE else c(FALSE, TRUE)) {
      model <- ResidualPhyloDAE(p_latent = 2L, hidden = 8L, k_eigen = 2L,
        cov_dim = 4L, n_gnn_layers = 1L, dropout = 0,
        use_transformer_blocks = transformer, use_attention = attention,
        n_heads = 2L, per_column_rs = FALSE, gate_cap = 0.6)
      model$eval()
      out <- model(x, coords, covs, adj)
      blend <- model(x, coords, covs, adj, baseline_mu = baseline)
      expect_equal(as.array(blend), as.array((1 - out$rs) * baseline + out$rs * out$delta), tolerance = 1e-6)
      expect_true(as.numeric(out$rs$item()) > 0 && as.numeric(out$rs$item()) < 0.6)
      blend$sum()$backward()
      expect_true(all(is.finite(as.array(model$res_raw$grad))))
      expect_equal(dim(as.array(blend)), c(4L, 2L))
    }
  }
})

test_that("covariates refine observations sharing a species", {
  skip_if_no_libtorch()
  torch::torch_manual_seed(7)
  model <- ResidualPhyloDAE(p_latent = 2L, hidden = 8L, k_eigen = 2L,
    cov_dim = 5L, n_user_cov = 1L, n_gnn_layers = 1L, dropout = 0,
    use_transformer_blocks = TRUE, n_heads = 2L)
  model$eval()
  torch::with_no_grad({
    for (layer in list(model$cov_encoder[[1]], model$cov_encoder[[4]],
                       model$obs_refine[[1]], model$obs_refine[[3]],
                       model$dec1, model$dec2)) {
      layer$weight$fill_(0.1)
      layer$bias$fill_(0.1)
    }
    model$cov_linear$weight$fill_(0.25)
    model$cov_linear$bias$zero_()
  })
  x <- torch::torch_zeros(c(4L, 2L))
  coords <- torch::torch_zeros(c(2L, 2L))
  covs <- torch::torch_tensor(cbind(matrix(0, 4, 4), c(0, 1, 0, 1)))
  species <- torch::torch_tensor(c(1L, 1L, 2L, 2L), dtype = torch::torch_long())
  out <- model(x, coords, covs, torch::torch_eye(2L), obs_to_species = species)
  expect_equal(dim(as.array(out$delta)), c(4L, 2L))
  expect_false(isTRUE(all.equal(as.array(out$delta[1, ]), as.array(out$delta[2, ]))))
  expect_equal(as.array(out$delta[1, ]), as.array(out$delta[3, ]), tolerance = 1e-6)
  expect_equal(as.array(out$fixed_effects[1, ]), as.array(out$fixed_effects[3, ]), tolerance = 1e-6)
})

test_that("trait attention and per-column gates return finite differentiable output", {
  skip_if_no_libtorch()
  torch::torch_manual_seed(10)
  model <- ResidualPhyloDAE(input_dim = 2L, hidden_dim = 8L, coord_dim = 2L,
    cov_dim = 4L, n_gnn_layers = 1L, dropout = 0,
    use_transformer_blocks = TRUE, n_heads = 2L, use_trait_attention = TRUE,
    n_trait_heads = 2L, trait_embed_dim = 8L)
  x <- torch::torch_randn(c(3L, 2L), requires_grad = TRUE)
  out <- model(x, torch::torch_zeros(c(3L, 2L)),
               torch::torch_zeros(c(3L, 4L)), torch::torch_eye(3L),
               D_sq = torch::torch_ones(c(3L, 3L)))
  expect_equal(dim(as.array(out$rs)), c(1L, 2L))
  expect_true(all(as.array(out$rs) > 0 & as.array(out$rs) < 0.8))
  out$delta$sum()$backward()
  expect_true(all(is.finite(as.array(x$grad))))
  expect_true(all(is.finite(as.array(out$delta))))
})

test_that("refinement validates arguments and calibrated blend", {
  skip_if_no_libtorch()
  model <- ResidualPhyloDAE(input_dim = 1L, hidden_dim = 8L, coord_dim = 2L,
    cov_dim = 2L, n_gnn_layers = 1L, dropout = 0, n_heads = 2L)
  model$eval()
  x <- torch::torch_zeros(c(3L, 1L))
  args <- list(model = model, X_init = x, t_coords = torch::torch_zeros(c(3L, 2L)),
    t_adj = torch::torch_eye(3L), t_obs_to_sp = NULL, t_D_sq = NULL,
    t_MU = torch::torch_ones(c(3L, 1L)), t_covariates = NULL, pin_mask = NULL,
    refine_steps = 1L, device = torch::torch_device("cpu"))
  expect_error(do.call(pigauto_refine_forward, modifyList(args, list(refine_steps = 0L))),
               "refine_steps must be >= 1")
  expect_error(do.call(pigauto_refine_forward, c(args, list(blend = "calibrated"))),
               "calibrated blend requires")
  w <- torch::torch_tensor(0)
  result <- do.call(pigauto_refine_forward, c(args, list(blend = "calibrated",
    t_w_bm = w + 1, t_w_gnn = w, t_w_mean = w, t_mean_baseline = w)))
  expect_equal(as.array(result$pred), matrix(1, 3, 1), tolerance = 1e-6)
  expect_equal(dim(result$delta_mat), c(3L, 1L))
})
