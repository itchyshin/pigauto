# B2 (2026-08-16): conformal_method = "mondrian" -- stratified conformal
# conditioned on phylogenetic sampling locality. Fixes the exchangeability
# failure diagnosed in
# docs/dev-log/2026-08-16-mechanism-coverage-results.md (split conformal
# calibrates on the observed complement, which sits closer -- in cophenetic
# distance -- to other observed cells than genuinely-missing cells do,
# undercovering worst under clade-structured / MAR_phylo missingness).

skip_if_no_libtorch()

# ---------------------------------------------------------------------------
# unit: locality statistic on a hand-computable 6-tip tree
# ---------------------------------------------------------------------------

test_that("mondrian_locality matches hand-computed cophenetic distances", {
  # (((A:1,B:1):1,(C:1,D:1):2):1,(E:1,F:1):5);  cophenetic distances:
  #    A B  C  D  E  F
  #  A 0 2  5  5  9  9
  #  B 2 0  5  5  9  9
  #  C 5 5  0  2 10 10
  #  D 5 5  2  0 10 10
  #  E 9 9 10 10  0  2
  #  F 9 9 10 10  2  0
  nwk <- "(((A:1,B:1):1,(C:1,D:1):2):1,(E:1,F:1):5);"
  tr  <- ape::read.tree(text = nwk)
  D   <- ape::cophenetic.phylo(tr)[tr$tip.label, tr$tip.label]
  D_sq <- D^2
  idx <- stats::setNames(seq_along(tr$tip.label), tr$tip.label)

  # k = 5, all of B..F observed: distances from A are 2,5,5,9,9 -> mean 6
  loc <- pigauto:::mondrian_locality(
    D_sq, obs_idx = idx[c("B", "C", "D", "E", "F")],
    target_idx = idx["A"], k = 5L
  )
  expect_equal(unname(loc), 6, tolerance = 1e-10)

  # k = 2: nearest two distances are 2, 5 -> mean 3.5
  loc2 <- pigauto:::mondrian_locality(
    D_sq, obs_idx = idx[c("B", "C", "D", "E", "F")],
    target_idx = idx["A"], k = 2L
  )
  expect_equal(unname(loc2), 3.5, tolerance = 1e-10)

  # fewer than k observed: only C, D -> use both, mean(5, 5) = 5
  loc3 <- pigauto:::mondrian_locality(
    D_sq, obs_idx = idx[c("C", "D")], target_idx = idx["A"], k = 5L
  )
  expect_equal(unname(loc3), 5, tolerance = 1e-10)

  # self-exclusion: obs_idx includes the target itself (A)
  loc4 <- pigauto:::mondrian_locality(
    D_sq, obs_idx = idx[c("A", "B", "C")], target_idx = idx["A"], k = 5L
  )
  expect_equal(unname(loc4), mean(c(2, 5)), tolerance = 1e-10)

  # no eligible observed species (only self) -> NA
  loc5 <- pigauto:::mondrian_locality(
    D_sq, obs_idx = idx["A"], target_idx = idx["A"], k = 5L
  )
  expect_true(is.na(loc5))

  # vectorised over multiple targets at once
  loc_vec <- pigauto:::mondrian_locality(
    D_sq, obs_idx = idx[c("B", "C", "D", "E", "F")],
    target_idx = idx[c("A", "B")], k = 5L
  )
  expect_equal(unname(loc_vec[1]), 6, tolerance = 1e-10)
  # obs_idx is still {B,C,D,E,F}; B excludes itself, leaving only 4
  # eligible neighbours: C,D,E,F at 5,5,9,9 -> mean 7.
  expect_equal(unname(loc_vec[2]), 7, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# stratification: per-stratum scores differ when residuals differ by stratum
# ---------------------------------------------------------------------------

test_that("compute_conformal_scores mondrian produces different near/far scores", {
  n_obs_species <- 20L
  n_val_species <- 40L
  n <- n_obs_species + n_val_species

  # Species laid out on a line 1..n; D_sq = squared linear distance. The
  # first 20 are the observed (calibration donor) pool; val species
  # 21:40 sit just past it (near), 41:60 sit far past it (far).
  D_sq <- outer(seq_len(n), seq_len(n), function(a, b) (a - b)^2)
  trait_map <- list(list(
    name = "t1", type = "continuous", latent_cols = 1L, mean = 0, sd = 1
  ))

  near_idx <- (n_obs_species + 1L):(n_obs_species + 20L)   # species 21:40
  far_idx  <- (n_obs_species + 21L):n                       # species 41:60

  mu_cal <- matrix(0, n, 1)
  mu_cal[near_idx, 1] <- -0.1   # |truth(0) - pred| = 0.1
  mu_cal[far_idx, 1]  <- -5     # |truth(0) - pred| = 5
  delta_cal <- mu_cal           # gate = 0, so delta_cal is unused
  X_truth_r <- matrix(0, n, 1)

  val_mask_mat <- matrix(FALSE, n, 1)
  val_mask_mat[c(near_idx, far_idx), 1] <- TRUE
  obs_mask_mat <- matrix(FALSE, n, 1)
  obs_mask_mat[seq_len(n_obs_species), 1] <- TRUE

  scores <- pigauto:::compute_conformal_scores(
    trait_map = trait_map,
    calibrated_gates = 0,
    mu_cal = mu_cal, delta_cal = delta_cal,
    X_truth_r = X_truth_r, val_mask_mat = val_mask_mat,
    method = "mondrian", D_sq = D_sq, obs_mask_mat = obs_mask_mat
  )

  mo <- attr(scores, "mondrian")[["t1"]]
  expect_false(mo$fallback)
  expect_equal(mo$near_score, 0.1, tolerance = 1e-10)
  expect_equal(mo$far_score, 5, tolerance = 1e-10)
  expect_true(mo$far_score > mo$near_score)
  # Base named vector carries a finite global score for callers that
  # ignore the "mondrian" attribute.
  expect_true(is.finite(unname(scores["t1"])))
})

# ---------------------------------------------------------------------------
# fallback: < 10 residuals per stratum -> global score, flag set
# ---------------------------------------------------------------------------

test_that("compute_conformal_scores mondrian falls back below the 19-per-stratum floor", {
  n_obs_species <- 20L
  n_val_species <- 12L  # 6 near + 6 far -- below the 19-cell floor
  # (19 = the smallest stratum size whose achievable coverage ceiling
  # n_s/(n_s+1) reaches 0.95; measured failure at 13-cell strata in the
  # 2026-08-16 mech_cov_mondrian verification)
  n <- n_obs_species + n_val_species

  D_sq <- outer(seq_len(n), seq_len(n), function(a, b) (a - b)^2)
  trait_map <- list(list(
    name = "t1", type = "continuous", latent_cols = 1L, mean = 0, sd = 1
  ))

  val_idx <- (n_obs_species + 1L):n
  set.seed(1)
  mu_cal <- matrix(0, n, 1)
  mu_cal[val_idx, 1] <- -stats::runif(length(val_idx), 0.1, 5)
  delta_cal <- mu_cal
  X_truth_r <- matrix(0, n, 1)

  val_mask_mat <- matrix(FALSE, n, 1); val_mask_mat[val_idx, 1] <- TRUE
  obs_mask_mat <- matrix(FALSE, n, 1)
  obs_mask_mat[seq_len(n_obs_species), 1] <- TRUE

  scores_mond <- pigauto:::compute_conformal_scores(
    trait_map = trait_map, calibrated_gates = 0,
    mu_cal = mu_cal, delta_cal = delta_cal,
    X_truth_r = X_truth_r, val_mask_mat = val_mask_mat,
    method = "mondrian", D_sq = D_sq, obs_mask_mat = obs_mask_mat
  )
  scores_split <- pigauto:::compute_conformal_scores(
    trait_map = trait_map, calibrated_gates = 0,
    mu_cal = mu_cal, delta_cal = delta_cal,
    X_truth_r = X_truth_r, val_mask_mat = val_mask_mat,
    method = "split"
  )

  mo <- attr(scores_mond, "mondrian")[["t1"]]
  expect_true(mo$fallback)
  expect_equal(mo$near_score, mo$far_score, tolerance = 1e-12)
  expect_true(is.na(mo$threshold))
  expect_equal(unname(mo$near_score), unname(scores_split["t1"]), tolerance = 1e-10)
  expect_equal(unname(scores_mond["t1"]), unname(scores_split["t1"]), tolerance = 1e-10)
})

test_that("compute_conformal_scores errors clearly when mondrian inputs are missing", {
  trait_map <- list(list(
    name = "t1", type = "continuous", latent_cols = 1L, mean = 0, sd = 1
  ))
  m <- matrix(0, 4, 1)
  expect_error(
    pigauto:::compute_conformal_scores(
      trait_map = trait_map, calibrated_gates = 0,
      mu_cal = m, delta_cal = m, X_truth_r = m,
      val_mask_mat = matrix(TRUE, 4, 1), method = "mondrian"
    ),
    "D_sq.*obs_mask_mat"
  )
})

# ---------------------------------------------------------------------------
# end-to-end smoke: wider intervals for an isolated long-branch clade
# ---------------------------------------------------------------------------

test_that("impute() with conformal_method = 'mondrian' widens intervals in an isolated clade", {
  set.seed(11)
  # 300 main tips so the per-trait val set (~45 cells at missing_frac 0.6,
  # val_frac 0.25) yields >= 19 residuals per stratum -- below that the
  # 19-per-stratum floor (see compute_conformal_scores) falls back to the
  # global score and near/far widths would be identical by design.
  main <- ape::rcoal(300)
  main$tip.label <- paste0("M", seq_len(300))
  iso <- ape::rcoal(6)
  iso$tip.label <- paste0("I", seq_len(6))
  tree <- ape::bind.tree(main, iso, where = 1L, position = 0)
  iso_mrca <- ape::getMRCA(tree, paste0("I", seq_len(6)))
  stem_edge <- which(tree$edge[, 2] == iso_mrca)
  # Very long stem branch: the isolated clade is genuinely far (in
  # cophenetic distance) from every observed species.
  tree$edge.length[stem_edge] <- tree$edge.length[stem_edge] + 80
  tree <- ape::reorder.phylo(tree)

  set.seed(11)
  y <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  df <- data.frame(row.names = names(y), trait1 = as.numeric(y))
  iso_sp <- paste0("I", seq_len(6))
  df[iso_sp, "trait1"] <- NA  # genuinely missing, never observed

  fit_split <- suppressWarnings(impute(
    df, tree, conformal_method = "split", epochs = 30L,
    missing_frac = 0.6, verbose = FALSE, seed = 5
  ))
  fit_mond <- suppressWarnings(impute(
    df, tree, conformal_method = "mondrian", epochs = 30L,
    missing_frac = 0.6, verbose = FALSE, seed = 5
  ))

  lo_s <- fit_split$prediction$conformal_lower[, "trait1"]
  hi_s <- fit_split$prediction$conformal_upper[, "trait1"]
  lo_m <- fit_mond$prediction$conformal_lower[, "trait1"]
  hi_m <- fit_mond$prediction$conformal_upper[, "trait1"]

  expect_true(all(is.finite(lo_s)) && all(is.finite(hi_s)))
  expect_true(all(is.finite(lo_m)) && all(is.finite(hi_m)))

  w_split <- hi_s - lo_s
  w_mond  <- hi_m - lo_m

  expect_true(all(w_mond[iso_sp] > w_split[iso_sp]))
})

# ---------------------------------------------------------------------------
# back-compat: default "split" results unchanged
# ---------------------------------------------------------------------------
# NOTE: an earlier version compared fit$conformal_scores to numbers recorded
# on one machine; torch RNG differs across platforms, so that failed on CI.
# The invariant is tested deterministically instead: the "split" branch of
# compute_conformal_scores() must still produce exactly the pre-mondrian
# formula, quantile(|truth - blend|, ceiling((1-alpha)(n+1))/n).

test_that("default conformal_method = 'split' reproduces the split-quantile formula", {
  set.seed(7)
  n <- 30L
  tm <- list(list(name = "x", type = "continuous", latent_cols = 1L,
                  mean = 0, sd = 1))
  mu    <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "x"))
  delta <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "x"))
  truth <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "x"))
  val   <- matrix(FALSE, n, 1L)
  val[1:25, 1L] <- TRUE
  w_bm <- 0.6; w_gnn <- 0.3; w_mean <- 0.1; mean_col <- 0.2

  scores <- pigauto:::compute_conformal_scores(
    trait_map = tm, calibrated_gates = c(x = w_gnn),
    mu_cal = mu, delta_cal = delta, X_truth_r = truth, val_mask_mat = val,
    method = "split",
    r_cal_bm = c(x = w_bm), r_cal_gnn = c(x = w_gnn), r_cal_mean = c(x = w_mean),
    mean_baseline_per_col = c(x = mean_col)
  )

  blend <- w_bm * mu[1:25, 1] + w_gnn * delta[1:25, 1] + w_mean * mean_col
  res <- abs(truth[1:25, 1] - blend)
  q_level <- min(ceiling(0.95 * 26) / 25, 1)
  expected <- as.numeric(stats::quantile(res, q_level))

  expect_equal(unname(scores["x"]), expected, tolerance = 1e-12)
  expect_null(attr(scores, "mondrian"))
})

test_that("a 'split' fit carries no mondrian state and records the method", {
  set.seed(42)
  tree <- ape::rtree(40)
  sp <- tree$tip.label
  set.seed(42)
  df <- data.frame(
    row.names = sp,
    tr1 = abs(stats::rnorm(40)) + 0.5,
    tr2 = abs(stats::rnorm(40)) + 0.5
  )
  pd  <- preprocess_traits(df, tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 1, trait_map = pd$trait_map)
  fit <- suppressWarnings(fit_pigauto(
    pd, tree, splits = spl, epochs = 20L, eval_every = 10L, patience = 5L,
    verbose = FALSE, seed = 1, conformal_method = "split"
  ))
  expect_true(all(is.finite(fit$conformal_scores)))
  expect_null(fit$conformal_mondrian)
  expect_identical(fit$conformal_method, "split")
})
