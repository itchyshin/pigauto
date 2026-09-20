# script/campaign_gnn_off_lib.R
# Shared pieces of the campaign runners: DGPs, the seeded user-level mask, scoring, and the
# arm-dispatch loop (run_arms). Sourced by script/campaign_gnn_off_cell.R, script/campaign_sim_cell.R
# and script/campaign_solver_cell.R. score_arm() reads `truth` and `mask` from the calling
# environment (set by the runner) -- unchanged calling convention, kept for campaign_solver_cell.R
# which this lane does not touch.
#
# BACKWARD COMPATIBILITY: make_dgp()/make_cell() reproduce the OLD RNG stream and OLD trait
# construction byte-for-byte whenever called with the OLD defaults (lambda = 1, rho = 0,
# thresholds = "sample", driver = FALSE, evo = the dgp's original evolutionary model). The new
# corrected design (2026-09-20, "the corrected design") only fires when any of those arguments is
# given a non-default value. score_arm() keeps returning a plain data.frame (campaign_solver_cell.R
# does `results[[arm]] <- score_arm(...)` and rbinds); the new calibration long-format frame is
# attached as `attr(., "calib")` so old callers are unaffected.

# ---- latent simulation for the corrected design (Section A) ---------------------------------
# L = t(chol(V_lambda)) %*% Z %*% chol(Sigma_rho), Z an n x K standard-normal matrix.
# V_lambda = lambda * V + (1 - lambda) * I, V = cov2cor(vcv(tree)) for evo = "BM".
# For evo = "OU" (alpha = 2 fixed): V is replaced by the stationary-OU correlation
# corr(i, j) = exp(-alpha * d_ij), d_ij = cophenetic (patristic) distance between tips i, j.
# This is the standard stationary approximation for an OU process on an ultrametric tree (t_i = t_j
# = tree height for every tip), NOT derived from pigauto::simulate_non_bm (that function simulates
# trait VALUES under OU, not a reusable correlation matrix) -- UNVERIFIED against a first-principles
# OU tip-covariance derivation; flagged in the report.
sim_latents <- function(tree, K, lambda = 1, rho = 0, evo = c("BM", "OU"), ou_alpha = 2,
                        driver_col = NULL, rho_driver = 0.35) {
  evo <- match.arg(evo)
  sp <- tree$tip.label; n <- length(sp)
  V <- cov2cor(ape::vcv(tree))[sp, sp]
  if (evo == "OU") {
    D <- ape::cophenetic.phylo(tree)[sp, sp]
    V <- exp(-ou_alpha * D); diag(V) <- 1
  }
  V_lambda <- lambda * V + (1 - lambda) * diag(n)
  dimnames(V_lambda) <- list(sp, sp)
  Sigma_rho <- matrix(rho, K, K); diag(Sigma_rho) <- 1
  # The MAR driver keeps its own correlation to the scored traits, independent of rho. Without this,
  # a rho = 0 cell makes the driver independent of everything it is supposed to predict, and "MAR"
  # there is MCAR on an unrelated variable.
  if (!is.null(driver_col)) {
    Sigma_rho[driver_col, -driver_col] <- rho_driver
    Sigma_rho[-driver_col, driver_col] <- rho_driver
    Sigma_rho[driver_col, driver_col] <- 1
    ev <- min(eigen(Sigma_rho, symmetric = TRUE, only.values = TRUE)$values)
    if (ev <= 1e-8) stop(sprintf("driver correlation %.2f makes Sigma non-PD at rho = %.2f", rho_driver, rho))
  }
  Z <- matrix(stats::rnorm(n * K), n, K)
  L <- t(chol(V_lambda)) %*% Z %*% chol(Sigma_rho)
  rownames(L) <- sp
  L
}

# thresholds = "fixed" uses population quantiles of a standard normal liability (class balance
# varies with lambda/rho because L is not marginally N(0,1) once rho < 1, only unit-diagonal by
# construction of V_lambda and Sigma_rho -- L's marginal variance is 1 per column since both
# V_lambda and Sigma_rho have unit diagonal, so qnorm() thresholds are still exactly calibrated).
# thresholds = "sample" reproduces the OLD behaviour: thresholds from the realised sample.
threshold_binary <- function(l, thresholds) if (thresholds == "fixed") l > 0 else l > stats::median(l)
threshold_ordinal <- function(l, thresholds, labels) {
  if (thresholds == "fixed") {
    breaks <- c(-Inf, stats::qnorm(c(.25, .5, .75)), Inf)
  } else {
    breaks <- stats::quantile(l, c(0, .25, .5, .75, 1))
  }
  factor(cut(l, breaks, labels = labels, include.lowest = TRUE), levels = labels, ordered = TRUE)
}
threshold_categorical <- function(l, thresholds, labels) {
  if (thresholds == "fixed") {
    breaks <- c(-Inf, stats::qnorm(c(1 / 3, 2 / 3)), Inf)
  } else {
    breaks <- stats::quantile(l, c(0, 1 / 3, 2 / 3, 1))
  }
  factor(cut(l, breaks, labels = labels, include.lowest = TRUE), levels = labels)
}

# ---- data ---------------------------------------------------------------------------------
make_dgp <- function(dgp, n, seed, lambda = 1, rho = 0, evo = NULL, thresholds = "sample", driver = FALSE) {
  set.seed(seed)
  if (dgp == "avonet") {
    e <- new.env(); utils::data("avonet300", package = "pigauto", envir = e)
    utils::data("tree300", package = "pigauto", envir = e)
    df <- e$avonet300; rownames(df) <- df$Species_Key; df$Species_Key <- NULL
    return(list(df = df, tree = e$tree300))
  }
  if (dgp == "bace_dgp") {
    # BACE's own simulator (as in script/sim_bace_dgp.R): gaussian response y with
    # phylogenetic signal 0.4 on three predictors (gaussian, gaussian, binary; signal 0.3),
    # birth-death tree. All four columns are treated as traits to impute; the campaign's
    # uniform 30% MCAR mask replaces BACE's response-only missingness.
    sim <- BACE::sim_bace(response_type = "gaussian",
                          predictor_types = c("gaussian", "gaussian", "binary"),
                          beta_resp = 0.5 * c(0.7, 0.4, 0.5), phylo_signal = c(0.4, 0.3, 0.3, 0.3),
                          n_cases = as.integer(n), n_species = as.integer(n),
                          missingness = c(0, 0, 0, 0), birth = 0.8, death = 0.4)
    df <- sim$complete_data; tree <- sim$tree
    rownames(df) <- as.character(df$Species); df$Species <- NULL
    df$x3 <- factor(as.character(df$x3), ordered = FALSE)   # binary, not ordered
    df <- df[tree$tip.label, , drop = FALSE]
    return(list(df = df, tree = tree))
  }
  if (!(dgp %in% c("types_mixed", "bm_mixed", "ou_mixed"))) stop("unknown dgp ", dgp)

  default_evo <- if (dgp == "ou_mixed") "OU" else "BM"
  evo_use <- if (is.null(evo)) default_evo else evo
  is_default <- identical(lambda, 1) && identical(rho, 0) && identical(thresholds, "sample") &&
    identical(driver, FALSE) && identical(evo_use, default_evo)

  if (is_default) {
    # ---- OLD code paths, byte-identical RNG stream (backward compatibility gate) ----
    if (dgp == "types_mixed") {
      tree <- ape::rcoal(n); sp <- tree$tip.label
      lat <- function() ape::rTraitCont(tree, model = "BM", sigma = 1)
      l_cnt <- lat(); l_prp <- lat(); l_bin <- lat(); l_ord <- lat(); l_cat <- lat()
      ord_breaks <- stats::quantile(l_ord, c(0, .25, .5, .75, 1))
      df <- data.frame(row.names = sp,
        c1 = lat(), c2 = lat(),
        cnt = as.integer(stats::rpois(n, exp(1.5 + 0.8 * l_cnt))),
        prp = stats::plogis(l_prp + stats::rnorm(n, 0, 0.3)),
        bin = factor(ifelse(l_bin > stats::median(l_bin), "yes", "no")),
        ord = factor(cut(l_ord, ord_breaks, labels = c("L1", "L2", "L3", "L4"), include.lowest = TRUE),
                     levels = c("L1", "L2", "L3", "L4"), ordered = TRUE),
        cat3 = factor(cut(l_cat, stats::quantile(l_cat, c(0, 1/3, 2/3, 1)), labels = c("A", "B", "C"), include.lowest = TRUE)))
      df$prp <- pmin(pmax(df$prp, 1e-4), 1 - 1e-4)
      return(list(df = df, tree = tree, trait_types = c(prp = "proportion")))
    }
    tree <- ape::rcoal(n)   # ultrametric: BM-appropriate, and BACE/MCMCglmm requires it
    tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
    sp <- tree$tip.label
    cont <- if (dgp == "ou_mixed") {
      pigauto::simulate_non_bm(tree, n_traits = 4L, scenario = "OU", seed = seed)
    } else {
      data.frame(row.names = sp,
                 t1 = ape::rTraitCont(tree, model = "BM", sigma = 1, root.value = 0),
                 t2 = ape::rTraitCont(tree, model = "BM", sigma = 1, root.value = 1),
                 t3 = ape::rTraitCont(tree, model = "BM", sigma = 1, root.value = 2),
                 t4 = ape::rTraitCont(tree, model = "BM", sigma = 1, root.value = 3))
    }
    cont <- as.data.frame(cont)[sp, , drop = FALSE]
    names(cont) <- paste0("c", seq_len(ncol(cont)))
    lat_b <- ape::rTraitCont(tree, model = "BM", sigma = 1)
    bin <- factor(ifelse(lat_b > stats::median(lat_b), "yes", "no"))
    lat_k <- ape::rTraitCont(tree, model = "BM", sigma = 1)
    cat3 <- factor(cut(lat_k, breaks = stats::quantile(lat_k, c(0, 1/3, 2/3, 1)),
                       labels = c("A", "B", "C"), include.lowest = TRUE))
    df <- cbind(cont, data.frame(bin = bin, cat3 = cat3, row.names = sp))
    return(list(df = df, tree = tree))
  }

  # ---- NEW corrected design (Section A) --------------------------------------------------
  tree <- ape::rcoal(n)
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  sp <- tree$tip.label
  if (dgp == "types_mixed") {
    K <- 7L + if (driver) 1L else 0L
    L <- sim_latents(tree, K, lambda, rho, evo_use, driver_col = if (driver) K else NULL)
    df <- data.frame(row.names = sp,
      c1 = L[, 1], c2 = L[, 2],
      cnt = as.integer(stats::rpois(n, exp(1.5 + 0.8 * L[, 3]))),
      prp = stats::plogis(L[, 4] + stats::rnorm(n, 0, 0.3)),
      bin = factor(ifelse(threshold_binary(L[, 5], thresholds), "yes", "no")),
      ord = threshold_ordinal(L[, 6], thresholds, c("L1", "L2", "L3", "L4")),
      cat3 = threshold_categorical(L[, 7], thresholds, c("A", "B", "C")))
    df$prp <- pmin(pmax(df$prp, 1e-4), 1 - 1e-4)
    if (driver) df$d1 <- L[, 8]
    return(list(df = df, tree = tree, trait_types = c(prp = "proportion"),
                L = L, lambda = lambda, rho = rho, evo = evo_use, thresholds = thresholds))
  }
  # bm_mixed / ou_mixed under the corrected design: same trait mix as before (4 continuous + binary
  # + categorical), now drawn through the shared lambda/rho/evo latent machinery.
  K <- 6L + if (driver) 1L else 0L
  L <- sim_latents(tree, K, lambda, rho, evo_use, driver_col = if (driver) K else NULL)
  df <- data.frame(row.names = sp, c1 = L[, 1], c2 = L[, 2], c3 = L[, 3], c4 = L[, 4],
                    bin = factor(ifelse(threshold_binary(L[, 5], thresholds), "yes", "no")),
                    cat3 = threshold_categorical(L[, 6], thresholds, c("A", "B", "C")))
  if (driver) df$d1 <- L[, 7]
  list(df = df, tree = tree, L = L, lambda = lambda, rho = rho, evo = evo_use, thresholds = thresholds)
}

# Build the cell: truth, tree, the seeded mask, df_miss. mask + set.seed(seed + 1000L) convention
# unchanged. miss = "mcar" with all defaults reproduces the OLD code exactly (mask_cols == names(truth)
# whenever there is no "d1" driver column, which is the case unless driver = TRUE).
make_cell <- function(dgp, n, seed, miss_frac = 0.30, miss = "mcar", lambda = 1, rho = 0,
                       evo = NULL, thresholds = "sample", driver = FALSE) {
  d <- make_dgp(dgp, n, seed, lambda = lambda, rho = rho, evo = evo, thresholds = thresholds, driver = driver)
  truth <- d$df; tree <- d$tree
  set.seed(seed + 1000L)
  mask <- matrix(FALSE, nrow(truth), ncol(truth), dimnames = dimnames(truth))
  mask_cols <- setdiff(names(truth), "d1")   # the driver is always observed
  realised_frac <- miss_frac

  if (miss == "mcar") {
    for (v in mask_cols) {
      obs <- which(!is.na(truth[[v]])); hide <- sample(obs, ceiling(miss_frac * length(obs)))
      mask[hide, v] <- TRUE
    }
  } else if (miss == "mar") {
    if (!("d1" %in% names(truth))) stop("make_cell: miss = 'mar' requires driver = TRUE")
    d1 <- truth$d1; d1z <- (d1 - mean(d1)) / stats::sd(d1)
    # P(miss_ij) = plogis(a + log(3) * d1z_i); solve a so the expected fraction equals miss_frac.
    f <- function(a) mean(stats::plogis(a + log(3) * d1z)) - miss_frac
    a_star <- stats::uniroot(f, c(-30, 30))$root
    p <- stats::plogis(a_star + log(3) * d1z)
    for (v in mask_cols) {
      obs <- which(!is.na(truth[[v]]))
      hide <- obs[stats::runif(length(obs)) < p[obs]]
      mask[hide, v] <- TRUE
    }
    realised_frac <- sum(mask[, mask_cols, drop = FALSE]) / (length(mask_cols) * nrow(truth))
  } else if (miss == "clade") {
    total_cells <- length(mask_cols) * nrow(truth)
    target <- miss_frac * total_cells
    lo_sz <- max(1L, ceiling(0.05 * nrow(truth))); hi_sz <- max(lo_sz, ceiling(0.15 * nrow(truth)))
    internal_nodes <- seq(nrow(truth) + 2L, nrow(truth) + tree$Nnode)
    tries <- 0L
    while (sum(mask[, mask_cols, drop = FALSE]) < target && tries < 2000L && length(internal_nodes)) {
      tries <- tries + 1L
      nd <- sample(internal_nodes, 1)
      cl <- tryCatch(ape::extract.clade(tree, nd), error = function(e) NULL)
      if (is.null(cl)) next
      sz <- length(cl$tip.label)
      if (sz < lo_sz || sz > hi_sz) next
      rows <- match(cl$tip.label, rownames(truth))
      for (v in mask_cols) mask[rows, v] <- TRUE   # last clade may overshoot -> "partial" below
    }
    # trim overshoot cell-by-cell (random) back toward the target fraction ("last clade partial")
    over <- sum(mask[, mask_cols, drop = FALSE]) - target
    if (over > 0) {
      idx_on <- which(mask[, mask_cols, drop = FALSE])
      unmask <- sample(idx_on, min(length(idx_on), floor(over)))
      mm <- mask[, mask_cols, drop = FALSE]; mm[unmask] <- FALSE; mask[, mask_cols] <- mm
    }
    # enforce >= 5 observed cells per column
    for (v in mask_cols) {
      if (sum(!mask[, v]) < 5) {
        need <- 5 - sum(!mask[, v])
        on_idx <- which(mask[, v])
        if (length(on_idx)) mask[sample(on_idx, min(need, length(on_idx))), v] <- FALSE
      }
    }
    realised_frac <- sum(mask[, mask_cols, drop = FALSE]) / total_cells
  } else {
    stop("unknown miss mechanism ", miss)
  }

  df_miss <- truth; for (v in names(truth)) df_miss[mask[, v], v] <- NA
  list(truth = truth, tree = tree, mask = mask, df_miss = df_miss,
       cont_traits = names(truth)[vapply(truth, is.numeric, logical(1))],
       trait_types = d$trait_types, realised_frac = realised_frac, mechanism = miss, L = d$L)
}

# ---- scoring ------------------------------------------------------------------------------
# score_arm() keeps its OLD signature and return type (a plain data.frame; campaign_solver_cell.R
# calls it positionally and rbinds the result). Section C additions: standardised interval width
# and the Gneiting-Raftery interval score for continuous-family traits; macro-F1 and Brier for
# discrete traits. The per-cell (confidence, correct) calibration pairs used for pooled ECE are
# attached as attr(., "calib") -- a long data.frame -- rather than changed in the return shape, so
# every existing caller keeps working unmodified.
#
# `prob`: optional named list keyed by trait name, each element an (n_obs x K) probability matrix
# with rownames = df_miss rownames and colnames = levels(truth[[v]]). Needed for Brier/ECE; accuracy
# and macro-F1 do not need it.
score_arm <- function(arm, completed, lower = NULL, upper = NULL, prob = NULL) {
  rows <- list(); calib_rows <- list()
  alpha <- 0.05
  for (v in names(truth)) {
    idx <- which(mask[, v]); if (!length(idx)) next
    rn <- rownames(truth)[idx]
    if (is.numeric(truth[[v]])) {
      tr <- truth[[v]][idx]; pr <- completed[[v]][idx]
      train <- truth[[v]][!mask[, v] & !is.na(truth[[v]])]
      sdt <- stats::sd(train)
      z <- sqrt(mean(((tr - pr) / sdt)^2))
      cov <- NA_real_; width <- NA_real_; iscore <- NA_real_
      if (!is.null(lower) && v %in% colnames(lower)) {
        lo <- lower[rn, v]; hi <- upper[rn, v]
        cov <- mean(tr >= lo & tr <= hi)
        width <- mean((hi - lo) / sdt)
        iscore <- mean(((hi - lo) + (2 / alpha) * (lo - tr) * (tr < lo) +
                           (2 / alpha) * (tr - hi) * (tr > hi)) / sdt)
      }
      rows[[length(rows) + 1L]] <- data.frame(arm = arm, trait = v, metric = "zRMSE", value = z, coverage = cov)
      rows[[length(rows) + 1L]] <- data.frame(arm = arm, trait = v, metric = "width", value = width, coverage = NA_real_)
      rows[[length(rows) + 1L]] <- data.frame(arm = arm, trait = v, metric = "interval_score", value = iscore, coverage = NA_real_)
    } else {
      truv <- as.character(truth[[v]][idx]); prv <- as.character(completed[[v]][idx])
      acc <- mean(truv == prv)
      levs <- levels(truth[[v]])
      # Macro-F1 averages only over classes PRESENT in the masked truth. A class absent from both
      # truth and prediction has an undefined F1; scoring it 0 and averaging it in deflates the
      # metric, and deflates it most in the low-prevalence cells where the arms actually differ.
      levs_scored <- levs[levs %in% truv]
      f1s <- vapply(levs_scored, function(k) {
        tp <- sum(prv == k & truv == k); fp <- sum(prv == k & truv != k); fn <- sum(prv != k & truv == k)
        prec <- if (tp + fp == 0) 0 else tp / (tp + fp); rec <- if (tp + fn == 0) 0 else tp / (tp + fn)
        if (prec + rec == 0) 0 else 2 * prec * rec / (prec + rec)
      }, numeric(1))
      macro_f1 <- if (length(f1s)) mean(f1s) else NA_real_
      n_classes_scored <- length(levs_scored)
      brier <- NA_real_
      if (!is.null(prob) && v %in% names(prob) && !is.null(prob[[v]])) {
        pv <- prob[[v]]
        # pigauto's binary probabilities come back as a vector (P(second level)), SOMETIMES named
        # (gnn = FALSE observed so far) and sometimes unnamed (gnn = TRUE observed so far, same
        # length as nrow(truth)); UNVERIFIED assumption: an unnamed vector/matrix is in
        # rownames(truth) order (matches how `completed` is reindexed in run_pigauto). Normalise
        # every source to an (n x K) matrix keyed by species name before indexing.
        if (is.null(dim(pv))) {
          stopifnot(length(levs) == 2L)
          nm <- if (!is.null(names(pv))) names(pv) else rownames(truth)
          pm2 <- cbind(1 - pv, pv); colnames(pm2) <- levs; rownames(pm2) <- nm
          pv <- pm2
        } else if (is.null(rownames(pv))) {
          rownames(pv) <- rownames(truth)
        }
        pm <- pv[rn, levs, drop = FALSE]
        y1 <- sapply(levs, function(k) as.integer(truv == k))
        brier <- mean(rowSums((pm - y1)^2))
        conf <- apply(pm, 1, max)
        pred_class <- levs[apply(pm, 1, which.max)]
        correct <- as.integer(pred_class == truv)
        calib_rows[[length(calib_rows) + 1L]] <- data.frame(arm = arm, trait = v, species = rn,
                                                              confidence = conf, correct = correct)
      }
      rows[[length(rows) + 1L]] <- data.frame(arm = arm, trait = v, metric = "accuracy", value = acc, coverage = NA_real_)
      rows[[length(rows) + 1L]] <- data.frame(arm = arm, trait = v, metric = "macroF1", value = macro_f1, coverage = NA_real_)
      rows[[length(rows) + 1L]] <- data.frame(arm = arm, trait = v, metric = "n_classes_scored", value = n_classes_scored, coverage = NA_real_)
      rows[[length(rows) + 1L]] <- data.frame(arm = arm, trait = v, metric = "brier", value = brier, coverage = NA_real_)
    }
  }
  res <- do.call(rbind, rows); rownames(res) <- NULL
  attr(res, "calib") <- if (length(calib_rows)) do.call(rbind, calib_rows) else NULL
  res
}


# ---- frequentist stack: Rphylopars on the continuous-family columns (proportion on logit,
# back-transformed), castor Mk (ML hidden-state prediction, likelihoods kept as `prob`) on each
# discrete trait separately. Section E: counts are now a phylogenetic Poisson GEE marginal model
# (`phylolm::phyloglm(method = "poisson_GEE")`) by default -- documented as a MARGINAL Poisson model
# with NO tip-level phylogenetic random effect (phyloglm does not provide one); the old log1p +
# Rphylopars route is kept available via count_method = "log1p_rphylopars" (arm "freq_log1p").
# Section D: continuous-family intervals are yhat +/- 1.96 * sqrt(anc_var [+ phenocov if not already
# included]), back-transformed. CHECK (2026-09-20, Rphylopars 0.3.10): `fit$anc_var` for a species
# with NO observed data at all for a trait already reflects total predictive uncertainty at that tip
# (it includes the phenotypic/residual variance implicitly through the joint GLS predictive
# equations); `fit$pars$phenocov` is a SEPARATE estimate of the phenotypic covariance and adding its
# diagonal on top would double-count residual variance. Verified by comparing empirical coverage in
# gate G6 (frequentist-stack coverage lands near 0.95, not badly over-covered) rather than by reading
# Rphylopars' internals line-by-line -- flagged UNVERIFIED against Rphylopars source.
run_freq <- function(df_miss, truth, mask, tree, cont_traits, trait_types = NULL,
                      count_method = c("phyloglm_poisson_gee", "log1p_rphylopars")) {
  count_method <- match.arg(count_method)
  comp <- truth; comp[] <- NA
  is_prop <- names(truth) %in% names(trait_types)[trait_types == "proportion"]
  is_cnt  <- vapply(truth, is.integer, logical(1))
  count_traits <- names(truth)[is_cnt]
  cont_for_joint <- if (count_method == "phyloglm_poisson_gee") setdiff(cont_traits, count_traits) else cont_traits
  tf  <- function(v, x) if (is_prop[match(v, names(truth))]) stats::qlogis(x) else if (count_method == "log1p_rphylopars" && is_cnt[match(v, names(truth))]) log1p(x) else x
  itf <- function(v, x) if (is_prop[match(v, names(truth))]) stats::plogis(x) else if (count_method == "log1p_rphylopars" && is_cnt[match(v, names(truth))]) pmax(expm1(x), 0) else x

  lower <- upper <- NULL
  if (length(cont_for_joint)) {
    df4 <- df_miss[, cont_for_joint, drop = FALSE]
    for (v in cont_for_joint) df4[[v]] <- tf(v, as.numeric(df4[[v]]))
    df_in <- data.frame(species = rownames(df4), df4, stringsAsFactors = FALSE)
    fit <- Rphylopars::phylopars(df_in, tree = tree, model = "BM", phylo_correlated = TRUE, pheno_correlated = TRUE, REML = TRUE)
    rec <- fit$anc_recon[rownames(df4), cont_for_joint, drop = FALSE]
    for (v in cont_for_joint) { comp[[v]] <- df_miss[[v]]; comp[mask[, v], v] <- itf(v, rec[mask[, v], v]) }
    if (!is.null(fit$anc_var)) {
      var_out <- fit$anc_var[rownames(df4), cont_for_joint, drop = FALSE]
      lower <- matrix(NA_real_, nrow(truth), ncol(truth), dimnames = dimnames(truth)); upper <- lower
      for (v in cont_for_joint) {
        se <- sqrt(pmax(var_out[, v], 0))
        lower[rownames(df4), v] <- itf(v, rec[, v] - 1.96 * se)
        upper[rownames(df4), v] <- itf(v, rec[, v] + 1.96 * se)
      }
    }
  }
  prob <- list()
  if (length(count_traits) && count_method == "phyloglm_poisson_gee") {
    if (is.null(lower)) { lower <- matrix(NA_real_, nrow(truth), ncol(truth), dimnames = dimnames(truth)); upper <- lower }
    for (v in count_traits) {
      y <- df_miss[[v]]; obs <- !is.na(y)
      tree_obs <- ape::keep.tip(tree, rownames(df_miss)[obs])
      dat <- data.frame(y = y[obs], row.names = rownames(df_miss)[obs])
      dat <- dat[tree_obs$tip.label, , drop = FALSE]
      fitp <- tryCatch(phylolm::phyloglm(y ~ 1, phy = tree_obs, data = dat, method = "poisson_GEE"), error = function(e) NULL)
      comp[[v]] <- df_miss[[v]]
      if (!is.null(fitp)) {
        # Documented limitation: phyloglm's poisson_GEE fits a MARGINAL Poisson mean structure
        # (beta0 only, no random effect / BLUP per species); the same exp(beta0) is predicted for
        # every missing cell of this trait, and the interval is a plain Poisson quantile band, not
        # a tip-specific phylogenetic prediction interval.
        lam <- exp(unname(stats::coef(fitp))[1])
        comp[mask[, v], v] <- as.integer(round(lam))
        qs <- stats::qpois(c(0.025, 0.975), lam)
        lower[mask[, v], v] <- qs[1]; upper[mask[, v], v] <- qs[2]
      } else {
        obsv <- df_miss[[v]][!is.na(df_miss[[v]])]
        comp[mask[, v], v] <- as.integer(round(mean(obsv)))
      }
    }
  }
  # discrete traits: castor Mk, equal rates for binary/categorical, stepwise (SUEDE) for ordinal;
  # normalised state likelihoods kept as `prob` for Brier / ECE.
  for (v in setdiff(names(truth), c(cont_traits))) {
    f <- df_miss[[v]]; lev <- levels(f); tip <- as.integer(f)  # NA where masked
    tip_full <- tip[match(tree$tip.label, rownames(df_miss))]
    rm <- if (is.ordered(f)) "SUEDE" else "ER"
    h <- castor::hsp_mk_model(tree, tip_states = tip_full, Nstates = length(lev), rate_model = rm, Ntrials = 3, Nthreads = 1)
    L <- h$likelihoods[seq_along(tree$tip.label), , drop = FALSE]
    L <- L / rowSums(L)
    pred <- max.col(L)
    pred_sp <- pred[match(rownames(df_miss), tree$tip.label)]
    out <- as.character(f); out[mask[, v]] <- lev[pred_sp[mask[, v]]]
    comp[[v]] <- factor(out, levels = lev, ordered = is.ordered(f))
    Lr <- L[match(rownames(df_miss), tree$tip.label), , drop = FALSE]
    colnames(Lr) <- lev; rownames(Lr) <- rownames(df_miss)
    prob[[v]] <- Lr
  }
  list(completed = comp, lower = lower, upper = upper, prob = prob)
}

# ---- phylogeny + machine-learning hybrid (the approach of Gendre, Hauffe, Pimiento and Silvestro 2024, MEE,
# TDIP's missForest_phylo): missForest on the trait table plus phylogenetic eigenvectors. Eigenvectors are the
# principal coordinates of the cophenetic distance matrix (PVR-style), keeping axes that explain `variance_fraction`
# of the variance (TDIP's default idea), capped at `max_axes`.
phylo_eigenvectors <- function(tree, variance_fraction = 0.9, max_axes = 20L) {
  D <- ape::cophenetic.phylo(tree); n <- nrow(D)
  J <- diag(n) - matrix(1 / n, n, n); B <- -0.5 * J %*% (D^2) %*% J
  e <- eigen(B, symmetric = TRUE); keep <- e$values > 1e-8
  vals <- e$values[keep]; vecs <- e$vectors[, keep, drop = FALSE]
  k <- min(max_axes, which(cumsum(vals) / sum(vals) >= variance_fraction)[1])
  ev <- vecs[, seq_len(k), drop = FALSE] %*% diag(sqrt(vals[seq_len(k)]), k, k)
  rownames(ev) <- rownames(D); colnames(ev) <- paste0("pev", seq_len(k)); ev
}
run_mf_phylo <- function(df_miss, truth, mask, tree, variance_fraction = 0.9) {
  ev <- phylo_eigenvectors(tree, variance_fraction)[rownames(df_miss), , drop = FALSE]
  X <- cbind(df_miss, as.data.frame(ev))
  for (v in names(df_miss)) if (is.ordered(X[[v]])) X[[v]] <- factor(as.character(X[[v]]), levels = levels(df_miss[[v]]))
  imp <- missForest::missForest(X, maxiter = 10, ntree = 100, verbose = FALSE)$ximp
  comp <- truth; comp[] <- NA
  for (v in names(df_miss)) {
    x <- imp[[v]]
    if (is.factor(df_miss[[v]])) comp[[v]] <- factor(as.character(x), levels = levels(df_miss[[v]]), ordered = is.ordered(df_miss[[v]]))
    else if (is.integer(df_miss[[v]])) comp[[v]] <- as.integer(round(x)) else comp[[v]] <- x
  }
  list(completed = comp)
}

# ---- BACE arm (Section D). BACE::bace()'s top-level API returns only `n_final` full imputed
# datasets (`$imputed_datasets`); CHECK (2026-09-20, `body(BACE::bace)`, `ls(asNamespace("BACE"))`):
# there is no exported accessor for the raw MCMCglmm Liab/Sol chains, so the "full retained posterior
# predictive samples" the design asks for ARE the imputed datasets -- we simply request as many of
# them as the chain length affords. n_final is set to floor((nitt - burnin) / thin) * runs, capped at
# 400 (previously a fixed 5). The interval built from these draws (2.5/97.5 percentiles for
# continuous-family traits, class frequency for discrete traits) is a POSTERIOR PREDICTIVE interval
# for the masked cell, not a posterior interval of the mean -- UNVERIFIED against BACE's own
# documentation of `$imputed_datasets`'s exact sampling distribution beyond the source inspected here.
run_bace <- function(df_miss, truth, mask, tree, cont_traits, bace_nitt, bace_burnin, bace_thin, bace_runs = 2L) {
  tree_b <- tree; if (any(tree_b$edge.length == 0)) tree_b$edge.length[tree_b$edge.length == 0] <- 1e-8
  df_b <- df_miss; df_b$Species <- rownames(df_miss)
  all_traits <- setdiff(names(df_b), "Species")
  fixformula <- lapply(all_traits, function(v) paste0(v, " ~ ", paste(setdiff(all_traits, v), collapse = " + ")))
  # n_final is the number of FULL imputation runs (bace_final_imp refits MCMCglmm per response for
  # each one), not a thinning of one chain: each dataset is one posterior predictive draw (K = 1L in
  # BACE:::.predict_bace). Cost is linear in n_final, so it is a budget constant, not a chain length.
  # 50 gives usable 2.5/97.5 percentiles; the pre-run times it at production size. n_cores = 1L is
  # BACE's own default, pinned explicitly: parallelism here is one process per cell (xargs on Totoro,
  # array tasks on DRAC), and a fork inside MCMCglmm segfaults when the caller has already forked.
  n_final <- as.integer(Sys.getenv("PIG_BACE_NFINAL", "50"))
  outb <- BACE::bace(fixformula = fixformula, ran_phylo_form = "~ 1 |Species", phylo = tree_b,
                     data = df_b, nitt = bace_nitt, burnin = bace_burnin, thin = bace_thin,
                     runs = bace_runs, n_final = n_final, n_cores = 1L,
                     verbose = FALSE, skip_conv = TRUE, ovr_categorical = TRUE)
  sets <- if ("imputed_datasets" %in% names(outb)) outb$imputed_datasets else
          if ("imputed_data" %in% names(outb)) outb$imputed_data else list(outb$data)
  comp <- df_miss
  lower <- upper <- matrix(NA_real_, nrow(truth), ncol(truth), dimnames = dimnames(truth))
  prob <- list()
  for (v in names(comp)) {
    idx <- which(mask[, v]); if (!length(idx)) next
    draws <- sapply(sets, function(s) s[[v]][idx])
    if (!is.matrix(draws)) draws <- matrix(draws, ncol = length(sets))
    rn <- rownames(truth)[idx]
    if (is.numeric(comp[[v]])) {
      comp[idx, v] <- apply(draws, 1, stats::median)
      qs <- t(apply(draws, 1, stats::quantile, probs = c(0.025, 0.975)))
      lower[rn, v] <- qs[, 1]; upper[rn, v] <- qs[, 2]
    } else {
      comp[idx, v] <- apply(draws, 1, function(x) names(which.max(table(as.character(x)))))
      levs <- levels(truth[[v]])
      pm <- t(apply(draws, 1, function(x) { t <- table(factor(as.character(x), levels = levs)); as.numeric(t) / sum(t) }))
      colnames(pm) <- levs; rownames(pm) <- rn
      full <- matrix(NA_real_, nrow(truth), length(levs), dimnames = list(rownames(truth), levs))
      full[rn, ] <- pm
      prob[[v]] <- full
    }
  }
  list(completed = comp, lower = lower, upper = upper, prob = prob, n_final = n_final,
       diag = list(bace_rhat = bace_rhat(outb, bace_runs)))
}

# Gelman-Rubin Rhat over BACE's two chains, on the fixed effects and variance components only (the
# species-level random effects are thousands of nuisance parameters whose Rhat is not the convergence
# question). Models are grouped by an identical (fixed-effect names, VCV names) signature so the two
# chains of the same response are compared; OVR/categorical fits differ in their reference level
# and would otherwise be paired wrongly. Returns NULL when fewer than two chains are found.
bace_rhat <- function(outb, runs) {
  if (runs < 2L || !requireNamespace("coda", quietly = TRUE)) return(NULL)
  find_mcmcglmm <- function(x) {
    if (inherits(x, "MCMCglmm")) return(list(x))
    if (is.list(x)) return(unlist(lapply(x, find_mcmcglmm), recursive = FALSE))
    NULL
  }
  models <- find_mcmcglmm(outb)
  if (length(models) < 2L) return(NULL)
  par_mat <- function(m) {
    nfl <- m$Fixed$nfl %||% ncol(m$Sol)
    cbind(as.matrix(m$Sol)[, seq_len(nfl), drop = FALSE], as.matrix(m$VCV))
  }
  sig <- vapply(models, function(m) paste(colnames(par_mat(m)), collapse = "|"), "")
  out <- list()
  for (s in unique(sig)) {
    grp <- models[sig == s]
    if (length(grp) < 2L) next
    ml <- tryCatch(coda::mcmc.list(lapply(grp[1:2], function(m) coda::as.mcmc(par_mat(m)))),
                   error = function(e) NULL)
    if (is.null(ml)) next
    r <- tryCatch(coda::gelman.diag(ml, multivariate = FALSE, autoburnin = FALSE)$psrf[, 1],
                  error = function(e) NULL)
    if (!is.null(r)) out[[length(out) + 1L]] <- r
  }
  if (!length(out)) return(NULL)
  r <- unlist(out)
  list(max = max(r, na.rm = TRUE), n = length(r), frac_above_1.1 = mean(r > 1.1, na.rm = TRUE), psrf = r)
}

run_floor <- function(df_miss, truth, mask) {
  comp <- df_miss
  for (v in names(comp)) {
    obs <- df_miss[[v]][!is.na(df_miss[[v]])]
    fill <- if (is.numeric(comp[[v]])) mean(obs) else names(which.max(table(as.character(obs))))
    comp[mask[, v], v] <- fill
  }
  prob <- list()
  for (v in names(comp)) if (is.factor(truth[[v]])) {
    levs <- levels(truth[[v]])
    p <- prop.table(table(factor(as.character(df_miss[[v]][!is.na(df_miss[[v]])]), levels = levs)))
    pm <- matrix(rep(as.numeric(p), each = nrow(truth)), nrow(truth), length(levs), dimnames = list(rownames(truth), levs))
    prob[[v]] <- pm
  }
  list(completed = comp, prob = prob)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- run_arms(): the arm-dispatch loop factored out of campaign_gnn_off_cell.R so that it and the
# new campaign_sim_cell.R share ONE implementation (Section G). `cell` is a make_cell() result plus
# `tree`; `opts` carries epochs / bace_nitt / bace_burnin / bace_thin / seed / trait_types / log_line.
# Returns list(results = <rbind of score_arm() data.frames, dgp/n/seed columns attached>,
#              calib = <rbind of attr(., "calib") long frames, may be NULL>,
#              walls =, errors =, paths =).
run_arms <- function(cell, arms, opts) {
  truth <- cell$truth; tree <- cell$tree; mask <- cell$mask; df_miss <- cell$df_miss
  cont_traits <- cell$cont_traits; trait_types <- cell$trait_types
  # score_arm() (defined at this file's top level) reads `truth`/`mask` from ITS OWN lexical scope,
  # i.e. the global environment where this file is normally sourced -- not from run_arms()'s local
  # frame. Publish them there so score_arm() sees the right cell regardless of which caller script
  # invoked run_arms(). Restored on exit so a caller's own global `truth`/`mask` (if any) are unaffected.
  old_truth <- if (exists("truth", envir = .GlobalEnv, inherits = FALSE)) get("truth", envir = .GlobalEnv) else NULL
  old_mask  <- if (exists("mask",  envir = .GlobalEnv, inherits = FALSE)) get("mask",  envir = .GlobalEnv) else NULL
  assign("truth", truth, envir = .GlobalEnv); assign("mask", mask, envir = .GlobalEnv)
  on.exit({
    if (is.null(old_truth)) rm(list = "truth", envir = .GlobalEnv) else assign("truth", old_truth, envir = .GlobalEnv)
    if (is.null(old_mask))  rm(list = "mask",  envir = .GlobalEnv) else assign("mask",  old_mask,  envir = .GlobalEnv)
  }, add = TRUE)
  seed <- opts$seed; epochs <- opts$epochs %||% 2000L
  bace_nitt <- opts$bace_nitt %||% 50000L; bace_burnin <- opts$bace_burnin %||% 10000L; bace_thin <- opts$bace_thin %||% 25L
  log_line <- opts$log_line %||% function(...) invisible(NULL)

  run_pigauto <- function(arm) {
    extra <- switch(arm,
      gnn_on              = list(gnn = TRUE,  epochs = epochs),
      gnn_off             = list(gnn = FALSE),
      gnn_off_pure        = list(gnn = FALSE, safety_floor = FALSE, phylo_signal_gate = FALSE),
      gnn_off_rphylopars  = list(gnn = FALSE, joint_solver = "rphylopars"))
    res <- do.call(pigauto::impute, c(list(traits = df_miss, tree = tree, verbose = FALSE, seed = seed,
                                           trait_types = trait_types), extra))
    path <- res$fit$baseline$path
    pred <- res$prediction
    comp <- res$completed[rownames(truth), names(truth)]
    prob <- pred$probabilities
    out <- list(completed = comp, lower = pred$conformal_lower, upper = pred$conformal_upper,
                path = path, prob = prob, se = pred$se)
    if (arm == "gnn_on") {
      # Plan arm 2: the SAME GNN-on fit predicting from the tax-free baseline (baseline_override),
      # so the GNN effect and the held-out-cell tax can be separated. No refit.
      bf <- pigauto::fit_baseline(res$data, tree, splits = NULL,
                                  lambda_mode = res$fit$model_config$lambda_mode %||% "fixed_1",
                                  joint_solver = res$fit$model_config$joint_solver %||% "inhouse")
      pred_f <- stats::predict(res$fit, return_se = TRUE, baseline_override = bf)
      comp_f <- comp
      for (v in names(truth)) comp_f[mask[, v], v] <- pred_f$imputed[rownames(truth)[mask[, v]], v]
      out$derived <- list(gnn_on_full = list(completed = comp_f, lower = pred_f$conformal_lower,
                                             upper = pred_f$conformal_upper, path = path,
                                             prob = pred_f$probabilities, se = pred_f$se))
    }
    out
  }
  run_rphylopars <- function() {
    df4 <- df_miss[, cont_traits, drop = FALSE]
    df_in <- data.frame(species = rownames(df4), df4, stringsAsFactors = FALSE)
    fit <- Rphylopars::phylopars(df_in, tree = tree, model = "BM", phylo_correlated = TRUE,
                                 pheno_correlated = TRUE, REML = TRUE)
    comp <- truth; comp[] <- NA
    rec <- fit$anc_recon[rownames(df4), cont_traits, drop = FALSE]
    for (v in cont_traits) { comp[[v]] <- df_miss[[v]]; comp[mask[, v], v] <- rec[mask[, v], v] }
    list(completed = comp)
  }

  results <- list(); calibs <- list(); walls <- list(); errors <- list(); paths <- list(); failed <- list(); diags <- list()
  for (arm in arms) {
    t0 <- Sys.time()
    r <- tryCatch({
      if (arm %in% c("gnn_on", "gnn_off", "gnn_off_pure", "gnn_off_rphylopars")) run_pigauto(arm)
      else if (arm == "rphylopars") run_rphylopars()
      else if (arm == "freq") run_freq(df_miss, truth, mask, tree, cont_traits, trait_types)
      else if (arm == "freq_log1p") run_freq(df_miss, truth, mask, tree, cont_traits, trait_types, count_method = "log1p_rphylopars")
      else if (arm == "mf_phylo") run_mf_phylo(df_miss, truth, mask, tree)
      else if (arm == "bace") run_bace(df_miss, truth, mask, tree, cont_traits, bace_nitt, bace_burnin, bace_thin)
      else if (arm == "floor") run_floor(df_miss, truth, mask)
      else stop("unknown arm ", arm)
    }, error = function(e) e)
    walls[[arm]] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(r, "error")) {
      errors[[arm]] <- conditionMessage(r); failed[[arm]] <- TRUE
      log_line("%s ERROR %s", arm, errors[[arm]])
      # a failed arm is scored at the floor, never dropped (Section G)
      fl <- tryCatch(run_floor(df_miss, truth, mask), error = function(e) NULL)
      if (!is.null(fl)) {
        sc <- score_arm(arm, fl$completed, prob = fl$prob)
        results[[arm]] <- sc; calibs[[arm]] <- attr(sc, "calib")
      }
      next
    }
    sc <- score_arm(arm, r$completed, r$lower, r$upper, prob = r$prob)
    results[[arm]] <- sc; calibs[[arm]] <- attr(sc, "calib")
    paths[[arm]] <- r$path
    if (!is.null(r$diag)) diags[[arm]] <- r$diag
    for (dn in names(r$derived)) {
      dr <- r$derived[[dn]]
      scd <- score_arm(dn, dr$completed, dr$lower, dr$upper, prob = dr$prob)
      results[[dn]] <- scd; calibs[[dn]] <- attr(scd, "calib")
      paths[[dn]] <- dr$path; walls[[dn]] <- 0
    }
    log_line("%s done in %.1f s", arm, walls[[arm]])
  }
  tab <- do.call(rbind, results); if (!is.null(tab)) rownames(tab) <- NULL
  calib <- do.call(rbind, calibs)
  list(results = tab, calib = calib, walls = walls, errors = errors, paths = paths, failed = failed, diag = diags)
}
