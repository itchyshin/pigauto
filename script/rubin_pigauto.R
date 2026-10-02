# script/rubin_pigauto.R
#
# pigauto arm for the Rubin study: multi_impute(draws_method = "posterior"), proper Bayesian MI from the
# multivariate phylogenetic mixed model (PR #189). It fits the SAME continuous-family block as the
# frequentist arms (default_block_traits(): c1, c2, prp and the always-observed driver d1; counts excluded),
# on the SAME scale (prp through qlogis, via transform_block()), and writes draws back with
# apply_block_draw(), so observed cells stay bit-identical to df_miss. Non-block columns are left as in
# df_miss, as in mi_freq_A().
#
# No conformal arm: conformal draws bias the downstream slope (-0.20 to -0.46, coverage 0 to 17 percent,
# arc/mi-gls-attenuation, 16 regimes); Shinichi, 2026-10-01.
#
# Requires rubin_freq.R to be sourced first (default_block_traits, transform_block, apply_block_draw).

mi_pig_post <- function(cell, M, seed, control = list()) {
  if (!requireNamespace("pigauto", quietly = TRUE)) stop("pigauto is not installed")
  if (!"posterior" %in% eval(formals(pigauto::multi_impute)$draws_method))
    stop("installed pigauto has no draws_method = \"posterior\" (needs PR #189)")
  df_miss <- cell$df_miss; tree <- cell$tree; trait_types <- cell$trait_types
  block_traits <- default_block_traits(cell)
  is_prp <- block_traits %in% names(trait_types)[trait_types == "proportion"]
  Yt <- transform_block(df_miss, block_traits, is_prp)          # numeric, prp on the logit scale
  Yt <- as.data.frame(Yt); rownames(Yt) <- rownames(df_miss)
  ctrl <- utils::modifyList(list(seed = seed), control)
  # log_transform = FALSE: pigauto would otherwise log any all-positive column, putting this arm on a different
  # scale from freqA (review R1). The block is already on the frequentist fitting scale.
  mi <- pigauto::multi_impute(Yt, tree, m = M, draws_method = "posterior", log_transform = FALSE,
                              posterior_control = ctrl)
  logged <- names(Filter(function(tm) isTRUE(tm$log_transform), mi$data$trait_map))
  if (length(logged)) stop("pigauto log-transformed ", paste(logged, collapse = ", "), " despite log_transform = FALSE")
  datasets <- lapply(mi$datasets, function(d) {
    raw <- as.matrix(d[rownames(df_miss), block_traits, drop = FALSE])
    for (i in seq_along(block_traits)) if (is_prp[i]) raw[, i] <- stats::plogis(raw[, i])
    apply_block_draw(df_miss, block_traits, raw)
  })
  dg <- mi$posterior$diagnostics
  list(datasets = datasets,
       diag = list(block_traits = block_traits,
                   converged = isTRUE(mi$posterior$converged),
                   rhat_max = suppressWarnings(max(dg$rhat, na.rm = TRUE)),
                   ess_min = suppressWarnings(min(dg$ess_bulk, na.rm = TRUE)),
                   n_extensions = attr(dg, "n_extensions"),
                   wall_s = mi$posterior$wall_s,
                   mi_workflow = mi$mi_workflow,
                   control = mi$posterior$control,
                   logged = logged,
                   pigauto_version = as.character(utils::packageVersion("pigauto")),
                   pigauto_sha = pigauto_install_sha()))
}

# The sha of the INSTALLED pigauto (review R2): install it with remotes::install_github("itchyshin/pigauto@<sha>")
# so DESCRIPTION carries RemoteSha. A local install has none, and the arm then records NA, which the
# provenance gate rejects.
pigauto_install_sha <- function() {
  d <- utils::packageDescription("pigauto")
  sha <- d$RemoteSha %||% d$GithubSHA1
  if (is.null(sha) || !nzchar(sha)) NA_character_ else sha
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# pig_conf: pigauto's default conformal MI draws on the same block, as an opt-in NEGATIVE CONTROL only. Conformal draws
# are known to bias the downstream slope (arc/mi-gls-attenuation); Shinichi (2026-10-01) does not want them in the
# campaign, so this arm exists to show that failure on the Rubin datasets if asked, not as a contender.
mi_pig_conf <- function(cell, M, seed) {
  if (!requireNamespace("pigauto", quietly = TRUE)) stop("pigauto is not installed")
  df_miss <- cell$df_miss; tree <- cell$tree; trait_types <- cell$trait_types
  block_traits <- default_block_traits(cell)
  is_prp <- block_traits %in% names(trait_types)[trait_types == "proportion"]
  Yt <- as.data.frame(transform_block(df_miss, block_traits, is_prp)); rownames(Yt) <- rownames(df_miss)
  t0 <- proc.time()[["elapsed"]]
  mi <- pigauto::multi_impute(Yt, tree, m = M, draws_method = "conformal", log_transform = FALSE, seed = seed,
                              verbose = FALSE)
  datasets <- lapply(mi$datasets, function(d) {
    raw <- as.matrix(d[rownames(df_miss), block_traits, drop = FALSE])
    for (i in seq_along(block_traits)) if (is_prp[i]) raw[, i] <- stats::plogis(raw[, i])
    apply_block_draw(df_miss, block_traits, raw)
  })
  list(datasets = datasets,
       diag = list(block_traits = block_traits, draws_method = "conformal",
                   wall_s = proc.time()[["elapsed"]] - t0,
                   pigauto_version = as.character(utils::packageVersion("pigauto")),
                   pigauto_sha = pigauto_install_sha()))
}
