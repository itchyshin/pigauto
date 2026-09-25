# Pre-registered downstream slope-check pairs for the mi-posterior
# real-data validation (G8, .unlazy/mi-posterior/GATES.md). The pair list
# was written 2026-09-24, BEFORE any multi_impute(draws_method =
# "posterior") result existed on real data, so the pairs are locked in
# before the numbers are seen.
#
# Design that runs (01_run.R), stated 2026-09-24. This replaces the earlier
# header text, which described an MI leg of nlme::gls(corPagel) plus
# pool_mi() on ALL rows (observed + imputed) and an unpooled phylolm
# check. Neither is run on real data: dense gls(corPagel) is not feasible
# at these n, and pool_mi() has no phylolm adapter.
#
# For each pair, both legs are PAIRED: same species, same analysis model,
# same analysis scale. The only difference between them is that masked
# cells are imputed in the MI leg.
#   species   = rows where BOTH y and x are ORIGINALLY observed, i.e.
#     non-missing in the Mondrian mask receipt's `truth` (pre-masking),
#     never `masked`, and finite after the analysis-scale transform below;
#     the tree is ape::keep.tip(tree, species).
#   reference = phylolm(y ~ x, model = "lambda") on `truth` for those
#     species (the complete-rows slope).
#   MI slope  = the same phylolm(y ~ x, model = "lambda") fitted to each
#     of the m = 20 multi_impute(draws_method = "posterior") completions of
#     the MASKED data, on the same species and pruned tree, pooled by hand
#     with Rubin's rules and the Barnard-Rubin (1999) small-sample df with
#     complete-data df nu_com = n - 2.
#
# Analysis scale, pre-registered per pair (added 2026-09-24, still before
# any real-data result exists): `y_log` / `x_log` say whether the pair
# analysis takes log() of the response / predictor. Both legs use these
# flags; pigauto's own trait_map$log_transform is NOT used for the
# analysis scale. (It logs every all-positive continuous column, so it
# would have fitted the already-logged PanTHERIA traits on a log(log())
# scale.) Imputation itself stays at pigauto defaults.
#   PanTHERIA: FALSE / FALSE (the columns are already log() values, see
#     00_prepare_realdata_input.R's `kind = "log"` traits).
#   AVONET size traits ~ Mass: TRUE / TRUE (raw measurements, log-log
#     allometry).
#   FishBase Weight ~ Length: TRUE / TRUE; Troph ~ Length: y_log = FALSE
#     (trophic level is an index on its own scale), x_log = TRUE;
#     DepthRangeDeep ~ Length: y_log = FALSE, x_log = TRUE (see below).
#
# Change on 2026-09-24, made after the AVONET smoke cell and BEFORE any
# FishBase result existed: DepthRangeDeep ~ Length moved from y_log = TRUE
# to y_log = FALSE. The masked FishBase data contain an observed
# DepthRangeDeep of 0, so pigauto's default auto-log (all observed values
# > 0) leaves that column on its raw scale and the posterior can impute
# values <= 0. log() of those is not finite, which would drop whole
# completions from the MI pool (selective pooling) or leave none. The
# analysis scale now matches the imputation scale; imputation settings
# stay at pigauto defaults. Every other log-scale trait in this list
# (AVONET Mass, Wing.Length, Tarsus.Length, Beak.Length_Culmen; FishBase
# Weight, Length) has all observed values > 0 in every mask receipt
# (checked 2026-09-24), so pigauto imputes it on the log scale and its
# completions stay positive.
#
# Column names below are the ACTUAL columns produced by
# script/mondrian_confirmation/00_prepare_realdata_input.R (verified
# against script/mondrian_confirmation/returned/*/mask_receipt.rds on
# arc/mondrian-realdata, 2026-09-24), not the shorthand in the task brief.

mi_realdata_pairs <- list(

  # -- PanTHERIA (mammals; body_mass_g, head_body_length_mm, gestation_d,
  #    max_longevity_m are already log() in 00_prepare_realdata_input.R) --
  list(
    dataset = "pantheria", response = "body_mass_g", predictor = "head_body_length_mm",
    y_log = FALSE, x_log = FALSE,
    rationale = paste(
      "Kleiber-adjacent size allometry: adult body mass scales with adult",
      "head-body length. Both traits are well populated in PanTHERIA and",
      "the relationship is one of the best-established in mammalian",
      "comparative biology -- a strong-signal anchor pair.")
  ),
  list(
    dataset = "pantheria", response = "gestation_d", predictor = "body_mass_g",
    y_log = FALSE, x_log = FALSE,
    rationale = paste(
      "Life-history allometry: gestation length increases with body mass",
      "along the mammalian slow-fast continuum. A second slope from the",
      "same predictor as pair 1, testing whether MI pooling holds up for a",
      "different response with a different missingness pattern.")
  ),
  list(
    dataset = "pantheria", response = "max_longevity_m", predictor = "body_mass_g",
    y_log = FALSE, x_log = FALSE,
    rationale = paste(
      "Life-history allometry: maximum recorded longevity increases with",
      "body mass. Noisier than gestation length (longevity is",
      "under-sampled and outlier-prone in PanTHERIA) -- a harder,",
      "more realistic test case.")
  ),

  # -- AVONET (birds; Mass, Wing.Length, Tarsus.Length,
  #    Beak.Length_Culmen are raw, not pre-logged by the prep script) --
  list(
    dataset = "avonet", response = "Wing.Length", predictor = "Mass",
    y_log = TRUE, x_log = TRUE,
    rationale = paste(
      "Flight-morphology allometry: wing length scales with body mass.",
      "AVONET's best-populated, strongest-signal morphology pair -- the",
      "anchor case for the bird dataset.")
  ),
  list(
    dataset = "avonet", response = "Tarsus.Length", predictor = "Mass",
    y_log = TRUE, x_log = TRUE,
    rationale = paste(
      "Skeletal allometry: tarsus length scales with body mass. A second,",
      "largely independent morphological slope from the same predictor as",
      "the wing-length pair.")
  ),
  list(
    dataset = "avonet", response = "Beak.Length_Culmen", predictor = "Mass",
    y_log = TRUE, x_log = TRUE,
    rationale = paste(
      "Beak allometry: culmen length scales with body mass, but more",
      "weakly and more noisily than wing or tarsus (beak shape is under",
      "strong ecological/dietary selection independent of size) -- a",
      "harder test case for the same dataset.")
  ),

  # -- FishBase (Length, Weight, DepthRangeDeep, Troph all raw numeric) --
  list(
    dataset = "fishbase", response = "Weight", predictor = "Length",
    y_log = TRUE, x_log = TRUE,
    rationale = paste(
      "Length-weight allometry: the standard fisheries cube-law scaling",
      "relationship, strongly phylogenetically structured -- the anchor",
      "case for the fish dataset.")
  ),
  list(
    dataset = "fishbase", response = "DepthRangeDeep", predictor = "Length",
    y_log = FALSE, x_log = TRUE,  # was TRUE / TRUE until 2026-09-24, see header
    rationale = paste(
      "Depth-size relationship: larger-bodied fish tend to reach greater",
      "maximum depths. Ecologically motivated but noisier and more",
      "heterogeneous across clades than length-weight -- a harder case.")
  ),
  list(
    dataset = "fishbase", response = "Troph", predictor = "Length",
    y_log = FALSE, x_log = TRUE,
    rationale = paste(
      "Trophic-level allometry: larger fish tend to occupy higher trophic",
      "positions. A weak, well-studied ecological slope near the",
      "detection limit -- a stress test for whether MI pooling preserves",
      "a small true effect rather than just not being obviously wrong on",
      "a strong one.")
  )
)

# Sanity: every pair's dataset is one of the three planned datasets, no
# pair names the same response twice within a dataset, and every pair
# carries a pre-registered analysis scale (a single TRUE/FALSE per leg).
is_flag <- function(v) is.logical(v) && length(v) == 1L && !is.na(v)
stopifnot(
  all(vapply(mi_realdata_pairs, function(p) p$dataset, character(1)) %in%
        c("pantheria", "avonet", "fishbase")),
  !anyDuplicated(vapply(mi_realdata_pairs, function(p) paste(p$dataset, p$response), character(1))),
  all(vapply(mi_realdata_pairs, function(p) is_flag(p$y_log) && is_flag(p$x_log), logical(1)))
)
rm(is_flag)
