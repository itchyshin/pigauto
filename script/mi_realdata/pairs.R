# Pre-registered downstream slope-check pairs for the mi-posterior
# real-data validation (G8, .unlazy/mi-posterior/GATES.md; design.md
# section 5). Written 2026-09-24, BEFORE any
# multi_impute(draws_method = "posterior") result exists on real data --
# these are locked in before the numbers are seen.
#
# For each pair:
#   reference slope = PGLS via phylolm(y ~ x, model = "lambda"), fit on the
#     rows where BOTH y and x are ORIGINALLY observed, i.e. in the Mondrian
#     mask receipt's `truth` (pre-masking), never `masked`.
#   MI slope        = pool_mi() over the m = 20
#     multi_impute(draws_method = "posterior") completions of the MASKED
#     data, each fitted with nlme::gls(y ~ x, correlation = corPagel(1,
#     tree)) on ALL rows (originally observed + imputed). gls (not
#     phylolm) is used for the MI leg because pool_mi() ships an automatic
#     coefficient/vcov adapter for nlme::gls fits; phylolm fits are not in
#     its documented adapter list (see R/pool_mi.R). phylolm(model =
#     "lambda") is still run once per cell as an unpooled check, per
#     design.md section 5 ("nlme::gls(corPagel) as a check").
#   scale: both legs use the SAME per-trait transform, decided by pigauto's
#     own auto-detected `trait_map$log_transform` flag at fit time (read
#     off the multi_impute() fit, not hardcoded here) -- so the reference
#     and MI models can never silently diverge in scale. For PanTHERIA the
#     traits below are already on the log scale by construction (see
#     00_prepare_realdata_input.R's `kind = "log"` traits), so
#     log_transform there is expected to read FALSE at the pigauto layer
#     (it would be re-logging an already-logged variable) -- this is
#     recorded per pair in the result, not assumed here.
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
    rationale = paste(
      "Kleiber-adjacent size allometry: adult body mass scales with adult",
      "head-body length. Both traits are well populated in PanTHERIA and",
      "the relationship is one of the best-established in mammalian",
      "comparative biology -- a strong-signal anchor pair.")
  ),
  list(
    dataset = "pantheria", response = "gestation_d", predictor = "body_mass_g",
    rationale = paste(
      "Life-history allometry: gestation length increases with body mass",
      "along the mammalian slow-fast continuum. A second slope from the",
      "same predictor as pair 1, testing whether MI pooling holds up for a",
      "different response with a different missingness pattern.")
  ),
  list(
    dataset = "pantheria", response = "max_longevity_m", predictor = "body_mass_g",
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
    rationale = paste(
      "Flight-morphology allometry: wing length scales with body mass.",
      "AVONET's best-populated, strongest-signal morphology pair -- the",
      "anchor case for the bird dataset.")
  ),
  list(
    dataset = "avonet", response = "Tarsus.Length", predictor = "Mass",
    rationale = paste(
      "Skeletal allometry: tarsus length scales with body mass. A second,",
      "largely independent morphological slope from the same predictor as",
      "the wing-length pair.")
  ),
  list(
    dataset = "avonet", response = "Beak.Length_Culmen", predictor = "Mass",
    rationale = paste(
      "Beak allometry: culmen length scales with body mass, but more",
      "weakly and more noisily than wing or tarsus (beak shape is under",
      "strong ecological/dietary selection independent of size) -- a",
      "harder test case for the same dataset.")
  ),

  # -- FishBase (Length, Weight, DepthRangeDeep, Troph all raw numeric) --
  list(
    dataset = "fishbase", response = "Weight", predictor = "Length",
    rationale = paste(
      "Length-weight allometry: the standard fisheries cube-law scaling",
      "relationship, strongly phylogenetically structured -- the anchor",
      "case for the fish dataset.")
  ),
  list(
    dataset = "fishbase", response = "DepthRangeDeep", predictor = "Length",
    rationale = paste(
      "Depth-size relationship: larger-bodied fish tend to reach greater",
      "maximum depths. Ecologically motivated but noisier and more",
      "heterogeneous across clades than length-weight -- a harder case.")
  ),
  list(
    dataset = "fishbase", response = "Troph", predictor = "Length",
    rationale = paste(
      "Trophic-level allometry: larger fish tend to occupy higher trophic",
      "positions. A weak, well-studied ecological slope near the",
      "detection limit -- a stress test for whether MI pooling preserves",
      "a small true effect rather than just not being obviously wrong on",
      "a strong one.")
  )
)

# Sanity: every pair's dataset is one of the three planned datasets, and no
# pair names the same response twice within a dataset.
stopifnot(
  all(vapply(mi_realdata_pairs, function(p) p$dataset, character(1)) %in%
        c("pantheria", "avonet", "fishbase")),
  !anyDuplicated(vapply(mi_realdata_pairs, function(p) paste(p$dataset, p$response), character(1)))
)
