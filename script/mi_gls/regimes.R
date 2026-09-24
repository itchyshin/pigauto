# script/mi_gls/regimes.R
#
# Regime grid for the MI-GLS-attenuation campaign (arc/mi-gls-attenuation).
# Background: multi_impute(draws_method = "conformal") halves the pooled
# PGLS slope of y ~ x under bivariate BM (see useful/mondrian-mi-se-
# justification.md "## Result" and
# script/mondrian_confirmation/13_mi_gls_attenuation_diag.R). This grid
# drives script/mi_gls/01_cell.R across the regimes where that attenuation
# should (or should not) show up.
#
# rho (trait correlation between x and y) is fixed at 0.7 throughout --
# it is not part of the grid.
#
# `lambda` is a Pagel lambda applied to the tree FOR SIMULATION ONLY (via
# the internal transform_tree_pagel(), R/pagel_lambda.R): traits are
# generated on the lambda-transformed tree, but every downstream step
# (imputation and both lm/gls fits) uses the ORIGINAL (lambda = 1) tree.
# At lambda = 0.5 this deliberately misspecifies the phylogenetic model
# used for imputation/inference relative to the true DGP.
#
# `missing`:
#   - "x_only": 30% of x cells missing (drawn per the `mechanism` below);
#     y always fully observed.
#   - "both": 30% of x cells AND (independently) 30% of y cells missing,
#     each drawn separately per `mechanism`.
#
# `mechanism`:
#   - "MCAR": each candidate cell missing independently with p = 0.30.
#   - "MAR_phylo": missingness concentrated in 2 randomly chosen clades
#     (15-35% of tips each), as in script/mondrian_confirmation/06_mi_se_sim.R
#     (itself following ~/pigauto_regime_map/mech_cell.R) -- clade tips get
#     7x the baseline miss probability, rescaled so the expected count is
#     30% of n. "both" draws a FRESH pair of clades independently for y.

regimes <- expand.grid(
  lambda    = c(1, 0.5),
  n         = c(300L, 1000L),
  mechanism = c("MCAR", "MAR_phylo"),
  missing   = c("x_only", "both"),
  KEEP.OUT.ATTRS  = FALSE,
  stringsAsFactors = FALSE
)
regimes$rho <- 0.7
regimes$regime_id <- seq_len(nrow(regimes))
regimes <- regimes[, c("regime_id", "lambda", "n", "mechanism", "missing", "rho")]

# Oracle proper-MI is only run when lambda == 1 (correctly-specified model)
# AND missing == "x_only" (the oracle formula in 01_cell.R conditions on
# ALL of y being observed; see 13_mi_gls_attenuation_diag.R's derivation).
regimes$run_oracle <- regimes$lambda == 1 & regimes$missing == "x_only"

if (sys.nframe() == 0L) {
  print(regimes)
}
