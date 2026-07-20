# pigauto capability surface

**Current development version:** 0.10.0.9000  
**Lifecycle:** experimental  
**Scope:** phylogenetic trait imputation before a user-chosen downstream analysis. This is a capability inventory, not a release certificate or a claim that every enabled component improves every dataset.

## What the package is for

`pigauto` completes missing cells in a species-by-trait matrix using a phylogenetic baseline and a gated graph-neural correction. The calibrated gate may close (`r_cal = 0`), leaving the baseline as the valid prediction. It supports a mixed latent representation rather than a separate package per trait type, and exposes point predictions, uncertainty outputs, multiple imputation, and tree-aware imputation workflows.

## Capability matrix

| Capability | Status | User-facing entry point | Evidence / boundary |
|---|---|---|---|
| Mixed-type trait preprocessing | covered | `preprocess_traits()`, `impute()` | One pipeline handles continuous, binary, categorical, ordinal, count, proportion, zero-inflated count, and multi-proportion inputs; multi-proportion groups require explicit declaration. |
| Phylogenetic baseline with graceful fallback | covered | `fit_baseline()`, `impute()` | Per-column BM / label propagation remains available; joint MVN, threshold-joint, and OVR paths activate only when their data and optional dependencies permit it. |
| Gate-safe graph correction | covered | `fit_pigauto()`, `impute()` | Predictions blend baseline and GNN with a per-latent-column calibrated gate; `r_cal = 0` is a valid fallback and the safety floor protects against degradation on validation data. |
| Optional covariates and multi-observation refinement | covered | `impute(..., covariates =, species_col =)` | Multi-observation inputs aggregate at species level for graph/baseline work, then may use observation-level covariate refinement. Discrete within-species aggregation remains lossy. |
| Point predictions and type-aware decoding | covered | `predict()`, `impute()` | Outputs are decoded to the original trait scale; categorical probabilities and multi-proportion simplex constraints are retained. |
| Conformal prediction intervals | partial | `predict.pigauto_fit()` | Continuous-family and explicit proportion traits receive conformal intervals when calibration scores are available. Discrete `se` values are uncertainty scores, not Gaussian standard errors; `conformal_split_val = TRUE` is the more defensible coverage option when the calibration split is small. |
| Multiple imputation and Rubin pooling | covered | `multi_impute()`, `with_imputations()`, `pool_mi()` | Supports conformal (default) and MC-dropout draws plus downstream pooling. Discrete draws are categorical/Bernoulli rather than Normal draws. |
| Posterior-tree uncertainty propagation | partial | `multi_impute_trees()` | Produces tree-indexed completed datasets for downstream refitting and Rubin pooling. `share_gnn = TRUE` is a speed-oriented approximation; use `FALSE` when an exact per-tree GNN refit is required. |
| Active-imputation guidance | covered | `suggest_next_observation()` | Ranks currently missing cells or species by expected variance/entropy reduction across all eight trait types; it guides data collection but does not itself collect observations. |
| Cross-trait attention | experimental | `fit_pigauto(..., use_trait_attention = TRUE)` | Opt-in and off by default. The 60-replicate BIEN-scale ablation did not show benefit where the joint-MVN baseline already captured cross-trait covariance. |
| BACE comparison baseline | partial | `fit_baseline_bace()` | Optional comparison/wrapper path only; BACE remains a separate package and is not modified or bundled as the main pigauto engine. |
| Benchmark-derived performance claims | experimental | benchmark scripts and pkgdown articles | Claims must name DGP, sample size, missingness, trait type, seeds, and comparator. No generic claim that the GNN beats the baseline or other methods is warranted. |

## Important limits

- pigauto is a preprocessing / multiple-imputation workflow, not an in-model FIML engine and not a downstream comparative-model fitter.
- A fitted gate or a successful run is not evidence of improved accuracy. The baseline may be the right answer, especially in high-phylogenetic-signal or small-validation regimes.
- Tree uncertainty is only fully propagated when the downstream model is refit with the matching tree for every completed dataset; pigauto supplies the tree-indexed imputations, not a universal downstream model.
- The package is R-only. There is no Julia twin or R-to-Julia parity claim.

## Evidence anchors

- Interface and trait-type contract: `README.md`, `R/impute.R`, `R/fit_baseline.R`, `R/predict_pigauto.R`.
- Regression coverage: `tests/testthat/test-mixed-types.R`, `test-fit-predict.R`, `test-multi-impute.R`, `test-multiobs-levelc.R`, `test-ovr-categorical.R`.
- Current caveats and evidence boundaries: `README.md` (recommendations and known limitations) and `NEWS.md` (prediction-path, tree-MI, conformal, and trait-attention changes).
