# Runtime findings kept separate from coverage tests

## Uncertainty plot ranks

In R/plot_pigauto.R, plot_uncertainty() sorts predictions and interval bounds but retains their original indices as x. For predictions 3, 1, 2, the plot data has x 2, 3, 1. geom_line then connects the original order, contrary to its axis label. A regression test checks ranks, line values, interval bounds and aligned truth. It failed three assertions before the one-line repair. Repair work uses codex/uncertainty-plot-order in the older completed-issue worktree; no package runtime changes enter the coverage branch.

## Single candidate clade

In R/mask_missing.R, generate_missing_probs() calls sample(candidates, n_clades). With one numeric candidate, R samples from 1:candidate instead of sampling that vector. A three-tip tree has one non-root internal node, so the wrong node may be selected and the intended clade may receive no elevated missingness. A separate repair must sample positions and index candidates, with a deterministic three-tip regression and empty-candidate control.

## Empty forest panel

The initial empty-data-frame probe errored in aggregate() before the guard. Independent review later found a real supported input: with na.action set to "na.pass", finite metrics whose trait labels are all missing aggregate successfully into zero groups. The public plot_comparison() call then returns invisible NULL through the guard. Preserve this guard. The dedicated real-input test is committed in coverage head 0942bf9c and included in the frozen 3015dbe1 profile now running; its line hits have not yet been measured.

These are pre-existing findings, not regressions introduced by PR 283. They do not change the CTmax study's findings or justify a new scientific recommendation.

## Missing uncertainty SE diagnostic

In R/plot_pigauto.R, the missing-SE guard follows matrix indexing. A NULL SE or absent trait column errors before reaching the intended message. The guard remains unresolved in the coverage inventory. This diagnostic repair is distinct from the sorted-rank fix in PR 284.

## Prepared repairs

PR 284 contains the sorted-rank fix and its regression. PR 285 contains the single-candidate clade fix and its regression, including the no-candidate control. Both remain drafts and require separate merge approval. The masking regression and existing masking suite passed 35 assertions with no warnings or skips.

## Method attribution in comparison legends

The default comparison plot labels improvement as GNN better/worse even though the pigauto method can use the default GNN-off fit. Coverage tests retain checks for the two actual method labels but do not endorse that component attribution. A separate repair should use the supplied method name in the improvement legend.

## Coverage measurement omission

The original coverage XML contains 48 R source classes, while the source tree contains 50 files. data.R contains six character-literal documentation anchors and no calls or function definitions. graph_transformer_block.R contains executable torch module methods, which covr's ordinary namespace traversal does not discover through the factory's module attribute. The supplemental collector exposes that actual R6 generator to covr's R6 traversal, runs the existing transformer tests, and adds the disjoint file to the denominator. The first probe covered 60 of 64 tracked lines, with four covariate-injection lines unexercised. This is measurement evidence, not a new model-validation claim.

## Proof-backed obsolete branches

Draft PR 294 removes only the invalid Bayesian weight-sum fallback in R/pagel_lambda.R. Once any finite likelihood exists, subtracting the finite maximum gives at least one exact zero log weight and hence a unit weight; every other exponential weight lies in [0, 1]. Their sum must be finite and positive for any addressable R vector. Sparse-data and empty-grid fallbacks remain. Independent review checked the numerical argument and preservation tests. The original removed location remains unresolved pending approval.

Draft PR 295 contains a second separate cleanup. In the GNN-off fit path, Mondrian requests are rejected before metadata extraction, so its Mondrian attribute arm cannot run. In pool_mi(), fewer than two retained fits are rejected before Rubin arithmetic, so the later one-fit between-imputation fallback cannot run. Independent review found no supported-input counterexample. All 18 real supported-method and two-fit arithmetic assertions pass before removal. The related GNN-off and pooling suites pass 182 assertions afterward, with 17 existing small-calibration warnings and no skips or failures. These are control-flow cleanups, not fixes to scientific validity, and their removed locations will remain unresolved until approved.
