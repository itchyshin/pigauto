# Runtime findings kept separate from coverage tests

## Uncertainty plot ranks

In R/plot_pigauto.R, plot_uncertainty() sorts predictions and interval bounds but retains their original indices as x. For predictions 3, 1, 2, the plot data has x 2, 3, 1. geom_line then connects the original order, contrary to its axis label. A regression test checks ranks, line values, interval bounds and aligned truth. It failed three assertions before the one-line repair. Repair work uses codex/uncertainty-plot-order in the older completed-issue worktree; no package runtime changes enter the coverage branch.

## Single candidate clade

In R/mask_missing.R, generate_missing_probs() calls sample(candidates, n_clades). With one numeric candidate, R samples from 1:candidate instead of sampling that vector. A three-tip tree has one non-root internal node, so the wrong node may be selected and the intended clade may receive no elevated missingness. A separate repair must sample positions and index candidates, with a deterministic three-tip regression and empty-candidate control.

## Empty forest panel

In R/plot.R, the empty-trait guard follows aggregate(), which errors on a zero-row input first. The guard is not safely reachable by passing an empty data frame through the existing implementation. Keep its uncovered lines unresolved. Do not remove or exclude them solely to raise coverage. Any change to the desired empty-input contract belongs in a separately reviewed repair.

These are pre-existing findings, not regressions introduced by PR 283. They do not change the CTmax study's findings or justify a new scientific recommendation.

## Missing uncertainty SE diagnostic

In R/plot_pigauto.R, the missing-SE guard follows matrix indexing. A NULL SE or absent trait column errors before reaching the intended message. The guard remains unresolved in the coverage inventory. This diagnostic repair is distinct from the sorted-rank fix in PR 284.

## Prepared repairs

PR 284 contains the sorted-rank fix and its regression. PR 285 contains the single-candidate clade fix and its regression, including the no-candidate control. Both remain drafts and require separate merge approval. The masking regression and existing masking suite passed 35 assertions with no warnings or skips.

## Method attribution in comparison legends

The default comparison plot labels improvement as GNN better/worse even though the pigauto method can use the default GNN-off fit. Coverage tests retain checks for the two actual method labels but do not endorse that component attribution. A separate repair should use the supplied method name in the improvement legend.

## Coverage measurement omission

The original coverage XML contains 48 R source classes, while the source tree contains 50 files. data.R contains six character-literal documentation anchors and no calls or function definitions. graph_transformer_block.R contains executable torch module methods, which covr's ordinary namespace traversal does not discover through the factory's module attribute. The supplemental collector exposes that actual R6 generator to covr's R6 traversal, runs the existing transformer tests, and adds the disjoint file to the denominator. The first probe covered 60 of 64 tracked lines, with four covariate-injection lines unexercised. This is measurement evidence, not a new model-validation claim.
