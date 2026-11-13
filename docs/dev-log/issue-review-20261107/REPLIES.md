# Posted issue replies

All ten comments were read back. The changes remain proposed in draft PR #262.

## #261

Thanks, Ayumi. You are right: CRAN 0.11.0 does not support `predict(fit, evaluation = TRUE)`, and matching version numbers hid the difference from GitHub.

[PR #262](https://github.com/itchyshin/pigauto/pull/262) adds the CRAN workaround, `evaluate(fit)`, to Getting Started Step 7 and the mixed-types article, and explains its long-table output. It also distinguishes the development build as 0.11.0.9000, rejects unused prediction arguments so typos cannot silently change an evaluation, and notes that training-history plots require the GNN. The existing GitHub evaluation fix remains in place.

These changes do not update CRAN: a patch release is still needed to bring the released API into line. Does the proposed documentation address the immediate problem? Please close the issue only if you consider it addressed; otherwise keep it open, including if you want it to track the CRAN release.

The changes are in a draft PR and are not yet merged or released on CRAN.

## #259

Thanks, Patrice. I reproduced the factor-versus-numeric discrepancy. [PR #262](https://github.com/itchyshin/pigauto/pull/262) selects independent columns from the observed covariate design for both lambda estimation and prediction, while retaining the full one-hot encoding for the GNN. The regression also checks redundant numeric columns. On your n = 200, seed = 1 example, the factor and treatment-coded covariates now both give lambda = 0.5015374; their completed values differ by less than 3e-10.

The patch also distinguishes an unavailable lambda estimate from the effective fallback value used for prediction. If a covariate effect cannot be identified because a factor level occurs only among missing outcomes, it warns about that limitation rather than promising coding-invariant predictions. Boundary estimates can also be legitimate; this patch does not warn for every fitted boundary value.

Could you check the patch against your amphibian workflow? Does it address the issue? If you are satisfied, please close it; otherwise please describe what still differs.

The changes are in a draft PR and are not yet merged or released on CRAN.

## #258

Thanks, Iwo. Your distinction is right: those rows are excluded from modeling but retained unchanged in `completed`. [PR #262](https://github.com/itchyshin/pigauto/pull/262) changes the shared warning to say that explicitly, so the input check and imputation no longer describe output rows as dropped. A smoke check confirmed that five data-only rows remain unchanged.

CSV files do not preserve factor classes or their scientific order. The patch adds `read_traits(ordinal_levels = list(Migration = c(...)))` so the order can be supplied explicitly; the example must use the order appropriate to the data. `check_pigauto()` continues to show the inferred trait types. The imbalance and small-validation warnings describe the actual fit and remain visible.

Does this address the warning and input-declaration concerns? If satisfied, please close the issue; otherwise please tell me what remains confusing.

The changes are in a draft PR and are not yet merged or released on CRAN.

## #257

Thanks, Szymon. [PR #262](https://github.com/itchyshin/pigauto/pull/262) adds a worked visualisation article and `plot_missingness()` with four views: per-trait missing fractions, an observation-by-trait matrix, a matrix ordered by phylogenetic traversal, and pairwise co-missing fractions. Trait labels include their types, expanded categorical traits count once, and repeated observations stay together in the tree-ordered view. The co-missing map includes numerical fractions.

The guide also demonstrates conformal-score plots and held-out evaluation. It explains that a missingness pattern can motivate investigation but cannot establish the missingness mechanism. Baseline-only HTML reports now show held-out metrics and conformal scores as offline SVG figures. GNN report charts still use the existing online chart library.

Do these views address the first set of missingness summaries you had in mind? If satisfied, please close the issue; otherwise please identify the view that still needs work.

The changes are in a draft PR and are not yet merged or released on CRAN.

## #254

Thanks, Dan. [PR #262](https://github.com/itchyshin/pigauto/pull/262) includes your accessor and bundled-data quick-start contributions from #252 and #253, preserving your commits and attribution. It also addresses the remaining review points:

1. The opening and a short walkthrough explain latent encoding, phylogenetic and cross-trait information, optional neural predictions, held-out calibration and decoding. Detailed output caveats come later.
2. `species_col` explicitly names a column in the traits frame. Repeated-observation instructions are shorter and explain species-level aggregation and the role of observation-level covariates with the GNN.
3. Baseline-only reports show held-out metric and conformal-score figures offline. GNN report charts retain their current internet dependency.
4. `read_traits(ordinal_levels = list(...))` accepts an explicit level order. For an existing data frame, use `ordered(..., levels = ...)`; a type override alone cannot infer the ecological order.
5. `?impute` explains that `n_imputations` controls diagnostic draws while the result still has one completed table. The MI documentation distinguishes diagnostic draws from software-eligible inferential routes and retains the validation limits.

Does this cover your checklist? If satisfied, please close the issue; otherwise please point to any item that remains incomplete.

The changes are in a draft PR and are not yet merged or released on CRAN.

## #135

The evidence has advanced since the July update. The [public validation article](https://itchyshin.github.io/pigauto/articles/validation-results.html#when-does-the-neural-component-add-predictive-value) now reports the paired neural-on/off results, gate-policy sensitivity, comparator results and reproducibility records.

The two confirmatory scenarios each used 500 independent datasets with 500 species, repeated continuous observations and ten covariates. The neural contribution reduced RMSE by 17.3% and 24.9% relative to the same fit with that contribution removed. These are within-fit gains, not evidence of universal superiority: phylogenetic GLS was better in the interactive scenario, and the no-covariate-signal control showed a small neural loss.

The package defaults to `gnn = FALSE`, so the baseline workflow requires no neural training. [PR #262](https://github.com/itchyshin/pigauto/pull/262) adds missingness guidance and useful baseline report figures. The wider external-comparator campaign across trait types and missingness settings remains incomplete, as does a matched test isolating graph message passing from the observation-level network. This issue should remain open for that remaining evidence rather than be closed on the narrower results.

## #250

Thank you for sharing this comparison. The reported agreement is encouraging for this two-trait mammal dataset and three trees; it does not establish held-out accuracy or validate the GNN. A worked comparison with in-model imputation could be useful, with assumptions and the two-step analysis clearly distinguished. A public example would need shareable inputs; AVONET may be an option.

[PR #262](https://github.com/itchyshin/pigauto/pull/262) proposes expanded help to describe `posterior$params`, `posterior$cell_interval`, `pooled_point` and `se`. For transformed traits, the point completion decodes the mean latent draw; the intervals and SD use decoded retained draws. A Gaussian model does not enforce a biological lower bound automatically.

The tree-1 ties remain unexplained. Could you share a minimal synthetic or anonymised example that retains the pattern, with its tree, settings and package version? The manuscript data need not be shared. Does this address the documentation and article questions? Please keep the tie investigation open until a reproducer is available.

The changes are in a draft PR and are not yet merged or released on CRAN.

## #251

[PR #262](https://github.com/itchyshin/pigauto/pull/262) proposes the documentation changes below.

Thank you for testing the two-step route. In your reported mammal analysis, posterior pooled estimates were close to the in-model estimates, while residual-correlation intervals were narrower; the reason remains uncertain. The single point fill changed several estimates in that example. This is a useful warning that treating one completion as observed ignores imputation uncertainty and can affect estimates as well as intervals. It does not establish the size or direction of bias in other analyses.

The narrower residual-correlation intervals cannot yet be attributed to finite Monte Carlo sampling or differences in priors. A useful next comparison would increase the number of imputations and repeat the pooling with independent runs, examining pooled total variance and Monte Carlo uncertainty while keeping the priors and model specifications consistent. This can help distinguish sampling variation from a persistent difference; it would not by itself prove a cause.

[PR #262](https://github.com/itchyshin/pigauto/pull/262) clarifies the single-fill distinction in the multiple-imputation article. The posterior route has eligibility conditions, and the registered validation has not established general inferential reliability. Does this address your questions and documentation suggestion? If you are satisfied, please close the issue; otherwise the interval-width question remains open.

The changes are in a draft PR and are not yet merged or released on CRAN.

## #256

[PR #262](https://github.com/itchyshin/pigauto/pull/262) expands the documentation on these questions.

Both scores are useful, but answer different questions. Comparing `pooled_point` with held-out truth assesses point prediction; it is not cheating. A single completed dataset includes draw variation, so its error also reflects that variation. To assess uncertainty, compare the predictive distribution with held-out values using interval coverage and width or a proper distributional score. To assess inference, evaluate the target estimate and its interval over replicated datasets.

Please check `mi$draws_method`: `auto` selects posterior only when eligible and otherwise falls back to conformal diagnostic draws. Conformal and MC-dropout draws are not supported for downstream Rubin pooling. The eligible posterior route can be passed through the software pooling workflow, but its inferential validity remains unestablished in the registered tests. Does this answer your question? If you are satisfied, please close the issue.

The changes are in a draft PR and are not yet merged or released on CRAN.

## #260

[PR #262](https://github.com/itchyshin/pigauto/pull/262) proposes the documentation changes below.

Thank you for the careful cross-validation. The reported GNN improvement and seed sensitivity are specific to these 375 held-out CTmax values and three seeds; the eastern newt contributed substantially to the gain. The seed average was explored after seeing the spread, so it needs a separate evaluation before recommending it as an improvement.

[PR #262](https://github.com/itchyshin/pigauto/pull/262) clarifies that the multi-observation baseline aggregates to species level and does not model within-species covariate differences when the GNN is off. Covariate use also depends on the selected baseline route. For untransformed continuous traits, a global conformal score gives a common width within a fit; other scoring methods or transformations can differ. The coverage in this dataset does not establish species-specific coverage. The factor-covariate lambda bug is addressed by the proposed fix discussed in #259.

Does this address your documentation question? If you are satisfied, please close this issue, with #259 retaining the separate bug discussion.

The changes are in a draft PR and are not yet merged or released on CRAN.
