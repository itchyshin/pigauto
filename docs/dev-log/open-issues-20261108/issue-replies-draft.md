# Draft issue replies

Drafts only. None has been posted. The accompanying PR adds regression coverage; the replies distinguish that work from changes already on main. Link that PR when posting.

## #254: README and documentation

The changes on main now include the function map, the distinction between `impute()` and `multi_impute()`, explicit ordinal levels and the revised report explanations. The focused ordinal and report checks passed again in this follow-up. Codecov also has a processed coverage report and a working percentage badge.

The follow-up PR adds regression coverage and an issue ledger; it does not introduce another imputation method. Your confirmation that the revised workflow is clear is still useful. The earlier installation probe did not identify what made your installation slow, so that part remains unexplained. Please keep any remaining documentation or installation concern in this thread.

## #257: visualisation documentation

The missingness views and worked visualisation guide are on main. The follow-up checks passed for trait grouping, repeated-observation ordering, single-trait plots and baseline-only report figures.

The phylogenetic view is still a matrix ordered by tree traversal, rather than a plot of the tree itself; trait types appear in labels rather than colours. If those distinctions leave part of your request unmet, please describe the view you need. We have not treated passing plot tests as confirmation that the visual design meets your needs.

## #258: data-only species warning

The warning fix is already on main. The follow-up PR adds a public-workflow regression for both in-memory traits and a CSV import: data-only species keep their rows unchanged, their missing values remain missing, and an in-tree control is filled. The test also checks that the misleading “-- dropped” wording is absent. Both routes passed.

This follow-up used a small continuous-trait fixture with the GNN off. It does not extend the claim to every input mode or trait type. Your original example remains a useful independent confirmation.

## #259: factor-covariate lambda

The rank-deficient covariate-design fix is already on main. The follow-up PR adds an end-to-end regression using your n = 200, seed = 1 example: factor and numeric treatment coding agree on lambda, completed values and prediction uncertainty. The regression passed; temporarily removing independent-column selection made three assertions fail, so the test detects that fault.

This verifies coding invariance for the supplied simulated example. It does not establish lambda recovery generally, and we have not rerun your amphibian data on current source. Please retain that distinction when checking the original application.

## #260: amphibian predictions and seed variation

The documentation on species-level aggregation, row-level covariates and interval widths is already on main. I also checked the public cached prediction attachment without refitting. After the same fold-5 deduplication, it contains 375 scored values from 77 species per seed. The median prediction range is 0.852 degrees C, and 30 values span more than 3 degrees C across seeds, matching your report.

That confirms the arithmetic in the attachment from the older source version. It does not show that seed sensitivity is resolved on current source. The post-hoc three-seed average also needs a separate held-out comparison before we recommend it. This issue should remain open for that validation, including a sensitivity analysis for the eastern newt's contribution.

## #250: mammal comparison and tied predictions

The posterior output fields and biological-bound limitations are documented on main. The repeated values on tree 1 remain unexplained. Identical covariance against observed species can produce identical univariate conditional means, but that is only a possible mechanism; it does not diagnose your joint-trait workflow.

Could you provide a synthetic or anonymised example retaining the ties, with the tree, exact settings/version, resolved per-trait baseline route and unrounded predictions? The manuscript data need not be made public. Those inputs would let us locate the ties before changing prediction or decoding. We have not added jitter or altered rounding to hide them.

## #251: narrower residual-correlation intervals

The warning against analysing a single point-filled dataset as complete observations is documented on main. The narrower residual-correlation intervals after the two-step analysis remain an open question; this follow-up does not attribute them to either finite imputation counts or model differences.

The next useful inputs are the exact model/prior specifications, scaling constants, retained imputation indices, per-imputation parameter draws and reference draws, or a small shareable synthetic counterpart. Repeated subsets of the existing fits could first assess finite-m variability without refitting. Increasing m can reduce Monte Carlo variation, but cannot establish compatibility between different imputation and analysis models. No interval adjustment is proposed yet.

## #261: CRAN evaluation API

The evaluation fix and unused-argument check are on main, and the affected evaluation tests passed in this follow-up. Getting Started and the mixed-types article retain the `evaluate(fit)` workaround for CRAN 0.11.0.

CRAN still listed 0.11.0 when checked on 8 November. This PR adds verification and tracking; it does not publish a corrected CRAN package. Please keep this issue open until the release containing the corrected API is available.

## #135: comparative evidence

This follow-up verifies specific engineering paths and reproduces the arithmetic in one public seed-prediction attachment. It does not replace the broader multi-seed, mixed-type and real-data comparisons requested here.

The earlier one-seed feasibility run remains timing and workflow evidence only. The wider comparison needs its own predeclared scenarios, paired evaluation cells, seed plan and compute budget. No additional campaign was launched, and no broader performance or calibration claim follows from these regression tests.
