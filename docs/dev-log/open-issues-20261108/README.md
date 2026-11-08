# Remaining issue work, 8 November 2026

Base source: `30ab1e1b50b34a732485f6cccb195e9fa75299fc`. Branch: `codex/open-issues-nov08`.
The user approved verification, reproducible fixes and separate data, validation and release tracks. No issue closure, message, merge or release is included.

## Verified engineering and documentation

| Issue | Evidence and disposition |
|---|---|
| [#254](https://github.com/itchyshin/pigauto/issues/254) | README function map, impute versus multi_impute distinction and ordinal_levels are already on main. Existing ordinal and report tests pass. No further behavioural change justified by this pass. Reporter confirmation and the cause of slow installation remain open. |
| [#257](https://github.com/itchyshin/pigauto/issues/257) | Existing visualisation tests exercise latent-to-trait grouping, repeated-observation tree ordering, single-trait plots and offline baseline report figures. A tree-ordered matrix remains distinct from a tree plot; agreement on the desired visual is still needed. |
| [#258](https://github.com/itchyshin/pigauto/issues/258) | New public-path regression checks retained data-only rows and a finite in-tree imputation, both in memory and after CSV import. It also rejects the old "-- dropped" warning. Passed. |
| [#259](https://github.com/itchyshin/pigauto/issues/259) | New public-path regression uses the reporter's n=200, seed=1 DGP. Factor and treatment-coded covariates agree on lambda, completed values and uncertainty. This is coding invariance in one regime, not recovery validation or a rerun of the amphibian data. |
| [#261](https://github.com/itchyshin/pigauto/issues/261) | GitHub has evaluation mode and rejects unused arguments; both vignette sources document the CRAN workaround evaluate(fit). CRAN's package page still lists 0.11.0. Keep the release task separate. |

The already merged fixes were retained. New work adds two end-to-end regression tests and a reproducible descriptive audit; it changes no package runtime behaviour.

## Data-dependent investigations

### #250: repeated point predictions on tree 1

The original manuscript data were not attached. Existing discussion demonstrates that missing sisters with identical covariance against observed species can have identical univariate conditional means; that is a possible mechanism, not a diagnosis. Joint prediction can also use the other trait.

Needed next: a synthetic or anonymised input retaining the ties, the corresponding tree, exact settings/version and the resolved per-trait baseline route. Include the unrounded point predictions before counting values rounded to eight decimals, and compare the same missing cells on all trees. Inspect whether ties arise in the production baseline, decoding or completed-data assembly. Do not add jitter or change rounding to conceal the pattern.

Posterior fields and biological-bound limitations are already documented. A public article needs shareable inputs and an explicit distinction between model agreement and accuracy against hidden truth.

### #251: residual-correlation intervals narrower after analysis

Needed next: analysis and imputation model specifications with exact priors, scaling constants, retained imputation indices, parameter draws per completed-data fit, and the in-model reference draws. A shareable small synthetic example can substitute for the manuscript inputs.

First reuse the stored fits to estimate variability across repeated subsets of m imputations, without refitting. Keep trees separate and retain complete-data uncertainty. Then distinguish finite-m variation from differing priors/models. Increasing m addresses Monte Carlo variability; it does not establish compatibility of the models. No cause is identified and no interval correction is proposed here. Private inputs must not enter a public issue or repository.

## Validation track

### #260: seed sensitivity

The public attachment is [gnn_predictions.csv](https://github.com/user-attachments/files/34432606/gnn_predictions.csv), supplied in [issue #260](https://github.com/itchyshin/pigauto/issues/260). It represents the reporter's old-source fits, not new fits on this branch.

Run:

```sh
python3 docs/dev-log/open-issues-20261108/summarize_seed_predictions.py /path/to/gnn_predictions.csv /tmp/seed-summary.json
```

The script checks matched truth/species/folds across seeds, finite predictions, valid intervals and duplicates, then follows the issue's fold-5 deduplication rule. The retained summary records the attachment SHA256. No raw prediction rows are committed.

The attachment has 1,281 rows before deduplication and 375 distinct scored observations across 77 species per seed. Median prediction range is 0.85234 degrees C; 30 observations span more than 3 degrees C. Per-seed RMSE is 1.78869, 1.74793 and 1.69146 degrees C; interval counts are 354, 352 and 359 of 375. These are descriptive summaries of cached predictions, with clustered observations and only three seeds. No independent inferential test or coverage guarantee follows.

A remedy needs a predeclared comparison on fresh held-out data, with fixed seed sets, paired cells, species-level summaries and a sensitivity analysis excluding the eastern newt. Compare any proposed ensemble using new validation data; the post-hoc three-seed average is insufficient to recommend a new default. No GNN refits or campaign launched.

### #135: broader comparative evidence

The one-seed feasibility run already discussed in the issue is timing/workflow evidence only. Multi-seed, mixed-type and real-data comparisons remain a separate compute campaign. Reuse the existing validation inventory and route serious runs to Totoro or DRAC. This pass neither launches that campaign nor substitutes the small public-path tests for it.

## Release track

[CRAN](https://cran.r-project.org/web/packages/pigauto/index.html) lists version 0.11.0 on this audit date. #261 remains open until a corrected package is available there. No release checks, submission or release settings were changed. Codecov is already working and requires no token change.

## Acceptance and remaining scope

Completed: verify engineering fixes, add regression coverage for public paths, reproduce descriptive seed arithmetic, state data requirements and retain separate scientific/release work. Pending: reporter confirmation, causal diagnosis of mammal discrepancies, current-source seed evaluation, broader method validation and CRAN publication. Issue states are unchanged.
