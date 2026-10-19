# GNN incremental-value simulation: frozen ADEMP contract

Status: frozen for the bounded pilot after independent statistical and source review. Results are conditional on this simulation design and do not establish universal GNN superiority.

## Aim

Estimate the incremental predictive value of pigauto's calibrated neural contribution for repeated observations with observation-level covariates, compared with the same fitted model after removing only the neural contribution. The primary policy uses the package's current `phylo_signal_gate = FALSE` default. A paired `TRUE` policy-sensitivity arm is also run because the user requested both settings. Source inspection found that the option is false by default. The TRUE arm follows the ordinary fit path and then can route a weak-signal trait to the Brownian baseline corner, or the mean corner when the safety floor is enabled. The current source also used internal validation/test outcomes to estimate signal despite documenting a training-only calculation. This study includes a tested source fix that masks those split cells before estimating signal; the `TRUE` arm therefore measures the corrected policy implementation in this study branch. It does not by itself validate historical releases.

## Data-generating mechanisms

Each replicate independently generates a 500-tip birth-death tree with `ape::rtree()`. Let `C = ape::vcv(tree)` and `R = cov2cor(C)`, ordered by tree tip labels. The species effect is `u ~ N(0, variance_multiplier * R)`. To avoid an O(n^3) Cholesky for every replicate, draw Brownian `u_raw` with `ape::rTraitCont(tree, sigma = sqrt(variance_multiplier), root.value = 0)` and divide each tip value by `sqrt(diag(C))`; this produces covariance `variance_multiplier * R`.

For each species, draw `N_s = min(max(Poisson(5), 1), 15)` observations. For observation `i`, draw ten independent `Z_ij ~ N(0,1)`. Define `f` by scenario and standardize it within the generated dataset to mean zero and sample SD one:

- nonlinear: `sin(Z1) exp(0.3 Z2) + 0.5 tanh(Z3) + 0.3 sin(2 Z4) - 0.2 cos(Z5)`;
- interactive: `Z1 Z2 + 0.5 Z3^2 + 0.3 Z1 Z4 - 0.2 Z2 Z5`;
- linear recovery: `0.5 Z1 + 0.3 Z2 + 0.2 Z3 + 0.1 Z4 - 0.1 Z5 + 0.15 Z6 - 0.2 Z7 + 0.05 Z8 + 0.1 Z9 - 0.15 Z10`;
- null: `f = 0`.

Generate `Y_is = u_s + beta_strength * f_i + epsilon_is`, where `epsilon_is ~ N(0, 0.5^2)` independently. The `variance_multiplier` is a Brownian variance component, not Pagel's lambda. The signal function and coefficients above are locked for all runs.

## Estimands and methods

The outer mask is generated before fitting. Select exactly half the species for complete withholding. Among remaining species, mask a target of 20% of their rows, subject to a per-species cap of `N_s - 1`; this keeps every species in the partial-observation stratum represented by at least one observed row. Scored cells are the union of fully withheld and partly observed strata. Their outcomes are never supplied to preprocessing or fitting. The normal internal validation/test split applies only to remaining observed cells and is recorded separately.

For each dataset and gate policy, fit `impute(..., gnn = TRUE)` with the fixed seed and identical data, tree, outer mask, covariates, and other settings. The frozen fit ceiling is 2,000 epochs, with the recorded early-stopping rule; the pilot uses the same limit. The primary paired contrast is the held-out prediction from that fit minus prediction from the same fit with its GNN gate set to zero and the GNN gate weight transferred to the BM weight. The mean weight, stored mean, model state, and fitted linear covariate effect are preserved. This targets within-fit incremental value, not a comparison of independently selected models.

Secondary comparators are: observed-training mean; observation-level linear OLS; observation-level OLS with all ten covariates, their ten squares, and all 45 pairwise interactions; and phylogenetic random-intercept GLS with the same basis and tree correlation `R`. These are specified reference methods, not a performance ceiling. The GLS fit uses ML for variance components and conditional Gaussian predictions. A nonconverged GLS optimizer is a failed attempt.

## Scenario manifest and replication

The five scenario cells and allocation are recorded in `scenario_manifest.csv`. The nonlinear and interactive cells at variance multiplier 0.3 are confirmatory for the current default policy, `phylo_signal_gate = FALSE`. The `TRUE` setting is a policy sensitivity, not a second confirmatory family: it can force the GNN weight to zero by design. The low-phylogenetic-signal nonlinear, linear-recovery, and null cells are controls. Each cell has 500 independently generated datasets. Both gate settings use paired data and seeds, for 5,000 pigauto fits total. Seed stream `s` is `1 + ((20261019 + scenario_index * 1000003 + replicate * 7919 + s * 104729) %% 2147483646)`; stream 1 is DGP and stream 2 is fit. The pilot may explicitly override the 500-replicate production allocation; the production audit requires exactly 500 per scenario. No random fit seed is drawn from ambient process state.

## Outcomes and statistical analysis

For each dataset, gate policy, method, and missingness stratum, retain row-level predictions, truth, squared error, RMSE, and status. The inferential primary scale is the dataset-level paired MSE difference (full fit minus same-fit neural-off); the relative RMSE reduction is descriptive and its point estimate must be at least 2% to meet the practical-effect criterion. Only the default `FALSE` policy supports confirmatory claims. Use Holm-adjusted paired tests over the two confirmatory scenario contrasts and a simultaneous Bonferroni family-wise 95% t interval for their mean MSE differences. The `TRUE` policy is reported as a sensitivity, including trigger rates, and cannot independently establish confirmatory benefit. Declare benefit supported only when the adjusted test rejects in the improvement direction, the interval is below zero, and the point relative RMSE reduction is at least 2%. Declare failed when the interval rules out a benefit of that practical size; otherwise indeterminate. Report MCSE across independent datasets, with a cluster bootstrap over datasets as a sensitivity check. Do not select or headline whichever gate policy happens to perform better.

Report gate weights, estimated phylogenetic signal, trigger frequency, failure, nonconvergence, and paired results separately by fully withheld versus partly observed species. Failed or nonconverged attempts remain in the requested denominator. The primary paired estimate is complete-pair; report `n_paired`, requested count, and paired fraction. If fewer than 95% of requested datasets form a valid pair in a confirmatory scenario, classify that scenario as indeterminate regardless of the complete-pair estimate. The Gaussian outcomes have unbounded support, so no finite worst-case error bound is asserted for failed fits. Failures are never silently dropped.

## Harness controls and reproducibility

Before production, deterministic checks must show: the generator is reproducible; its species covariance matches `variance_multiplier * R`; the two outer-mask strata are disjoint and satisfy their row-retention definitions over many seeds; a named-tree GLS prediction matches a dense conditional-covariance calculation; a deliberately broken covariance control is detected; the pigauto ablation preserves mean weight and direct linear covariate contribution; predictions align to original input rows; and outer scoring outcomes are absent from all fit inputs. Every run record binds scenario, replicate, gate policy, seeds, source commit when available, source-bundle label, per-file SHA-256 source hashes, full R session/platform/external-software identity, settings, elapsed time, status, and failure reason. The TRUE arm requires `phytools`; if signal is unestimable for a replicate, that arm is explicitly marked as not exercised for that replicate.
