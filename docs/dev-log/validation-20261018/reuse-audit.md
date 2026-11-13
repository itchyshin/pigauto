# Prior evidence source comparison

Current base: b2e6054685e6bdedb073ba46e058af2c0a619a18

## Posterior sep

Prior source cd3a4cb. Changed paths:

```text
DESCRIPTION
NAMESPACE
R/build_phylo_graph.R
R/check_pigauto.R
R/cross_validate.R
R/data.R
R/evaluate.R
R/fit_baseline.R
R/fit_helpers.R
R/fit_pigauto.R
R/impute.R
R/joint_mvn_solver.R
R/joint_threshold_baseline.R
R/mi_posterior.R
R/multi_impute.R
R/multi_impute_trees.R
R/ovr_categorical.R
R/phylo_signal.R
R/plot.R
R/plot_pigauto.R
R/pool_mi.R
R/predict_pigauto.R
R/preprocess_traits.R
R/read_tree.R
R/report.R
R/result_accessors.R
R/save_load.R
R/simulate_benchmark.R
R/utils.R
script/mi_gls/01_cell_v2.R
```

Verdict: no blanket equivalence established. Retain historical evidence with its source identity; rerun affected current cases. Dependencies still need comparison and sentinel replay before any narrower reuse.

## GNN defaults

Prior source b565cad. Changed paths:

```text
DESCRIPTION
NAMESPACE
R/build_phylo_graph.R
R/check_pigauto.R
R/cross_validate.R
R/data.R
R/evaluate.R
R/fit_baseline.R
R/fit_helpers.R
R/fit_pigauto.R
R/impute.R
R/joint_mvn_solver.R
R/joint_threshold_baseline.R
R/mi_posterior.R
R/multi_impute.R
R/multi_impute_trees.R
R/ovr_categorical.R
R/phylo_signal.R
R/plot.R
R/plot_pigauto.R
R/pool_mi.R
R/predict_pigauto.R
R/preprocess_traits.R
R/read_tree.R
R/report.R
R/result_accessors.R
R/save_load.R
R/simulate_benchmark.R
R/utils.R
script/campaign_gnn_off_lib.R
```

Verdict: no blanket equivalence established. Retain historical evidence with its source identity; rerun affected current cases. Dependencies still need comparison and sentinel replay before any narrower reuse.
