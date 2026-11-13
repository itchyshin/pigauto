# Point-prediction tie diagnostic

This no-refit helper is for the single-observation numeric trait in issue #250. Save a local bundle containing the original `traits`, matching `tree`, saved `result` and `trait = "log_ls"`, then run:

```sh
Rscript script/issue_250_diagnostics/diagnose_ties.R bundle.rds report.rds
```

The report counts exact ties and ties after rounding to eight decimal places, and checks up to 20 exact-tied pairs. For each pair it reports the maximum difference between their phylogenetic correlation rows against species with observed target values. It also returns the saved production baseline route and lambda slots when available. It does not fit a model or alter any predictions. The covariance calculation is dense and uses memory proportional to the square of the number of tree tips.

The supplied tree must be identical to the tree saved in `result$tree`, and observed trait values must agree with the saved completion. This prevents a different tree or dataset from producing a misleading comparison under the same species labels. Use the saved tree rather than reconstructing an equivalent tree with different attributes.

Identical covariance rows can produce identical conditional means under a univariate Brownian baseline. They do not establish why a joint-trait prediction tied: another observed trait can distinguish the species. A zero difference is therefore a lead for the selected route, not proof of a defect or its cause. Nonzero differences likewise do not locate a decoding or numerical problem by themselves.

The helper emits input row indices rather than species names in its pair table. The saved route/lambda objects and any other fields should still be reviewed before sharing. Keep the original manuscript bundle outside Git and public issue attachments. A synthetic or anonymised reproducer retaining the pattern is still needed for independent investigation.
