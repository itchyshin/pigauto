# Goal: test pigauto's neural contribution

Determine whether the current gated neural correction reduces paired held-out observation-level RMSE beyond the same fitted model with the neural gate set to zero, while preserving the fitted mean-baseline mixture and linear covariate contribution. Compare separately with an observation-level phylogenetic random-intercept GLS using the same covariate basis.

Deliver a reproducible current-source simulation, retained replicate-level results and provenance, a qualified scientific report, and a concise pkgdown results summary. No public API, method, default, deployment, or publication changes are in scope.

Success means the planned scenarios and attempts are fully accounted for, the harness and ablation are validated, summaries reproduce, and limitations are clear. A positive predictive-value claim requires at least 2% relative RMSE reduction and a Holm-adjusted paired interval below zero in a confirmatory scenario. Otherwise report failed or indeterminate evidence honestly.
