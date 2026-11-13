# Residual-correlation interval subset diagnostic

This script is ready for shareable stored draws from issue #251. It fits no model and cannot diagnose the reported discrepancy without those inputs.

Supply an RDS containing a uniquely named list of trees. Each tree entry contains:

- `completed_draws`: a uniquely named list, one numeric residual-correlation draw vector per completed-data fit;
- `reference_draws`: numeric residual-correlation draws from the corresponding in-model fit on that same tree.

Use the raw correlation scale, within [-1, 1]. Completed-data draw vectors must have equal lengths so stacking gives each imputation equal weight. Do not supply summaries, Fisher-z values, phylogenetic correlations in place of residual correlations, or fits from different trees under one tree entry. Retain model code, priors, scaling constants and convergence diagnostics alongside the inputs. Private manuscript material should remain outside Git and public issue attachments.

```r
source("script/issue_251_diagnostics/subset_intervals.R")
x <- readRDS("shareable_rho_E_draws.rds")
out <- subset_intervals(x, m_grid = c(5L, 10L), repetitions = 200L)
write.csv(out, "subset_widths.csv", row.names = FALSE)
```

The command-line form uses the default m grid 5, 10 and 20. Request only counts available for every tree; unavailable levels cause an error rather than being silently omitted. The function samples subsets without replacement within each repetition, retaining selected imputation IDs. Across repetitions, the same subset can recur. When m equals every available fit, it evaluates that single complete set once.

The reported ratio is the equal-weight stacked-draw 95% interval width divided by the corresponding reference width, using R's type-7 quantiles. It also records the pooled mean and interval endpoints. This is descriptive sensitivity to selecting from the stored fits. The repetitions are not independent datasets or independent MCMC runs. The distribution of ratios is not a confidence interval for the true ratio and does not establish frequentist coverage. The diagnostic cannot extrapolate to more imputations than were fitted or establish compatibility between different model specifications or priors.

Synthetic controls exercised a valid complete grid and rejected invalid correlations, a zero-width reference interval, unavailable m levels, nonfinite and missing m values. No manuscript data or actual residual-correlation draws were available for this pass.
