# Methodology Claim Audit - 2026-05-11

Scope: active pkgdown Methodology pages only, as listed in `_pkgdown.yml`.
No benchmark drivers were rerun. HTML pages were regenerated only from existing
`script/*.rds` artifacts after generator wording changes.

## Active Pages Audited

The active Methodology menu currently exposes seven per-trait pages and six
cross-dataset / simulation pages (`_pkgdown.yml:69-94`):

- Per-trait: continuous, binary, ordinal, count, categorical, proportion,
  zero-inflated counts.
- Cross-dataset / simulation: AVONET missingness, missingness mechanisms,
  tree uncertainty, Delhey, multi-observation, covariate effectiveness.

## Fixes Applied

| Page | Problem | Fix |
|---|---|---|
| `bench_continuous` | Claimed the GNN broadly "adds value" under all non-BM models and implied preserved ordering. Existing results include mixed deltas. | Reworded verdict and interpretation to report scenario-specific deltas and state that BM-vs-pigauto ordering is scenario- and trait-dependent (`pkgdown/assets/dev/bench_continuous.html:44-45`, `pkgdown/assets/dev/bench_continuous.html:241-242`). |
| `bench_binary` | Claimed pigauto "matches or slightly improves" across signal levels, while table rows include trailing label propagation. | Reworded to say pigauto stays close but can trail LP, and that GNN contribution is not uniformly positive in this run (`pkgdown/assets/dev/bench_binary.html:235`). |
| `bench_ordinal` | Claimed more ordinal levels give the GNN room to improve. Current rendered page cannot support a general improvement claim. | Reworded to "finer-grained variation" and explicitly avoids a general improvement claim (`pkgdown/assets/dev/bench_ordinal.html:201`). |
| `bench_count` | Dense-count language implied GNN room to improve even when dense-count RMSE moved from 0.484 to 0.485. | Reworded dense-count and overdispersion interpretation as scenario-dependent, with ordering read from the table (`pkgdown/assets/dev/bench_count.html:44`, `pkgdown/assets/dev/bench_count.html:258`). |
| `bench_categorical` | Rendered HTML was stale relative to the current RDS / markdown summary and claimed pigauto matched K = 3 LP accuracy. | Regenerated from existing RDS and reworded interpretation: K = 3 now shows LP 0.746 vs pigauto 0.690, and prose says pigauto can trail LP (`pkgdown/assets/dev/bench_categorical.html:43-55`, `pkgdown/assets/dev/bench_categorical.html:167`). |
| `bench_proportion` | Claimed pigauto matches or improves across all signal levels and no degradation. | Reworded to "usually close" with modest average gains and possible trait-level ties or slight trails (`pkgdown/assets/dev/bench_proportion.html:226`). |
| `bench_missingness_mechanism` | Claimed the GNN maintains relative advantage across mechanisms and ranked mechanisms too broadly. | Reworded MCAR as the reference setting and pigauto as mostly tracking the phylogenetic baseline in this run (`pkgdown/assets/dev/bench_missingness_mechanism.html:44`, `pkgdown/assets/dev/bench_missingness_mechanism.html:228`). |
| `bench_multi_obs` | Claimed covariates outperform and do not degrade too absolutely. | Reworded to "usually improves" in this simulation and tells readers to use the beta = 0 rows as evidence (`pkgdown/assets/dev/bench_multi_obs.html:347-355`). |
| `bench_covariate_sim` | Marked all covariate rows visually as "best" and said gating prevented noise injection. | Removed unconditional "best" styling and reworded the beta = 0 interpretation around the exact RMSE ratio (`pkgdown/assets/dev/bench_covariate_sim.html:150`). |

## Reviewed Without Patch

- `bench_zi_count` uses scenario-specific values in its top verdict and does
  not make a broad dominance/no-degradation claim in the active rendered page.
- `bench_avonet_missingness` was already tightened in the prior cleanup pass;
  the remaining "ordering is preserved" sentence is tied to explicit 80%
  missingness numbers on the same page.
- `bench_tree_uncertainty` and `bench_delhey` did not surface high-risk
  dominance/no-degradation wording in this pass.

## Remaining Work

1. Audit pending-but-hidden methodology assets before re-adding them to the
   navbar: `bench_multi_proportion`, `bench_scaling_v090`, Phase 8 pages, and
   BACE head-to-head pages.
2. Decide whether the generated HTML pages should keep wall-clock "Run on"
   timestamps from regeneration when the underlying benchmark RDS was not
   rerun. The current generators do this automatically.
3. Consider adding a tiny helper for all benchmark generators: a standard
   `claim_scope_note()` block that prints species, reps, missingness, DGP, and
   compared methods above each verdict.
