# TabPFN phylo-feature pilot note

Generated from `script/bench_tabpfn_phylo_features.rds` on 2026-06-10.

This is a branch-local benchmark note, not package documentation. It records
the first clean TabPFN pilot after installing TabPFN and accepting the local
inference license.

## Scope

- Data: bundled AVONET continuous traits only (`Mass`,
  `Beak.Length_Culmen`, `Tarsus.Length`, `Wing.Length`).
- Scales: `n = 50, 75, 300`.
- Replicates: 3 per scale.
- Missingness: pigauto `make_missing_splits()` MCAR holdout.
- Methods: BM baseline, pigauto, `tabpfn_plain`, `tabpfn_lappe`,
  `tabpfn_lappe_nfa`, `tabpfn_knn`.
- Comparison target: test-set latent-scale RMSE, Pearson r, and split-conformal
  coverage where available.

## Result

All 216 result rows completed with `status == "ok"`.

Best TabPFN variant beat pigauto on RMSE in all 12 scale-by-trait comparisons.
The mean best-TabPFN / pigauto RMSE ratio was 0.742; the median ratio was 0.758.

| n | trait | best TabPFN | TabPFN RMSE | pigauto RMSE | RMSE ratio | coverage |
|---:|---|---|---:|---:|---:|---:|
| 50 | Beak.Length_Culmen | tabpfn_plain | 0.4830 | 0.6153 | 0.785 | 0.874 |
| 50 | Mass | tabpfn_lappe_nfa | 0.2962 | 0.3937 | 0.752 | 0.790 |
| 50 | Tarsus.Length | tabpfn_plain | 0.3754 | 0.4429 | 0.848 | 0.917 |
| 50 | Wing.Length | tabpfn_knn | 0.3185 | 0.4602 | 0.692 | 0.861 |
| 75 | Beak.Length_Culmen | tabpfn_lappe | 0.6148 | 0.6441 | 0.955 | 0.745 |
| 75 | Mass | tabpfn_knn | 0.2394 | 0.3185 | 0.752 | 1.000 |
| 75 | Tarsus.Length | tabpfn_lappe | 0.4026 | 0.5266 | 0.765 | 0.939 |
| 75 | Wing.Length | tabpfn_plain | 0.3063 | 0.3962 | 0.773 | 0.852 |
| 300 | Beak.Length_Culmen | tabpfn_lappe | 0.4777 | 0.5725 | 0.835 | 0.983 |
| 300 | Mass | tabpfn_knn | 0.2004 | 0.3719 | 0.539 | 0.965 |
| 300 | Tarsus.Length | tabpfn_lappe | 0.3882 | 0.5218 | 0.744 | 0.964 |
| 300 | Wing.Length | tabpfn_lappe_nfa | 0.2311 | 0.4908 | 0.471 | 0.973 |

Average RMSE across scales and continuous traits:

| method | average RMSE |
|---|---:|
| tabpfn_lappe | 0.3874 |
| tabpfn_plain | 0.3897 |
| tabpfn_lappe_nfa | 0.4296 |
| tabpfn_knn | 0.4462 |
| pigauto | 0.4795 |

## Interpretation

This is encouraging evidence for a TabPFN benchmark lane. It is not evidence
that TabPFN should replace pigauto's core model. The measured regime is narrow:
continuous AVONET traits, MCAR holdout, small-to-mid `n`, no mixed-type cells,
no discrete traits, no MNAR, no multi-tree Rubin pooling, no multi-observation
covariate refinement, and no active-imputation guidance.

Honest current wording:

> TabPFN with cross-trait and phylogenetic features outperformed pigauto on
> continuous AVONET trait RMSE in this MCAR benchmark at n = 50, 75, and 300.
> This does not yet cover pigauto's mixed-type, uncertainty, multiple
> imputation, multi-tree, or active-sampling workflows.

## Larger check

Completed one larger continuous-trait check at `n = 2000`, `rep = 1`, using
`plain`, `lappe`, and `lappe_nfa`.

Command shape:

```sh
PIGAUTO_TABPFN_OUT_STEM=bench_tabpfn_phylo_features_n2000 \
PIGAUTO_TABPFN_SCALES=2000 \
PIGAUTO_TABPFN_REPS=1 \
PIGAUTO_TABPFN_VARIANTS=plain,lappe,lappe_nfa \
PIGAUTO_TABPFN_RUN_PIGAUTO=true \
Rscript script/bench_tabpfn_phylo_features.R
```

All 20 result rows completed with `status == "ok"`. Best TabPFN variant beat
pigauto on RMSE in all four continuous-trait comparisons. The mean
best-TabPFN / pigauto RMSE ratio was 0.803; the median ratio was 0.791.

| n | trait | best TabPFN | TabPFN RMSE | pigauto RMSE | RMSE ratio | coverage |
|---:|---|---|---:|---:|---:|---:|
| 2000 | Beak.Length_Culmen | tabpfn_lappe | 0.2844 | 0.4048 | 0.703 | 0.950 |
| 2000 | Mass | tabpfn_lappe | 0.2289 | 0.3428 | 0.668 | 0.957 |
| 2000 | Tarsus.Length | tabpfn_lappe_nfa | 0.2872 | 0.3267 | 0.879 | 0.959 |
| 2000 | Wing.Length | tabpfn_lappe_nfa | 0.8417 | 0.8750 | 0.962 | 0.938 |

This supports continuing the benchmark lane. The claim is still limited:
continuous AVONET traits, MCAR holdout, single `n = 2000` replicate, and no
mixed-type/discrete/multi-tree/downstream-inference workflows.

## Discrete AVONET check

Completed a local categorical/ordinal AVONET check on 2026-06-10, using
`script/bench_tabpfn_discrete_avonet.R`.

Command shape:

```sh
PIGAUTO_TABPFN_OUT_STEM=bench_tabpfn_discrete_avonet_local \
PIGAUTO_TABPFN_SCALES=50,75,300 \
PIGAUTO_TABPFN_REPS=3 \
PIGAUTO_TABPFN_VARIANTS=plain,lappe,lappe_nfa,knn \
PIGAUTO_TABPFN_RUN_PIGAUTO=true \
Rscript script/bench_tabpfn_discrete_avonet.R
```

Scope:

- Data: bundled AVONET discrete traits (`Trophic.Level`,
  `Primary.Lifestyle`, `Migration`).
- Scales: `n = 50, 75, 300`.
- Replicates: 3 per scale.
- Missingness: pigauto `make_missing_splits()` MCAR holdout.
- Methods: majority baseline, pigauto's phylogenetic baseline, pigauto,
  `tabpfn_plain`, `tabpfn_lappe`, `tabpfn_lappe_nfa`, `tabpfn_knn`.
- Comparison target: test-set accuracy and balanced accuracy only. This check
  does not implement classification prediction sets.

All 189 result rows completed with `status == "ok"`.

This did not reproduce the continuous-trait TabPFN win. Using the most
favorable TabPFN reading, where the best TabPFN variant is selected separately
for each scale, replicate, and trait, best TabPFN beat pigauto on balanced
accuracy in 9/27 comparisons and beat the phylogenetic baseline in 5/27
comparisons. It tied pigauto in 6/27 comparisons and tied the baseline in 7/27
comparisons. Mean balanced accuracy across all scale-trait-replicate cells was:

| method | mean balanced accuracy | mean accuracy |
|---|---:|---:|
| baseline | 0.5916 | 0.7382 |
| pigauto | 0.5272 | 0.6926 |
| best TabPFN variant | 0.5216 | 0.7180 |
| majority | 0.3944 | 0.6176 |

Scale-by-trait balanced accuracy means:

| n | trait | baseline | pigauto | best TabPFN | majority |
|---:|---|---:|---:|---:|---:|
| 50 | Migration | 0.7381 | 0.7381 | 0.8333 | 0.8333 |
| 50 | Primary.Lifestyle | 0.5417 | 0.5417 | 0.4861 | 0.3889 |
| 50 | Trophic.Level | 0.5926 | 0.4630 | 0.5648 | 0.3889 |
| 75 | Migration | 0.6111 | 0.4444 | 0.4583 | 0.4444 |
| 75 | Primary.Lifestyle | 0.5833 | 0.3750 | 0.4167 | 0.2500 |
| 75 | Trophic.Level | 0.5602 | 0.4861 | 0.5671 | 0.3889 |
| 300 | Migration | 0.4764 | 0.4764 | 0.3510 | 0.3333 |
| 300 | Primary.Lifestyle | 0.5530 | 0.5526 | 0.5108 | 0.2167 |
| 300 | Trophic.Level | 0.6677 | 0.6677 | 0.5061 | 0.3056 |

Honest current wording:

> TabPFN with cross-trait and phylogenetic features looks promising for
> continuous AVONET imputation, including a single `n = 2000` check, but the
> local categorical/ordinal AVONET check is mixed and generally favors
> pigauto's phylogenetic baseline on balanced accuracy. This branch should be
> treated as a benchmark lane, not as evidence for a broad replacement of
> pigauto's mixed-type model.

## Simulated Mixed Scalar-Type Check

Completed a local mixed scalar-type simulation check on 2026-06-10, using
`script/bench_tabpfn_sim_types.R`.

Command shape:

```sh
PIGAUTO_TABPFN_OUT_STEM=bench_tabpfn_sim_types_local \
PIGAUTO_TABPFN_SCALES=75,150 \
PIGAUTO_TABPFN_REPS=2 \
PIGAUTO_TABPFN_SCENARIOS=mixed_moderate,mixed_high_phylo,mixed_sparse_imbalanced \
PIGAUTO_TABPFN_VARIANTS=plain,lappe,lappe_nfa \
PIGAUTO_TABPFN_RUN_PIGAUTO=true \
PIGAUTO_TABPFN_EPOCHS=200 \
Rscript script/bench_tabpfn_sim_types.R
```

Scope:

- Data: simulated mixed scalar datasets with two continuous traits and one
  trait each for count, proportion, binary, ordinal, and categorical data.
- Scales: `n = 75, 150`.
- Scenarios: `mixed_moderate`, `mixed_high_phylo`,
  `mixed_sparse_imbalanced`.
- Replicates: 2 per scale-scenario cell.
- Missingness: pigauto `make_missing_splits()` MCAR holdout.
- Methods: mean/mode baseline, pigauto's phylogenetic baseline, pigauto,
  `tabpfn_plain`, `tabpfn_lappe`, `tabpfn_lappe_nfa`.
- Comparison targets: latent-scale RMSE for continuous/count/proportion, and
  balanced accuracy for binary/ordinal/categorical.

The run produced 504 result rows. All baseline, mean/mode, and pigauto rows
completed with `status == "ok"`. TabPFN completed 249/252 rows; all three
skips were the same `n = 75`, `mixed_high_phylo`, `rep = 2`, `binary_1`
target where too few held-out cells were available.

Using an optimistic "best TabPFN variant per scale-scenario-replicate-trait"
reading, best TabPFN beat pigauto in 34/83 comparable target-cells and beat
the phylogenetic baseline in 32/83 comparable target-cells. This is the
favorable reading for TabPFN because it selects the best feature variant
separately for each target.

Overall means by type:

| type | metric | baseline | pigauto | best TabPFN | mean/mode |
|---|---|---:|---:|---:|---:|
| continuous | RMSE | 0.5521 | 0.5798 | 0.6852 | 1.0369 |
| count | RMSE | 0.7528 | 0.7910 | 0.7912 | 1.0164 |
| proportion | RMSE | 0.8010 | 0.7841 | 0.7543 | 0.9772 |
| binary | balanced accuracy | 0.7359 | 0.7099 | 0.6832 | 0.5417 |
| categorical | balanced accuracy | 0.6147 | 0.6008 | 0.6345 | 0.4167 |
| ordinal | balanced accuracy | 0.3758 | 0.3364 | 0.4283 | 0.1500 |

Best TabPFN beat pigauto most often for proportion (8/12), categorical
(6/12), and ordinal (8/12), but not for continuous (5/24), count (4/12), or
binary (3/11 comparable cells). In the `n = 150`, `mixed_high_phylo` scenario,
the phylogenetic baseline/pigauto were stronger than best TabPFN across
continuous, count, proportion, binary, and categorical targets.

Honest current wording:

> The simulated mixed scalar-type check supports keeping TabPFN as an
> experimental benchmark lane and suggests possible promise for proportion,
> categorical, and ordinal targets. It does not support promoting TabPFN as
> pigauto's replacement GNN/ML model. Pigauto's phylogenetic baseline remains
> stronger for continuous, count, and binary targets in this local grid, and
> the high-phylogenetic-signal regime favors the existing baseline/pigauto
> path. Classification prediction sets, multi-tree pooling,
> multi-observation covariates, and active-imputation guidance remain untested
> by this branch.

### Larger Simulated Scalar Check

Completed one additional local `n = 300` mixed scalar check on 2026-06-10:

```sh
PIGAUTO_TABPFN_OUT_STEM=bench_tabpfn_sim_types_n300 \
PIGAUTO_TABPFN_SCALES=300 \
PIGAUTO_TABPFN_REPS=1 \
PIGAUTO_TABPFN_SCENARIOS=mixed_moderate,mixed_high_phylo,mixed_sparse_imbalanced \
PIGAUTO_TABPFN_VARIANTS=plain,lappe,lappe_nfa \
PIGAUTO_TABPFN_RUN_PIGAUTO=true \
PIGAUTO_TABPFN_EPOCHS=200 \
Rscript script/bench_tabpfn_sim_types.R
```

All 126 result rows completed with `status == "ok"`.

Using the same favorable best-TabPFN reading, best TabPFN beat pigauto in 6/21
target-cells and beat the phylogenetic baseline in 4/21 target-cells.

`n = 300` means by type:

| type | metric | baseline | pigauto | best TabPFN | mean/mode |
|---|---|---:|---:|---:|---:|
| continuous | RMSE | 0.5144 | 0.5144 | 0.6818 | 0.9385 |
| count | RMSE | 0.6733 | 0.6892 | 0.7790 | 0.9941 |
| proportion | RMSE | 0.7446 | 0.7502 | 0.7850 | 1.0377 |
| binary | balanced accuracy | 0.6680 | 0.6680 | 0.6577 | 0.5000 |
| categorical | balanced accuracy | 0.6459 | 0.5683 | 0.5830 | 0.3333 |
| ordinal | balanced accuracy | 0.3279 | 0.3279 | 0.4331 | 0.2000 |

This larger check strengthens the conclusion above: in this local simulation,
TabPFN is most interesting for ordinal classification, but it is not ready to
be promoted as pigauto's replacement ML/GNN model.

## Complex-Type Check

Completed a local `zi_count` and `multi_proportion` check on 2026-06-10, using
`script/bench_tabpfn_complex_types.R`.

Command shape:

```sh
PIGAUTO_TABPFN_OUT_STEM=bench_tabpfn_complex_types_local \
PIGAUTO_TABPFN_SCALES=150,300 \
PIGAUTO_TABPFN_REPS=1 \
PIGAUTO_TABPFN_SCENARIOS=zi_moderate,zi_sparse,zi_many_zeros,multi_moderate,multi_high_phylo,multi_K8 \
PIGAUTO_TABPFN_VARIANTS=lappe,lappe_nfa \
PIGAUTO_TABPFN_RUN_PIGAUTO=true \
PIGAUTO_TABPFN_EPOCHS=200 \
Rscript script/bench_tabpfn_complex_types.R
```

Scope:

- Data: simulated `zi_count` datasets with two zero-inflated count traits, and
  simulated `multi_proportion` datasets with one grouped compositional trait.
- Scales: `n = 150, 300`.
- Scenarios: `zi_moderate`, `zi_sparse`, `zi_many_zeros`,
  `multi_moderate`, `multi_high_phylo`, `multi_K8`.
- Replicates: 1 per scale-scenario cell.
- Methods: mean/mode baseline, pigauto's phylogenetic baseline, pigauto,
  `tabpfn_lappe`, and `tabpfn_lappe_nfa`.
- TabPFN wrapper: `zi_count` uses one TabPFN classifier for the non-zero gate
  and one TabPFN regressor for the conditional non-zero magnitude;
  `multi_proportion` uses independent TabPFN regressions on the z-scored CLR
  latent columns.
- Comparison targets: expected-count RMSE, zero accuracy, and Brier score for
  `zi_count`; Aitchison distance, CLR RMSE, simplex MAE, and
  dominant-component accuracy for `multi_proportion`.

All 90 result rows completed with `status == "ok"`.

`zi_count` means across both scales and all three ZI scenarios:

| method | RMSE | zero accuracy | Brier |
|---|---:|---:|---:|
| baseline | 13.6395 | 0.7281 | 0.1921 |
| mean/mode | 14.2382 | 1.0000 | 0.1765 |
| pigauto | 14.7313 | 0.6998 | 0.2187 |
| tabpfn_lappe | 14.0938 | 0.5628 | 0.3006 |
| tabpfn_lappe_nfa | 17.2431 | 0.4867 | 0.4396 |

Metric-specific optimistic TabPFN wins, selecting the better of `lappe` and
`lappe_nfa` separately for each scale-scenario-replicate-trait row:

| metric | wins vs pigauto | wins vs baseline |
|---|---:|---:|
| expected-count RMSE | 6/12 | 3/12 |
| Brier score | 2/12 | 1/12 |
| zero accuracy | 2/12 | 2/12 |

`multi_proportion` means across both scales and all three composition
scenarios:

| method | Aitchison | CLR RMSE | simplex MAE | accuracy |
|---|---:|---:|---:|---:|
| baseline | 4.4104 | 0.8087 | 0.1353 | 0.5395 |
| mean/mode | 5.3588 | 1.0128 | 0.1794 | 0.3377 |
| pigauto | 4.4104 | 0.8087 | 0.1353 | 0.5395 |
| tabpfn_lappe | 4.4932 | 0.8443 | 0.1508 | 0.4649 |
| tabpfn_lappe_nfa | 5.6242 | 1.0546 | 0.1826 | 0.3728 |

Metric-specific optimistic TabPFN wins:

| metric | wins vs pigauto | wins vs baseline |
|---|---:|---:|
| Aitchison distance | 3/6 | 3/6 |
| CLR RMSE | 3/6 | 3/6 |
| simplex MAE | 1/6 | 1/6 |
| dominant-component accuracy | 3/6 | 3/6 |

Honest current wording:

> The complex-type check does not support promoting TabPFN as pigauto's
> replacement ML/GNN model. The ZI wrapper has occasional expected-count RMSE
> wins, but its non-zero gate calibration is worse on average by Brier score
> and zero accuracy. The multi-proportion wrapper is close in some cells but is
> worse than pigauto/baseline on the averaged Aitchison, CLR RMSE, simplex MAE,
> and dominant-component accuracy metrics. TabPFN should remain an
> experimental benchmark lane unless a future integration preserves pigauto's
> grouped mixed-type semantics, conformal uncertainty, multi-tree pooling,
> multi-observation refinement, and active-imputation guidance.

## Russell-Style Reconciliation Check

Completed a local reconciliation ablation on 2026-06-10 / 2026-06-11, using
`script/bench_tabpfn_reconciliation.R`.

Command shape:

```sh
PIGAUTO_TABPFN_OUT_STEM=bench_tabpfn_reconciliation_local \
PIGAUTO_TABPFN_SCALES=50,75,300 \
PIGAUTO_TABPFN_REPS=3 \
PIGAUTO_TABPFN_REGIMES=same_row,phylo_only,same_row_lappe,shuffled_same_row_lappe,same_row_lappe_nfa \
PIGAUTO_TABPFN_RUN_PIGAUTO=true \
PIGAUTO_TABPFN_EPOCHS=500 \
Rscript script/bench_tabpfn_reconciliation.R
```

Scope:

- Data: bundled AVONET continuous traits (`Mass`,
  `Beak.Length_Culmen`, `Tarsus.Length`, `Wing.Length`).
- Scales: `n = 50, 75, 300`.
- Replicates: 3 per scale.
- Missingness: pigauto `make_missing_splits()` MCAR holdout.
- Methods: pigauto's phylogenetic baseline, pigauto, and five TabPFN feature
  regimes.
- Aim: determine whether Russell-style TabPFN gains are driven by same-row
  cross-trait information, phylogenetic features, or their combination.

Feature regimes:

| regime | features |
|---|---|
| `same_row` | other observed same-row traits only |
| `phylo_only` | Laplacian phylogenetic eigenvectors only |
| `same_row_lappe` | same-row traits plus Laplacian eigenvectors |
| `shuffled_same_row_lappe` | same features as `same_row_lappe`, but other-trait rows shuffled across species |
| `same_row_lappe_nfa` | same-row traits, Laplacian eigenvectors, and nearest-training-target aggregates |

All 252 result rows completed with `status == "ok"`.

Mean RMSE across all scales, replicates, and continuous traits:

| method | mean RMSE |
|---|---:|
| tabpfn_same_row_lappe | 0.3871 |
| tabpfn_same_row | 0.4006 |
| tabpfn_same_row_lappe_nfa | 0.4714 |
| pigauto | 0.5299 |
| baseline | 0.6532 |
| tabpfn_phylo_only | 0.8493 |
| tabpfn_shuffled_same_row_lappe | 0.8737 |

Mean RMSE by scale:

| n | pigauto | same-row | same-row + LapPE | phylo-only | shuffled same-row + LapPE |
|---:|---:|---:|---:|---:|---:|
| 50 | 0.5364 | 0.4007 | 0.3920 | 0.7275 | 0.7653 |
| 75 | 0.5274 | 0.4401 | 0.4347 | 1.0122 | 1.0507 |
| 300 | 0.4713 | 0.3611 | 0.3345 | 0.8082 | 0.8052 |

Mean RMSE by trait:

| trait | pigauto | same-row | same-row + LapPE | phylo-only | shuffled same-row + LapPE |
|---|---:|---:|---:|---:|---:|
| Beak.Length_Culmen | 0.6843 | 0.5874 | 0.5566 | 0.9811 | 0.9972 |
| Mass | 0.3982 | 0.2678 | 0.2627 | 0.8670 | 0.8659 |
| Tarsus.Length | 0.5572 | 0.4303 | 0.4204 | 0.7642 | 0.7865 |
| Wing.Length | 0.4798 | 0.3170 | 0.3087 | 0.7849 | 0.8454 |

Pigauto gate audit:

| n | trait | pigauto RMSE | r_cal_bm | r_cal_gnn | r_cal_mean |
|---:|---|---:|---:|---:|---:|
| 50 | Beak.Length_Culmen | 0.7305 | 0.5833 | 0.2103 | 0.2063 |
| 50 | Mass | 0.4949 | 0.4500 | 0.5500 | 0.0000 |
| 50 | Tarsus.Length | 0.6083 | 0.4833 | 0.5167 | 0.0000 |
| 50 | Wing.Length | 0.3118 | 0.2833 | 0.7167 | 0.0000 |
| 75 | Beak.Length_Culmen | 0.6642 | 0.3500 | 0.5759 | 0.0741 |
| 75 | Mass | 0.3950 | 0.0000 | 0.9833 | 0.0167 |
| 75 | Tarsus.Length | 0.5734 | 0.4670 | 0.5330 | 0.0000 |
| 75 | Wing.Length | 0.6951 | 0.3833 | 0.2833 | 0.3333 |
| 300 | Beak.Length_Culmen | 0.6583 | 0.7500 | 0.2500 | 0.0000 |
| 300 | Mass | 0.3046 | 0.4000 | 0.6000 | 0.0000 |
| 300 | Tarsus.Length | 0.4899 | 0.8140 | 0.1860 | 0.0000 |
| 300 | Wing.Length | 0.4324 | 0.2500 | 0.6807 | 0.0693 |

Interpretation:

> Russell's stronger TabPFN result is compatible with these local results.
> The large continuous-trait gain is real in this AVONET setup, but the
> reconciliation ablation says it is driven mostly by aligned same-row
> cross-trait information. When same-row trait rows are shuffled, TabPFN
> collapses to the phylogeny-only regime. Therefore the next useful pigauto
> development direction is not "replace the GNN with TabPFN"; it is to examine
> why pigauto's current gated model does not exploit cross-trait information as
> strongly as same-row TabPFN in these continuous AVONET cells, while preserving
> pigauto's mixed-type semantics and uncertainty workflow.

## Pigauto Cross-Trait Diagnostic

Completed a focused pigauto-side diagnostic on 2026-06-11, using
`script/bench_pigauto_cross_trait.R`.

Command shape for the expensive pigauto configs:

```sh
PIGAUTO_XTRAIT_OUT_STEM=bench_pigauto_cross_trait_n300 \
PIGAUTO_XTRAIT_SCALES=300 \
PIGAUTO_XTRAIT_REPS=3 \
PIGAUTO_XTRAIT_CONFIGS=default,trait_attention,relaxed_gate \
PIGAUTO_XTRAIT_EPOCHS=500 \
Rscript script/bench_pigauto_cross_trait.R
```

A cheap ridge-only rerun used the same seeds with
`PIGAUTO_XTRAIT_RUN_PIGAUTO=false`. The ridge comparator drops near-zero
training-variance features (`PIGAUTO_XTRAIT_RIDGE_MIN_SD=1e-6`) because some
LapPE columns can be effectively zero in the training rows but non-zero in the
held-out rows.

Scope:

- Data: bundled AVONET continuous traits at `n = 300`.
- Replicates: 3.
- Missingness: pigauto `make_missing_splits()` MCAR holdout.
- Pigauto configs: current default, trait attention, and a relaxed-gate
  diagnostic (`gate_cap = 1`, `lambda_gate = 0`, `phylo_signal_gate = FALSE`).
- Simple comparators: ridge same-row traits, and ridge same-row traits plus
  LapPE.

All 72 combined rows completed with `status == "ok"`.

Mean RMSE at `n = 300`:

| method | mean RMSE | mean Pearson r |
|---|---:|---:|
| pigauto_relaxed_gate | 0.4682 | 0.8914 |
| pigauto_trait_attention | 0.4932 | 0.8765 |
| pigauto_default | 0.5023 | 0.8774 |
| ridge_same_row | 0.5052 | 0.8651 |
| ridge_same_row_lappe | 0.5378 | 0.8508 |
| baseline | 0.5452 | 0.8581 |
| tabpfn_same_row_lappe | 0.3345 | -- |
| tabpfn_same_row | 0.3611 | -- |

The TabPFN rows are copied from the reconciliation check above for the same
`n = 300` scale, so they are not new rows from this pigauto diagnostic.

Mean pigauto RMSE by trait:

| trait | default | trait attention | relaxed gate | TabPFN same-row + LapPE |
|---|---:|---:|---:|---:|
| Beak.Length_Culmen | 0.6689 | 0.6949 | 0.6629 | 0.4772 |
| Mass | 0.3201 | 0.3365 | 0.3226 | 0.2203 |
| Tarsus.Length | 0.5065 | 0.5117 | 0.4982 | 0.3424 |
| Wing.Length | 0.5138 | 0.4296 | 0.3893 | 0.2980 |

Gate audit:

| trait | relaxed-gate RMSE | mean BM gate | mean GNN gate | mean mean-gate |
|---|---:|---:|---:|---:|
| Beak.Length_Culmen | 0.6629 | 0.6667 | 0.3333 | 0.0000 |
| Mass | 0.3226 | 0.4833 | 0.5167 | 0.0167 |
| Tarsus.Length | 0.4982 | 0.8500 | 0.1500 | 0.0000 |
| Wing.Length | 0.3893 | 0.1000 | 0.8649 | 0.0351 |

Interpretation:

> More local testing supports the same direction as the reconciliation check.
> The relaxed-gate diagnostic improves pigauto at `n = 300`, especially for
> `Wing.Length`, but it does not close the gap to same-row TabPFN. Trait
> attention by itself is not enough in this AVONET continuous setup. A simple
> ridge same-row model is close to current pigauto but far behind TabPFN,
> suggesting TabPFN is exploiting more than a plain linear allometry signal.

Actionable next step:

> Keep TabPFN as an experimental benchmark/integration lane, and open a
> separate pigauto-side development issue for cross-trait learning/gate
> calibration. The issue should start from this bounded claim: in AVONET
> continuous MCAR cells, same-row TabPFN still outperforms default pigauto and
> the relaxed-gate diagnostic, while the current mixed-type and uncertainty
> workflows remain pigauto advantages that this test does not replace.
