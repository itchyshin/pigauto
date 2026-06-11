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
