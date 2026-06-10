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
