# Real-data check: lambda fixed at 1 against estimated lambda

Script: `script/bench_lambda_datasets.R`. Seven trait datasets (bundled AVONET 300, and the AVONET,
AmphiBIO, BIEN, GlobTherm, LepTraits and PanTHERIA snapshots), each as a random 300-species subset and a
2,000-species subset (GlobTherm: all 1,969). For each case: 5 seeds, 20% of observed cells of every
continuous trait masked MCAR, `impute(..., gnn = FALSE)` under `lambda_mode` fixed_1, estimate, cv and
bayes. Score: z-RMSE on the masked cells (log scale for positive skewed traits). Run on Totoro
2026-09-23, 260 of 260 jobs, 0 failures. Full tables: `lambda_datasets_report.txt`.

| dataset | species | fixed at 1 | estimated | change |
|---|---|---|---|---|
| PanTHERIA | 300 | 0.583 | 0.496 | -14.8% |
| PanTHERIA | 2000 | 0.415 | 0.358 | -13.7% |
| GlobTherm | 300 | 0.712 | 0.665 | -6.5% |
| GlobTherm | 1969 | 0.655 | 0.622 | -5.1% |
| AmphiBIO | 2000 | 0.685 | 0.665 | -2.8% |
| AVONET 300 (bundled) | 300 | 0.569 | 0.560 | -1.7% |
| LepTraits | 300 | 0.947 | 0.935 | -1.3% |
| AmphiBIO | 300 | 0.884 | 0.874 | -1.1% |
| BIEN | 2000 | 0.813 | 0.807 | -0.7% |
| AVONET | 300 | 0.629 | 0.626 | -0.3% |
| BIEN | 300 | 1.054 | 1.053 | -0.2% |
| AVONET | 2000 | 0.419 | 0.419 | +0.1% |
| LepTraits | 2000 | 0.893 | 0.901 | +0.9% |

Estimated lambda helps or ties in 12 of 13 cases and is never worse by more than 1%. Binary and
categorical accuracy (habitat, trophic level, terrestriality, primary lifestyle) is identical under all
four modes. The ordinal trait migration moved slightly on AVONET at 2,000 species (0.820 at lambda = 1,
0.823 otherwise): an ordinal code path was picking up the lambda setting, which the final review found and
which was fixed after this run, so ordinal traits now stay at lambda = 1 as designed.
Median fit time under "estimate" is within 0 to 12% of fixed_1; "cv" and "bayes" cost 10 to 50% more for
little extra gain, which supports "estimate" as the default.

Largest per-trait gains come from weak-signal life-history traits: PanTHERIA gestation (0.370 to 0.263)
and maximum longevity (0.672 to 0.589) at 2,000 species; BIEN specific leaf area (1.006 to 0.963), the
spec's motivating trait. Strong-signal traits (AVONET morphology, lambda 0.98 to 0.995) do not move.

Known limitation: LepTraits flight duration at 2,000 species. Lambda is estimated at about 0.35 and the
prediction is worse than at lambda = 1 (1.000 against 0.968); cross-validated lambda does best (0.972).
This is consistent with the downward bias of the profile-REML lambda estimate measured in the recovery
tests (about -0.05 at true lambda 0.3 and 0.7, n = 300).
