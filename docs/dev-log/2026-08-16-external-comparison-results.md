# External comparison — the first working pigauto-vs-BACE head-to-head

2026-08-16 · Claude lane · branch `arc/bace-comparators`.

## What was actually wrong

Every committed BACE head-to-head in this repo reports **"BACE skipped (not installed or
failed)"**, and `bace_ran = FALSE` in every rds. The cause was **not** a missing install —
BACE 0.0.0.9000 and MCMCglmm have been installed locally since 2026-08-09. The bench
scripts called a signature that does not exist:

```r
BACE::bace(data = df_miss, tree = tree_sub, n_iter = 2000L, ovr = TRUE)   # no such args
```

against the real `bace(fixformula, ran_phylo_form, phylo, data, nitt, ...)`. The `tryCatch`
swallowed `"argument fixformula is missing"` and the md writer printed a **fixed string**
that conflated "not installed" with "errored" — so the failure was invisible for months.
`script/bench_avonet_bace.R:17–22` had already diagnosed this in a comment and contained a
*correct* call, but was never re-run after the install.

Fixed on `arc/bace-comparators`: correct signatures in
`bench_pantheria_bace_head_to_head.R` and `bench_bace_avonet_head_to_head.R`, and all three
scripts now surface the actual `conditionMessage` instead of the fixed string. A second
blocker surfaced locally: **pigauto was not installed in the R 4.6 user library at all**,
so the first re-run died at `library(pigauto)` — worth remembering when a local bench
"skips".

## The result: BACE wins the continuous traits, pigauto wins the categorical

AVONET300, n = 300 × 7 traits, seed 2026, 30% MCAR, 90 held-out cells/trait, m = 20 for
pigauto. **Single seed — a direction, not an interval.**

| trait | metric | pigauto | BACE | winner |
|---|---|---:|---:|---|
| Mass | RMSE | 2434.9 | **1921.1** | BACE (−21%) |
| Beak.Length_Culmen | RMSE | 22.27 | **14.91** | BACE (−33%) |
| Tarsus.Length | RMSE | 23.30 | **11.84** | BACE (−49%) |
| Wing.Length | RMSE | 75.59 | **48.78** | BACE (−35%) |
| Trophic.Level | accuracy | **0.789** | 0.500 | pigauto (+29 pp) |
| Primary.Lifestyle | accuracy | **0.689** | 0.389 | pigauto (+30 pp) |
| Migration | accuracy | 0.578 | 0.578 | tie |
| — | wall | 362 s | 59 s | BACE 6× faster |

**The continuous result is if anything understated.** BACE ran with short chains
(`nitt = 2000`, `burnin = 500`, `n_final = 2`) against pigauto's 20 imputations, so its
point estimate averaged only two draws and still won by 21–49%.

pigauto beats the mean/mode floor everywhere (Mass 2435 vs 2622; Trophic.Level 0.789 vs
0.567), so this is not a broken pipeline — it is a genuine continuous-trait deficit against
a well-tuned Bayesian competitor, alongside a genuine categorical-trait advantage.

## Why, and what to do about it

BACE is MCMCglmm with a phylogenetic random effect: its posterior over the
phylo/residual variance ratio is **Bayesian model averaging over what pigauto calls Pagel's
λ**. pigauto's default baseline is fixed at λ = 1.

The λ-attribution campaign run the same day
(`2026-08-16-lambda-attribution-results.md`) independently measured that
`lambda_mode = "bayes"` gives the best baseline repair under λ-misspecification and
eliminates the ML boundary collapse that `estimate` suffers (P(λ̂<0.05) = 0.47 → 0.00 at
n=100). **The two findings point at the same mechanism**, from simulation and from real
data respectively.

Follow-up run: `script/bench_avonet_lambda_modes.R` — same data, seed, and mask, comparing
`lambda_mode = "fixed_1"` against `"bayes"`, so the delta is directly comparable to BACE's
column above. Results appended below when it completes.

## Regime and honesty notes

- ONE dataset, ONE seed, n = 300, 30% MCAR, single-obs. No MCSE. This is a direction to
  investigate, **not** a publishable head-to-head — that needs multi-seed and ≥2 datasets.
- BACE's own published Trophic.Level figure (~72%) is **BACE's number on BACE's setup**,
  not reproduced here; this run gives BACE 0.500 on that trait with short chains and no OVR.
  Do not cite either number as the other's.
- BACE remains `Suggests:`-only and is 404 on CRAN, so this comparison cannot ship in a
  CRAN build; it is dev-log evidence.
- Not yet run: the Rphylopars/phylolm standalone arm (`bench_external_comparators.R`) —
  its agent hit the session limit mid-build. That comparator is still owed.
