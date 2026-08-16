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
`lambda_mode = "fixed_1"` against `"bayes"`. **Result (single seed, m=20):**

| trait | fixed_1 | bayes | BACE (ref) | reading |
|---|---:|---:|---:|---|
| Mass RMSE | 2625.3 | **2301.7** | 1921.1 | bayes closes ~46% of the gap; Pearson 0.15 → 0.83 (the fixed_1 Mass fit was pathological this run) |
| Beak RMSE | 22.19 | 21.27 | 14.91 | small lift |
| Tarsus RMSE | 27.06 | 24.42 | 11.84 | partial |
| Wing RMSE | 67.77 | 67.89 | 48.78 | no change |
| Trophic.Level acc | **0.789** | 0.600 | 0.500 | **bayes COSTS categorical** |
| Primary.Lifestyle acc | **0.667** | 0.578 | 0.389 | same |
| Migration acc | 0.622 | 0.756 | 0.578 | noise-level swing |

Three honest readings:

1. **Direction confirmed, magnitude partial.** `bayes` improves every continuous trait it
   can move (Mass dramatically — the λ=1 kriging was actively mis-extrapolating there) but
   closes only part of the BACE gap. BACE's remaining edge plausibly includes its joint
   MCMC over the full covariance, not just λ.
2. **The trade-off nobody would guess from the docs:** λ modes set `force_per_column`,
   which switches off the joint/threshold/OVR baselines — so **discrete traits silently
   drop to label propagation** and Trophic.Level loses 19 pp. `lambda_mode = "bayes"` is
   currently a continuous-trait tool with a categorical price on mixed data. The clean
   fix is per-type dispatch (λ per-column for continuous, joint for discrete) — filed as
   the natural next code slice.
3. **Single-seed noise is material:** the two fixed_1 runs on the *same mask* differ
   (Mass 2434.9 in the BACE bench vs 2625.3 here — different call settings and GNN
   stochasticity). Nothing in this section is publishable without the multi-seed version.

## Regime and honesty notes

- ONE dataset, ONE seed, n = 300, 30% MCAR, single-obs. No MCSE. This is a direction to
  investigate, **not** a publishable head-to-head — that needs multi-seed and ≥2 datasets.
- BACE's own published Trophic.Level figure (~72%) is **BACE's number on BACE's setup**,
  not reproduced here; this run gives BACE 0.500 on that trait with short chains and no OVR.
  Do not cite either number as the other's.
- BACE remains `Suggests:`-only and is 404 on CRAN, so this comparison cannot ship in a
  CRAN build; it is dev-log evidence.
- ~~Not yet run: the Rphylopars/phylolm standalone arm~~ — **run, 5 seeds, below.**

## Rphylopars + phylolm standalone (multi-seed, MCSEs) — the uncomfortable one

`script/bench_external_comparators.R` (arc/bace-comparators `c7f6bf0`): AVONET300, same
30% MCAR masks across methods, 5 reps (seeds 2027–2031), scored on the 4 continuous
traits. pigauto fits all 7 traits (its selling point); rivals fit the continuous 4.
z-RMSE from training-portion mean/sd, no leakage. mean ± MCSE:

| trait | pigauto | rphylopars | phylolm(λ) | col-mean |
|---|---:|---:|---:|---:|
| Mass | 1.594 ± 0.508 | **1.360 ± 0.341** | 1.629 ± 0.473 | 1.734 |
| Beak.Length | 0.912 ± 0.162 | **0.445 ± 0.086** | 0.672 ± 0.083 | 1.204 |
| Tarsus.Length | 1.220 ± 0.270 | **0.639 ± 0.120** | 1.046 ± 0.213 | 1.516 |
| Wing.Length | 0.688 ± 0.083 | **0.409 ± 0.066** | 0.675 ± 0.105 | 1.120 |

**Raw Rphylopars beats pigauto's delivered predictions on all four traits** (clearly on
3, within noise on Mass), and phylolm(λ) beats pigauto on Beak/Tarsus. Consistent with
the BACE result measured independently the same day — two external Bayesian/joint
solvers now outperform pigauto's continuous-trait output on this dataset.

**What this does and does not mean.** Rphylopars is pigauto's own internal joint-baseline
solver, so this is NOT "pigauto loses to an unrelated method" — it measures everything
pigauto's pipeline layers on top of (and around) the raw joint fit, and that layer is
currently **net negative for continuous traits on AVONET300**. Candidate mechanisms, in
testable order:
1. **Calibration tax**: pigauto's baseline sees ~49% of cells (70% observed × val/test
   held out) vs the standalone's 70% — the price of gate calibration and conformal
   scores. Quantifiable by refitting the baseline on all observed cells post-calibration.
2. **Which path fired is unverified**: with ordinal/categorical traits present the
   threshold-joint (liability) path dispatches, not the pure continuous joint — its
   continuous columns may be degraded by the liability machinery, or the dispatch may
   even have fallen back per-column. One diagnostic run answers this.
3. The λ=1 default (partially — `bayes` closed some of it in the section above).

Mixed-type remains pigauto's win: no rival here can touch the categorical traits at all,
and BACE loses them by ~30 pp. But the continuous gap is now a measured, replicated,
multi-seed finding and belongs on the next arc's front page, not in a footnote.

Regime: AVONET300 only, n=300, 30% MCAR, single-obs, 5 reps. pigauto wall ~12–16
min/rep single-threaded.
