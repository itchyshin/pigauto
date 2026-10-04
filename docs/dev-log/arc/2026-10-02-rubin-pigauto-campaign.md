# Rubin study: pigauto posterior MI arm (pig_post), campaign results

Lane `claude:pigauto-mi-posterior` · branch `arc/rubin-freq-bace` · 2026-10-02 · Claude Code (Opus 5.5)

Results page (private): https://claude.ai/artifact/Y644K1sLxbrWsAmTQvZspC

## Run

- Approved by Shinichi 2026-10-01 ("go Totoro, option b"). Driver `script/rubin_pig_campaign_totoro.sh`.
- Totoro, 110 cores, 2026-10-01 20:38 to 2026-10-02 07:59 (11 h 21 min). 3,600 of 3,600 cell files.
- pigauto b565cad (#189 merged), private library. `multi_impute(draws_method = "posterior", m = 20, log_transform = FALSE)`
  on freqA's block (c1, c2, logit prp, driver d1). Option (b): `max_extend = 6` at n = 100, pigauto's default 3 elsewhere.
- Same datasets as the stored freq and BACE campaign: the complete-data estimates agree with the freq files to within
  6.2e-10 in all 7,200 pig/freq pairs (`sanity.txt`; corrected 2026-10-04, see below).
- Raw files: `~/pigauto_rubin_pool/pig/totoro/` (169 MB with logs; never delete). Tables:
  `script/rubin_study/data/agg_pig/` from `script/rubin_campaign_aggregate.R` (extended; the original tables are
  reproduced byte for byte when no `pig/` pool exists).

## Correction (2026-10-02): the lambda = 1 numbers above are pooled over rho

The 0.85 to 0.86 coverage and the -0.021 bias quoted above for n = 1000, lambda = 1 average rho = 0 (where pig_post
covers 0.97) with rho = 0.5. The problem is entirely in the correlated case, it is larger, and it already shows at
n = 300 (found by the study-report pass; checked against `script/rubin_study/data/agg_pig/down.csv`):

| lambda = 1, rho = 0.5 | n = 100 | n = 300 | n = 1000 |
|---|---|---|---|
| pig_post slope coverage | 0.955 | 0.875 | 0.750 |
| pig_post correlation coverage | 0.940 | 0.890 | 0.735 |
| freqA slope coverage | 0.955 | 0.920 | 0.935 |
| slope bias, pig_post / freqA | -0.061 / -0.021 | -0.050 / -0.026 | -0.045 / -0.022 |

At lambda 0.3 and 0.7 with rho = 0.5, pig_post covers 0.93 to 0.98 with bias under 0.01; at rho = 0 it is fine at
every lambda. The pigauto caveat merged in #201 quoted the pooled 0.86; PR #202 corrects it to these numbers.

## Failures and convergence

| n | lambda | sampler errors | unconverged |
|---|---|---|---|
| 100 | 0.3 | 10 | 56 |
| 100 | 0.7 | 1 | 51 |
| 300 | 0.3 | 1 | 0 |
| all other cells | | 0 | 0 |

The 12 errors are a Cholesky failure inside the sampler ("leading principal minor of order ~8n is not positive"),
mid-fit, all at lambda 0.3 or 0.7. Each aborts the whole fit; it is a pigauto robustness bug for the package, not a
property of the study. Failed fits have no estimand rows and are counted here. Unconverged fits are kept and flagged.
At n = 100 they cover at least as well as converged fits (0.96 to 0.98 against 0.91 to 0.96), so they do not explain
the n = 100 results below.

## Results (200 datasets per cell; BACE 100 at n = 1000; MC SE in `down_l.csv`, `paired.csv`)

Per-cell imputation (c1, c2 pooled): pig_post equals freqA at n = 300 and 1000 (coverage 0.949 and 0.954 against
0.949 and 0.953; same zRMSE and width). At n = 100 it is slightly worse: coverage 0.934 against 0.941, zRMSE 0.620
against 0.603; paired coverage difference -0.014 and -0.021 (SE 0.002) at lambda 0.3 and 0.7, +0.013 at lambda 1.

Downstream estimands (Rubin-pooled PGLS slope of c2 on c1, and phylogenetic correlation):

| n | lambda | estimand | pig_post coverage | freqA | chained BACE | pig_post minus freqA (paired, SE) |
|---|---|---|---|---|---|---|
| 1000 | 1 | slope | 0.860 | 0.963 | 0.968 | -0.102 (0.016) |
| 1000 | 1 | cor | 0.853 | 0.948 | 0.937 | -0.095 (0.016) |
| 1000 | 0.3, 0.7 | slope | 0.963 to 0.970 | 0.970 | 0.975 to 0.980 | -0.008 to 0.000 |
| 300 | all | slope | 0.952 overall | 0.952 | 0.964 | -0.007 to +0.008 |
| 100 | 0.7 | cor | | | | -0.030 (0.012) |

- Finding of record: pig_post matches proper frequentist MI (freqA) everywhere except n = 1000 at lambda = 1, where it
  under-covers (0.85 to 0.86) with a bias of -0.021 on both estimands (freqA -0.013 to -0.014; similar interval
  width, 0.162 against 0.160). The #189 sweep already found a small negative bias at lambda = 1 (-0.004 to -0.014);
  at n = 1000 the intervals are narrow enough for it to break coverage.
- Tested and ruled out (2026-10-02): lambda is not under-estimated. Sixteen campaign fits were refitted exactly
  (`script/rubin_pig_lambda_diag.R`, summary `script/rubin_pig_lambda_diag_summary.R`; n = 1000, lambda = 1, rho 0 and
  0.5, seeds 1 to 6, plus lambda = 0.7, rho = 0.5, seeds 1 to 4) with the posterior draws saved. At lambda = 1 the
  posterior mean lambda is 1.000 for c1, c2 and d1 in all 12 fits, every draw above 0.99 (only prp sits lower, at
  0.93, and prp is not in the slope); Sigma_E for c1 is about 0; the posterior phylogenetic correlation of c1 and c2
  is 0.47 to 0.58 at a true rho of 0.5. freqA also estimates lambda at the boundary (0.9996).
- Where the error is: at lambda = 1 and rho = 0.5 the mean slope error against complete data is -0.034 for pig_post
  and -0.013 for freqA (pig_post worse in 5 of 6 datasets); at rho = 0 they agree (+0.006 and +0.004); at
  lambda = 0.7, rho = 0.5 pig_post is unbiased (-0.000 against +0.006). Cause found (below).

## Cause of the lambda = 1 slope bias (2026-10-02)

Eight campaign datasets were refitted keeping the completed datasets (`script/rubin_pig_cond_diag.R`; n = 1000,
rho = 0.5, lambda = 1 seeds 1 to 6, lambda = 0.7 seeds 1 to 2) and compared with the exact conditional under the true
model (`script/rubin_pig_cond_summary.R`, `script/rubin_pig_oracle.R`, `script/rubin_pig_oracle_plus_e.R`).

1. Cell level, pig_post is right. For missing c2 with c1 observed, at lambda = 1: mean difference from the exact
   conditional 0.000, shrinkage slope 0.999, SD ratio 1.02, RMSE against truth 0.037 against the exact 0.036. The
   same holds for c1, for rows with both missing, and at lambda = 0.7. freqA also matches (SD ratio 0.94 to 0.97).
2. Oracle MI (20 joint draws from the exact conditional, scored with the campaign's estimator and Rubin pooling):
   mean slope error against complete data -0.004 at lambda = 1, against freqA -0.013 and pig_post -0.034. pig_post is
   0.023 to 0.039 below the oracle in every one of the 6 datasets, so the bias is pigauto's and systematic.
3. Mechanism. Adding pig_post's own posterior Sigma_E (from the lambda refits of the same datasets) as noise on the
   oracle's imputed cells moves the oracle from -0.004 to -0.028, about 80 percent of pig_post's error, downward in
   all 6 datasets. That noise is tiny (variance about 1e-4 per trait on the data scale) and nearly uncorrelated
   between c1 and c2 (0.05 to 0.17 against a true rho of 0.5). At lambda = 1 the true Sigma_E is 0; the prior
   Sigma_E ~ IW(K + 1, 0.01 diag(observed latent variances)) keeps it above 0 and pulls its correlation towards 0.
   With n = 1000 and lambda = 1 the PGLS slope rests on sister-pair differences of a similar size (imputation RMSE
   0.04), so this residual noise attenuates the slope. At rho = 0 there is no correlation to attenuate, and at
   lambda < 1 Sigma_E is real and estimated, which matches where the bias appears.

What would fix it is a pigauto package change, not a study change, and the choice is Shinichi's: a smaller or
data-adaptive Sigma_E prior scale, a prior that lets Sigma_E reach 0 (parameter expansion on Sigma_E as on Sigma_P),
or a parameterisation in which the residual correlation follows the phylogenetic one. Any of these needs re-running
the #189 acceptance sweep and this cell (n = 1000, lambda = 1).
- n = 100: pig_post is close to freqA (slope coverage 0.941 against 0.949 overall), lower at lambda 0.7 on cor.
- On the downstream estimands pig_post covers better than freqB (improper MI) in every cell except n = 1000 at lambda = 1 (paired +0.005 to +0.046); per-cell coverage at n = 100 goes the other way (-0.011 and -0.018 at lambda 0.3 and 0.7).

## Prior trial (2026-10-02)

`script/rubin_pig_prior_trial.R` patches `.mip_fit()` in memory only (no package edit) and refits 8 datasets
(n = 1000, rho = 0.5; lambda = 1 seeds 1 to 6, lambda = 0.7 seeds 1 to 2) under three Sigma_E priors; scored by
`script/rubin_pig_prior_trial_summary.R` against the oracle on the same datasets.

| variant | Sigma_E prior scale | lambda = 1: slope error minus oracle | covered rho | per-cell coverage (20-draw) | lambda = 0.7: minus oracle |
|---|---|---|---|---|---|
| base (current) | diag(0.01 obs var) | -0.030 | 5 of 6 | 0.898 | -0.009 |
| A | diag(1e-4 obs var) | -0.006 | 6 of 6 | 0.876 | -0.010 |
| B | 0.01 x observed covariance | -0.019 | 5 of 6 | 0.893 | -0.011 |
| oracle | | 0 | | 0.867 | |

Variant A removes about 80 percent of the gap to the oracle at lambda = 1, leaves lambda = 0.7 unchanged, and brings
per-cell calibration closer to the oracle's (the current prior's intervals are too wide by the same spurious noise).
It is also about 8 percent faster. Before it ships it must pass, on the package side: the #189 acceptance sweep (G3,
G6, G7), the low-signal and small-n cells where the sampler failures and non-convergence occurred (lambda 0.3,
n = 100), and this cell. Eight datasets is a screen, not evidence.

### Small-n check of variant A (2026-10-02): not a clean win

Campaign datasets at n = 100, rho = 0.5, lambda 0.3 and 0.7, seeds 1 to 10, current prior against A (same scripts).

| cell | prior | slope error minus oracle | converged | per-cell coverage (20-draw) |
|---|---|---|---|---|
| lambda 0.3 | current | +0.025 | 8 of 10 | 0.865 |
| lambda 0.3 | A | +0.020 | 8 of 10 | 0.847 |
| lambda 0.7 | current | +0.017 | 8 of 10 | 0.838 |
| lambda 0.7 | A | +0.024 | 4 of 9 | 0.830 |

A also crashed once (lambda 0.7, seed 9, the CHOLMOD "not positive definite" error) where the current prior did not.
With real residual variance at small n, the tiny prior scale lets Sigma_E wander towards 0 and mixing suffers. A flat
1e-4 trades the large-n lambda = 1 bias for worse small-n convergence, so it is not shipped. Candidates that might
avoid the trade-off, untested: a scale that shrinks with n (for example 0.1 / n, which equals A at n = 1000 and
1e-3 at n = 100), or parameter expansion on Sigma_E as on Sigma_P. Branch `fix/mi-posterior-sigma-e-prior` holds the
flat-1e-4 change as work in progress (the pinned first-run test is deliberately not re-pinned).

### Variant C, S_E scale 0.1 / n (2026-10-02): does not avoid the trade-off

Same datasets and scoring. At n = 1000 C equals A by construction (1e-4) and gives identical results (lambda = 1 gap
to the oracle -0.006 against -0.030 for the current prior). At n = 100 (scale 1e-3):

| cell | current | A (1e-4) | C (0.1 / n) |
|---|---|---|---|
| lambda 0.3: converged, gap to oracle | 8 of 10, +0.025 | 8 of 10, +0.020 | 10 of 10, +0.021 |
| lambda 0.7: converged, gap to oracle | 8 of 10, +0.017 | 4 of 9, +0.024 | 4 of 9, +0.043 |

C also crashed on lambda 0.7, seed 9, as A did. Any smaller S_E scale, fixed or shrinking with n, worsens mixing
when the residual variance is real at small n. Ten seeds per cell is a screen (4 of 9 against 8 of 10 is suggestive,
not decisive), but A and C agree. The remaining options are parameter expansion on Sigma_E (as on Sigma_P), which
targets mixing near Sigma_E = 0 without moving the prior's scale, or keeping the current prior and documenting the
large-n, lambda = 1 caveat.

## Separation-strategy residual prior (2026-10-02): passes the screen

pigauto branch `research/mi-posterior-sep-prior` (opt-in `posterior_control$residual_prior = "sep"`, default "iw"):
Sigma_E = diag(s) R diag(s), s_k ~ half-Cauchy(0, observed SD of trait k), R ~ LKJ(1), as a density on Sigma_E's
entries (Jacobian 2^K prod s_k^K) so the existing Metropolis moves stay valid; the conjugate IW draw is skipped. A
first version crashed in all 6 n = 1000, lambda = 1 fits (the half-Cauchy let Sigma_E collapse and the precision of
the current state could not be factorised); a floor of 1e-6 on the eigenvalues of the standardised Sigma_E fixed it.
Screen with the floor (variant SEPF; same datasets and oracle as above; 28 fits, 0 crashes):

| cell | current prior | flat 1e-4 (A) | SEPF |
|---|---|---|---|
| n = 1000, lambda = 1: gap to oracle | -0.030 | -0.006 | -0.002 |
| n = 1000, lambda = 0.7: gap to oracle | -0.009 | -0.010 | -0.019 |
| n = 100, lambda = 0.3: converged, gap | 8 of 10, +0.025 | 8 of 10, +0.020 | 10 of 10, -0.002 |
| n = 100, lambda = 0.7: converged, gap | 8 of 10, +0.017 | 4 of 9, +0.024 | 10 of 10, -0.016 |
| per-cell coverage n = 100 (oracle about 0.87) | 0.865, 0.838 | 0.847, 0.830 | 0.867, 0.880 |
| mean time per fit, n = 100 | about 470 s | about 570 s | about 320 s |

SEPF removes the lambda = 1 bias without the small-n convergence cost of the scale fixes, and is faster at n = 100.
The n = 1000, lambda = 0.7 gap (-0.019, 2 datasets) needs the full run to settle. This is a screen (2 to 10 datasets
per cell), not evidence; making "sep" the default needs the #189 acceptance sweep and the Rubin campaign re-run under it.

## Full re-run under "sep" (pig_sep, 2026-10-03)

Approved by Shinichi 2026-10-03 ("go validation"). Same driver (`ARMS=pig_sep`), 100 cores on Totoro, 05:16 to 18:17
(13 h, run beside the 40-core acceptance sweep). pigauto cd3a4cb (`research/mi-posterior-sep-prior`, PR #204), private
library `~/pigauto_rubin/rlib_sepf`; all 3,600 fits record that sha and `residual_prior = "sep"`. Same seed offset (909)
and datasets as pig_post: the complete-data estimates agree with freq to within 6.2e-10 on all 7,200 pairs
(`sanity.txt`). Raw files
`~/pigauto_rubin_pool/pig_sep/totoro/`; tables `script/rubin_study/data/agg_pig_sep/`.

**Failures and convergence.** 3,600 of 3,600 fits, 0 sampler errors, 0 unconverged (pig_post: 12 errors and 107
unconverged, all at n = 100 or 300, lambda 0.3 or 0.7). 6 fits needed one chain extension. Median max R-hat 1.004,
median min ESS 1,270 to 1,450.

**The lambda = 1, rho = 0.5 bias is gone.** Slope coverage (bias against complete data), 200 datasets per cell; the MC
SE of a coverage near 0.95 is about 0.015:

| lambda = 1, rho = 0.5 | n = 100 | n = 300 | n = 1000 |
|---|---|---|---|
| pig_post | 0.955 (-0.061) | 0.875 (-0.046) | 0.750 (-0.048) |
| pig_sep | 0.955 (-0.037) | 0.940 (-0.015) | 0.965 (-0.009) |
| freqA | 0.955 (-0.021) | 0.920 (-0.023) | 0.935 (-0.025) |
| bace_chain | 0.945 (-0.092) | 0.953 (-0.048) | 0.968 (-0.031) |
| correlation coverage, pig_post / pig_sep | 0.940 / 0.965 | 0.890 / 0.925 | 0.735 / 0.945 |

Paired on the same datasets (pooled over rho), pig_sep minus pig_post slope coverage at n = 1000, lambda = 1 is +0.103
(SE 0.018). The gain comes from removing the bias, not from wider intervals: pig_sep's slope CI widths are within 3.5% (0.968 to 1.033 times) of
pig_post's in every cell.

**Elsewhere.** In the other 15 cells pig_sep covers 0.940 to 0.975 for both slope and correlation (pig_post 0.915 to
0.980). No downstream paired contrast against pig_post is below zero by more than 1.3 SE. At n = 100, lambda 0.3 and
0.7, pig_sep lifts pig_post's under-coverage (0.915 to 0.934 become 0.945 to 0.960). One change to watch: at n = 100,
rho = 0.5, lambda 0.3 and 0.7, pig_sep's slope bias is -0.016 and -0.021 (pig_post +0.005, +0.007; freqA +0.007,
+0.005), with coverage 0.960. This looks like small-n shrinkage from the proper prior; BACE's is larger (-0.030, -0.037).

**Per-cell imputation intervals** (c1 and c2 pooled; `cells.csv`) move toward 0.95 everywhere: pig_sep 0.942 to 0.959
in all 18 cells, pig_post 0.914 to 0.963. pig_post's n = 100, lambda < 1 cells were 0.914 to 0.925 (now 0.942 to
0.947); at lambda = 1 it over-covered at 0.960 to 0.963 (now 0.951 to 0.959).

**Cost.** 1,294 core-hours against pig_post's 1,238. Median time per fit is 0 to 18% longer, most at n = 1000,
lambda < 1. The two runs had different loads (110 cores alone; 140 with the sweep), so this is not a clean timing
comparison, but the screen's "faster at n = 100" does not hold: 420 s against 380 s at lambda < 1.

## #189 acceptance sweep under "sep" (2026-10-03)

The 8,000-cell sweep that gated #189 (40 regimes x 200 reps, same seeds and masks), re-run at ebcd466 with
`MI_POST_RESIDUAL_PRIOR=sep`, 40 cores on Totoro, 05:21 to 18:38: 8,000 of 8,000 cells, 0 failed. Evidence on the #204
branch: `docs/dev-log/mi-posterior/sep_validation/` (summaries, gate logs, iw against sep comparison).

| | iw (of record) | sep |
|---|---|---|
| G6 (`04_acceptance.R`) | pass | pass |
| pooled relative SE ratio, band [0.95, 1.10] | 1.059 | 1.077 |
| gated rows outside [0.90, 1.15] (reported) | 3 | 3 |
| gated mean abs paired bias | 0.0052 | 0.0051 |
| converged fits | 15,998 of 16,000 | 16,000 of 16,000 |
| G7 (`05_cell_coverage.R`), gated per-cell coverage | pass, 0.933 to 0.955 | pass, 0.945 to 0.956 |

The sweep's own lambda = 1 caveat (small negative bias in the twin rows) is about the same under "sep" (-0.002 to
-0.016 against -0.004 to -0.014). That effect is a few thousandths and the sweep never showed the Rubin study's large
lambda = 1 bias, so the sweep could not have caught it: the two checks answer different questions.

## Verdict on making "sep" the default (#204)

Evidence supports it. On the Rubin datasets it removes the one real failure of posterior MI (lambda = 1, rho = 0.5:
slope coverage 0.750 to 0.965 at n = 1000), removes all 12 sampler errors and 107 unconverged fits, and moves per-cell
coverage toward 0.95 in every cell, with no cell made worse beyond Monte Carlo noise. On the #189 sweep it passes both
gates with numbers close to the current default. Costs and caveats, for the decision:

1. A small negative slope bias at n = 100, rho = 0.5 (about -0.02; coverage still 0.96), absent under "iw".
2. G6's pooled SE ratio moves from 1.059 to 1.077, closer to the 1.10 edge (MI SEs slightly conservative relative to
   complete data).
3. Not faster: about 5% more core-hours in the Rubin run.
4. Continuous traits only, simulated data, MCAR and the sweep's MAR/MNAR regimes; no real-data re-run under "sep" yet
   (the #189 G8 real-data cells ran under "iw").

The default was Shinichi's decision: "sep" became the default on 2026-10-04 (#204 merged, 45011fc), with the docs,
NEWS and tests updated in the same PR.

## Real data under the new default (2026-10-04)

9 of the 10 #189 real-data cells (PanTHERIA 4,027 species x 6 masks, AVONET 1,500 species x 3 masks; all but
FishBase) re-run at pigauto main 45011fc, same masks and harness (`script/mi_realdata/01_run.R`, launcher and tables in
`script/rubin_study/data/realdata_sep/`). Totoro, 9 cores, 04:31 to 06:01 (estimate 1.5 to 2 h). 9 of 9 cells ok.

| | "iw" (#189, 69670d4) | "sep" (45011fc) |
|---|---|---|
| converged cells | 7 of 9 | 9 of 9 |
| model coverage, PanTHERIA mean (range), 24 trait-cells | 0.931 (0.869 to 0.955) | 0.932 (0.873 to 0.970) |
| model coverage, AVONET mean (range), 12 trait-cells | 0.957 (0.913 to 0.987) | 0.955 (0.917 to 0.980) |
| interval width, sep / iw | | 0.993 to 1.011 |
| pooled slopes within 5% of the complete-row reference | 20 of 27 | 20 of 27 |
| largest slope shift, in reference SEs | 4.05 | 3.92 |
| wall time per cell | 1.0 to 1.5 h | 0.6 to 1.5 h |

The real-data results are essentially unchanged. The PanTHERIA slope shifts recorded under "iw" (body mass ~ head-body
length positive under MCAR masks, negative under structured masks; longevity ~ body mass) are the same under "sep", so
the residual prior does not cause them.

**FishBase** (10,484 species, 5 traits, clade-structured mask; approved by Shinichi 2026-10-04, "run FishBase"): one
core, 06:38 to 15:13, 8.6 h against 21.1 h under "iw" (estimate was 21 to 25 h). Converged (max R-hat 1.006, min ESS
689).

| Trait | model coverage, iw | model coverage, sep | split conformal | width, sep / iw |
|---|---|---|---|---|
| DepthRangeDeep | 0.957 | 0.957 | 0.966 | 1.00 |
| Length | 0.942 | 0.942 | 0.939 | 1.00 |
| Troph | 0.933 | 0.937 | 0.939 | 1.00 |
| Vulnerability | 0.950 | 0.950 | 0.954 | 1.00 |
| Weight | 0.939 | 0.942 | 0.950 | 0.94 |

Pooled slopes against the complete-row reference: Weight ~ Length -0.6% (iw -0.7%), DepthRangeDeep ~ Length +0.8%
(iw +1.1%), Troph ~ Length +2.6% (iw +2.9%).

**G8 on all 10 cells under "sep"** (`03_acceptance.R`, `g8.log`): `REALDATA_COMPLETE`; every fit converged (under "iw"
8 of 10; the two PanTHERIA cells short on ESS now reach 459 and 655); 23 of 30 pair-cells within 5% of the reference,
as under "iw".

## Correction (2026-10-04): the pig/freq identity check

Until 2026-10-04 `script/rubin_campaign_aggregate.R` compared the pigauto files' complete-data estimates with the freq
files after removing duplicate rows, which kept the BACE or freq copy of each dataset, so the reported "7,200 pairs, max
|diff| 0" compared freq with itself. The check now reads each pigauto pool's own rows: pig and pig_sep both agree with
freq to within 6.2e-10 on all 7,200 pairs (floating-point differences between machines). The datasets are the same;
"exactly" was wrong. No table changed.

## Does NOT cover

Discrete traits (posterior MI is continuous-only); MAR or clade missingness in the Rubin study; real trees (the
simulation); the real-data cells are one mask set per dataset, not a coverage study.

## Next

1. Done 2026-10-04: Shinichi chose "sep" as the default; #204 merged (45011fc) with local tests (0 failures),
   `--as-cran` (0 errors, 0 warnings) and CI green.
2. Done 2026-10-04: all 10 real-data cells under "sep"; G8 `REALDATA_COMPLETE`, results unchanged (above).
3. Done 2026-10-04: `study.qmd` reports pig_post and pig_sep (db2d32d).
