# Rubin study: pigauto posterior MI arm (pig_post), campaign results

Lane `claude:pigauto-mi-posterior` · branch `arc/rubin-freq-bace` · 2026-10-02 · Claude Code (Opus 5.5)

Results page (private): https://claude.ai/artifact/Y644K1sLxbrWsAmTQvZspC

## Run

- Approved by Shinichi 2026-10-01 ("go Totoro, option b"). Driver `script/rubin_pig_campaign_totoro.sh`.
- Totoro, 110 cores, 2026-10-01 20:38 to 2026-10-02 07:59 (11 h 21 min). 3,600 of 3,600 cell files.
- pigauto b565cad (#189 merged), private library. `multi_impute(draws_method = "posterior", m = 20, log_transform = FALSE)`
  on freqA's block (c1, c2, logit prp, driver d1). Option (b): `max_extend = 6` at n = 100, pigauto's default 3 elsewhere.
- Same datasets as the stored freq and BACE campaign: the complete-data estimates agree exactly in all 7,200 pig/freq
  pairs (`sanity.txt`).
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

## Does NOT cover

Discrete traits (posterior MI is continuous-only); MAR or clade missingness; real trees; the 12 sampler failures (a pigauto fix); `study.qmd` and the published report pages are not
yet updated with the pig_post arm.

## Next

1. pigauto: fix the Sigma_E prior behaviour at lambda = 1 (cause above) and the Cholesky failures; both are package
   work for a separate lane.
2. Study: add pig_post to `study.qmd` and the report pages once (1) says whether the lambda = 1 result is a bug.
