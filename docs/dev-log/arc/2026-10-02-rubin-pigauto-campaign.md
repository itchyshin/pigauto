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
  lambda = 0.7, rho = 0.5 pig_post is unbiased (-0.000 against +0.006). The cause is open: lambda, Sigma_E and the
  Sigma_P correlation all look right, so the next step is to compare pig_post's conditional imputation of c2 given
  an observed c1 with the exact conditional under the true model on these datasets.
- n = 100: pig_post is close to freqA (slope coverage 0.941 against 0.949 overall), lower at lambda 0.7 on cor.
- On the downstream estimands pig_post covers better than freqB (improper MI) in every cell except n = 1000 at lambda = 1 (paired +0.005 to +0.046); per-cell coverage at n = 100 goes the other way (-0.011 and -0.018 at lambda 0.3 and 0.7).

## Does NOT cover

Discrete traits (posterior MI is continuous-only); MAR or clade missingness; real trees; the cause of the lambda = 1
bias (inference above); the 12 sampler failures (a pigauto fix); `study.qmd` and the published report pages are not
yet updated with the pig_post arm.

## Next

1. pigauto: diagnose the lambda = 1 bias (posterior lambda at n = 1000, lambda = 1) and the Cholesky failures; both
   are package work for a separate lane.
2. Study: add pig_post to `study.qmd` and the report pages once (1) says whether the lambda = 1 result is a bug.
