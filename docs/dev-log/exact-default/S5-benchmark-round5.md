# S5 round 5: auto after the second review's fixes (split by trait row; no split when routes agree)

Code: commit 2be557f (branch feat/exact-default), installed on Totoro as `~/R/lib-exact8`. Same design as
round 4: 18 core cells, GNN off, 200 seeds, 3,600 jobs at 120 in parallel, about 20 minutes. Summary:
`core_lambda_auto5_off_agg_summary.csv`. Baseline for comparison: `core_lambda_main_off_agg_summary.csv`.

What changed from round 4 (see `rose-review-2.md`, B1 and B2): the validation split now samples held-out
trait rows (species in multi-obs data) rather than latent cells, so the K cells of a categorical row stay in
one half; and a trait whose exact and per-column fits agree is not split.

## Change against main

| true lambda | z-RMSE (c1, c2, cnt, prp), n 100 / 300 / 1000 | coverage (c1, c2) | accuracy (bin, cat3, ord) |
|---|---|---|---|
| 0.3 | -3.1% / -4.6% / -6.6% | +0.003 / +0.004 / +0.009 | +0.009 / +0.014 / +0.021 |
| 0.7 | -6.1% / -8.5% / -9.4% | +0.008 / +0.005 / +0.010 | +0.013 / +0.027 / +0.032 |
| 1 | -1.2% / -0.3% / +0.4% | -0.007 / -0.003 / +0.002 | +0.000 / +0.000 / -0.001 |

z-RMSE and coverage match round 4 to the digits shown. That is expected: for a trait with one latent column
and differing fits, sampling rows draws the same cells as sampling cells did. Accuracy moved only through the
categorical trait, whose rows no longer leak between halves: the lambda = 1 losses of round 4 (-0.002) are
gone, and the largest single-scenario accuracy fall is now -0.004 (cat3, lambda 0.7, n 1000, unpaired
z = -0.3).

Coverage below main by more than 0.005 in a single scenario (rho 0 or 0.5): all at lambda 1, n 100 (c1, c2,
cnt, prp: -0.006 to -0.018, |z| <= 1.2, where main covers 0.86 to 0.92), plus cnt at lambda 0.7, n 100
(-0.006) and two at lambda 1, n 300 (prp -0.008, c1 -0.005).

## Real data

Not re-run. The 13 real-data cases have continuous traits only in single-observation data, where the
row split draws the same cells as round 4's cell split (inference from the code; round 5's unchanged
continuous results agree), so the round-4 table in `S5-benchmark-round4.md` stands.

Checks on this code: suite FAIL 0 / PASS 2735 (`S5e-suite.log`); R CMD check in `S5e-check.log`.
