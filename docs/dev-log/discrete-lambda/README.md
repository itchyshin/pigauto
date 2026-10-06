# Discrete-trait lambda prototype: screen against BACE (2026-10-05)

Requested by Shinichi ("prototype λ for discrete traits"; aim: catch up with BACE on categorical accuracy).
Branch `feat/discrete-lambda`. Prototype switch, default off: `options(pigauto.discrete_lambda = "estimate")` adds
the binary and ordinal liability columns to the joint solver's `lambda_cols` in `fit_joint_threshold_baseline()`;
categorical traits get it through the one-vs-rest fits. Before this, every discrete trait was imputed at lambda = 1
(the "section 7 B iii cut" of D-278).

## Design

- Datasets: the four-arm study's core slice (branch `arc/imputation-sim`, `script/campaign_sim_cell.R`,
  `types_mixed`, BM, MCAR 30%, fixed thresholds, always-observed driver), n = 100 and 300, lambda = 0.3, 0.7, 1,
  rho = 0 and 0.5, seeds 1 to 100: the datasets BACE ran on. BACE's stored per-seed results
  (Totoro `~/pigauto_sim/results/core_bace/`) are compared dataset by dataset; 1,197 of 1,200 paired (BACE files
  missing for 3). The re-run mode floor matched the stored floor exactly on the smoke dataset.
- pigauto: the prototype (30eeedb) in a private Totoro library, `impute(..., gnn = FALSE)` as in the study, with the
  safety machinery on (`gnn_off`, the default) and off (`gnn_off_pure`: `safety_floor = FALSE,
  phylo_signal_gate = FALSE`), each with the switch off and on. Totoro, 80 cores, 2 x 1,200 jobs, about 4 min each.
- Screen 2 also tried giving the ordinal Brownian-motion candidate an estimated lambda (a4718d8); it lowered ordinal
  accuracy by 0.009 to 0.018 and was reverted (8540543). The tables below use screen 2's other columns, which are
  unaffected; `dlam_screen1.txt` has the ordinal figures without that change.

## Result (`dlam_screen2.txt`)

Discrete accuracy pooled over binary, ordinal and three-class traits and rho:

| n | lambda | mode floor | current default (gate on, lambda 1) | gate off, lambda estimated | BACE |
|---|---|---|---|---|---|
| 100 | 0.3 | 0.459 | 0.462 | 0.526 | 0.526 |
| 100 | 0.7 | 0.545 | 0.617 | 0.664 | 0.640 |
| 100 | 1 | 0.626 | 0.900 | 0.923 | 0.768 |
| 300 | 0.3 | 0.469 | 0.476 | 0.546 | 0.566 |
| 300 | 0.7 | 0.560 | 0.647 | 0.685 | 0.687 |
| 300 | 1 | 0.641 | 0.952 | 0.953 | 0.792 |

(screen 1 figures for the gate-off, lambda-estimated column, i.e. without the reverted ordinal change.)

- Paired against BACE, the gap at lambda = 0.3 shrinks from -0.064 / -0.090 (n = 100 / 300) to 0.000 / -0.020.
- Binary and categorical: at or above BACE in every cell (n = 300, lambda = 0.3: 0.712 vs 0.710 and 0.548 vs 0.544).
- Ordinal is the remaining gap at low signal: 0.356 / 0.379 vs BACE 0.410 / 0.443 at lambda = 0.3.
- Both parts are needed. Lambda alone (gate on) moves little, because the phylo-signal gate replaces weak-signal
  discrete predictions with the mode; gate off alone (lambda 1) leaves categorical below the mode floor at
  n = 300, lambda = 0.3 (0.448 vs 0.450).
- Discrete Brier score: gate off, lambda estimated is the lowest of all arms in 4 of 6 cells (at lambda = 1 the lambda-1 arms are 0.006 to 0.009 lower), and below BACE in all 6.

Turning the gate off also changes continuous traits (zRMSE pooled over c1, c2, cnt, prp; lower is better):

| n | lambda | gate on | gate off | BACE |
|---|---|---|---|---|
| 100 | 0.3 | 0.971 | 0.892 | 1.013 |
| 300 | 0.3 | 0.930 | 0.857 | 0.921 |
| 300 | 1 | 0.389 | 0.386 | 0.738 |

Gate off is better in every cell, and conformal coverage is unchanged or slightly higher (n = 300: 0.963 to 0.966;
n = 100: 0.880 to 0.893). The September GNN-off campaign found the opposite on its low-signal `bace_dgp` (the
safety machinery helped), when the continuous baseline still ran at lambda = 1; since D-278 lambda is estimated, which
may make the gate redundant. That has not been tested on `bace_dgp`, OU or AVONET.

## Screen 3: which part of the safety machinery blocks discrete lambda (`screen3_attribution.txt`)

Same 1,197 paired datasets; new arm `gnn_off_nogate` (`phylo_signal_gate = FALSE`, safety floor on), added to a copy
of the harness on Totoro (`~/pigauto_sim/dlam/script/campaign_gnn_off_lib.R`). With the floor off the gate cannot
impose the mode (it falls back to the BM baseline), so "floor off, gate on" equals "both off".

| n = 300, lambda = 0.3 | discrete accuracy | vs BACE (paired) | continuous zRMSE |
|---|---|---|---|
| gate + floor, lambda 1 (default) | 0.476 | -0.090 (0.004) | 0.930 |
| floor only, lambda estimated | 0.529 | -0.036 (0.003) | 0.861 |
| neither, lambda estimated | 0.546 | -0.020 (0.003) | 0.857 |
| BACE | 0.566 | | 0.921 |

The gate is the larger cost (most of both the discrete and the continuous gain); the floor costs more on top,
largest at n = 100 (0.497 vs 0.526) and at lambda = 1 (binary 0.924 vs 0.959 at n = 100). Continuous coverage at
n = 100: 0.874 to 0.880 with the floor, 0.883 without; at n = 300 all about 0.96.

## Screen 5: cumulative ordinal decomposition (`screen5_ordinal.txt`, branch `feat/discrete-lambda-ordinal`, 8e99670)

`options(pigauto.ordinal_method = "cumulative")`: K-1 "class >= k" binary liabilities through the threshold/one-vs-rest
machinery, each with its own lambda, monotone-combined, mode taken. Paired against the current ordinal route (both
with discrete lambda estimated, safety machinery off), ordinal accuracy:

| | n = 100 | n = 300 |
|---|---|---|
| lambda = 0.3 | +0.066 (0.011) | +0.032 (0.006) |
| lambda = 0.7 | -0.036 (0.010) | -0.067 (0.007) |
| lambda = 1 | -0.005 (0.004) | 0.000 (0.001) |

Macro-F1 is lower everywhere (it collapses toward the middle classes). Not adopted: a negative result. Ordinal stays
on the current route; against BACE that is -0.053 / -0.064 at lambda = 0.3 and +0.016 / -0.025 at 0.7 (n = 100 / 300).

## Screen 4 and the gate verdict (`screen4_crossdgp.txt`, `fisher-review.md`)

Cross-DGP screen: 56 factorial cells (types_mixed, BM and OU, MAR 30%, clade 30%, MCAR 10%; lambda 0.3 and 1;
n = 100 with 100 seeds, n = 1000 with up to 30, paired with stored BACE) plus bace_dgp, bm_mixed and AVONET300
(20 seeds each). 3,785 of 3,785 job pairs usable. Continuous zRMSE with gate and floor off is never worse (best
-0.094); discrete accuracy at lambda = 0.3 gains 0.043 to 0.085 in all 14 types_mixed cells; costs at lambda = 1,
n = 1000 (accuracy -0.002 to -0.005, Brier +0.011 to +0.016) and on AVONET (Brier +0.017).

Correction: the `ou_mixed` runs duplicated `bm_mixed` because the jobs passed `--evo BM`, overriding that DGP's
OU default. OU evidence comes from the factorial OU cells only; `ou_mixed` is re-run with `--evo OU`.

Independent statistical review (`fisher-review.md`, 2026-10-05): **not supported as a default change yet;
supported as an opt-in.** Main points: the AVONET Brier loss removes about 90% of pigauto's probability skill over
the mode on the only real dataset, and it comes from discrete lambda (floor,est also +0.017; none,l1 +0.001); the
large-n Brier cost grows with n; the regime the floor and gate were built for (real weak-signal data, lambda near
0, all-discrete data on the LP path) was not tested; "not flagged" in the small cells is a power problem, not
evidence of no harm. Recommended next: (i) types_mixed at lambda = 0 and 0.1, (ii) BIEN plants, (iii) AVONET
per-trait Brier and an n = 3000 cell; change defaults only if non-inferior on all three.

VERDICT_WRITTEN (gate GA3)

## Option-C tests (`screen6_optionC.txt`)

types_mixed at lambda = 0 and 0.1 (n = 300 and 1000, 120 datasets per cell), all-discrete data (bin, ord, cat3 only:
the label-propagation path; lambda 0.1 to 1, n = 300), `ou_mixed` with OU evolution (replacing the duplicated run),
and one n = 3000, lambda = 1 cell (10 datasets). At lambda = 0 and 0.1, the regime the gate and floor were built
for, the default is WORSE than the mean or mode (continuous zRMSE 1.008 to 1.020 vs 0.999 to 1.007; discrete
accuracy 0.359 to 0.414 vs mode 0.363 to 0.426); gate and floor off with lambda estimated gives zRMSE 0.891 to 0.911
and accuracy 0.470 to 0.509. All-discrete: +0.038 to +0.040 accuracy at lambda 0.1 to 0.3, +0.003 at lambda 1.
ou_mixed (OU): +0.005 to +0.073 accuracy, continuous unchanged or better. n = 3000, lambda = 1: accuracy -0.003
(SE 0.001), Brier +0.014 (SE 0.009).

## Screen 7: "auto" (validation-chosen lambda per discrete trait; branch `feat/discrete-lambda-auto`, 6969a20)

Per binary/categorical trait, fit lambda = 1 and lambda estimated and keep the lower validation Brier. Result
(`screen7_auto.txt`, 2,635 datasets): accuracy never worse than the default; the AVONET Brier cost disappears (0.517
vs default 0.516); but most of the low-signal Brier gain is lost (lambda 0.3, n = 1000: 0.518 vs estimate 0.483) and
on bace_dgp it is worse than both candidates (0.531 vs 0.491 and 0.513). `auto_diag.R` shows why: (1) a bug, the
recorded choice and the applied prediction can disagree (lambda = 1, seed 1: "fixed_1" recorded, estimated
prediction applied); (2) the choice itself is noisy (it picked lambda = 1 for the binary trait in 3 of 3 datasets at
lambda = 0.3), because 20 to 40 validation cells per trait cannot separate the two fits. Parked as a negative result.
"floor + auto" matches the default at lambda = 1 and on AVONET, so the floor alone protects those cases.

## Incident: thread cap

A hand-written timing smoke on Totoro (2026-10-05 ~20:14 to 20:25) set `OPENBLAS_NUM_THREADS=1` but not
`OMP_NUM_THREADS`, so 5 R processes ran about 72 threads each (about 360 cores, over the 150-core cap) for about
10 minutes before I killed them. The screens themselves use `run_one.sh`, which sets both; their timings (about 8 s per
job) confirm one core each. Re-run with both caps: 4 to 11 s per job at n = 300, 55 s at n = 1000.

ORD_SCREEN 1197 (gate GB2) · PAIRED 1197 (gate GA1)

## Does NOT cover

One DGP family (the four-arm core: BM, MCAR, fixed thresholds); n = 1000; multi-observation data; the LP path
(categorical traits with no continuous trait alongside); real data; the gate change on other DGPs; posterior MI
(continuous-only).
