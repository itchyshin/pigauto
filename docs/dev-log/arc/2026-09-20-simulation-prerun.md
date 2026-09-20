# Pre-run: the four-arm imputation simulation, measured before launch

2026-09-20, Totoro, pigauto 0.11.0, branch `arc/imputation-sim` (from `origin/main` fd0b513).
Runner `script/campaign_sim_cell.R`; design table `script/campaign_sim_design.R`.

**This note exists so Shinichi can approve a launch against measured numbers rather than an estimate.
Nothing beyond the pre-run has been launched.**

## What was run

16 design cells, one replicate each: the 12 core cells at n = 100 and n = 1000 (lambda 0.3, 0.7, 1.0
x rho 0, 0.5; BM; MCAR 0.30), plus 4 sentinels at n = 1000 that sample the regimes the factorial adds
(OU lambda 0.3; MAR 0.30; clade-biased 0.30; MCAR 0.10). Every arm on every cell. Separately, a
convergence probe: BACE at n = 100 with runs in {3, 5, 10}, two seeds each.

A first attempt duplicated the Bayesian arm across both waves through a wrong environment variable and
was killed after an hour; the walls below come from the corrected run.

## Walls per replicate (seconds, 4 threads, one process per cell)

| arm | n = 100 | n = 1000 |
|---|---:|---:|
| 1 frequentist stack (Rphylopars + castor + phyloglm) | 0.35 | 2.0 |
| 3a pigauto GNN off, in-house solver | 0.7 | 38 |
| 3b pigauto GNN off, Rphylopars solver | 27 | 260 |
| 4 pigauto GNN on | 134 | 319 |
| 2 BACE (runs 2, n_final 20) | 835 | > 4,500 (killed at 75 min) |
| floor | 0.0 | 0.0 |

Two things follow. **Arm 3b is affordable**: it was the unmeasured risk in the plan, and at n = 1000 it
costs less than the GNN arm (260 s against 319 s), so both solvers can run as planned. **BACE sets the
budget**: it is roughly 2,400 times the frequentist stack at n = 100, and it is the only arm whose cost
at n = 1000 had to be inferred rather than measured, because the run was stopped.

## BACE convergence: the pre-run's main correction

Every cell first reported `converged = FALSE`. The cause was ours. BACE judges convergence with
`assess_convergence()` over its **initial imputation iterations** (`runs`) on the imputed values
themselves, using autocorrelation, percent change, trend and Geweke; that function has
`min_iterations = 3`, so the `runs = 2` we passed could not be assessed at all. Its vignette uses
`runs = 5` for demonstration and `runs = 15` with `nitt = 100000` for a real analysis.

Two related mistakes were corrected at the same time. `n_final` counts **full imputation runs**, each
refitting MCMCglmm per response, not posterior draws; a chain-length formula had set it to 400, about
eighty times the intended cost per cell. And Gelman-Rubin across `runs` is not a diagnostic here, since
successive runs condition on different completed datasets by construction: it read 5.7 on a fit BACE
itself called converged, and 20 when the species-level random effects were pooled in.

Probe, n = 100, `nitt` 50000 / burnin 10000 / thin 25, `n_final` 20:

| runs | converged | wall (s) | median ESS |
|---:|---|---:|---:|
| 3 | 0 of 2 | 739, 762 | 203, 77 |
| 5 | **2 of 2** | 798, 828 | 201, 70 |
| 10 | 2 of 2 | 968, 1001 | 225, 87 |

`runs = 5` converges at n = 100 and costs 18% less than `runs = 10`. It is the proposed setting, with
the caveat that this was measured at n = 100 only: convergence at n = 1000 is not yet established, and
the campaign records BACE's verdict per cell so the **convergence rate is reported as a result** rather
than assumed. `skip_conv` stays `TRUE`; the alternative retries up to `max_attempts` and makes per-cell
cost unbounded across thousands of cells.

Effective sample size is reported, not gated on its minimum. Median ESS per cell ran 42 to 1,457, while
the worst single parameter sat near 1 — expected for MCMCglmm threshold and categorical models, and a
disclosure item for the methods section.

## Failures

Two of twelve BACE cells failed outright with "mixed model equations singular", one of them the
clade-masked n = 1000 sentinel. Whole-clade masking leaves the mixed model equations singular, which is
a result about phylogenetically biased missingness rather than a defect. Failed fits are scored at the
mean/mode floor and their failure rate is reported beside every metric; nothing is dropped.

## Re-derived budget

Measured walls above, applied to the design, with n = 300 interpolated log-linearly and BACE at
n = 1000 extrapolated from its per-fit cost (38 s per fit at n = 100, times the 9.2-fold n = 100 to
n = 1000 ratio measured in the 2026-09-19 campaign). BACE at `runs = 5` plus `n_final = 20` is 25 fits
per replicate.

Per replicate: fast arms 162 s (n = 100), 307 s (n = 300), 619 s (n = 1000).
BACE 950 s, 2,739 s, **8,740 s (2.4 h)**.

| stage | fast arms | BACE |
|---|---:|---:|
| core slice, 18 cells | 362 h | 2,071 h |
| factorial, 56 cells | 1,214 h | 7,537 h |
| AVONET300 case study | 17 h (both) | |
| covariate sensitivity | 362 h | |
| **total** | | **11,564 slot-hours, BACE 83%** |

Against the plan's estimate of about 3,600 slot-hours, this is 3.2 times larger, entirely because the
plan costed BACE at `n_final = 5` without the imputation runs its own convergence check requires.

## Proposed allocation

Shinichi's instruction was that BACE is expected to be slow and that resources should be allocated
accordingly rather than the arm trimmed. So:

| work | machine | shape |
|---|---|---|
| BACE at n = 1000 (8,264 h, the bulk) | nibi + rorqual + fir arrays | one replicate per task, `--time 03:30:00`, `%400` each |
| BACE at n = 100 and n = 300 | nibi + rorqual | 5 replicates per task |
| fast arms, core + covariate sensitivity | Totoro, 62 concurrent cells (250-core allowance) | `xargs`, resume on existing rds |
| fast arms, factorial | fir + whichever cluster drains first | per-n arrays |
| AVONET300 | Totoro after the core slice | |

Wall time: **9 h if the three clusters fill (1,200 slots), about 17 h at half fill.** Totoro alone
would take 187 h, so the clusters are the plan, not a fallback. All three are bootstrapped and
smoke-proven today; rorqual's queue wait measured 2 min 39 s. Every finished (cell, seed) is skipped on
resubmission, so a cluster that stalls can hand its remainder to another without losing work.

One operational note: rorqual's `/project` allocation is at its **file-count quota** (about 499,000 of
500,000 inodes), so its R library and torch home were placed on `/home`. Results there must be
aggregated and pulled rather than left as per-cell files.

## What is still not measured

- BACE at n = 1000 end to end: the only inferred number in the budget, and the one the budget is most
  sensitive to. The first core-slice task at n = 1000 will measure it, and the estimate is re-stated
  then; if it exceeds 3.5 h the array `--time` is wrong and the tasks requeue.
- Whether `runs = 5` converges at n = 1000.
- G6, the coverage gate for arms 1 and 2, has not returned a number yet.

## The decision

Approve the launch with `runs = 5`, `n_final = 20`, the allocation above, and the 11,564 slot-hour
budget; or trim the design first. Nothing runs until then.
