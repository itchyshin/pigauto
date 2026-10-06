# After-task: the four-arm imputation simulation lane

2026-09-20, Claude Code, branch `arc/imputation-sim` off `origin/main` fd0b513.
Status at writing: the campaign is **running**, not finished. This report covers the lane that built
and launched it; the results close-out follows when the cells land.

## 1. Goal

Run the ADEMP simulation study comparing four phylogenetic trait-imputation arms to completion, then
produce a private results board, a methods write-up for the BACE paper, and an unlisted pkgdown
article. Approved plan: `/Users/z3437171/.claude/plans/distributed-snacking-turtle.md`.

## 2. Implemented

- Runner extended to the corrected design: Pagel-lambda latents with cross-trait correlation, OU as a
  sensitivity, fixed population thresholds, an always-observed MAR driver, three missingness
  mechanisms including clade-biased, and per-arm interval, calibration and interval-score scoring.
- A new cell wrapper (`campaign_sim_cell.R`), a design table (`campaign_sim_design.R`), gate checks
  (`campaign_sim_checks.R`), cross-machine pooling (`campaign_sim_pool.sh`), a pre-run summariser,
  and the results-board build chain (`campaign_sim_page_data.R`, `_page_build.sh`, `.template.html`).
- Drivers for Totoro (`xargs`, resume, thread caps) and for DRAC arrays (nibi, rorqual, fir), with
  arm, seed-range and per-host overrides.
- Environments bootstrapped and smoke-proven on Mac, Totoro, nibi, rorqual and fir.
- Pre-run executed and its note published; campaign launched after G0.
- Deliverables (b) and (c) drafted; deliverable (a) published privately, marked partial.

## 3a. Decisions and Rejected Alternatives

- **Both arm-3 solvers run**, not one. The pre-run showed the Rphylopars solver costs 260 s at
  n = 1000, under the GNN arm, so the arc C question can be settled with data rather than assumed.
- **BACE kept at full scope** rather than trimmed, on Shinichi's instruction that it is expected to be
  slow and should be resourced accordingly. Rejected: cutting its replicates or its n = 1000 cells.
- **`skip_conv = TRUE` with the verdict recorded**, rather than letting BACE retry. Retrying makes
  per-cell cost unbounded across thousands of cells; the convergence rate becomes a reported result.
- **Median rather than minimum effective sample size** as the gate. MCMCglmm threshold and categorical
  models mix slowly by construction, so the minimum fails on a healthy fit.
- **Factorial trimmed** of n = 300 and lambda = 0.7, which remain in the core slice. Rejected: the full
  108-cell design, which spends compute on filling a grid rather than on the contrasts.
- **rorqual retired** from the campaign rather than repaired: its project allocation is at its
  file-count quota and its nodes are slower than fir's.

## 4. Files Touched

`script/campaign_gnn_off_lib.R`, `campaign_gnn_off_cell.R`, `campaign_gnn_off_aggregate.R`,
`campaign_sim_cell.R`, `campaign_sim_checks.R`, `campaign_sim_design.R`, `campaign_sim_pool.sh`,
`campaign_sim_totoro.sh`, `campaign_sim_nibi_array.sh`, `campaign_sim_prerun_summary.R`,
`campaign_sim_page_data.R`, `campaign_sim_page_build.sh`, `campaign_sim_page.template.html`;
`vignettes/articles/simulation-study.Rmd`; `.Rbuildignore`;
`docs/dev-log/arc/2026-09-20-simulation-prerun.md`, `2026-09-20-simulation-methods-bace.md`.
No change to `R/`, `BACE/`, tests, or another lane's files.

## 5. Checks Run

Gates G1a/G1b/G1c (three environments), G2 (six-arm smoke), G3 (mask integrity and rownames), G4
(Pagel lambda recovered to 0.26 / 0.68 / 1.00 against targets 0.3 / 0.7 / 1.0 over 20 replicates at
n = 1000), G5 (realised missing fraction within 0.001 of target for MAR and clade), G7 and G8 (Totoro
and nibi smokes on compute nodes), Gold (byte-level no-regression against a pre-change reference,
max abs diff 3.9e-11), G9a and G9c (pre-run complete, arm 3b measured, budget re-derived), G12
(aggregate carries MCSE and regime columns), G10 exercised and correctly refusing on partial data.

## 6. Tests of the Tests

Three gates were found to pass or fail for the wrong reason and were rewritten:

- G3 originally asserted that scored cells equal the mask, which `score_arm()` makes true by
  construction. It now asserts the NA pattern of the masked frame against the mask, and each arm's
  rownames against the truth.
- G9b originally computed Gelman-Rubin across BACE's `runs`. Those are sequential imputation
  iterations, not parallel chains, so it read 5.7 on a fit BACE itself called converged. It now reads
  BACE's own verdict plus effective sample size.
- G6 was scoped to n = 1000 with a forked worker; MCMCglmm segfaults inside a fork. Now sequential
  and at n = 300, with its regime printed alongside the number.

G10 was checked by running it against an incomplete directory and confirming it refuses.

## 7a. Issue Ledger

An independent Opus review raised four blocking findings, all fixed before compute: BACE `n_final`
set from chain length (about 80x the intended cost), macro-F1 deflated by absent classes, the MAR
driver degenerate at rho = 0, and Gelman-Rubin as the wrong convergence diagnostic. A scope review
raised six more, of which the material ones (pre-run not sampling the factorial regimes, the article
breaking `R CMD check`, tautological gates, no AVONET slice) are fixed or scheduled.

## 8. Consistency Audit

The design in the plan, the runner code, the methods note and the article agree on the corrected
lambda parameterisation, the fixed thresholds, the driver loading, the mechanisms and the arm list. A
symbolic-alignment table mapping each design term to its implementing line was produced and reviewed.

## 9. What Did Not Go Smoothly

- The first pre-run wave duplicated the Bayesian arm through a wrong environment variable and burned
  an hour before it was killed.
- rorqual lost 379 tasks in six seconds to a disk quota and 44 more to a wall-clock limit.
- Pooling with `rsync --ignore-existing` silently discarded every BACE cell, because the fast arms and
  the Bayesian arm write the same filename for the same cell. Caught by checking the pooled count.
- Two agents converged on one file when a reviewer's fixes landed while a builder was still working;
  resolved by standing the builder down.

## 10. Known Residuals

- G6, the coverage gate for arms 1 and 2, has not returned a number.
- Whether `runs = 5` converges at n = 1000 is unmeasured; only n = 100 was probed.
- The factorial BACE half B runs with two replicates per task, a shape not yet timed end to end.
- rorqual holds nine salvaged cells that must be included when pooling.

## 11. Team Learning

Read the tool's own documentation before inventing a diagnostic: BACE defines convergence over its
imputation iterations and exposes a verdict, and an hour of Gelman-Rubin work was wasted before the
vignette was consulted. Second, measure the task shape before scaling it: `seff` showed 25% CPU
efficiency and 200 out-of-memory kills, and fixing both quadrupled throughput per core-hour.

## 12. Cross-Product Coverage

pigauto only. The BACE package is used, never modified. No other repository is touched.

## Does NOT cover

The results themselves, the completeness gates against pooled data, AVONET300, the covariate
sensitivity, publication of the article, and the merge. Those follow the campaign.
