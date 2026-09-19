# Handover to Claude: the phylogenetic imputation simulation study (new lane, new account)

You are Claude, picking up in the pigauto repository on a fresh account. You have no access to Shinichi's second-brain vault, to his Mac, or to the Totoro / DRAC compute servers, and you inherit no chat. This file, `AGENTS.md`, and the repository at `origin/main` are your whole context. Written 2026-09-19 by Claude Code (Shinichi's Mac session) after `main` reached d72ad06.

## Critical context

- **Mission for this lane.** Run the simulation study that compares four phylogenetic trait-imputation arms on accuracy, probability calibration and interval coverage: (1) the frequentist stack (Rphylopars for continuous-family traits, castor's Mk model for discrete traits), (2) BACE (Bayesian, joint), (3) pigauto with the GNN off, (4) pigauto with the GNN on. Two papers from one set of runs: the BACE paper reports arms 1 and 2; the pigauto paper (later) reports 1 to 4.
- **The plan is written and approved as a plan, not as a launch.** `docs/dev-log/arc/2026-09-19-imputation-simulation-plan.md` (ADEMP structure, Williams et al. 2024 self-audit, compute estimate, self-contained environment section). Shinichi shared it with Szymek for review. Do not launch the full study until the four open decisions below are answered and the pre-run has been shown.
- **Rules that bind here** (from `AGENTS.md`, which you must read first): state the regime of every number you cite (DGP, n, seeds, arms); no superlatives; a run over 30 minutes needs an explicit estimate, a pre-run with its output shown, and Shinichi's approval before the full submit; never merge a PR yourself; `r_cal = 0` (the GNN off) must stay a valid fallback; do not modify `BACE/` or package code for this lane, the runners live in `script/`.

## What was accomplished (all merged to `main`)

| PR | what | evidence |
|---|---|---|
| #180 | `gnn = FALSE` on `impute()`, `fit_pigauto()`, `multi_impute()`, `multi_impute_trees()`: the phylogenetic baseline through the whole pipeline with zero torch calls; `fit$baseline_full` (no held-out cells) for production predictions while every scorer stays on the held-out `fit$baseline`; `fit_baseline()` returns `path` | `docs/dev-log/after-task/2026-09-18-gnn-off.md`; 15 tests in `tests/testthat/test-gnn-off.R`; suite green, `R CMD check` clean |
| #181 | 200-cell with/without-GNN campaign vs raw Rphylopars and BACE (BM, OU, BACE DGP, AVONET300; n 100/300/1000; 20 seeds) | `docs/dev-log/arc/2026-09-19-campaign-gnn-off-results.md`, `script/campaign_gnn_off_results/` |
| #182 | Arc C: the AVONET300 continuous-trait gap is pigauto's in-house joint solver, not the mixed-type path; `joint_solver = "rphylopars"` inside pigauto beats raw Rphylopars by 0.11 z-RMSE (paired, 20 seeds) at 50x the fit time. The simulation plan. The per-type feasibility test (every arm on every trait type, 0 errors). castor kept as the frequentist discrete arm after beating the missForest + eigenvector hybrid | `docs/dev-log/arc/2026-09-19-avonet-gap-results.md`, `...-avonet-gap-decision-map.md`, `...-imputation-simulation-plan.md`, `script/campaign_types_results/`, `script/campaign_solver_results/` |

Headline findings, each with its regime in the linked note: on BM/OU simulations the GNN adds nothing and the shipped GNN-on loss is the held-out-cell cost of calibration; on the low-signal BACE DGP the default GNN-off safety machinery wins and raw Rphylopars falls below the mean floor; on AVONET300 raw Rphylopars beats pigauto's default baseline by about 30% because of the in-house solver; castor matches or beats the phylogeny + missForest hybrid on discrete traits; Rphylopars on log1p counts is not a count model (below the floor at n = 100).

## Current working state

- **Working:** `script/campaign_gnn_off_lib.R` (DGPs `bm_mixed`, `ou_mixed`, `bace_dgp`, `types_mixed`, `avonet`; the seeded user-level mask; scoring; the frequentist-stack arm `run_freq`; the hybrid arm `run_mf_phylo`), `script/campaign_gnn_off_cell.R` (one cell, every arm, one rds; `--smoke`), `script/campaign_gnn_off_aggregate.R`, `script/campaign_gnn_off_tables.R`, `script/campaign_gnn_off_figures.R`, `script/campaign_solver_cell.R`, `script/campaign_solver_aggregate.R`.
- **In progress:** nothing. No simulation is running anywhere.
- **Blocked on decisions (Shinichi and Szymek):** see Open questions.

## Key decisions and rationale

- `gnn = FALSE` keeps the GNN arm's safety floor and phylo-signal gate so the two arms differ only in the GNN term; the pure traditional-statistics arm is `gnn = FALSE, safety_floor = FALSE, phylo_signal_gate = FALSE`. Chosen after a Fable plan review showed a pure-BM default would lose low-signal traits to the mean corner, not to the GNN; the BACE-DGP results confirmed it.
- Two baselines rather than one refit: `fit$baseline` (held-out) feeds every scorer; `fit$baseline_full` is read only by production `predict()`. Routing a scorer to `baseline_full` is test-cell leakage.
- The default `joint_solver` stays `"inhouse"` until Shinichi decides; the plan runs both as arms 3a/3b if the decision is still open.
- castor, not TDIP's corHMM or ensemble, is the frequentist discrete arm: same model class as TDIP's phylogenetic imputer, seconds per trait, and it beat the hybrid on simulated discrete traits.
- The masks are derived from `seed + 1000` inside `make_cell()`, so every arm in a cell sees the same missing cells; keep that convention.

## Files created or modified (this lane, all on `main` at d72ad06)

R: `R/fit_pigauto.R`, `R/fit_baseline.R`, `R/predict_pigauto.R`, `R/plot_pigauto.R`, `R/impute.R`, `R/multi_impute.R`, `R/multi_impute_trees.R`, `R/evaluate.R`, `R/check_pigauto.R`, `R/report.R`; tests: `tests/testthat/test-gnn-off.R`; docs: `NEWS.md`, `AGENTS.md`, `man/*.Rd`; scripts: `script/campaign_gnn_off_{lib,cell,aggregate,tables,figures}.R`, `script/campaign_solver_{cell,aggregate}.R`, results under `script/campaign_gnn_off_results/`, `script/campaign_gnn_off_prerun/`, `script/campaign_solver_results/`, `script/campaign_types_results/`; dev-log: `docs/dev-log/arc/2026-09-18-gnn-off-contract.md`, `...-gnn-off-verify.md`, `...-campaign-prerun.md`, `...-brain-draft.md`, `2026-09-19-campaign-gnn-off-results.md`, `...-avonet-gap-decision-map.md`, `...-avonet-gap-results.md`, `...-imputation-simulation-plan.md`, `docs/dev-log/after-task/2026-09-18-gnn-off.md`, `docs/dev-log/plan-actual/2026-09-18-gnn-off.md`, `docs/dev-log/handover/2026-09-19-gnn-off-morning.md`, and this file.

## Next immediate steps (OWED, in order)

1. Run lane preflight if the tool exists in your environment (`tools/lane_preflight.sh` lives in Shinichi's vault, not in this repo; if you do not have it, read `AGENTS.md`, `git status -sb`, `git branch -a` and PR #175's description instead) and classify this handover's items against the current `origin/main`.
2. Set up the environment exactly as the plan's "Environment and reproduction" section says (R 4.4+, pigauto from `main`, torch, Rphylopars, castor, BACE from github.com/daniel1noble/BACE, MCMCglmm, missForest). Prove it with `Rscript script/campaign_gnn_off_cell.R --dgp types_mixed --n 60 --seed 1 --out /tmp/smoke --arms gnn_on,gnn_off,gnn_off_pure,freq,bace,floor --smoke` (seconds; proves the invocation, never results), then one real cell at n = 100.
3. Get the four decisions from Shinichi (or "use your judgment", in which case: 3a/3b both solvers; covariates as sensitivity only; phyloglm Poisson for counts in the frequentist stack; lambda on the tree covariance as written).
4. Add to the runner what the plan still lacks: the lambda and rho factors in `make_dgp` (currently the simulated DGPs are lambda = 1, rho = 0), a MAR mechanism and the phylogenetically biased missingness mechanism in `make_cell`, calibration metrics (Brier, ECE) and interval coverage for every arm in `score_arm` (pigauto's conformal interval is already scored; add Rphylopars' `anc_var` interval and BACE's draw percentiles, which needs `n_final = 20` and keeping the draws), and a phyloglm Poisson route for counts in `run_freq`. Each addition gets a `--smoke` run before anything else.
5. Pre-run (under 30 minutes): one replicate per core cell at n = 100 and n = 1000, all arms, with Gelman-Rubin diagnostics on BACE's two chains. Write the output and the re-stated estimate to `docs/dev-log/arc/<date>-simulation-prerun.md`. Stop and show Shinichi.
6. Only after his approval: the core slice, then the factorial. Compute is on Shinichi's side (Totoro or a DRAC job array; you cannot reach them); your part is the runner, the driver script pattern in `script/campaign_gnn_off_cell.R`'s header comment, and the aggregation.

## Blockers and open questions

- Arm 3 solver default (in-house vs Rphylopars inside pigauto): Shinichi.
- Covariates: in the core slice or sensitivity only; if in, arm 1 becomes phylolm / phyloglm / castor and arm 3 needs a covariate-aware non-GNN route: Shinichi and Szymek.
- Count traits in the frequentist stack: phyloglm Poisson or the hybrid.
- Phylogenetic signal parameterisation to agree with Szymek.
- BACE chain length: fix by convergence on the pre-run, then re-estimate compute.
- TDIP's own ensemble and GAIN (github.com/Matgend/TDIP): untested; include only if the package installs cleanly.
- Compute access: you have none. The full study (about 3,600 core-hours) runs from Shinichi's side.

## Gotchas and failed approaches

- `sigma_method = "fisher_ml"` on the in-house solver fell back to single-pass in 200 of 200 cells: it does not repair the AVONET gap.
- Rphylopars on log1p-transformed counts, back-transformed, fell below the mean floor at n = 100; do not present it as a count model.
- BACE refuses non-ultrametric trees: simulated trees come from `ape::rcoal`, never `rtree`.
- BACE's `bace()` signature is `fixformula, ran_phylo_form, phylo, data, nitt, burnin, thin, runs, n_final, ovr_categorical`; older bench scripts in `script/` used a wrong signature and silently skipped BACE for months.
- `impute(gnn = FALSE)` ignores `covariates` with a warning; under `gnn = TRUE` they enter only through the GNN.
- The default `gnn = FALSE` arm spends its time in `calibrate_gates()` (117 s at n = 1000); the pure arm takes 2 s.
- `skip_on_cran()` in the tests fires under plain `Rscript`; set `NOT_CRAN=true` when running test files directly.
- Never run `xargs -P` on a shared machine without a thread cap per process (`OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS` explicit).
- A detached driver launched through an ssh heredoc can leave the wrapper alive with the script text in its command line; a `pgrep -f` waiter keyed on that text then never fires.

## Landing state

All work of this lane is merged to `origin/main` (d72ad06). `CARRIED-OVER`: none from this lane. The repository also carries many older branches with unpushed commits (`codex/*`, `experiment/*`, `feature/*`, `analysis/*`, `arc/pertype-benches`, `chore/dedupe-agent-files`); they belong to earlier lanes and are `PROTECTED`: not this lane's, not touched, not to be pushed or rebased by you. PR #175 (`evidence/gnn-sentinel-prerun`) is another live lane's; do not edit its files.

FINDING-OF-RECORD: `gnn = FALSE` semantics and the two-baseline contract; the campaign result (GNN adds nothing on BM/OU simulations; the shipped loss is the held-out-cell cost); the arc C result (in-house solver owns the AVONET gap; Rphylopars solver inside pigauto beats raw Rphylopars). vault-note: not yet written (pigauto's rule is to stage and propose); the staged draft is `docs/dev-log/arc/2026-09-18-brain-draft.md`, and the findings themselves are on `main` in the notes listed above.

## Mission control

| repo | main | CI | shipped this lane | plan by leverage |
|---|---|---|---|---|
| itchyshin/pigauto | d72ad06 | R-CMD-check green on 92edc73 (last package change); pkgdown red since 2026-08-28, pre-existing | `gnn = FALSE` + `baseline_full`; 200-cell campaign; arc C solver diagnosis; simulation plan; per-type feasibility | 1. decisions and pre-run; 2. core slice; 3. factorial; 4. arm-3 solver default as a package change (separate lane) |

## How to resume

Environment: a clone of github.com/itchyshin/pigauto at `main`, R 4.4 or later, packages per the plan's environment table; `NOT_CRAN=true` for tests; no credentials needed for anything in this lane. Safe verification: `Rscript -e 'devtools::test(filter = "gnn-off")'` and the `--smoke` invocation above. Do not stage `.unlazy/`, `LOOP/`, `results*/`, `logs*/`, or anything under `docs/` that is a build output (the `docs/` directory is git-ignored except force-added dev-log files; add dev-log files with `git add -f`).

One-command resume, run from the repo root in your own terminal:

```text
Read AGENTS.md and docs/dev-log/handover/2026-09-19-claude-handover.md. Run the handover rehydration steps, reconcile them with the current git state, then continue only the OWED Next Immediate Steps.
```
