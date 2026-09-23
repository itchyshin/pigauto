# Session Handoff: imputation simulation v1 closed; v2 redo is next

Meta: 2026-09-23 · from Claude Code (Opus 5.5, session 84cc2f15) · to Claude · branch `arc/imputation-sim`
at `601628b` plus this commit · draft PR https://github.com/itchyshin/pigauto/pull/184

You are Claude, picking up a finished simulation arc whose redo (v2) has been ordered but not planned.
Nothing from this arc is running. Read this document, then the v1 summary, before touching anything.

## Critical Context

1. **v1 is closed; do not re-run or extend it.** Shinichi stopped it on 2026-09-23 to redo the study with
   longer MCMCglmm settings. Its results stand as the record: every cell present in both stages (G10, G11
   PASS), acceptance ledger all met (26 of 26), publication decisions recorded (G13c). The authoritative
   summary and the v2 requirements are in `docs/dev-log/arc/2026-09-23-simulation-v1-summary.md`.
2. **v2 depends on another lane that is running now.** Shinichi's own session holds the lambda-default lane
   (brain D-278) in `../pigauto-lambda-default`, branch `feat/joint-lambda-default`: pigauto's joint baseline
   now estimates Pagel's lambda by default (commit `b8eb384`), and its core benchmark is running on Totoro
   (about 140 R processes writing `~/pigauto_sim/results/core_lambda_core_fast`). **Do not touch that
   worktree, its branch, or those Totoro processes.** v2's pigauto arms must use that baseline once it lands.

## Goals and plan

- **Mission (repo AGENTS.md / CLAUDE.md):** measure where pigauto stands against the frequentist stack and
  BACE on mixed-type phylogenetic imputation, with every number carrying its regime and MCSE. Two consumers:
  Dan Noble's BACE paper (frequentist vs BACE) and pigauto's own article (all arms).
- **v2 plan (ordered, not yet written):** (a) wait for the lambda-default lane to land; (b) pre-run BACE on
  low-signal cells over a grid of `runs` 5/10/15 by `nitt` 50k/100k, measuring BACE's own convergence-pass
  rate and wall time; (c) choose settings, re-derive the budget, and present both to Shinichi before any
  campaign (D-139: a run over 30 minutes needs a shown pre-run and his approval); (d) rerun with
  `freq_lambda` in from the start.

## What Was Accomplished (v1)

- Core slice (18 cells x 200 replicates, BACE 100) and factorial (56 cells x 200, BACE 100 at n = 100 and
  30 to 100 at n = 1000), AVONET300 case study, covariate sensitivity. Six arms plus a floor.
- Aggregation with paired MCSE, cross-host dedupe, pooled ECE, and a divergence rule (finite-but-absurd fits
  floored and counted as `n_divergent`).
- Deliverables: private results board (https://claude.ai/artifact/M5HtGRnNGfwsK2Se4gMX24, v10, all seven
  decisions marked); BACE-paper methods note; pkgdown article (unlisted, unpublished); after-task report;
  plan-vs-actual reconcile; the v1 summary.
- Headline results, in the v1 summary: the ranking inverts with phylogenetic signal; Rphylopars'
  evolutionary model (BM default vs `lambda`) matters more than the method; BACE's real advantage is
  discrete traits at low signal, at a cost of 16 to 63% failures and a 15 to 20% convergence-pass rate at
  lambda 0.3; only pigauto's conformal intervals reach 95% coverage, and only at n >= 300.

## Current Working State

- **Working:** everything on `arc/imputation-sim` is committed and pushed; the article renders from
  `script/campaign_sim_results/*.csv`.
- **In progress:** nothing in this arc. The lambda-default lane (Shinichi's) is running its benchmark.
- **Blocked / waiting on humans:** Shinichi reads the board in full, then says whether #184 goes ready for
  review. Szymek signs off on the corrected Pagel-lambda parameterisation, then the article goes public.

## Key Decisions and Rationale

Recorded in the ledger (`.unlazy/imputation-sim/gates/leaf-results.md`, G13c), the after-task (section 3a)
and brain `memory/DECISIONS.md` D-278:

- The BACE paper reports both frequentist specifications with the gap as a finding, and leads with the
  discrete result with failure and convergence costs beside it (Shinichi, 2026-09-23).
- The article stays unlisted until Szymek signs off (Shinichi, 2026-09-23).
- PR #184 stays a draft until Shinichi has read the board; **never merge it yourself**.
- D-278: the joint baseline estimates lambda by default (Shinichi, 2026-09-22).
- BACE at n = 1000 in the factorial used at least 30 replicates, not 100 (Shinichi, 2026-09-21, budget).

Standing constraints: never edit `R/`, `BACE/`, PR #175 files, or the dirty `handover/2026-08-09-cursor`
main checkout; never trigger Duo (use the `~/.ssh/cm-*` ControlMaster sockets); Totoro is shared (cap 250
cores for snakagaw, and read per-user RSS before launching: another user held 930 GB on 2026-09-22); never
compute on a DRAC login node; `ScheduleWakeup` is forbidden in this project; no brain-vault write without
Shinichi's approval.

## Landing State

`handoff_gate.sh` run 2026-09-23: **"GATE INCONCLUSIVE -- 5 acceptance ledger(s) exist but could not be run"**,
because the gate's shell could not find `node`; the same ledger read by `gate-check.mjs --status` from an
interactive shell returns **"ALL MET (26 met, 6 never executed)"**. The six never-executed gates are
evidence-only or exceed the checker's 120 s cap and were run directly; their evidence lines say so.

| Artifact / branch | Committed | Pushed | PR | State |
|---|---|---|---|---|
| `pigauto` `arc/imputation-sim` (worktree `../pigauto-imputation-sim`), HEAD = this commit | y | y | #184 draft | LANDED; not merged (Shinichi's call after he reads the board) |
| `.unlazy/imputation-sim/**` (goal, progress, ledger) | y (tracked on this branch) | y | #184 | LANDED; goal file marked COMPLETE |
| Results board source `scratchpad/sim-results.html` (session scratchpad) | n | published v10 | none | CARRIED-OVER: the live page is the durable copy; to change it, `Artifact` read the URL, edit, republish with `url` |
| Raw rds: Totoro `~/pigauto_sim/results/*`, nibi `~/projects/def-snakagaw/snakagaw/pigauto_sim/results/{core,factorial}`, fir `~/pigauto_sim/results/*`; Mac working pools `/tmp/pig_pool4`, `/tmp/pig_pool6`, `/tmp/pig_pool8` | n (data, never git) | n/a | none | CARRIED-OVER: keepers on Totoro and nibi; the Mac `/tmp` pools are disposable |
| `feat/joint-lambda-default` (worktree `../pigauto-lambda-default`) | y | y | none yet | **PROTECTED**: Shinichi's lane; not yours |

FINDING-OF-RECORD: pigauto's joint baseline assumed lambda = 1 and so sat at the floor at low signal, where BACE and Rphylopars at model = "lambda" did not  vault-note: [[DECISIONS#D-278|D-278]]

The remaining v1 results live on branch `arc/imputation-sim` (PR #184, canonical source:
`docs/dev-log/arc/2026-09-23-simulation-v1-summary.md`) and the private board. A fuller vault distillation of
them has not been approved.

## Files Created or Modified (session diff vs `origin/main`)

- `script/campaign_gnn_off_lib.R`, `campaign_gnn_off_cell.R`, `campaign_gnn_off_aggregate.R`, `campaign_sim_cell.R`,
  `campaign_sim_design.R`, `campaign_sim_checks.R`, `campaign_sim_prerun_summary.R`
- `script/campaign_sim_totoro.sh`, `campaign_sim_nibi_array.sh`, `campaign_sim_pool.sh`
- `script/campaign_sim_page_data.R`, `campaign_sim_page_build.sh`, `campaign_sim_page.template.html`
- `script/campaign_sim_results/{summary,paired,per_trait,failures,ece}.csv`; `script/campaign_gnn_off_prerun/bm_mixed_n100_s1.rds`
- `docs/dev-log/arc/2026-09-20-simulation-prerun.md`, `2026-09-20-simulation-methods-bace.md`,
  `2026-09-22-imputation-sim-after-task.md`, `2026-09-22-joint-lambda-default-plan.md`, `2026-09-23-simulation-v1-summary.md`
- `docs/dev-log/handover/2026-09-20-imputation-sim-handover.md`, `2026-09-22-claude-handover.md`, and this file
- `docs/dev-log/plan-actual/2026-09-20-imputation-sim-reconcile.md`
- `vignettes/articles/simulation-study.Rmd`, `.Rbuildignore`
- `.unlazy/imputation-sim/{GOAL,GATES,PROGRESS}.md`, `.unlazy/imputation-sim/gates/leaf-*.md`

No snapshot pointer was edited: `AGENTS.md`/`CLAUDE.md` on main carry no Live Phase Snapshot, and four lanes
are live in this repo, so a single pointer would orphan the others.

## Other live lanes (do not collide)

| lane | where | owner | what |
|---|---|---|---|
| lambda-default (D-278) | `../pigauto-lambda-default`, `feat/joint-lambda-default` | Shinichi's Claude session | joint baseline estimates lambda by default; benchmark on Totoro |
| mondrian-realdata | lease `claude:pigauto-mondrian-realdata` on `R/fit_helpers.R`, `R/predict_pigauto.R`, `R/multi_impute.R`, `R/fit_pigauto.R`, tests, `useful/`, `docs/dev-log/`, `NEWS.md` | another Claude lane | Mondrian conformal on real data |

Both hold paths in `R/`, which this arc never edits. There is no coordination board in the repo
(`lane_preflight.sh` reports none), so read `lane_lease.sh --list pigauto` before claiming anything.

## Next Immediate Steps

1. Run `~/shinichi-brain/tools/lane_preflight.sh "$PWD"` and `~/shinichi-brain/tools/lane_lease.sh --list pigauto`;
   state which lane you take. Classify each item below as OWED, DONE, RETRACTED or PROTECTED against the
   current git state.
2. **OWED:** read `docs/dev-log/arc/2026-09-23-simulation-v1-summary.md` ("Why v2 is needed") and the
   lambda-default lane's latest commits (read-only: `git -C ../pigauto-lambda-default log --oneline -5`).
3. **OWED, after the lambda-default lane has landed (ask Shinichi if unsure):** write the v2 BACE settings
   pre-run plan as `docs/dev-log/arc/<date>-simulation-v2-prerun-plan.md`: cells (low-signal core cells at
   n = 100 and 300, plus one clade-masked n = 1000), grid `runs` 5/10/15 x `nitt` 50k/100k (burnin and thin
   scaled with `nitt`), metrics (BACE `converged` rate, drift, median ESS, wall, z-RMSE), a time estimate,
   and the compute target. **Stop and show it to Shinichi; do not launch.**
4. **PROTECTED:** PR #184 (do not merge or mark ready without Shinichi's word), the lambda-default worktree,
   the Totoro `core_lambda_core_fast` processes, the `R/` files held by other lanes.

## Blockers and Open Questions

- Whether v2 waits for the lambda-default lane to merge, or starts with the BACE pre-run in parallel (only
  the BACE arm is independent of it). A drafted question for Shinichi: *"Start the BACE settings pre-run now,
  in parallel with the lambda lane, since BACE does not depend on pigauto's baseline? yes / no"*.
- Szymek's sign-off on the corrected Pagel-lambda form.

## Gotchas and Failed Approaches (all in the v1 summary's "Operational lessons")

- Resume-skip is per host: seed a host's results directory from the pool before any recovery array, or it
  recomputes the campaign.
- Never aggregate a pool while an rsync into it is running; pools use per-host subdirectories; dedupe on
  (filename, arm set).
- The Mac's `rsync` is 2.6.9 (no `--info=stats1`); `git add -f` next to `git push` in one command trips the
  destructive-command hook (use `--force` on `git add`, or split the commands).
- `gate-check.mjs --approve` has a 120 s per-check cap and unchecks slow gates; run G10/G11 directly.
- Clade-masked BACE at n = 1000 needs a 12 h Slurm limit; nibi caps a user at 1,000 submitted array tasks;
  some BACE replicates segfault or never return; the runner cannot floor a crash (three hand records in v1).

## How to Resume

Working directory: `/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim` (branch `arc/imputation-sim`).
Toolchain: R 4.6 on the Mac (`Rscript`), `rmarkdown` to render the article; Totoro, nibi and fir through
`ssh -o ControlPath=~/.ssh/cm-snakagaw@<host>:22 -o ControlMaster=no <host>`. Safe verification command:
`Rscript -e 'rmarkdown::render("vignettes/articles/simulation-study.Rmd", output_file = "/tmp/check.html", quiet = TRUE)'`.
Never stage anything in the main checkout `../pigauto` (dirty, another lane's).

```text
Read AGENTS.md and docs/dev-log/handover/2026-09-23-claude-handover.md. Run the handover rehydration steps, reconcile them with the current git state, then continue only the OWED Next Immediate Steps.
```

## Mission control

| repo | branch | CI | what shipped | next by leverage |
|---|---|---|---|---|
| pigauto | `arc/imputation-sim`, PR #184 draft | local checks only (render OK; ledger ALL MET) | v1 simulation: both stages complete, board v10, methods note, article, after-task, v1 summary | 1. v2 BACE pre-run plan (after or beside the lambda lane) · 2. Shinichi reads board, #184 ready · 3. Szymek sign-off, article public |
| pigauto | `feat/joint-lambda-default` (Shinichi's lane) | benchmark running on Totoro | lambda-estimating joint baseline as default | PROTECTED; v2's pigauto arms wait on it |
