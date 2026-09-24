# Session Handoff: freq vs BACE under proper multiple imputation and Rubin's rules

Meta: 2026-09-24 · from Claude Code (session 84cc2f15, which ran simulation v1) · to a fresh Claude session
· branch `arc/rubin-freq-bace` (pushed) · worktree `/Users/z3437171/Dropbox/Github Local/pigauto-rubin-freq-bace`

You are Claude, opening a new lane. Nothing in it has been built or run yet: it holds a plan, a goal file and
the v1 runner it was branched from. Read this document first; it stands alone.

## Where this lane lives

| what | value |
|---|---|
| folder | `/Users/z3437171/Dropbox/Github Local/pigauto-rubin-freq-bace` (a git worktree of the pigauto repo) |
| branch | `arc/rubin-freq-bace`, on GitHub as `origin/arc/rubin-freq-bace`; no PR yet |
| based on | `arc/imputation-sim` at `14b919d` (simulation v1: runner, design table, aggregator, drivers) |
| plan | `docs/dev-log/arc/2026-09-24-rubin-freq-bace-plan.md` |
| goal file | `.unlazy/rubin-freq-bace/GOAL.md` |
| v1 record | `docs/dev-log/arc/2026-09-23-simulation-v1-summary.md` (results, corrections, operational lessons) |

If the desktop app does not list this folder, open it directly: File, Open Folder, choose the path above, or
in a terminal `cd "/Users/z3437171/Dropbox/Github Local/pigauto-rubin-freq-bace" && claude`.

## Critical Context

1. **Why this lane exists.** Simulation v1 scored per-cell prediction intervals, not multiple imputation, and
   measured no downstream estimand. Two faults surfaced after it closed:
   - BACE's v1 interval was the 2.5 and 97.5 percentile of 20 imputed datasets. That construction covers a new
     draw only 0.872 of the time under a correct model; the Rubin interval from the same 20 draws,
     mean ± `qt(0.975, M - 1) * sqrt((1 + 1/M) * B)`, covers 0.950 (200,000 simulated sets).
   - **BACE's continuous imputations are posterior means, not predictive draws.** In `BACE/R/model_functions.R`,
     `.predict_bace()` takes `pred_prob[,1]` ("Extract posterior mean") from `.pred_cont()`, the posterior of the
     fitted values `X beta + Z u`; no residual is drawn. So BACE's multiple imputation is improper for continuous
     traits. The direction (under-coverage, worse at low lambda) is inferred; the size is not measured.
2. **Scope.** Frequentist stack vs BACE only, for Dan Noble's BACE paper. No pigauto arm: Shinichi is building
   a pigauto posterior in another lane, which can plug into these estimands later. **Never edit `R/` or `BACE/`.**

## Decisions already taken (Shinichi, 2026-09-24, "yes yes yes")

- Q1: add a **BACE + residual draw** arm, built in our runner from BACE's returned model objects, beside BACE as shipped.
- Q2: frequentist MI = **parametric bootstrap then draw** as the headline (proper), with **draws from one fit** as the
  improper contrast.
- Q3: downstream estimands = **phylogenetic-GLS slope of c2 on c1** and **their correlation**, on the v1 core cells.
- Compute: **Totoro + DRAC** (Totoro for the frequentist arms and the BACE settings pre-run; nibi arrays for the BACE
  campaign, 12 h limit at n = 1000).

## Next Immediate Steps

1. Run `~/shinichi-brain/tools/lane_preflight.sh "$PWD"` and claim the lane:
   `LANE_ID='claude:pigauto:rubin' ~/shinichi-brain/tools/lane_lease.sh --claim pigauto --paths script/rubin_,docs/dev-log/arc/2026-09-24-rubin,.unlazy/rubin-freq-bace/`
   (this session released its lease, so the claim is free). Classify each step below OWED / DONE / RETRACTED / PROTECTED.
2. **OWED, step 1 of the plan (Mac, minutes):** write `script/rubin_lib.R` + `script/rubin_cell.R`, reusing
   `make_cell()`, `sim_latents()` and `score_arm()` from `script/campaign_gnn_off_lib.R`. Smoke one cell at
   n = 60 with M = 20: freq A, freq B, BACE as shipped, BACE + residual draw; per-cell Rubin intervals, PGLS slope,
   correlation; check the Rubin pooling against a hand computation. Tests or a smoke script ship with the code.
   Implementation hint for freq A: from each bootstrap refit's parameters (phylogenetic and residual covariance,
   lambda, means), build the joint covariance `kron(Sigma_p, C_lambda) + kron(Sigma_e, I)` and draw the missing
   cells of the original data jointly from their conditional normal; that keeps cross-trait and cross-species
   dependence, which the downstream estimands need.
3. **OWED, step 2:** write the BACE settings pre-run as a ready-to-launch plan (grid `runs` {5, 10, 15} x
   `nitt` {50k, 100k}, core cells at lambda 0.3 and 0.7, n = 100 and 300, 5 seeds; about 720 fits, roughly 2 h on
   Totoro) with a time estimate. **Stop and show Shinichi before launching anything over 30 minutes (D-139).**
4. **PROTECTED:** `R/`, `BACE/`, PR #184, the `arc/imputation-sim` and `feat/joint-lambda-default` worktrees,
   and any Totoro processes that are not this lane's.

## Other live lanes (do not collide)

| lane | where | notes |
|---|---|---|
| simulation v1 (closed) | `../pigauto-imputation-sim`, `arc/imputation-sim`, PR #184 draft | finished; do not re-run; merge is Shinichi's call |
| lambda-default (D-278) | `../pigauto-lambda-default`, `feat/joint-lambda-default` | Shinichi's session; its Totoro benchmark used about 140 cores on 2026-09-23 |
| pigauto posterior | another lane of Shinichi's | owns `R/`; this lane never touches it |
| mondrian-realdata | lease on several `R/` files | another Claude lane |

Read `~/shinichi-brain/tools/lane_lease.sh --list pigauto` before claiming.

## Landing State

| Artifact / branch | Committed | Pushed | PR | State |
|---|---|---|---|---|
| `arc/rubin-freq-bace` at `92ca468` plus this handover commit | y | y | none | LANDED; no code yet, plan and goal only |
| `.unlazy/rubin-freq-bace/GOAL.md` | y | y | none | LANDED |

FINDINGS-OF-RECORD: BACE's continuous imputations are posterior means of the fitted values with no residual draw,
so its multiple imputation is improper for continuous traits  vault-note: none yet; canonical source is
`docs/dev-log/arc/2026-09-24-rubin-freq-bace-plan.md` on this branch (a vault note needs Shinichi's approval).

## Gotchas carried from v1

- Totoro is shared: read `free -g` and per-user RSS before every launch; cap at 150 cores (D-143).
- nibi caps a user at 1,000 submitted array tasks; clade-masked BACE at n = 1000 needs a 12 h limit.
- Resume-skip is per host: seed a host's results directory from the pool before any recovery array.
- Pools use per-host subdirectories and (filename, arm set) dedupe; never aggregate while an rsync runs.
- A segfault or timeout leaves no rds; the runner cannot floor it (record by hand and label it).
- Attach to hosts only through the `~/.ssh/cm-*` ControlMaster sockets; never trigger Duo; never compute on a
  DRAC login node. `ScheduleWakeup` is forbidden in this project. Mac `rsync` is 2.6.9 (no `--info=stats1`).
- `git add -f` beside `git push` in one command trips the destructive-command hook; use `--force` or split them.

## How to Resume

```text
Read AGENTS.md and docs/dev-log/handover/2026-09-24-claude-handover.md. Run the handover rehydration steps, reconcile them with the current git state, then continue only the OWED Next Immediate Steps.
```

## Mission control

| repo | branch | CI | what shipped | next by leverage |
|---|---|---|---|---|
| pigauto | `arc/rubin-freq-bace` | none yet | plan, decisions, compute plan | 1. step 1 smoke on the Mac · 2. step 2 pre-run plan, shown to Shinichi · 3. campaign after approval |
