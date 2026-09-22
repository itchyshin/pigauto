# Session Handoff: four-arm imputation simulation, closing the campaign

Meta: 2026-09-22 (Mac clock ~10:20 MDT) · from Claude Code (Fable 5.1, session 84cc2f15) · context healthy
· platform: Claude Code only; Codex/Cursor not involved in this arc.

## Critical Context

1. **Everything durable lives on disk, not in chat.** Re-read, in this order:
   `.unlazy/imputation-sim/GOAL.md` (arc state lines), `.unlazy/imputation-sim/PROGRESS.md` (append-only
   run log; the last entries are the resume point), `.unlazy/imputation-sim/gates/leaf-*.md` (acceptance
   ledger with CHECK/EXPECT/EVIDENCE), then `docs/dev-log/arc/2026-09-22-imputation-sim-after-task.md`.
   `.unlazy/` is git-ignored by design and exists only in this worktree.
2. **Resume-skip is per host.** A recovery array on any host must be seeded with every finished rds from the
   pool first, or it recomputes the campaign (nibi array 22471569 was cancelled for exactly this on 09-22).
3. **The aggregator is the single authority for every reported number** (`script/campaign_gnn_off_aggregate.R`,
   divergence rule included). The article, the methods note and the board all read the same committed csv.
   Do not report a number from a per-replicate script unless it has been checked against the aggregator.
4. **Publication is gated on Shinichi reading the board** (G13c). The pkgdown article is unlisted and
   unpublished; nothing goes public until he records the six decisions on the board's last tab.

## What Was Accomplished

- Core slice (18 x 200, every arm), factorial (56 x 200, every fast arm; BACE 98.4%), AVONET300 (20 seeds,
  every arm) computed on the corrected design; covariate sensitivity at 3,599 of 3,600 plus its freq_lambda wave.
- Aggregation with paired MCSE, cross-host dedupe, the missing-fraction key, coverage rows, pooled ECE, and
  the divergence rule; committed csv under `script/campaign_sim_results/`.
- (a) Private results Artifact v3, https://claude.ai/artifact/M5HtGRnNGfwsK2Se4gMX24; (b) BACE-paper methods
  note with core and factorial results; (c) four-arm pkgdown article rendering both; after-task draft;
  plan-vs-actual Reconcile 2; draft PR #184 current.
- Gates: G14 PASS x3, G12 PASS, G13a/b/d PASS. G10/G11 fail only on BACE tails now computing; G6d on covsens.

## Current Working State

- Working: all three deliverables render/publish from the committed aggregate; PR #184 (draft) is current.
- In progress (machines): Totoro BACE tails (24 core + 57 factorial seeds, pgids 1840413 / 1840886, ~3 h from
  10:18); nibi seeding then three recovery arrays as backup; fir arrays PD with no start estimate; the two
  patched aggregations' ECE stage (failures.csv with n_divergent, ece.csv) still running on the Mac; fir
  covsens being pulled so its last seed can run on Totoro.
- Not working / blocked: fir scheduler (3,000+ nodes drained); the in-app browser cannot show me the board.

## Key Decisions & Rationale

See after-task section 3a and the DECISION RECEIPT in
`docs/dev-log/plan-actual/2026-09-20-imputation-sim-reconcile.md`. Headline ones: report the frequentist
stack at both `model = "BM"` and `model = "lambda"`; divergence rule (3x floor or interval score > 1e3) scored
at the floor and counted; BACE 30 replicates at n = 1000 in the factorial (Shinichi 09-21); in-house solver
stays pigauto's default. Standing constraints: never edit `R/`, `BACE/`, PR #175 files or the dirty
`handover/2026-08-09-cursor` checkout; never merge; never trigger Duo; Totoro <= 250 cores for snakagaw;
`ScheduleWakeup` forbidden in this project; no publish before G13c.

## Landing State

Run `~/shinichi-brain/tools/handoff_gate.sh pigauto` before trusting this table; refresh it at close.

| Artifact / branch | Committed | Pushed | PR | State |
|---|---|---|---|---|
| `pigauto` `arc/imputation-sim` `077b4de` (worktree `../pigauto-imputation-sim`) | y | y | #184 draft | LANDED on the branch; NOT merged (merge is Shinichi's call after G13c) |
| `.unlazy/imputation-sim/**` (goal, progress, ledger) | n (git-ignored by design) | n | none | CARRIED-OVER: lives only in the worktree; resume = read the files above |
| Results board source `scratchpad/sim-results.html` | n (scratchpad) | published v3 | none | CARRIED-OVER: republish the same path to update; source also mirrored in PROGRESS notes |
| Raw rds pools `/tmp/pig_pool4`, `/tmp/pig_pool6`, `/tmp/pig_pool7` on the Mac; `~/pigauto_sim/results/*` on Totoro; `results/` on nibi and fir | n (data, never git) | n/a | none | CARRIED-OVER: keepers are on Totoro and on the clusters' /project; the Mac copies are working pools |
| `docs/dev-log/handover/2026-09-22-claude-handover.md` (this file) | pending | pending | #184 | to be committed at close |

FINDINGS-OF-RECORD: the simulation results themselves are on the branch and in the private board; a vault
distillation (3-8 durable bullets, WHAT-WORKS candidates in after-task section 11) is proposed and
awaits Shinichi's approval per the brain-write boundary. Until then: **FINDINGS-OF-RECORD: on branch
`arc/imputation-sim` and PR #184, not yet in the vault.**

## Next Immediate Steps

1. When Totoro's BACE tails finish (watcher biu95j9dm; logs/core.log DONE >= 3 and logs/factorial.log DONE
   >= 4): rsync `results/core_bace` and `results/factorial_bace` into `/tmp/pig_pool4/core/totoro_bace`
   and `/tmp/pig_pool6/factorial/totoro_bace`; re-run G10 and G11 (direct `Rscript`, not through
   gate-check's 120 s timeout); cancel the now-redundant nibi and fir recovery arrays.
2. When the ECE stage lands: copy `agg_failures.csv`/`agg_ece.csv` (core, pool7) and `fact_*` (factorial,
   pool6/agg) into `script/campaign_sim_results/` as `failures.csv` (concatenate) and `ece.csv`; re-render
   the article; commit.
3. Covsens: seed Totoro `results/covsens` from `/tmp/pig_pool6/covsens/fir`, run the covsens stage (fast
   arms) to fill the last seed, pull `covsens` and `covsens_fl` to `/tmp/pig_pool6/covsens/{totoro,totoro_fl}`,
   aggregate with the patched script, add one paragraph each to the methods note and the article, board v4,
   G6d.
4. Re-aggregate core and factorial once the BACE tails are pooled so the committed csv are final; re-render;
   commit; refresh PR #184 body.
5. Finish after-task sections 5 and 10 with the final gate lines; run
   `python3 ~/shinichi-brain/tools/closeout.py check <abs path>` from the worktree; run
   `gate-check.mjs --reverify --approve` on every leaf; commit.
6. Refresh this handover's Landing State from `handoff_gate.sh`; commit; push.
7. Wait for Shinichi: read the board, record the six decisions (G13c), decide article visibility and merge.

## Blockers / Open Questions

- G13c: Shinichi has not yet read the board or recorded the publication decisions.
- Szymek's sign-off on the corrected Pagel-lambda parameterisation.
- fir's scheduler state (informational; Totoro carries the tails).

## Gotchas & Failed Approaches

- Mac `rsync` 2.6.9 has no `--info=stats1`; a grep hid the usage error and an empty pull looked complete.
- The aggregator drops any path containing `_smoke`; name smoke pools without it.
- `gate-check.mjs` runs CHECK lines with CWD = the gates directory and a 120 s timeout; CHECK lines are now
  anchored with `cd '<worktree>' &&`, and G11 (39k rds) must be run directly.
- `closeout.py check` needs an absolute report path and scans `./.unlazy` of its CWD; run it from the worktree.
- G12 lists a flat directory; point it at a host subdirectory of a pool, not the pool root.
- Never aggregate a pool while an rsync into it is running (the 08:37 race on 09-21).
- BACE at n = 1000 with clade masking needs a 12 h Slurm limit, not 5 h.
- Totoro clock and the Mac clock disagree by minutes; PROGRESS labels are ordered, not synchronised.

## How to Resume

```
cd "/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim"
~/shinichi-brain/tools/lane_preflight.sh "$PWD"
cat .unlazy/imputation-sim/GOAL.md; tail -80 .unlazy/imputation-sim/PROGRESS.md
ssh totoro 'grep -c DONE ~/pigauto_sim/logs/core.log; grep -c DONE ~/pigauto_sim/logs/factorial.log; ls ~/pigauto_sim/results/core_bace | wc -l; ls ~/pigauto_sim/results/factorial_bace | wc -l'
```
Then follow "Next Immediate Steps" from wherever the counts say the machines are.
