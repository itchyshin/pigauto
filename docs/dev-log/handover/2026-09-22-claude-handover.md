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

Validators at close (2026-09-23 02:25), quoted rather than summarised:

- `gate-check.mjs --status` over the five leaves: **"UNMET: 1 (met: 25, 6 never executed)"**, the one being
  `leaf-results:G13c` (Shinichi's six publication decisions). The "6 never executed" are the evidence-only gates
  (G0, G6c, G6d, G9c, Galign, G13c) plus checks whose CHECK exceeds the checker's 120 s cap (G11, 39k rds) and were
  run directly; their evidence lines say so.
- `handoff_gate.sh`: **"GATE FAIL -- 5 acceptance ledger(s) have UNMET gates"**. It counts the never-executed
  gates above as unmet and G13c is genuinely open, so this is the expected reading until Shinichi records his
  decisions; every unlanded item is declared in the table below.
- `closeout.py check`: structure PASS ("after-task structure check passed"); its acceptance-ledger stage halts on
  unmet gates in *other* projects' ledgers under the vault's `.unlazy/` (book-format, brain-campaign) and on
  G13c here. Not this arc's to fix.

`handoff_gate.sh` first run 2026-09-22 14:35: PR #184 OPEN; branch current; acceptance ledger 25 of 26 gates met after the tails landed on 2026-09-23 (G10, G11 PASS); the one open gate is
G13c (Shinichi's six publication decisions). G6 and G9b were resolved by the campaign's
own measurements and are marked so in the ledger. Two unrelated branches (`shannon-install`, `spec/vulcan-gpu-avonet9993`)
carry unpushed commits from other lanes and are not this arc's.

| Artifact / branch | Committed | Pushed | PR | State |
|---|---|---|---|---|
| `pigauto` `arc/imputation-sim` (worktree `../pigauto-imputation-sim`), HEAD at close | y | y | #184 draft | LANDED on the branch; NOT merged (merge is Shinichi's call after G13c) |
| `.unlazy/imputation-sim/**` (goal, progress, ledger) | n (git-ignored by design) | n | none | CARRIED-OVER: lives only in the worktree; resume = read the files above |
| Results board source `scratchpad/sim-results.html` | n (scratchpad) | published v3 | none | CARRIED-OVER: republish the same path to update; source also mirrored in PROGRESS notes |
| Raw rds pools `/tmp/pig_pool4`, `/tmp/pig_pool6`, `/tmp/pig_pool7` on the Mac; `~/pigauto_sim/results/*` on Totoro; `results/` on nibi and fir | n (data, never git) | n/a | none | CARRIED-OVER: keepers are on Totoro and on the clusters' /project; the Mac copies are working pools |
| `docs/dev-log/handover/2026-09-22-claude-handover.md` (this file) | pending | pending | #184 | to be committed at close |

FINDINGS-OF-RECORD: the simulation results themselves are on the branch and in the private board; a vault
distillation (3-8 durable bullets, WHAT-WORKS candidates in after-task section 11) is proposed and
awaits Shinichi's approval per the brain-write boundary. Until then: **FINDINGS-OF-RECORD: on branch
`arc/imputation-sim` and PR #184, not yet in the vault.**

## Next Immediate Steps

Items 2, 3 and 5 of the earlier list are DONE (bootstrap csv committed 4e61cc7; covsens complete and reported;
after-task finalised). Remaining:

1. DONE 2026-09-23: both BACE tails landed (factorial on nibi, core on nibi after Totoro was ruled out by
   another user's 930 GB process), G10 and G11 PASS, both stages re-aggregated, csv rebuilt, article
   re-rendered, board republished, PR #184 body refreshed. Keepers: nibi `results/core` and `results/factorial`
   hold every BACE rds including the three hand-assembled failure records.
2. Wait for Shinichi: read the board, record the six decisions (G13c), decide article visibility and merge.
3. Then run `closeout.py check` from the worktree (it passes only once G13c is closed) and the
   `handoff_gate.sh`; commit.

## Ordered follow-on lane (not this arc)

D-278 (2026-09-22): estimate Pagel's lambda in the joint baseline and make it the default; plan stub
`docs/dev-log/arc/2026-09-22-joint-lambda-default-plan.md`; opens after PR #184 merges, own worktree and branch.

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
