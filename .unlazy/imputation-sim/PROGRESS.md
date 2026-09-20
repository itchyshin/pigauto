# PROGRESS RECORD — pigauto four-arm imputation simulation study

**This file is the single source of truth for resuming. Read it first, update it after every
meaningful step, and never start a second copy of work already running.**

Status: **IN PROGRESS** (not complete). Last updated 2026-09-20 by the Claude session that started the lane.

---

## 1. The goal (approved by Shinichi 2026-09-20)

Run to completion the ADEMP simulation study comparing four phylogenetic trait-imputation arms, then
produce three deliverables. Full approved plan: `/Users/z3437171/.claude/plans/distributed-snacking-turtle.md`
(read it; it carries the corrected design, the slice table and the acceptance ledger).

Arms: **1** frequentist stack (Rphylopars joint BM + castor Mk + phyloglm Poisson for counts) ·
**2** BACE (Bayesian, MCMCglmm) · **3a** pigauto GNN off, in-house solver · **3b** pigauto GNN off,
`joint_solver = "rphylopars"` · **4** pigauto GNN on · plus a mean/mode floor.

Primary contrast, pre-registered before any result: **BACE vs the frequentist stack on z-RMSE and
95% coverage (interval score beside it), core slice, pooled over trait types.** Everything else is
secondary.

## 2. Completion criteria (all must hold)

- [ ] Core slice complete: 18 cells x 200 replicates (BACE on seeds 1..100), failure rate recorded per cell.
- [ ] Trimmed factorial complete: 56 cells, same conditions.
- [ ] AVONET300 case study, all arms, 20 seeds.
- [ ] Covariate sensitivity on the 18 core cells.
- [ ] Aggregated with paired MCSE on every number; regime columns filled.
- [ ] **S7a** results Artifact shown to Shinichi (Dan's freq-vs-BACE view + the four-arm view); his
      publication decisions recorded here.
- [ ] **S7b** BACE-paper methods write-up (`docs/dev-log/arc/<date>-simulation-methods-bace.md`).
- [ ] **S7c** unlisted pkgdown article (`vignettes/articles/simulation-study.Rmd`), published only
      AFTER Shinichi reads the Artifact.
- [ ] After-task report, Melissa plan-vs-actual reconcile, handover, PR opened (never merged by an agent).

## 3. Where everything lives

| what | where |
|---|---|
| Worktree / branch | `/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim`, branch `arc/imputation-sim` (from `origin/main` fd0b513) |
| Approved plan | `/Users/z3437171/.claude/plans/distributed-snacking-turtle.md` |
| Acceptance ledger | `.unlazy/imputation-sim/gates/leaf-{env,runner,prerun,campaign,results}.md` |
| Goal card | `.unlazy/imputation-sim/GOAL.md` |
| Reviews that shaped the design | `/private/tmp/claude-503/-Users-z3437171-Dropbox-Github-Local-pigauto/84cc2f15-470d-4aaf-82bc-9553e6d0c894/scratchpad/review-{rose,emmy,fisher}.md`, `s3-verify.md` (may be purged; findings are folded into the plan and the commits) |
| Runner | `script/campaign_gnn_off_lib.R`, `script/campaign_sim_cell.R`, `script/campaign_sim_checks.R`, `script/campaign_sim_design.R` |
| Drivers | `script/campaign_sim_totoro.sh`, `script/campaign_sim_nibi_array.sh` (works on nibi, rorqual, fir; `HALF=A|B` splits cells across two clusters) |
| Summariser | `script/campaign_sim_prerun_summary.R` |
| Totoro workspace | `snakagaw@totoro.biology.ualberta.ca:~/pigauto_sim/` (script, prerun_fast, prerun_bace, results, logs) |
| nibi workspace | `~/projects/def-snakagaw/snakagaw/pigauto_sim/` (R lib on /project) |
| rorqual workspace | `/project/def-snakagaw/snakagaw/pigauto_sim/`, **R lib and torch home on `/home`** (its /project is at its ~500k inode quota) |
| fir workspace | `/home/snakagaw/pigauto_sim` (bootstrap PASS 2026-09-20: env + compute-node smoke job 60661032; its `/project` is at 500K/500K files so ROOT, R lib and torch home all live on `/home`). Always `export PIG_SIM_ROOT=/home/snakagaw/pigauto_sim`. |

## 4. Compute rules (binding)

- Connect ONLY via the existing ControlMaster sockets: `ssh -o BatchMode=yes -o ConnectTimeout=15 snakagaw@<host>`.
  A fresh interactive login triggers Duo 2FA, which is forbidden. Totoro needs no Duo.
- **Totoro: up to 250 cores for snakagaw** (Shinichi, 2026-09-20, raising this lane above the standing
  150 default) = 62 concurrent cells at 4 threads each.
- DRAC: `--account=def-snakagaw_cpu`. Never compute on a login node. Size `--time` from `seff`, not guesses.
- One process per cell; BACE's own `n_cores` stays 1. Never wrap MCMCglmm in `mclapply` (it segfaults).
- `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`; torch threads via
  `torch_set_num_threads()`, not an env var.
- Resume is free: a (cell, seed) whose rds exists is skipped, so re-submitting an array only runs what is missing.

## 5. Completed so far

1. **S1** worktree, branch, lane lease (`claude:pigauto:84704`), acceptance ledger written before dispatch.
2. **S2** environments green and smoke-proven on **Mac, Totoro, nibi, rorqual** (gates G1a, G1b, G1c, G7, G8,
   plus rorqual's own pair). pigauto 0.11.0 everywhere. **fir PASS too** (login-node install clean, smoke job
   60661032 COMPLETED exit 0 on compute node fc30557, 6/6 arms). All five hosts are green.
3. **S3** runner rebuilt to the corrected design (commits b58980b, 4ccc6ac, 80dc686, 3b62c78, 9268c16).
   Gates **G2, G3, G4, G5, Gold PASS**. G12 passed on a small real directory.
4. **S3-verify** independent Opus review; four blocking findings all fixed (see section 7).
5. **S4 pre-run** on Totoro: n = 100 complete for every arm; n = 1000 complete for the fast arms.

## 6. Measured numbers (regime: types_mixed, MCAR 0.30, one replicate, 4 threads, seconds per fit)

| arm | n = 100 | n = 1000 |
|---|---:|---:|
| frequentist stack | 0.35 | 2.0 |
| pigauto GNN off (3a) | 0.7 | 38 |
| pigauto GNN off + Rphylopars solver (3b) | 27 | **260** |
| pigauto GNN on (4) | 134 | 319 |
| BACE (runs = 2, n_final = 20) | 835 | **> 4,500, killed at 75 min** |
| floor | 0.0 | 0.0 |

Arm 3b at n = 1000 was the headline unknown: **260 s, cheaper than the GNN arm** — it is affordable.
BACE dominates the budget (roughly 84% of it).

## 7. Findings that changed the design (do not re-derive)

- BACE `n_final` counts **full imputation runs**, each refitting MCMCglmm per response — not chain
  draws. A chain-length formula set it to 400, about 80x the intended cost. Now a budget constant
  `PIG_BACE_NFINAL` (default 50; the pre-run used 20).
- BACE `runs` are **sequential imputation iterations** for its own `assess_convergence()`, not parallel
  chains. Gelman-Rubin across them is meaningless (it read 5.7 on a converged fit, 20 with random
  effects pooled in). `assess_convergence()` has **min_iterations = 3**, so `runs = 2` can never be
  assessed — the cause of `converged = FALSE` on all 12 pre-run cells. Default is now
  `PIG_BACE_RUNS = 10` (its vignette uses 5 for demos, 15 with nitt = 100000 for a real analysis).
  `skip_conv` stays TRUE to keep per-cell cost bounded; the verdict is recorded and the convergence
  RATE is reported as a result.
- Cells store BACE's own verdict, the imputed-mean drift over the last two iterations, and effective
  sample size (median gated, minimum reported: threshold/categorical models mix slowly by design).
- macro-F1 now averages only over classes present in the masked truth, and records `n_classes_scored`.
- The MAR driver `d1` carries its own loading (0.35, the largest keeping Sigma positive definite with
  seven traits) so MAR stays MAR at rho = 0. Verified: masked and observed truth differ by 0.38 SD.
- BACE fails on some cells ("mixed model equations singular"), notably clade-masked n = 1000. Failures
  are scored at the floor and REPORTED, never dropped.
- Never wrap MCMCglmm in `mclapply` (SIGSEGV). Parallelism is one process per cell.

## 8. NEXT UNFINISHED STEPS (in order)

**CURRENT STATE (checked 2026-09-20 13:50 MDT): G0 APPROVED. EVERY HOST IS SATURATED AND HEALTHY.
DO NOT RELAUNCH ANYTHING. Measured this check:**

| host | job | tasks | state |
|---|---|---:|---|
| Totoro | core fast arms, 62 slots | 3600 | 139 rds landed, 127 procs, load 144 of the 250-core cap, **~29 results/min, ETA ~2.1 h** |
| nibi | core BACE n=100 (22339703), n=300 (22339823), factorial BACE n=100 (22340187) | 964 | 468 RUNNING, 496 PENDING |
| rorqual | core BACE n=1000 (21475732) | 600 | all PENDING, not yet started |
| fir | factorial BACE n=1000 HALF=A (60661338) | 1356 | 194 RUNNING (started 13:47), 1162 PENDING |

A landed core cell was opened and checked: 5 arms, 144 result rows, `failed` and `errors` both empty,
walls recorded (gnn_on 173 s, gnn_off_rphylopars 36 s at n = 100, both within ~30% of the pre-run under
62-way contention). No error strings anywhere in the Totoro logs.

**The n = 1000 BACE wall check is still BLOCKED: nothing n = 1000 has finished.** rorqual has not started
a single task and fir's array is three minutes old. That measurement gates the remaining submissions,
so the next run should check it first.

**Submission capacity is the binding constraint, not approval.** DRAC MaxSubmit is ~1000 array tasks per
user per cluster (fir accepted 1400, so its limit is higher). Right now nibi is at 964 and rorqual at 600.
The queued item 1 below needs 1400 tasks at BLOCK=1, which fits nowhere; at BLOCK=2 it is 700 tasks and
would fit rorqual only after its current 600 drain. Do not force it in: a refusal costs nothing but a
half-submitted array does.

**WAS: at the G0 gate. The pre-run note is written, committed and shown to Shinichi
(`docs/dev-log/arc/2026-09-20-simulation-prerun.md`, commit ffd6d85). The pre-run is complete:
16/16 fast cells and 6 BACE cells at n = 100 on Totoro, plus the runs 3/5/10 convergence probe.
Re-derived budget 11,564 slot-hours, BACE 83%. Proposed settings: runs = 5, n_final = 20, with
BACE n = 1000 on nibi/rorqual/fir arrays and fast arms on Totoro. NOTHING LAUNCHES until Shinichi
approves. If he has not answered, do not re-run the pre-run and do not launch; just say it awaits him.**

1. ~~Finish the pre-run note~~ DONE 2026-09-20. (Former step 1 text kept below for reference.)
   ~~**Finish the pre-run note**~~ `docs/dev-log/arc/2026-09-20-simulation-prerun.md`: the table in section 6,
   BACE convergence-vs-runs probe results (`~/pigauto_sim/logs/probe_runs.log` on Totoro: runs 3, 5, 10 x
   2 seeds), failures, and a **re-derived budget** from measured walls x the design
   (`Rscript script/campaign_sim_design.R --stage core|factorial`). Force-add it (`git add -f`, docs/ is ignored).
2. **G0: STOP and show Shinichi that note.** Do not launch the core slice or factorial without his approval.
   He has already said BACE is expected to be slow and to allocate resources accordingly, so the note should
   propose an allocation rather than a trim.
3. After G0 only: **S6a** core slice on Totoro (62 slots), **S6b** factorial split across nibi / rorqual / fir
   (`HALF=A|B`, per-n arrays, `--time` from `seff` + 30%), **S6c** AVONET300, **S6d** covariate sensitivity.
4. MECHANICAL-VERIFY (cell counts, failures, cross-host reproducibility), then **S7** aggregate.
5. **S7a** Artifact -> Shinichi decides -> **S7b** methods note, **S7c** pkgdown article.
6. **S8** after-task, Melissa reconcile, handover, PR.

## 8b. PENDING CLUSTER SUBMISSIONS (submit as capacity frees; DRAC MaxSubmit = 1000 job/array tasks per user per cluster)

Submitted and running as of 2026-09-20 13:50:
- Totoro: core fast arms, 3600 jobs, 62 slots -> `~/pigauto_sim/results/core`
- nibi: core BACE n=100 (22339703), n=300 (22339823), factorial BACE n=100 (22340187)
- rorqual: core BACE n=1000 (21475732, BLOCK=1, --time 03:30:00)
- fir: factorial BACE n=1000 HALF=A (60661338, BLOCK=1, --time 03:30:00)

STILL TO SUBMIT (both were refused with AssocMaxSubmitJobLimit; retry when that cluster's
`squeue -u snakagaw -h | wc -l` drops well below 1000):
1. factorial BACE n=1000 **HALF=B** -> whichever of rorqual / fir / nibi has room:
   `export PIG_SIM_ROOT=<that cluster's root>; cd $PIG_SIM_ROOT && HALF=B ARMS_BACE=bace ARMS=bace SEEDS=bace BLOCK=1 THROTTLE=200 bash script/campaign_sim_nibi_array.sh factorial 1000 03:30:00`
2. factorial FAST arms, n=100 and n=1000 (arms gnn_on,gnn_off,gnn_off_rphylopars,freq,floor; no SEEDS cap):
   `ARMS_BACE=gnn_on,gnn_off,gnn_off_rphylopars,freq,floor ARMS=$ARMS_BACE BLOCK=5 THROTTLE=200 bash script/campaign_sim_nibi_array.sh factorial 100 01:00:00` and the same for `1000 02:00:00`.
   These can also run on Totoro once its core fast arms finish.
3. AVONET300 case study (stage `avonet`) and the covariate sensitivity: Totoro, after the core slice.

Cluster roots: nibi `~/projects/def-snakagaw/snakagaw/pigauto_sim`; rorqual
`/project/def-snakagaw/snakagaw/pigauto_sim`; **fir `/home/snakagaw/pigauto_sim`** (its /project is
at 500K/500K files, so everything there lives on /home). Always `export PIG_SIM_ROOT=` on rorqual and fir.

## 9. Pauses that require Shinichi (never proceed past these alone)

- **G0** after the pre-run note (launching the campaign is the irreversible compute commitment).
- **S7a** after the Artifact, before anything is published.
- Never merge a PR. Never edit `R/`, `BACE/`, PR #175 files, or the dirty `handover/2026-08-09-cursor` checkout.

## 10. Run log (append one line per scheduled run)

- 2026-09-20 ~13:40 — record created; pre-run measured; awaiting probe results then G0 note.
- 2026-09-20 ~14:10 — pre-run COMPLETE (16/16 fast, 6 BACE at n=100, convergence probe runs 3/5/10).
  Note written and committed (ffd6d85). Budget re-derived: 11,564 slot-hours. **AT G0, awaiting Shinichi.**
- 2026-09-20 13:42 — **G0 APPROVED ("Go ahead"). Core slice LAUNCHED**: Totoro fast arms (3600 jobs),
  nibi BACE n=100/300 (22339703, 22339823), rorqual BACE n=1000 (21475732). runs=5, n_final=20.
- 2026-09-20 13:50 — factorial BACE n=1000 HALF=A on fir (60661338); factorial BACE n=100 on nibi
  (22340187). HALF=B refused on nibi and rorqual (MaxSubmit 1000); queued in section 8b for retry.
- 2026-09-20 13:50 — scheduled check, no new compute launched (all four hosts saturated).
  Verified Totoro core slice healthy (139 rds, 5 arms, zero failures, ETA ~2.1 h) and confirmed fir's
  bootstrap PASS. n=1000 BACE wall still unmeasurable: rorqual 600 tasks all PENDING, fir array 3 min old.
  Section 8b items 1 and 2 stay queued; no cluster has room.
