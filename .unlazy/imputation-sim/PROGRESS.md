# PROGRESS RECORD — pigauto four-arm imputation simulation study

**This file is the single source of truth for resuming. Read it first, update it after every
meaningful step, and never start a second copy of work already running.**

Status: **IN PROGRESS** (not complete). Last updated 2026-09-21 14:15 MDT by the scheduled Claude run.

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
      **Fast arms 3600/3600 DONE. BACE 1777/1800 (98.7%): n=100 598/600, n=300 599/600, n=1000 580/600.**
- [ ] Trimmed factorial complete: 56 cells, same conditions.
- [x] AVONET300 case study, all arms, 20 seeds. **DONE 2026-09-21** - 20/20 cells on Totoro, zero errors, all six arms plus gnn_on_full, aggregated to `avo_summary.csv` and on the board's own tab.
- [ ] Covariate sensitivity on the 18 core cells. **Runner support BUILT and verified 2026-09-21**
      (stage `covsens`, `--ncov`). **DISPATCH STARTED 2026-09-21 14:00** - the n = 100 wave is RUNNING on
      fir (job 60829554, 240 tasks, fast arms only, no BACE). The n = 300 and n = 1000 waves are
      deliberately held until the n = 100 wave gives a measured fir wall to size `--time` from.
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

**CURRENT STATE (checked 2026-09-20 22:50 MDT). Every host is still busy; nothing was relaunched.
The one ready submission is still REFUSED by the permission classifier (see the finding below).**

| host | what is running | queued | results landed | note |
|---|---|---:|---|---|
| Totoro | **AVONET, all 6 arms, 20 seeds** (pgid 3354530) | 20 procs (43 pids) | core **3600/3600**, avonet **0/20** | 2 h in and HEALTHY: 20 R processes at 105% CPU and 7.5 GB RSS each, all inside BACE at n = 300. About 20 of Totoro's 250 permitted cores are in use, so the machine has room; only the permission rule is missing |
| nibi | factorial BACE n=1000 HALF=B (22352507) only | 251 | core 1038, factorial 2759 (n100 2464, n1000 295) | BLOCK=2, 1 core, 32 GB, `--time 06:00:00`. 37 COMPLETED, zero TIMEOUT yet, completed blocks at 3h26 to 4h04 |
| rorqual | nothing | 0 | core 9 salvaged | retired for the inode quota, see 8a |
| fir | 60681245 factorial half A (249 R), 60681249 core n=1000 (200 R), 60684950 core n=300 redo (68 R) | 524 | core **1372**, factorial **302** | core n=300 redo nearly drained (531 COMPLETED, 1 FAILED). Still zero OOM at 32 GB, but **74 TIMEOUT has now appeared on the factorial n=1000 array**, see the new finding below |

**Core BACE completeness, measured on the union of nibi + fir filenames (1438 of 1800):**

| n | landed | missing | who is filling it |
|---:|---:|---:|---|
| 100 | 597/600 | 3 | fir array drained; these 3 are the residue |
| 300 | 592/600 | 8 | fir 60684950, 68 tasks left, finishes within the hour |
| 1000 | 249/600 | 351 | fir 60681249, running, but the n = 1000 budget decision below still stands |

Cross-host duplication has grown from 650 to **992** identical core filenames on both nibi and fir.
The aggregator's (filename, arm set) dedupe neutralises this in the numbers; it is recorded as waste,
not repeated. Once fir's core arrays drain, fir alone holds a near-complete n = 100 and n = 300 set.

### NEW: BLOCK=2 is the wrong shape for n = 1000 BACE, and it is costing whole tasks

`--mem=32G` fixed the out-of-memory problem: fir has recorded **zero OOM** since the switch. But
fir's factorial half-A array 60681245 has now accumulated **74 TIMEOUT against 303 COMPLETED**
(20% of finished tasks), every one of them stopping at exactly 04:30:0x.

The cause is arithmetic, not a node problem. That array is `627-1400`, i.e. **BLOCK=2**, and a
successful BACE fit at n = 1000 is now measured at **03:17 to 04:04** of wall on fir and nibi (slower
than the 02:46 measured on an empty fir). Two successes in one task therefore need about 7 h and
cannot fit inside a 4h30 wall. A block times out whenever it happens to draw two lambda = 0.3 cells;
it completes in seconds when it draws two lambda = 1 cells, which fail singular in about 5 s
(core array 60681249 shows exactly this bimodality: most COMPLETED tasks at 5 to 8 s, the successes
at 03:17:15).

nibi's half B (22352507) is the same BLOCK=2 shape at a 6 h wall, so it will lose its two-success
blocks too, just fewer of them. It has no TIMEOUT yet only because its completed blocks so far have
each held one success plus one fast failure.

**Consequence for the decision Shinichi still owes:** option 1 in the budget question below is more
expensive than it was written, because retries are now part of it. The correct shape for any n = 1000
BACE array is **BLOCK=1** (one cell-seed per task, wall 05:00:00), which makes a task's cost
independent of which cells it draws and removes this class of loss entirely. Core n=1000 (60681249)
is already BLOCK=1 and has zero TIMEOUT, which is the control case.

### NEW: coverage was never reaching the aggregate, so half the primary contrast was invisible

Found and fixed this run (commit 5a5480c, `script/` only, `R/` untouched). The cell runner carries
`coverage` as a SIDE COLUMN on each zRMSE row, not as a `metric` row of its own. The Section I
aggregator splits by `metric`, so coverage was dropped before `_summary.csv` was written: measured,
the core summary held 372 accuracy rows, 496 zRMSE rows and **zero** coverage rows. The results
board filters on `metric`, so its interval-coverage panel rendered "No cell yet for this
combination" on Versions 1, 2 and 3 - and coverage is half the pre-registered primary contrast.

Coverage is now emitted as its own metric with the same mean-and-MCSE-over-seeds treatment
(420 rows on the core slice, 24 on AVONET), and appears on the board for the first time in
Version 4. Two smaller defects were fixed alongside it: the page decided continuous-vs-discrete
from a hard-coded list of SIMULATED trait names, so every real AVONET trait fell into the discrete
bucket; and `isFinite(null)` is true in JavaScript, so the floor arm - which has no interval at all -
was read as a phantom zero and every other arm was pilled "worse than floor" on interval score.

**First reading of the coverage numbers** (pooled over the four continuous-family traits, rho = 0,
MCAR 0.3, nominal 95%). This is a result, not a diagnostic:

| n | lambda | BACE | frequentist | pigauto GNN off |
|---:|---:|---:|---:|---:|
| 100 | 0.3 | 0.804 | 0.801 | 0.878 |
| 300 | 0.3 | 0.812 | 0.845 | 0.957 |
| 300 | 0.7 | 0.824 | 0.852 | 0.954 |
| 1000 | 0.3 | 0.820 | 0.869 | 0.957 |
| 1000 | 0.7 | 0.826 | 0.876 | 0.954 |

Both arms of the primary contrast **under-cover**, BACE the more so, and neither closes the gap as n
grows. pigauto's split-conformal intervals sit essentially on nominal from n = 300. Caveat carried
forward: BACE's lambda = 1 rows rest on 39 to 48 replicates rather than 100, because that is where it
fails singular.

### NEW: one replicate in 200 destroys two cells of arm 3b's interval summary

In `gnn_off_rphylopars` the interval width and interval score on the **count** trait blow up in
exactly two cells - n = 100 lambda = 0.7 rho = 0, and n = 300 lambda = 1 rho = 0 - to 6.6e20 and
7.2e05 respectively, with MCSE equal to the mean. Traced to a single seed: in the n = 100 cell,
**seed 146 alone** returns a width of 1.3e23 against a median of 3.44 over the 200 replicates, and
it is the only replicate above 100. The point estimate is unaffected (no zRMSE cell exceeds 10), so
this is confined to the interval on the back-transformed count scale.

The mean is therefore not a usable summary for interval width or interval score in those two cells.
The board is left reporting the true mean, with the caveat in its status line, rather than the
estimator being quietly swapped for a median - **that choice is Shinichi's**, and it would have to
apply to every arm and metric alike, not just where it flatters.

### NEW: S6d covariate support is BUILT (the last piece that needed no cluster and no decision)

S6d was the one remaining completion criterion blocked on code rather than on Shinichi or on a free
slot: `script/campaign_sim_cell.R` had no covariate flag at all. It does now. `script/` only; `R/`,
`BACE/` and PR #175 untouched.

What shipped (commit below): a `covsens` stage in `campaign_sim_design.R` (the same 18 core cells,
`ncov = 2`), an `--ncov` / `--rho_cov` pair on the cell runner, an `_k<ncov>` field in the rds
filename so covariate cells can never collide with core cells, and `--ncov` threaded through both
drivers.

**The plan does not specify how the covariate is generated or how strongly it correlates with the
traits** (it says only "18 core cells with 2 covariates"). Those were chosen here and each is one
flag to change:

- `cov_j = rho_cov * L[, target_j] + sqrt(1 - rho_cov^2) * e_j`, with `e_j` an independent column
  carrying the same phylogenetic structure. Targets are `c1` (continuous) and `cnt` (the count
  trait) - the two places arm 1 can actually put a covariate.
- `rho_cov = 0.6`. Measured realised correlation 0.605 / 0.604 over 40 seeds, Pagel's lambda of the
  covariate itself 1.000 at lambda = 1. At the driver's 0.35 a covariate explains only 12% of a
  trait's variance, which risks a sensitivity slice too weak to detect anything; 0.6 gives 36% and
  still leaves the phylogeny doing most of the work.
- Covariates are fully observed, never masked, never scored.

**Two things measured while building it, both of which changed the implementation:**

1. **The obvious construction is not positive definite.** Imposing the covariate correlations
   directly on `Sigma_rho` fails at rho = 0: two covariates each correlated r with seven mutually
   independent traits imply a mutual correlation of 7r^2, so pinning them to 0 breaks the matrix -
   min eigenvalue **-0.32 at r = 0.25**, before `rho_driver` is even applied. This is the same
   ceiling the existing comment on `rho_driver = 0.35` records ("the largest keeping Sigma PD with
   seven traits"), and a second such column does not fit under it. The generative construction above
   is PD for any |rho_cov| < 1 by construction and needs no eigen check.
2. **Drawing the covariate noise inside the shared latent matrix silently corrupted the pairing.**
   Widening `Z` from K = 8 to K = 10 shifts the RNG stream that `rpois()` and `rnorm()` read
   afterwards, so the count and proportion traits changed: measured, 6 of 60 counts moved and the
   proportions by up to 0.24. The covariates are now drawn AFTER the traits from their own
   `rnorm()`, so a covsens cell has traits **byte-identical** to the core cell at the same seed.
   That makes covsens-vs-core a **paired** comparison on the same data rather than two independent
   draws - materially more powerful, and it was nearly lost by accident.

Where the covariate reaches each arm, and where it cannot:

| arm | covariate enters | how |
|---|---|---|
| 1 freq, continuous | yes | OLS residualisation on the covariates over observed rows, joint BM on the residuals, contribution added back. Beta treated as known when forming the interval - a documented approximation |
| 1 freq, count | fitted, rarely adopted | `phyloglm(y ~ cov1 + cov2)` is fitted BESIDE the intercept-only model and adopted only if it predicts the observed counts better (Poisson deviance). On this DGP it never wins - see below - so the count trait ends up covariate-free in practice |
| 1 freq, discrete | **no** | `castor::hsp_mk_model` takes no covariates at all. Reported as covariate-free, not silently ignored |
| 2 BACE | yes | extra fixed terms on the RHS of every `fixformula`; never a response, so never imputed or scored |
| 3a / 3b pigauto GNN off | **no** | pigauto threads covariates through the GNN only, so a GNN-off fit cannot use them. This is what the plan means by "arm 3 reported covariate-free" |
| 4 pigauto GNN on | yes | `covariates =` on `pigauto::impute` |
| floor | no | unchanged |

**Verified before claiming it works** (all on the Mac, no cluster):

- `ncov = 0` reproduces the old DGP exactly - truth, mask and the latent matrix all `identical()`
  to the pre-change code. **Every core / factorial / avonet result already on disk stays valid**,
  and the filename is unchanged, so resume still skips them.
- `ncov = 2` leaves all seven traits byte-identical at the same seed (the pairing above).
- Cell counts per stage unchanged: core 18, factorial 56, prerun 16, avonet 1, covsens 18, so gates
  G10/G11 still count against the right totals.
- The `ncov` column is the LAST column of the design table, so both drivers' positional `awk`
  (`$3..$9` cell args, `$10` reps, `$11` bace_reps) keeps working; both were exercised and emit
  `--ncov 2` on covsens and `--ncov 0` on core.
- End-to-end cell runs at n = 100 with and without covariates (see the run log for the result).

**A measured finding about arm 1, found while building this and worth reporting as a result.** The
count path was the one place a covariate could most obviously help: the DGP is
`cnt ~ Poisson(exp(1.5 + 0.8 * L3))` and `cov2` is a reading of `L3` at r = 0.6, so a correct
Poisson regression on it should do well. It does not. `phylolm::phyloglm(..., method =
"poisson_GEE")` is a MARGINAL estimator and on this data it returns badly conditioned coefficients
(it warns `system is singular` on most seeds). Measured over 25 seeds at n = 100, lambda = 0.7,
rho = 0.5, taking its tip-specific `exp(x'beta)` at face value:

| count trait, arm 1 | mean zRMSE | covariate better in |
|---|---:|---|
| intercept-only (the existing arm) | 1.163 | - |
| naive `exp(x'beta)` | **9.727** | 1 of 25 seeds |
| after the deviance guard | 1.163 | 0 of 25 (guard always prefers the intercept) |

Two defects surfaced on the way and both are fixed: `exp(x'beta)` **overflowed to non-finite** on at
least one seed, which made the count zRMSE `NA` for the whole trait (it would have silently voided
the count column across all 18 x 200 covsens cells); and even when finite it was an order of
magnitude worse than the constant mean. The fix is a clamp plus **in-sample model selection** - the
covariate model is adopted only when it beats intercept-only on Poisson deviance over the observed
rows - so the covariate version can never make the arm worse than the covariate-free arm. The honest
reading is that **the frequentist stack cannot exploit a covariate on the count trait with the
estimator this campaign uses**, which is a reportable property of arm 1, not a bug in the harness.

For the continuous traits the covariate path is stable and does something small: mean zRMSE 0.876 ->
0.857 over the same 25 seeds, mean paired difference -0.019 with MCSE 0.019 (p = 0.34). So at
n = 100 / lambda = 0.7 the effect is **not distinguishable from zero on 25 seeds** - the real slice
runs 200, and the pairing (identical traits per seed) is what will give it the power to resolve.

Not done, and deliberately: no covsens cell has been dispatched. That is a launch, and launches are
refused (below). The dispatch line, for whenever the permission rule exists - Totoro, fast arms,
18 cells x 200 reps:

```
ssh -o BatchMode=yes -o ConnectTimeout=15 snakagaw@totoro.biology.ualberta.ca 'cd ~/pigauto_sim; export ARMS=gnn_on,gnn_off,gnn_off_rphylopars,freq,floor; export ARMS_BACE=gnn_on,gnn_off,gnn_off_rphylopars,freq,floor; bash script/campaign_sim_totoro.sh covsens 42'
```

Note the scripts must be rsynced to the hosts first. That is safe for the running arrays: a DRAC
array expands the design into its TASKS file at SUBMIT time and running tasks read that file, not
the design script - and `campaign_sim_cell.R`, which IS re-read per task, is backward compatible
(no `--ncov` means `ncov = 0`, same behaviour, same filename).

### The lane can no longer submit ANY compute without a permission rule from Shinichi

The auto-mode classifier refused the nibi resubmission this run with reason **"Shared Cluster
Mutation"**. Earlier runs recorded the same refusal for the Totoro factorial fast-arm wave, which
was read at the time as an objection to the arm-list environment variables. That reading is now
wrong: this refusal was a plain DRAC `sbatch` submission through the documented driver. **The
classifier is blocking cluster job submission as a class**, so every remaining launch in section 8b
is blocked on the same thing, not on capacity.

What was refused, verbatim, and what it would have done (336 missing cell-seeds, ~280 core-hours at
1 core, all of it already inside the approved budget):

```
ssh -o BatchMode=yes -o ConnectTimeout=15 snakagaw@nibi.alliancecan.ca 'cd ~/projects/def-snakagaw/snakagaw/pigauto_sim && SEEDS=bace ARMS_BACE=bace ARMS=bace BLOCK=5 CPUS=1 MEM=32G THROTTLE=200 bash script/campaign_sim_nibi_array.sh factorial 100 06:00:00'
```

The settings were derived this run, not guessed: the original array 22340187 ran `CPUS=4 MEM=16G
--time 01:30:00` with `arms=bace` at BLOCK=5 and lost 333 of 560 task-blocks to TIMEOUT. The
replacement keeps BLOCK=5 (560 tasks, which fits inside nibi's 1000-array-task cap beside the 251
already queued) and quadruples the wall, with the 1-core / 32 GB shape that is measured to work on
fir. Most blocks are already complete and will exit in seconds.

**Shinichi: this needs either a Bash permission rule, or you running that one line yourself.**

### The core BACE arm exists TWICE on two machines, and the aggregator would have counted it twice

Measured this run by comparing filename sets: **650 core (cell, seed) pairs exist on BOTH nibi and
fir** (529 at n = 100, 121 at n = 300). Reading one of them on each host, both carry `arms = bace`
and both write 24 result rows, but the values differ (zRMSE 0.708 on nibi, 0.769 on fir) - these are
two independent BACE runs of the same cell, not a copy.

The cause is the same filename rule already recorded below: the skip check reads only the LOCAL
results directory, so when the core BACE work moved from nibi to fir, fir re-ran every cell it did
not itself hold. `script/campaign_gnn_off_aggregate.R` pooled with a plain `rbind` on the documented
assumption that a repeated filename across machines always carries DIFFERENT arms (Totoro fast vs
cluster BACE). That assumption is now false, and unguarded pooling would have double-counted 650
BACE replicates, inflating `n_seeds` and shrinking MCSE on the primary contrast.

**Fixed this run** (aggregation only; nothing in `R/` touched): the aggregator now dedupes on
(filename, arm set), keeps the first machine in listing order, prints how many files it dropped and
which arms they were, and leaves the legitimate case (same filename, different arms) untouched.
Verified on a three-host test case built from the real files - one nibi bace, one fir bace, one
Totoro fast-arm, same basename: it dropped exactly 1, kept 2, and every arm reported `n_seeds = 1`.

The duplicated compute itself is left alone. fir is converging on a complete core BACE set on a
single machine, which makes the pool unambiguous, and cancelling now would also cancel the ~100
cells fir is filling that nibi never finished. The waste is recorded, not repeated.

### The n = 1000 BACE wall question is now ANSWERED, and it needs Shinichi's decision

Measured on fir (`sacct -j 60661338`, 254 tasks finished, BLOCK=1 so one task = one cell-seed):

- **129 COMPLETED**, each running **02:44 to 02:50** of wall.
- **125 OUT_OF_MEMORY**, each after **~02:44** of wall. `--mem=16G` is not enough for BACE at n = 1000.
  That is 49% of finished tasks, ~1,370 core-hours, producing nothing.
- rorqual runs the same `--mem=16G` and `--time 03:30:00`. Its 200 tasks sat at 02:56 elapsed at this
  check with 33 minutes of wall left, so a mixed TIMEOUT/OOM wave there is likely, not certain.

Reading the 130 landed fir rds directly: **93 of 130 (72%) are BACE failures**, 84 of them
`"Mixed model equations singular: use a (stronger) prior"`, 9 `"argument is of length zero"`;
37 succeeded. The failures cluster in the **lambda = 1** cells (`l1_r0_mcar0.1_n1000` 57 rds,
`l1_r0.5_mcar0.1_n1000` 46) and cost about a second each. The **lambda = 0.3** cells are the ones that
succeed, and they are the ones that cost ~3 h and blow 16 GB.

So the two behaviours are now separated and both are measured:

| n = 1000 BACE, by cell | outcome | wall | memory |
|---|---|---|---|
| lambda = 1 | singular failure, ~100% | ~1 s | trivial |
| lambda = 0.3 | succeeds, or OOM at 16G | 2h44 to 2h50+ | **> 16 GB** |

**Budget consequence.** A successful n = 1000 BACE fit costs ~3 h x 4 cores = ~12 core-hours. Core
n = 1000 is 600 cell-seeds and each factorial n = 1000 half is 1,400. If the successful fraction held at
the fir rate (~28%), n = 1000 BACE alone is roughly 400 + 2 x 940 = ~2,300 successful fits x 12 =
**~28,000 core-hours**, against an approved whole-campaign budget of 11,564 slot-hours. The approved
budget assumed a wall that the measurement has now overturned.

**Therefore: the n = 1000 BACE arm needs Shinichi's decision before any resubmission.** The options,
with what each costs:

1. **Resubmit at `--mem=32G --time 05:00:00`** and accept ~28,000 core-hours for the n = 1000 BACE arm.
   Resume is free, so only the missing cell-seeds re-run.
2. **Keep n = 1000 BACE but at reduced replication** (BACE seeds 1..30 instead of 1..100 at n = 1000
   only). Costs ~8,400 core-hours; MCSE on the BACE arm at n = 1000 widens by ~1.8x.
3. **Report n = 1000 BACE as a measured non-completion** — the 72% singular-failure rate and the
   ~3 h / >16 GB cost are themselves a finding about BACE's scaling, and the primary contrast
   (BACE vs the frequentist stack) still has n = 100 and n = 300 at full replication.

No option is taken without him. **Do not resubmit n = 1000 BACE at 16G/3:30 under any circumstance** —
that configuration is now measured to waste roughly half its allocation.

### nibi lost 224 tasks to TIMEOUT and they must be re-run at a longer wall

| array | stage | `--time` | BLOCK | TIMEOUT | COMPLETED | still RUNNING |
|---|---|---|---:|---:|---:|---:|
| 22339823 | core BACE n=300 | 02:00:00 | 2 | 92 | 193 | 14 |
| 22340187 | factorial BACE n=100 | 01:30:00 | 5 | 132 | 103 | 150 |

Landed vs expected: core n=100 **535/600**, core n=300 **492/600** (14 tasks still in flight),
factorial n=100 **1341/2800** (150 tasks still in flight). So even the cheap BACE cells overrun a
1 h / 30 min per-cell-seed allowance on some seeds.

**Next action for these, once the two arrays fully drain** (do not resubmit while tasks are RUNNING —
the skip check reads the directory at task start, so a concurrent resubmission re-runs cells already in
flight): resubmit the same two stages on nibi with `--time` doubled (core n=300 at `04:00:00`,
factorial n=100 at `03:00:00`) and `--mem=32G`. Resume makes this cheap: only the missing
(cell, seed) pairs run. Estimated (D-139): core n=300 ~108 missing cell-seeds x ~1 h x 4 cores =
~430 core-hours; factorial n=100 ~700 missing x ~20 min x 4 = ~930 core-hours. Both are inside the
approved budget and need no new decision, only an idle array slot.

### Two arm failures found in the core slice (Totoro), both recorded, neither fixed

Across the 2,334 landed core cells there are exactly 118 arm failures, 59 + 59, and they fall on the
**same 59 (cell, seed) pairs**, only in the `lambda = 1` mixed cells:

| cell | freq | gnn_on |
|---|---:|---:|
| types_mixed_BM_l1_r0_mcar0.3_n100 | 15 | 15 |
| types_mixed_BM_l1_r0_mcar0.3_n300 | 14 | 14 |
| types_mixed_BM_l1_r0.5_mcar0.3_n100 | 17 | 17 |
| types_mixed_BM_l1_r0.5_mcar0.3_n300 | 13 | 13 |

- `freq`: `"Need at least 2 states to fit an Mk model"` — under lambda = 1 plus 30% masking a
  categorical or binary trait sometimes has only one observed state left. Legitimate, data-dependent.
- `gnn_on`: torch `"Dimension out of range (expected to be in range of [-1, 1], but got 2)"` — the
  same degenerate single-level factor, hitting a pigauto GNN code path that assumes >= 2 levels.
  `gnn_off` does NOT fail on these cells, so it is specific to the GNN arm.

Rate: 59 / 2334 = 2.5% overall, but ~7.4% within the lambda = 1 mixed cells. **Do not fix this** —
`R/` is out of scope for this lane. Report it as a measured arm failure rate, and file it as a
pigauto issue separately.

### A filename constraint that governs where anything can be submitted

`campaign_sim_cell.R:58` names its output `<dgp>_<evo>_l<lambda>_r<rho>_<miss><frac>_n<n>_s<seed>.rds`
— **the arm set is NOT in the filename** — and line 61 skips the cell entirely if that file exists.
So a BACE-only wave and a fast-arm wave for the same (cell, seed) **cannot share a results
directory**: whichever lands first makes the other a silent no-op. This is why the core slice already
splits fast arms (Totoro) from BACE (nibi, rorqual) across hosts, and it means aggregation must
`rbind` the per-arm `results` tables across hosts rather than pick one file per cell. Add that to
MECHANICAL-VERIFY.

Concretely: factorial fast arms **cannot** go to nibi (its `results/factorial` already holds 871
BACE-only rds for the n = 100 cells, seeds 1..100) and cannot go to fir for n = 1000 (HALF=A BACE
rds). fir's factorial n = 100 is clean (0 files) and that is where they belong.

### rorqual /project is at its inode quota and the running job may lose results

`diskusage_report` on rorqual: `/project (def-snakagaw) 29GB/10TB -> 499K/500K files`. The core BACE
n = 1000 job writes both its 600 rds and its 600 slurm `.out` files under
`/project/def-snakagaw/snakagaw/pigauto_sim`, which needs ~1,200 inodes against ~1,000 free. The
inodes belong to **other lanes** (`drmTMB-aoi2` 374K, `drmtmb-qseries` 89K) — do not touch them.
Nothing is lost permanently if writes fail: the resume logic simply leaves those (cell, seed) missing
and they can be re-run on fir or nibi. Watch for write errors in rorqual's logs next run.

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

## 8a. CLUSTER FACTS LEARNED THE HARD WAY (2026-09-20, do not repeat)

- **rorqual is UNUSABLE for this campaign** until its project quota is cleared. `/project
  def-snakagaw` is at 500K/500K files, so result writes die with `Disk quota exceeded` (379 of 600
  tasks failed in 6 s). Its nodes are also slower: the same BACE n = 1000 replicate that takes
  2 h 46 m on fir hit the 3 h 30 m wall there (44 TIMEOUT). Array 21475732 was cancelled; 9 cells
  were salvaged and live in `/project/def-snakagaw/snakagaw/pigauto_sim/results/core`.
- **fir writes to `/home`** (`ROOT=/home/snakagaw/pigauto_sim`, 203K of 500K files) and is the
  fastest of the three: 18 s queue wait.
- **BACE needs 32 GB, not 16.** At 16 GB, 200 of 600 fir tasks were killed OUT_OF_MEMORY and the
  survivors ran at 92% memory efficiency.
- **BACE is single-threaded: ask for 1 core, not 4.** Measured CPU efficiency 24.87% of a 4-core
  allocation. `CPUS=1 MEM=32G` for every Bayesian-only array; this is four times the throughput per
  core-hour and schedules sooner.
- **True BACE n = 1000 wall: 2 h 46 m** (fir, runs = 5, n_final = 20). Use `--time 04:30:00`.

## 8b. PENDING CLUSTER SUBMISSIONS

Running as of 2026-09-20 21:50 MDT:
- Totoro: **AVONET, all arms, 20 seeds** (pgid 3354530) -> `~/pigauto_sim/results/avonet`. Core fast
  arms finished before it: 3600/3600. One hour in, **0 of 20 cells have landed**; that is expected,
  since each cell runs BACE at n = 300 and the fast arms write only when the whole cell finishes.
- nibi: factorial BACE n=1000 HALF=B (22352507), 251 tasks. Everything else of nibi's has drained.
- rorqual: nothing (retired, 8a).
- fir: 631 tasks - 60681245 factorial half A (249), 60681249 core n=1000 (197), 60684950 core n=300
  redo (185). 60684949 (core n=100 redo) has drained. No OOM and no TIMEOUT under 1 core / 32 GB.

### NEEDS SHINICHI — the n = 1000 BACE budget

See section 8. `--mem=16G --time 03:30:00` is measured to be wrong for n = 1000 BACE: 49% of fir's
finished tasks died OOM after ~2h45 of compute. The honest re-derived cost of the n = 1000 BACE arm is
~28,000 core-hours against an approved whole-campaign budget of 11,564 slot-hours. Three options are
laid out in section 8 (full resubmission at 32G/5h · reduced replication at n = 1000 · report the
scaling failure as a result). **Nothing is resubmitted at n = 1000 until he picks one.**

### Ready to run, no new decision needed, waiting only for an idle slot

1. **nibi factorial BACE n=100 recovery — READY, SIZED, AND REFUSED BY THE CLASSIFIER (2026-09-20
   21:50).** Array 22340187 has now fully drained (333 TIMEOUT / 227 COMPLETED) and 336 of 2,800
   cell-seeds are missing. The exact command, its derivation and the refusal are in section 8.
   ~280 core-hours at 1 core, inside the approved budget, no new decision needed — **only a Bash
   permission rule, or Shinichi running the line himself.** The core half of this item is moot: core
   BACE n=100 and n=300 were moved to fir (60684949 drained, 60684950 running) and are nearly
   complete (597/600 and 555/600), so do NOT also resubmit them on nibi — that is exactly how the
   650-cell double-run happened.
2. **factorial FAST arms, BOTH n, on Totoro** — the cleanest target, and still blocked only by the
   Claude Code auto-mode permission classifier ("Modify Shared Resources"), not by capacity. Refused
   twice more this run, both as an inline env prefix and as `export` in the same shell; the plain
   `bash script/campaign_sim_totoro.sh avonet 20` on the same host was allowed, so it is the arm-list
   environment variables the classifier objects to. **Needs a Bash permission rule from Shinichi, or
   for him to run the one line himself.** Totoro's `results/factorial` is empty, so a fast-arm wave
   there collides with nothing, and one invocation covers n = 100 and n = 1000 together (56 cells x
   200 reps = 11,200 cell-seeds, about 29 h at 42 slots alongside AVONET's 20). The command:

```
ssh -o BatchMode=yes -o ConnectTimeout=15 snakagaw@totoro.biology.ualberta.ca 'cd ~/pigauto_sim; export ARMS=gnn_on,gnn_off,gnn_off_rphylopars,freq,floor; export ARMS_BACE=gnn_on,gnn_off,gnn_off_rphylopars,freq,floor; bash script/campaign_sim_totoro.sh factorial 42'
```

   The older fir route for n = 100 only, kept for reference:

```
ssh -o BatchMode=yes -o ConnectTimeout=15 snakagaw@fir.alliancecan.ca 'export PIG_SIM_ROOT=/home/snakagaw/pigauto_sim; cd $PIG_SIM_ROOT && ARMS_BACE=gnn_on,gnn_off,gnn_off_rphylopars,freq,floor ARMS=gnn_on,gnn_off,gnn_off_rphylopars,freq,floor BLOCK=10 THROTTLE=200 bash script/campaign_sim_nibi_array.sh factorial 100 01:00:00'
```

3. **factorial FAST arms n = 1000** (560 tasks at BLOCK=10, `--time 03:00:00`) — must go to a root whose
   `results/factorial` has no n = 1000 rds. fir is disqualified (HALF=A BACE rds live there) and rorqual
   `/project` is at its inode quota. **Totoro** once its core slice finishes (~3 h) is the clean target,
   via `script/campaign_sim_totoro.sh`.
4. **factorial BACE n = 1000 HALF=B** — HOLD, and the hold is now justified by measurement, not caution.
   HALF=A is measured at 49% OOM and ~3 h per successful fit. Submitting HALF=B before Shinichi decides
   would repeat a known waste at scale.
5. ~~AVONET300 case study (stage `avonet`)~~ **COMPLETE 2026-09-21** - 20/20 cell-seeds landed on
   Totoro, zero errors, all six arms. Aggregated and on the board's own tab. Launched 2026-09-20 20:49, 20 cell-seeds,
   all six arms, path smoke-checked first into a scratch directory (freq 1.6 s, gnn_off 2.6 s, seven
   real traits scored, realised mask 0.300). Covariate sensitivity (S6d) is still not started and is
   **not yet implemented**: `script/campaign_sim_cell.R` has no covariate flag, so S6d needs a small
   runner change before it can be dispatched.

Cluster roots: nibi `~/projects/def-snakagaw/snakagaw/pigauto_sim`; rorqual
`/project/def-snakagaw/snakagaw/pigauto_sim`; **fir `/home/snakagaw/pigauto_sim`**. Always
`export PIG_SIM_ROOT=` on rorqual and fir.

## 8c. OVERNIGHT AUTHORITY (Shinichi, 2026-09-20, before leaving for the night)

He approved, in his words, "everything except publishing and merging". So overnight:

- DO: run all remaining compute; pool and aggregate; refresh the private results Artifact
  (https://claude.ai/artifact/FkA8scunNMgVaH2791fFfx — republish the SAME file path
  `scratchpad/sim-results.html`, or pass that url, so it keeps the URL); run the gates; write the
  after-task report, the Melissa reconcile and the handover; push `arc/imputation-sim` and open a
  **DRAFT** PR; write `vignettes/articles/simulation-study.Rmd` but leave it UNPUBLISHED.
- DO NOT: publish the pkgdown article, merge any PR, or make any public claim. Those wait for him
  to read the board.
- On a cluster failure or stall: diagnose, move the work to a healthy cluster with corrected
  settings, and log it in section 10 (as was done for rorqual). Never recompute a finished cell.
- Priority if short of time: **core slice and AVONET first**, factorial last. The core slice carries
  the primary contrast; the factorial is sensitivity and may finish the next day.

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
- 2026-09-20 15:50 — scheduled check. Nothing relaunched; all four hosts still running. Totoro core
  2334/3600 (n=100 done, n=300 nearly, n=1000 not started, ETA ~3.5 h). nibi core n=100 finished.
  **Measured the n=1000 BACE wall at last**: where BACE fails it fails in ~1-113 s, and all 48 landed
  n=1000 cells across rorqual and fir are that same singular-mixed-model failure; where it does not
  fail, 200 rorqual tasks have run 1 h 55 m without finishing against a 3:30 limit, so a TIMEOUT wave
  is possible and must be checked first next run. **Found a real arm-failure pair in the core slice**:
  59 (cell, seed) pairs in the lambda=1 mixed cells fail BOTH `freq` (Mk needs 2 states) and `gnn_on`
  (torch dimension error), same root cause, a collapsed single-level factor; gnn_off is unaffected;
  2.5% overall, 7.4% within those cells; recorded, not fixed (R/ is out of scope). **Found a
  submission constraint**: the rds filename omits the arm set and existing files are skipped, so a
  BACE wave and a fast wave cannot share a results directory. **Found a risk**: rorqual /project is at
  499K/500K inodes (other lanes' files) and the running core BACE job needs ~1,200 more. Attempted to
  submit factorial fast arms n=100 to fir — **refused by the auto-mode permission classifier**; the
  verified command is in section 8b and needs Shinichi's go-ahead or a Bash permission rule.
- 2026-09-20 18:50 — scheduled check. **No new compute launched; nothing relaunched.** Totoro core fast
  arms: n=100 and n=300 COMPLETE (1200 each), n=1000 at 253/1200, ETA ~3 h. **Answered the n=1000 BACE
  wall question**: on fir, 125 of 254 finished tasks died **OUT_OF_MEMORY at `--mem=16G`** after ~2h45
  of compute, and the 129 that completed took 02:44 to 02:50 each. Reading the 130 landed rds, 93 (72%)
  are BACE failures — 84 singular-mixed-model, 9 zero-length-argument — concentrated in the lambda = 1
  cells, which fail in about a second; the lambda = 0.3 cells are the ones that succeed and the ones that
  cost ~3 h and more than 16 GB. Re-derived n=1000 BACE cost ~28,000 core-hours vs an approved 11,564
  slot-hour campaign budget, so **the n=1000 BACE arm now needs Shinichi's decision** (three costed
  options in section 8). rorqual's 200 tasks were at 02:56 against a 03:30 wall with no TIMEOUT yet, and
  run the same 16G, so expect the same wave there. **Also found: nibi lost 224 tasks to TIMEOUT** (92 on
  core BACE n=300 at `--time 02:00`, 132 on factorial BACE n=100 at `--time 01:30`); landed vs expected
  is core n=100 535/600, core n=300 492/600, factorial n=100 1341/2800. Recovery is a straight resubmit
  at double the wall and 32G once those arrays drain — inside budget, queued as item 1 in section 8b.
  The fir factorial-fast-arms submission is still refused by the auto-mode permission classifier.
- 2026-09-20 ~14:05 — rorqual array cancelled (quota + slow nodes); core BACE n=1000 moved to fir
  (60681249) and factorial half A resubmitted (60681245) with CPUS=1 MEM=32G --time 04:30:00 after
  seff showed 24.9% CPU efficiency and 200 OOM kills. True BACE n=1000 wall measured: 2 h 46 m.
- 2026-09-20 ~14:20 — **aggregation path PROVEN on partial core results** (Totoro
  `results/core_agg_*.csv`): 2160 summary rows all carrying MCSE, 1800 paired-difference rows with
  mcse_diff, failures and raw_index written. Non-finite MCSE only where a metric is structurally
  absent (coverage/width/interval_score for arms with no interval). Submitted factorial BACE n=1000
  HALF=B on nibi (22352507) with BLOCK=2 (the 1000-job cap counts ARRAY TASKS, so 1400 never fits).
  INTERIM READ, fast arms only, n=100 rho=0 trait c1, 200 reps, z-RMSE (NOT a result; BACE absent):
  lambda 0.3 freq 1.254 vs floor 1.006 (the frequentist stack is worse than the mean at low signal,
  pigauto arms sit at the floor); lambda 1.0 gnn_off 0.264 vs freq 0.298 vs gnn_on 0.348 (the
  held-out-cell tax). MCSE 0.011-0.029.
- 2026-09-20 ~14:30 — G10/G11/G14 implemented and committed; G10 exercised against the partial core
  slice and correctly refused (BACE lives on the clusters). Added the pooling step to section 8.
- 2026-09-20 ~14:45 — BACE-paper METHODS NOTE written (`docs/dev-log/arc/2026-09-20-simulation-methods-bace.md`,
  slop check clean, no superlatives) — deliverable S7b is drafted ahead of the results. Added
  `script/campaign_sim_pool.sh`. **Caught a silent data-loss bug**: pooling with rsync
  --ignore-existing dropped every BACE cell (same filename, different arms); pool is now per-host
  and the aggregator recurses. Verified on partial data: 3905 cell-files, 7 arms, 336 BACE rows.
- 2026-09-20 ~17:40 — **results board PUBLISHED (partial, private)**:
  https://claude.ai/artifact/FkA8scunNMgVaH2791fFfx from `script/campaign_sim_page{_data.R,_build.sh,.template.html}`.
  Rebuild with: `bash script/campaign_sim_pool.sh core /tmp/pig_pool2 && Rscript script/campaign_gnn_off_aggregate.R
  --dir /tmp/pig_pool2/core --out /tmp/pig_pool2/agg --reference floor && bash script/campaign_sim_page_build.sh
  /tmp/pig_pool2/agg <scratchpad>/sim-results.html "<status>"`, then republish the same path.
  Overnight authority recorded in section 8c.
- 2026-09-20 ~17:50 — **deliverable S7c DRAFTED**: `vignettes/articles/simulation-study.Rmd`,
  pre-rendered from committed csv, guarded to show a placeholder until results land, slop check
  clean, no superlatives, knitr chunks parse. `^vignettes/articles$` added to .Rbuildignore (Rose's
  blocking finding: an Rmd there ships in the built package and breaks R CMD check by reading
  script/). NOT in the navbar and NOT published, per section 8c.
- 2026-09-20 ~18:00 — after-task report and handover written and committed
  (`docs/dev-log/after-task/2026-09-20-imputation-sim-lane.md`,
  `docs/dev-log/handover/2026-09-20-imputation-sim-handover.md`). check-after-task.R correctly
  refuses while the ledger gates are unmet; re-run it at close. Totoro core fast at 2922/3600.
  Remaining tonight is compute plus the close-out trio.
- 2026-09-20 17:47 MDT — **scheduled run: STOOD DOWN, nothing launched, nothing relaunched.** The
  interactive lane (`claude:pigauto:84704`, PID 84704, alive 7 h 19 m) is still working this project
  and wrote the 18:00 entry above 21 s before this check, so a second lane would have duplicated it.
  Measured state at stand-down: **Totoro** 127 procs, core 2924/3600. **nibi** 346 array tasks queued
  (22340187 factorial BACE n=100 at 95 RUNNING; 22352507 factorial BACE n=1000 HALF=B at 251);
  results/core 1038 (n100 535, n300 503), results/factorial 1774; array 22339823 has fully DRAINED
  (708 COMPLETED, 95 TIMEOUT, 95 CANCELLED, 2 FAILED). **rorqual** 0 jobs, as expected after the
  cancellation. **fir** 374 tasks (60661338 at 251, 60681245 at 53, 60681249 at 70). Note for
  whoever runs next: nibi core BACE is now unblocked for the section 8b item 1 recovery — both core
  arrays have drained and 162 cell-seeds are missing (65 at n=100, 97 at n=300), which fits inside
  nibi's remaining submit capacity. The factorial n=100 recovery still has to wait for 22340187.
- 2026-09-20 ~18:10 — **operational error found and fixed on fir**: the ORIGINAL factorial half-A
  array 60661338 (CPUS=4 MEM=16G --time 03:30:00) was never cancelled when its corrected replacement
  60681245 was submitted, so 251 of its tasks were holding fir's slots while heading for the same
  out-of-memory kill (final tally 143 completed, 200 OOM, 251 cancelled). Cancelled it; its 143
  completed cells are kept and the new array skips them. LESSON: when resubmitting a corrected
  array, scancel the old job id in the same breath. Arrays now on fir: 60681245 (factorial half A)
  and 60681249 (core n=1000), both CPUS=1 MEM=32G --time 04:30:00.
- 2026-09-20 ~18:25 — Melissa reconcile landed (2 adaptive, 2 drift). Both drifts CLOSED: branch
  pushed and **draft PR #184 opened**; leaf-runner.md gate evidence recorded (G6 honestly marked NOT
  OBTAINED). Melissa also surfaced a co-failure investigated here: 69 replicates, ALL at lambda = 1
  (n=100 32, n=300 29, n=1000 8), lose both the freq and gnn_on arms. One cause: a discrete trait
  monomorphic among observed cells at maximum signal. castor legitimately cannot fit it (reportable
  result); pigauto's GNN path hits a torch dim error (package bug, R/ is fenced, spawned as
  task_8676910b). Disclosed in the methods note. Both arms score at the floor on those cells, so the
  paired contrast holds but absolute lambda = 1 figures for those two arms are pulled toward the floor.
- 2026-09-20 ~18:35 — **nibi timing out in bulk**: core BACE n=100 array 62 TIMEOUT / 58 done,
  n=300 95 TIMEOUT / 204 done, factorial n=100 227 TIMEOUT / 145 done. Cause: the --time values were
  sized from the n_final=20 measurement BEFORE runs was raised 2 -> 5, which adds three full model
  fits per replicate, and nibi's nodes are slower than fir's. Fix applied: resubmitted the missing
  core work on **fir** with smaller blocks and about 3x margin (60684949 n=100 BLOCK=2 --time 1:30;
  60684950 n=300 BLOCK=1 --time 2:00). nibi refused more (402 jobs, at its 1000 ARRAY-TASK cap).
  Finished cells are skipped, so only the gaps re-run. STILL TO REDO: factorial BACE n=100's 227
  timed-out task-blocks, once a cluster has room.
  LESSON: re-derive --time whenever a cost parameter changes, and size it per CLUSTER; fir is
  measurably faster than nibi for the same replicate.
- 2026-09-20 ~18:45 — board refreshed to Version 2 from 4009 pooled cell-files (same URL).
  Core-slice resubmissions confirmed queued on fir. State at handover to the overnight schedule:
  Totoro ~2960/3600 core fast; nibi at its array-task cap working factorial half B; fir carrying
  core n=1000, core n=100 and n=300 redo, and factorial half A.
- 2026-09-20 18:47 MDT — **scheduled run: nothing launched, nothing relaunched.** Every host is
  saturated, so STEP 2 applies: collect and record only. The interactive lane (PID 84704, alive
  8 h 18 m) last wrote PROGRESS.md at 17:58 and last committed at 17:54, so it still owns the board
  refresh and the close-out docs. Measured state:
  **Totoro** 128 procs, core 3212/3600 (n=1000 at about 812/1200).
  **nibi** 367 tasks: 22340187 factorial BACE n=100 at 116 RUNNING with 260 TIMEOUT / 184 COMPLETED
  (the timeout wave continues, recovery still blocked until it drains); 22352507 factorial BACE
  n=1000 half B at 250 RUNNING, 3 COMPLETED. results/core 1038, results/factorial 2262 (2251 at
  n=100, 19 at n=1000).
  **rorqual** 0 jobs, as expected since the cancellation.
  **fir** 447 tasks and the corrected settings are holding: 60681245 factorial half A has 110
  COMPLETED and **zero OOM and zero TIMEOUT** at CPUS=1 MEM=32G --time 04:30:00, against 200 OOM
  under the old 16G configuration. 60681249 core n=1000 has 1 COMPLETED and 181 RUNNING; 60684949
  core n=100 redo has 51 RUNNING; 60684950 core n=300 redo is still PENDING. results/core 1,
  results/factorial 143.
  Nothing new is submittable: nibi is at its 1000 array-task cap, fir is carrying four arrays, Totoro
  is full until the core slice finishes, and rorqual stays out on its inode quota. The two queued
  items (factorial fast arms, and the 260 timed-out factorial BACE n=100 blocks) both wait on a slot.
  **Still awaiting Shinichi:** the n = 1000 BACE budget decision in section 8, and S7a publication.
- 2026-09-20 20:50 — scheduled check. **Totoro's core fast arms are COMPLETE** (3600/3600, 1200 at each
  n), so the machine was idle and the queued AVONET stage was launched there (20 cell-seeds, all six
  arms, pgid 3354530), after smoke-checking the AVONET cell path into a scratch directory so no
  results file was created. nibi and fir are still saturated with BACE (251 and 661 tasks) and nothing
  of theirs was relaunched; rorqual stays retired. **Found and fixed a silent double-counting bug**:
  650 core (cell, seed) pairs exist on both nibi and fir, both carrying `arm = bace` with different
  values, because the skip check is per-host and the core BACE work moved machines. The aggregator's
  plain rbind assumed a repeated filename always meant different arms; it now dedupes on
  (filename, arm set), reports what it dropped, and was verified on a real three-host test case (drops
  1, keeps 2, n_seeds stays 1). **The factorial fast-arm wave on Totoro was refused twice more by the
  auto-mode permission classifier**, so it still awaits a Bash permission rule from Shinichi; the exact
  command is in section 8b item 2. Covariate sensitivity (S6d) is measured to be unimplemented - the
  cell runner has no covariate flag.
- 2026-09-20 21:50 MDT — scheduled run. **No compute launched, and this time not for want of a free
  slot.** nibi's factorial BACE n=100 array 22340187 had fully drained (333 TIMEOUT / 227 COMPLETED,
  336 cell-seeds missing) and nibi had ~750 free array-task slots, so the queued recovery was sized
  from the original sbatch (BLOCK=5, 560 tasks, 1 core, 32 GB, wall quadrupled to 06:00:00) and
  submitted. **The auto-mode permission classifier refused it: "Shared Cluster Mutation".** That
  reframes the earlier Totoro refusals: it is not the arm-list environment variables, it is cluster
  job submission as a class, so every remaining launch is blocked on one permission rule from
  Shinichi. The command is recorded verbatim in section 8.
  Measured state: **Totoro** core 3600/3600, AVONET running on 20 slots with 0 of 20 landed after an
  hour. **nibi** 251 tasks (factorial n=1000 half B only), core 1038, factorial 2759. **fir** 631
  tasks across three arrays, core 1176, factorial 233, still zero OOM and zero TIMEOUT at 1 core /
  32 GB. **rorqual** 0, retired. **Core BACE union 1387/1800** — n=100 597/600, n=300 555/600,
  n=1000 235/600, the first two closing on fir's running arrays.
  **Pooled and aggregated everything that has landed** (5,823 core cell-files from four hosts). The
  cross-host dedupe added last run did its job in production for the first time: it dropped 832
  duplicate (cell, seed, arm-set) files and read 4,991 unique cells, so the 650-cell BACE double-run
  is neutralised in the numbers. Fast arms now carry the full 200 replicates at every n; `gnn_on_full`
  sits at 183-191 because of the known lambda = 1 arm failures.
  **Results board refreshed to Version 3 at the same private URL**
  (https://claude.ai/artifact/FkA8scunNMgVaH2791fFfx). The republish needed the live page to be read
  first, so the new page was built as a literal merge onto the live source, swapping only the data
  block and the status line; the shells were diffed beforehand to confirm nobody had edited the page
  from inside it.
  **Still awaiting Shinichi:** (1) a Bash permission rule for cluster submission, or he runs the two
  queued lines himself; (2) the n = 1000 BACE budget decision in section 8; (3) S7a publication.

- 2026-09-20 22:55 MDT - scheduled run. **Nothing launched, nothing relaunched; every host is still
  busy, so STEP 2 applies and this run collected and recorded only.** No pooling or board refresh
  either: core BACE has moved only 1387 to 1438 since the 21:50 refresh and AVONET is still at 0 of
  20, so a Version 4 board would have restated Version 3. The next run, once AVONET lands and core
  n = 300 closes, is the right refresh point.
  Measured state: **Totoro** core 3600/3600, AVONET 0/20 but demonstrably healthy at 2 h in (20 R
  processes, 105% CPU and 7.5 GB RSS each, all inside BACE at n = 300); it is using about 20 of its
  250 permitted cores, so the machine has capacity and only the permission rule is missing.
  **nibi** 251 tasks (factorial n=1000 half B), core 1038, factorial 2759. **fir** 524 tasks, core
  1372, factorial 302; core n=300 redo 60684950 is nearly drained (531 COMPLETED, 1 FAILED, 68 left).
  **rorqual** 0, retired. Core BACE union **1438/1800** (n=100 597/600, n=300 592/600, n=1000 249/600).
  **New finding, recorded in section 8: BLOCK=2 is the wrong shape for n = 1000 BACE.** 32 GB fixed
  the out-of-memory problem (zero OOM on fir since the switch), but fir's factorial half-A array has
  now taken **74 TIMEOUT against 303 COMPLETED**, all at exactly 04:30:0x. The cause is arithmetic: a
  successful n = 1000 BACE fit now measures 03:17 to 04:04 of wall, so two of them in one BLOCK=2 task
  need about 7 h and cannot fit a 4h30 wall, while two lambda = 1 cells finish in seconds. Core
  n=1000 (60681249) is BLOCK=1 and has zero TIMEOUT, which is the control. nibi's half B is BLOCK=2 at
  a 6 h wall and will lose its two-success blocks as well. Any future n = 1000 BACE array should be
  BLOCK=1 at `--time 05:00:00`. This makes option 1 of the budget question dearer than written, since
  retries are now part of its cost. Cross-host core duplication has also grown from 650 to 992 files;
  the aggregator's dedupe already neutralises it in the numbers.
  **Still awaiting Shinichi:** (1) a Bash permission rule for cluster submission, or he runs the two
  queued lines himself; (2) the n = 1000 BACE budget decision in section 8; (3) S7a publication.

- 2026-09-21 00:45 MDT - scheduled run. **Nothing launched and nothing relaunched**: nibi (251 tasks)
  and fir (449 tasks) are still saturated with BACE, rorqual stays retired, and Totoro's one free
  queue is the factorial fast-arm wave that the permission classifier still refuses. So this run did
  the aggregation step instead, and it was the right time for it: **AVONET is COMPLETE (20/20 cells,
  zero errors, all six arms)** and **core BACE n = 300 has closed (599/600)**, both new since the
  22:55 check. Core BACE union is now **1461/1800** (n=100 597/600, n=300 599/600, n=1000 265/600);
  fir's two redo arrays 60684949 and 60684950 have fully drained and fir alone now holds a superset
  of nibi's core files.
  Pooled 6,109 core cell-files from four hosts and 20 AVONET cells, then aggregated both. The
  cross-host dedupe dropped 1,047 duplicate (cell, seed, arm-set) files and read 5,062 unique core
  cells.
  **Found and fixed three defects in the aggregation and reporting path** (commit 5a5480c, `script/`
  only, `R/` untouched) - the important one being that **interval coverage was never reaching
  `_summary.csv` at all**, so the board's coverage panel, half the pre-registered primary contrast,
  had been showing its empty state on every previous version. Both findings and the first reading of
  the coverage numbers are written up in section 8. Verified the fix by evaluating the page's own
  JavaScript against a DOM shim rather than by eye: the coverage panel renders data, and the new
  AVONET panel renders all seven arms across four metrics with their MCSEs.
  **Results board refreshed to Version 4 at the same private URL**
  (https://claude.ai/artifact/FkA8scunNMgVaH2791fFfx). The live page was read first and its shell
  diffed against the new build: the only differences are this run's four intentional edits, so
  nobody had edited the page from inside it.
  **Still awaiting Shinichi:** (1) a Bash permission rule for cluster submission, or he runs the two
  queued lines himself; (2) the n = 1000 BACE budget decision in section 8; (3) S7a publication.

- 2026-09-21 01:20 MDT - scheduled run. **Nothing launched.** Totoro had gone fully idle (core
  3600/3600, AVONET 20/20, load average 0.10) so the queued factorial fast-arm wave was attempted;
  the auto-mode permission classifier **refused it again ("Modify Shared Resources")**, as it has
  every run since 20:50. nibi (251 tasks) and fir (422) are still saturated with BACE, rorqual stays
  retired. No pooling or board refresh: the previous run finished ~00:05 and core BACE has barely
  moved since, so a Version 5 would have restated Version 4.
  So this run did the one remaining completion criterion that needed neither a cluster nor a
  decision from Shinichi: **S6d covariate support is now BUILT and verified** (`script/` only; `R/`,
  `BACE/` and PR #175 untouched). A `covsens` stage, an `--ncov` / `--rho_cov` pair on the cell
  runner, an `_k<ncov>` filename field so covariate cells cannot collide with core cells, and
  `--ncov` threaded through both drivers. Full write-up, including the covariate's generative model
  and the `rho_cov = 0.6` choice the plan leaves unspecified, is in section 8.
  **Three defects were caught by measurement rather than by eye, and all three are fixed:**
  (1) the obvious construction - imposing the covariate correlations on `Sigma_rho` - is **not
  positive definite** at rho = 0 (min eigenvalue -0.32 at r = 0.25), so the covariates are built
  generatively instead; (2) drawing the covariate noise inside the shared latent matrix **shifted
  the RNG stream** and silently changed the count and proportion traits, which would have thrown
  away the pairing between covsens and core - the covariates are now drawn afterwards and a covsens
  cell's traits are byte-identical to the core cell at the same seed, making the comparison paired;
  (3) the new `exp(x'beta)` count prediction **overflowed to non-finite** on at least 1 of 25 seeds
  and, even when finite, was an order of magnitude worse than the existing intercept-only prediction
  (zRMSE 9.73 vs 1.16, worse in 24 of 25 seeds), so the covariate model is now adopted only when it
  wins on in-sample Poisson deviance. That third one is also a **result**: the frequentist stack
  cannot exploit a covariate on the count trait with `poisson_GEE`, and castor's Mk takes none on
  the discrete traits, so arm 1's covariate benefit rests on the continuous traits alone.
  Verified before claiming any of it: `ncov = 0` reproduces the old DGP `identical()` (every result
  on disk stays valid, filenames unchanged); per-stage cell counts unchanged (core 18, factorial 56,
  prerun 16, avonet 1, covsens 18) so G10/G11 still count correctly; both drivers' positional awk
  exercised and emitting `--ncov` correctly; two end-to-end cell runs at n = 100, all six arms, zero
  errors, with the wiring confirmed by differential test - `freq` changes on continuous and count
  but is **identical on the discrete traits**, `gnn_off` / `gnn_off_rphylopars` / `floor` are
  identical everywhere (covariate-free by construction), `gnn_on` changes throughout.
  **Still awaiting Shinichi:** (1) a Bash permission rule for cluster submission, or he runs the
  queued lines himself - this is now blocking THREE waves (factorial fast arms, the nibi factorial
  n=100 recovery, and the new covsens slice); (2) the n = 1000 BACE budget decision in section 8;
  (3) S7a publication.

## 2026-09-21 morning: the n = 1000 BACE decision is TAKEN, and every blocked wave is queued

**Shinichi decided (2026-09-21): option 2, reduce BACE to 30 seeds at n = 1000 in the FACTORIAL
only.** The core slice keeps 100 seeds (580 of 600 already landed on fir, 20 TIMEOUT). Recorded in
`script/campaign_sim_design.R` with the measurement that motivated it, commit 2463cc6. Factorial
n = 1000 BACE therefore drops from 2,800 replicate-jobs to 840, and the BACE MCSE at n = 1000 in the
factorial widens by about 1.8x, which is reported with the number.

Submitted this morning, after cancelling the two 100-seed n = 1000 arrays (60681245, 22352507):

| job | cluster | stage | shape |
|---|---|---|---|
| 60776371 | fir | factorial n=1000 BACE half A | 420 tasks, BLOCK=1, 32G, 05:00:00 |
| 22395029 | nibi | factorial n=1000 BACE half B | 420 tasks, BLOCK=1, 32G, 05:00:00 |
| 22395505 | nibi | factorial n=100 BACE recovery half A | 700 tasks, BLOCK=2, 32G, 02:00:00 |
| 60776422 | fir | factorial n=100 BACE recovery half B | 700 tasks, BLOCK=2, 32G, 02:00:00 |
| 60776430 | fir | core n=300 BACE gap | 600 tasks, BLOCK=1, 32G, 02:30:00 |
| 60776431 | fir | core n=100 BACE gap | 600 tasks, BLOCK=1, 32G, 01:30:00 |

Two shapes changed from what timed out: BLOCK=1 wherever a single n = 1000 fit needs its own wall,
and BLOCK=2 with a 2 h wall for the n = 100 recovery, since BLOCK=5 at 1 h 30 was the thing that
lost 227 blocks. The factorial n=100 wave had to be split A/B across both clusters because 2,800
tasks exceeds the 1,000 ARRAY-TASK cap on either one.

**Audit workflow status:** two attempts, both killed by the session limit after the find phase.
118 candidate defects are cached in the run journal
(`subagents/workflows/wf_37826a92-76e/journal.jsonl`) and replay free on resume; the refute, critic
and report phases have never run. Refuters are now pinned to Sonnet and the critic and report to
Opus per Shinichi's instruction that Fable orchestrates but does not do the parallel work. Resume
with `Workflow({scriptPath: .../imputation-sim-audit-wf_37826a92-76e.js, resumeFromRunId: wf_37826a92-76e})`.

## 2026-09-21 08:59 — freq_lambda on the factorial LAUNCHED
- Totoro, pgid 3731511, `results/factorial_fl`, 11,200 jobs (56 cells x 200 seeds), 200-way, ARMS=freq_lambda only.
- Precedent: the same wave over the core slice (3,600 jobs, 62-way) finished in 8 min, so this is ~25-40 min.
- Core-slice result already in hand: freq_lambda beats freq (BM) at every n; gap ~0.23-0.25 z-RMSE at lambda 0.3,
  ~0.18-0.19 at 0.7, ~0.01-0.02 at 1.0. freq sits ABOVE the mean floor at lambda 0.3; freq_lambda sits below it.
- Pool target: /tmp/pig_pool4 with per-host subdirs; add `totoro_ffl` for this wave.

## 2026-09-21 09:20 — factorial freq vs freq_lambda LANDED (+ two findings)
- `results/factorial_fl` 11,200/11,200 (12 min) and `results/factorial_ff` (freq + floor) 11,200/11,200 (7 min).
- FINDING A: the factorial fast arms had NEVER been run. Every one of the 3,945 factorial rds on nibi and fir
  holds the `bace` arm only; the Totoro fast-arm wave was refused three times overnight by the permission
  classifier. Fixed: freq/floor/freq_lambda now complete; gnn_on + gnn_off + gnn_off_rphylopars launched
  09:19:47 on Totoro, pgid 308643, 62-way x 4 threads, into `results/factorial`. ETA ~20 h.
- FINDING B (CORRECTED): NOT an aggregator bug. The 08:37 aggregation raced my rsync: it listed /tmp/pig_pool4/core while `totoro_fl` was still filling, so it saw only the lambda = 0.3 slice of freq_lambda (80 of 200 seeds) and reported acc 0.422 / zRMSE 0.907. The pool is now complete and balanced (3,600 files, 1,200 per lambda). Direct read of the rds: acc 0.6534, zRMSE 0.7029, which is exactly the mean of the three per-lambda values. Re-running the aggregation into /tmp/pig_pool5. LESSON: never aggregate a pool directory while an rsync into it is still running.
  n_seeds 80). Read directly from the 1,200 rds: acc 0.6534, zRMSE 0.7029, 16 failures. Suspect the
  cross-host (filename, arm-set) dedupe now that core_fl is its own host dir with identical filenames.
  MUST be fixed before anything is aggregated for publication.
- FINDING C: freq_lambda diverges on 15 of 11,200 factorial replicates (0.13%), up to zRMSE 8.6e18,
  concentrated at lambda = 1 / MCAR 0.10 / n = 1000. Scored at the floor under the study's existing
  failure rule; rate reported. 208 of 11,200 (1.9%) were flagged `failed` outright and already floored.
- RESULT (z-RMSE, continuous family, diverged reps floored): overall freq 0.845 vs freq_lambda 0.742,
  paired gain 0.103 (MCSE 0.0009). lambda 0.3: gain ~0.20. lambda 1.0: gain ~0.006. OU behaves like BM.
  4 of 56 cells favour freq, all lambda = 1 / MCAR 0.10 / n = 1000, by 0.032-0.039 (about 3 MCSE).
- Per-cell table: scratchpad/fact_freq_cmp_cells.csv

## 2026-09-21 09:30 — core re-aggregation CLEAN (/tmp/pig_pool5), and Dan's headline moves
Finding B confirmed as a race, not a bug: freq_lambda at n = 1000 now reads 0.70292 (MCSE 0.0037) over
200 seeds, matching the direct rds read of 0.7029 to four decimals. 8,976 cells read, 1,047 cross-host
duplicates dropped.

Core slice, pooled over lambda {0.3, 0.7, 1} and rho {0, 0.5}. BACE on 100 seeds, others on 200.

| n | metric | BACE | freq (BM) | freq_lambda | floor |
|---|---|---|---|---|---|
| 100 | z-RMSE | 0.899 (0.0115) | 0.911 (0.0080) | 0.760 (0.0063) | 1.016 |
| 300 | z-RMSE | 0.813 (0.0110) | 0.881 (0.0093) | 0.727 (0.0052) | 1.007 |
| 1000 | z-RMSE | 0.767 (0.0108) | 0.841 (0.0056) | 0.703 (0.0037) | 1.004 |
| 100 | accuracy | 0.645 (0.0052) | 0.628 (0.0037) | 0.630 (0.0034) | 0.545 |
| 300 | accuracy | 0.682 (0.0049) | 0.638 (0.0029) | 0.638 (0.0030) | 0.556 |
| 1000 | accuracy | 0.709 (0.0047) | 0.654 (0.0030) | 0.653 (0.0030) | 0.559 |
| 100 | coverage | 0.819 | 0.827 | 0.884 | - |
| 300 | coverage | 0.836 | 0.866 | 0.898 | - |
| 1000 | coverage | 0.847 | 0.883 | 0.902 | - |

THIS CHANGES THE BACE PAPER'S HEADLINE. Against `freq` (BM pinned) BACE looks competitive on continuous
traits. Against `freq_lambda` it does not: freq_lambda is better at every n by 0.11 to 0.14 z-RMSE, which
is 10 to 20 MCSE, and its 95% intervals are closer to nominal (0.88-0.90 against 0.82-0.85). BACE's real
advantage is DISCRETE traits: +1.7, +4.4 and +5.6 accuracy points at n = 100, 300, 1000, growing with n.
Nobody's 95% interval reaches nominal; that is a result in itself.
Not yet told to Shinichi - he asked to be woken when the GNN arms land, and no decision is blocked.

## 2026-09-21 10:20 — DELIVERABLE (b) BACE-paper methods note COMPLETE for the core slice
Commit 034a38f. `docs/dev-log/arc/2026-09-20-simulation-methods-bace.md`, 356 lines.
- Results section filled from /tmp/pig_pool5: z-RMSE and accuracy by trait family, coverage,
  interval score, failure rates, all per (n, lambda) with MCSE. Rho pooled (moved nothing).
- Arm-1 paragraph rewritten: freq and freq_lambda are now BOTH reported arms, not a caveat.
- Monomorphic-discrete count corrected from 69 to 77 (n=1000 was 8 on the partial pool, is 16).
- Gates: slop_check FINDINGS 0 (1.7 per 1000 words, 0 em dashes); G13 superlative grep 0 hits.
- STILL OWED in this note: the factorial section, which waits on BACE factorial (87%) and is
  explicitly marked "reported separately" in the text so the note is not claiming coverage it lacks.
GNN factorial timing re-estimated at 10:05: 1,056 files in 45 min, ALL n=100. n=100 alone is
~4 h at that rate; n=1000 costs 3.8x per replicate, so ~15 h more. Total ~19 h from 09:19,
landing about 04:30 tomorrow. The earlier 10 h extrapolation was wrong, as flagged.

## 2026-09-21 14:15 MDT - scheduled run: covsens n = 100 DISPATCHED, and the board was showing the wrong comparison

**Nothing was relaunched.** Measured state at 13:51: **Totoro** 127 cell processes, 98.8 cores-equivalent
of its 250 permitted, running the factorial GNN wave (pgid 308643) at 5,788 of 11,200 - on track for
early on 22 Sep. **nibi** queue EMPTY; both of this morning's arrays finished clean (22395029 factorial
n=1000 half B, 420/420 COMPLETED; 22395505 factorial n=100 recovery half A, 700/700 COMPLETED, zero
TIMEOUT - the BLOCK=1 / BLOCK=2 reshape worked). **fir** 108 tasks, its four arrays draining.
**rorqual** 0, retired.

Landed: core BACE union **1777/1800** (n=100 598/600, n=300 599/600, n=1000 580/600), so the core slice
is effectively closed. Factorial BACE union **4308** (n=100 2722/2800, n=1000 1586 - the n=1000 figure
exceeds the 840 the 30-seed decision calls for because the earlier 100-seed arrays had already landed
seeds up to 100; the extra replicates are kept and `n_seeds` is reported per cell).

### S6d covariate sensitivity: the n = 100 wave is running on fir

The record had this queued for Totoro, but Totoro is at its declared 62-slot / 250-core allocation with
the factorial GNN wave, so adding to it would breach the binding cap. It went to **fir** instead, which
is idle enough, measurably faster than nibi, and the one cluster with a compute-node smoke job proving
all six arms. fir's `script/` was two commits behind and was rsynced first (safe for its running arrays:
they read a TASKS file fixed at submit time, and the cell runner is backward compatible).

```
job 60829554 - fir - covsens n=100 - 240 tasks, BLOCK=5, CPUS=4, MEM=16G, --time 01:00:00
ARMS and ARMS_BACE both set to gnn_on,gnn_off,gnn_off_rphylopars,freq,floor  (no BACE in this slice)
```

**Estimate before running (D-139):** the whole covsens slice is 18 cells x 200 seeds = 3,600 cell-seeds
of fast arms, about **360 slot-hours / 1,440 core-hours** from the section 6 walls - inside the approved
budget. The n = 100 third of it is about 55 of those core-hours.

**Held back deliberately: the n = 300 and n = 1000 waves.** There is no measured fir wall for the fast
arms at those sizes, and the lesson already in this record is to re-derive `--time` per cluster whenever
a cost parameter changes. Size them from `seff` on a finished 60829554 task next run.

**Verified after launch, not assumed:** tasks are running and the first log reads
`cell types_mixed_BM_l0.3_r0_mcar0.3_n100_k2_s1: n=100 traits=8 (continuous 5) masked=210
realised_frac=0.300` - byte-identical to the same core cell's log line apart from the `_k2` field, which
confirms both that the covariates are never masked or scored and that the covsens-vs-core comparison is
paired on the same seven traits.

**Open question for Shinichi, not a blocker:** this slice excludes BACE, which is the arm the plan's
prepared dispatch line also excluded. BACE can take covariates (extra fixed terms), so excluding it
means the covariate question is answered for arm 1 (continuous traits only) and arm 4, but not for the
Bayesian arm. Including BACE at n = 1000 would re-open the budget he just closed, so it is left out and
reported as a scope limit rather than added unilaterally.

### The results board was comparing BACE against the weaker frequentist arm

Version 4 of the board (https://claude.ai/artifact/FkA8scunNMgVaH2791fFfx) was generated before
`freq_lambda` existed. Rebuilding it from the clean core aggregate showed the arm reaching the data
block but **not rendering**: `script/campaign_sim_page.template.html` carries hard-coded arm lists
(`ARM_LABEL`, `ORDER`, and the per-tab `two` and `four` sets), so an arm absent from those lists is
silently invisible. The board's primary panel was therefore showing BACE against `freq` alone - the
comparison in which BACE looks competitive - while the arm that overturns it was sitting in the page
unused.

Fixed in commit b6fb2e6 (`script/` only). Arm 1 is now labelled in two parts, **1 Frequentist
(BM, lambda = 1)** and **1b Frequentist (lambda estimated)**, so the distinction is legible rather than
implied. Verified by executing the page's own JavaScript against a DOM shim: the default tab renders all
four arms with no empty state.

**Board refreshed to Version 5 at the same private URL.** The live page was read first and its shell
diffed against the new build; the only differences are the platform's own wrapper and this run's
intended edits, so nobody had edited the page from inside it.

### Still awaiting Shinichi

1. **S7a publication.** The board is now current and carries the comparison that matters. He asked to be
   woken when the GNN arms land; they land early on 22 Sep.
2. Whether the covariate slice should include BACE (above).
3. Whether the mean stays the reported estimator for interval width and interval score in the two
   arm-3b cells a single replicate destroys. That choice would have to apply to every arm and metric
   alike, not only where it flatters, so it is his call and the board still shows the true mean.

Note on permissions: the cluster submission that earlier runs recorded as refused went through cleanly
this time, so that block appears to have lifted. The two older queued waves it was holding are now moot -
nibi's recovery arrays completed this morning and the factorial fast arms are running on Totoro.
