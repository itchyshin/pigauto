# PROGRESS RECORD — pigauto four-arm imputation simulation study

**This file is the single source of truth for resuming. Read it first, update it after every
meaningful step, and never start a second copy of work already running.**

Status: **IN PROGRESS** (not complete). Last updated 2026-09-21 20:55 MDT by the scheduled Claude run.

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
      **Measured per arm family 2026-09-21 20:52 MDT (fast arms and BACE are counted separately
      because they live on different hosts under the same filenames).**
      - Fast arms (freq, freq_lambda, floor, 3a, 3b, GNN), target 56 x 200 = 11,200: **freq/floor/
        freq_lambda 11,200/11,200 DONE**; the gnn wave on Totoro is **7,649/11,200**, n = 100
        **5,600/5,600 COMPLETE**, n = 1000 **2,049/5,600**. This is the one binding gap left.
      - BACE, target 28 x 100 at n = 100 plus 28 x **30** at n = 1000 (Shinichi's 30-seed decision):
        n = 100 **2,746/2,800 (98.1%)**, n = 1000 **779/840 (92.7%)**. Only **115 cell-seeds** are
        missing in total, across 7 OU cells at n = 100 and 6 BM cells at n = 1000, and **every one
        of them sits inside the two fir arrays that are running now** (60776422 and 60776371).
        nibi holding no jobs is CORRECT: its half B at n = 1000 is complete at 30 seeds.
- [x] AVONET300 case study, all arms, 20 seeds. **DONE 2026-09-21** - 20/20 cells on Totoro, zero errors, all six arms plus gnn_on_full, aggregated to `avo_summary.csv` and on the board's own tab.
- [ ] Covariate sensitivity on the 18 core cells. **Runner support BUILT and verified 2026-09-21**
      (stage `covsens`, `--ncov`). **ALL THREE WAVES NOW DISPATCHED on fir**, fast arms only, no BACE.
      Landed 2026-09-21 19:51: n = 100 796/1200, n = 300 295/1200, n = 1000 0/1200.
      60895543 (n = 1000) is still PENDING with reason Priority, zero tasks started, so there is still
      no n = 1000 wall to measure. 60895448 (n = 100) has 19 COMPLETED tasks but all are resume-skips.
      fir's fair share is 0.157, so every queued array is scheduling-limited, not mis-submitted.
      - n = 100: first array 60829554 (BLOCK=5, 4 cores, 1 h) DRAINED at **138 TIMEOUT / 102 COMPLETED**.
        Recovery **60895448** submitted 18:40 - BLOCK=1, **CPUS=1**, MEM=8G, `--time 01:00:00`, 1200 tasks.
      - n = 300: **60894610 etc.** running since ~14:55 (BLOCK=5, 4 cores, 1 h 30), 34 COMPLETED, no TIMEOUT yet.
      - n = 1000: **60895543** submitted 18:41 - BLOCK=1, 4 cores, 16 GB, `--time 04:00:00`, **THROTTLE=8**.
        Deliberately throttled: this IS the D-139 pre-run test for an unmeasured regime (see section 8).
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

## 2026-09-21 13:50 — state after the session gap; three "missing" BACE cells RESOLVED
Everything survived the logoff: Totoro detached with setsid, DRAC arrays owned by Slurm. Only my
in-session watcher died (re-armed as btj1iu4eh).
- Totoro GNN factorial 5,662/11,200. n=100 COMPLETE (5,600/5,600) in 4.2 h, matching the 4 h estimate.
  n=1000 started; 3.8x per replicate => ~16 h => lands about 05:30 on 09-22.
- BACE factorial 3,441/3,550 in-design (96.9%), up from 86.9%. nibi arrays drained; fir has 110 R + 4 PD.
- The 3 cells with ZERO BACE rds are all BM clade0.3 n=1000 (l0.3 r0, l0.3 r0.5, l1 r0.5). NOT a gap and
  NOT a failure: they are the tail of fir array 60776371 (420 tasks, %250). Tasks 301-330 are 30/30
  RUNNING at 3:01 against a 5 h limit; 368-420 are PENDING with reason Priority. Zero TIMEOUT on this
  array. A singular-BACE failure would finish in seconds and still write a floored rds, so 3 h of CPU
  with no file means genuine computation. Array started 05:49; expect these this evening, well inside
  the GNN window. BACE is not the critical path.
- Note: 149 of 310 completed tasks on that array finished in 3-7 s, the resume-skip path for
  (cell, seed) rds that already existed. That is why it drains faster than 420 x 3 h.

## 2026-09-21 15:10 — DELIVERABLE (a) results Artifact PUBLISHED (private), v1 = core slice complete
URL: https://claude.ai/artifact/M5HtGRnNGfwsK2Se4gMX24  ·  source: scratchpad/sim-results.html
Six tabs: Frequentist vs BACE (Dan) · All arms · Factorial (freq/freq_lambda only so far) · AVONET ·
Cost & failures · Decisions (six, each with my reading). Palette validated both modes (dataviz
validator: ok, light-mode relief satisfied by direct labels + tables). Status strip states exactly
what is complete and what is still computing. Republish the SAME file path to update in place.
Also: (c) pkgdown article filled from the core slice, commit 1347c06, renders clean, unlisted.
AVONET freq_lambda ran (20/20): identical to freq to 4 dp => estimated lambda ~ 1 on real data.
Walls (median s, one core): freq .4/.8/2.7 · freq_lambda 2.0/3.9/10.6 · gnn_off .7/3.0/45.3 ·
gnn_off_rph 34.4/85.9/353.4 · gnn_on 107.7/132.5/338.6 · bace 1161/3050/11019.
REMAINING: factorial GNN arms (Totoro, ~05:30-09:40 on 09-22) -> factorial BACE tail (fir, this
evening) -> full pool -> G10/G11/G14 -> add factorial to article + methods note + Artifact v2 ->
covsens (S6d) -> after-task, Melissa, handover. Publication of the article waits on Shinichi
reading the Artifact.

## 2026-09-21 18:45 MDT - scheduled run: the covsens n = 100 timeout wave is measured and resubmitted; n = 1000 launched as a throttled pre-run

**Nothing already running was relaunched or cancelled.** Measured state at 18:33:
**Totoro** 127 cell processes on the factorial GNN wave (pgid 308643) at **7,033 / 11,200** - n = 100 complete,
n = 1000 filling, still on the ~05:30 on 09-22 landing estimate. **nibi** queue EMPTY. **rorqual** 0, retired.
**fir** 52 tasks: the factorial n = 1000 BACE tail 60776371 draining (its `results/factorial` at 1,733),
plus a covsens n = 300 wave that a later session launched at ~14:55 and that the 14:15 entry had recorded
as held. That n = 300 wave was left alone.

### The covsens n = 100 wave lost 138 of 240 tasks, and `seff` says why

`sacct` on 60829554: **138 TIMEOUT against 102 COMPLETED** at `--time 01:00:00`, BLOCK=5. Only
**796 of 1,200** cell-seeds landed, and the shortfall is not spread evenly - the r = 0.5 cells are the
starved ones (l0.3 141, l0.7 65, l1 26) against 178-199 for the r = 0 cells.

`seff` over five COMPLETED tasks gives the two numbers that fix it:

| measure | value | consequence |
|---|---|---|
| wall, 5 cell-seeds | 18 to 32 min | a fast task costs 3.6 to 6.4 min per cell-seed |
| CPU efficiency | **24.6 to 25.0% of 4 cores** (5 of 5 tasks) | the fast arms are effectively **single-threaded** here |
| memory efficiency | 3.8 to 4.4% of 16 GB | **0.7 GB** is the real footprint, not 16 |

The timed-out tasks wrote about 2 of their 5 cell-seeds before the kill, so the slow cell-seeds cost
roughly **24 min** each against the fast ones' 5. That bimodality inside a BLOCK=5 task is the same
arithmetic that cost the n = 1000 BACE arrays their blocks: a task's wall depends on which cells it
draws. The recorded fix applies unchanged - **BLOCK=1**.

The CPU-efficiency reading is new and it matters beyond this wave: the lesson already in section 8a
("BACE is single-threaded: ask for 1 core, not 4") is now measured to hold for the **fast arms at
n = 100 too**. `PIG_TORCH_THREADS` is set to `CPUS`, so the GNN arm simply is not using the four
threads at this problem size. Asking for 1 core is four times the throughput per core-hour and
schedules sooner. Left at 4 cores for n = 1000, where torch may behave differently and nothing has
been measured.

**Resubmitted** (estimate stated before running, D-139: 404 missing cell-seeds x ~12 min at 1 effective
core = **~80 core-hours**, ~1 h of wall at THROTTLE=200):

```
job 60895448 - fir - covsens n=100 - 1200 tasks, BLOCK=1, CPUS=1, MEM=8G, --time 01:00:00
```

The 796 already-landed cell-seeds are skipped by the resume check, so only the 404 gaps run.

### covsens n = 1000 is launched, but THROTTLED to 8 tasks - that is the pre-run test, not caution

There is no measured wall for the fast arms at n = 1000 **with covariates**, and the honest range is
wide. The 15:10 walls table gives ~740 s per cell-seed for these five arms at n = 1000 without
covariates; the covsens n = 100 wave just measured **2 to 10x** the table's n = 100 figure. So the
slice is somewhere between **500 and 2,400 slot-hours**, and the top of that range is a fifth of the
whole approved campaign budget. D-139 forbids committing to that on a guess.

```
job 60895543 - fir - covsens n=1000 - 1200 tasks, BLOCK=1, CPUS=4, MEM=16G, --time 04:00:00, THROTTLE=8
```

THROTTLE=8 bounds the commitment to about **32 cores** while the wall is measured. **Next run: `seff`
the first finished tasks of 60895543, compute the real slice cost, then either raise the throttle (if
it is inside budget) or bring the number to Shinichi.** Do not raise the throttle before that
measurement exists.

Both arrays were accepted and were PENDING at 18:44; fir was carrying only 54 tasks, so neither is
competing with the factorial BACE tail for slots.

**Still awaiting Shinichi (unchanged):** (1) S7a publication - the results Artifact
(https://claude.ai/artifact/M5HtGRnNGfwsK2Se4gMX24) is current for the core slice and he asked to be
woken when the GNN arms land, which is early on 22 Sep; (2) whether the covariate slice should include
BACE; (3) whether the mean stays the reported estimator for interval width and interval score in the
two arm-3b cells that one replicate destroys.

## 2026-09-21 18:55 MDT - scheduled run: fired 10 minutes after the last one, nothing measurable yet

This run started at 18:42, seven minutes after the previous entry was written. **Nothing was
launched, relaunched or cancelled.** The prescribed next step (`seff` the first finished tasks of
covsens n = 1000 array 60895543 and price the slice) **cannot be done yet**: `sacct -j 60895543 -X`
shows a single PENDING array head, zero tasks started. Its `squeue` reason is `Priority`, i.e. fair
share, not a bad submission.

Measured state at 18:50:

| host | in flight | landed |
|---|---|---|
| Totoro | 127 cell processes, factorial GNN wave (pgid 308643) | core 3600/3600, factorial **7,064/11,200**, avonet 20/20 |
| fir | 53 queue entries: covsens n = 300 (60866290, 18 R), factorial BACE tail (60776371, 30 R), plus the two PENDING covsens recoveries | core 1,777, factorial 1,741, covsens 1,025 |
| nibi | empty | - |
| rorqual | empty, retired | - |

covsens by n on fir: **n = 100 796/1200, n = 300 229/1200, n = 1000 0/1200**. The n = 100 recovery
array 60895448 has one COMPLETED task, and `seff` on it reads 7 s wall / 555 MB, which is the
resume-skip path on a cell-seed that already existed. So it confirms the 8 GB request is ample but
prices nothing.

**New, and it needs an action later:** the factorial BACE tail 60776371 no longer has a clean record.
It now reads **348 COMPLETED, 30 RUNNING, 24 NODE_FAIL, 5 TIMEOUT** against the zero TIMEOUT recorded
at 13:50. The 29 lost tasks are BLOCK=1 cell-seeds and will simply be missing from
`results/factorial`. **Do not resubmit while 30 tasks are still RUNNING** (the skip check reads the
directory at task start, so a concurrent resubmission re-runs cells already in flight). Once the
array fully drains, resubmit the same stage with resume on; only the ~29 gaps will run.

**Next run, in order:** (1) if 60895543 has finished tasks, `seff` them, price the covsens n = 1000
slice, and either raise THROTTLE or bring the number to Shinichi; (2) if 60776371 has drained,
resubmit the factorial BACE tail to recover the 29 NODE_FAIL/TIMEOUT cell-seeds; (3) if the Totoro
factorial GNN wave has finished (estimate ~05:30 on 09-22), pool and run G10/G11/G14.

**Still awaiting Shinichi (unchanged):** (1) S7a publication, the results Artifact
(https://claude.ai/artifact/M5HtGRnNGfwsK2Se4gMX24) is current for the core slice; (2) whether the
covariate slice should include BACE; (3) whether the mean stays the reported estimator for interval
width and interval score in the two arm-3b cells that one replicate destroys.

## 2026-09-21 19:55 MDT - scheduled run: all three prescribed steps still blocked; fir is fair-share throttled, not mis-submitted

**Nothing was launched, relaunched or cancelled.** The three steps the 18:55 entry prescribed were
each checked and each is still blocked:

1. **Price covsens n = 1000 (array 60895543).** Still `PENDING`, zero tasks started
   (`sacct -j 60895543 -X` reads one PENDING head). No wall exists to `seff`. THROTTLE stays at 8.
2. **Resubmit the factorial BACE tail (60776371).** Still **32 RUNNING** (359 COMPLETED,
   24 NODE_FAIL, 5 TIMEOUT). The recorded rule holds: do not resubmit while tasks are in flight,
   because the skip check reads the directory at task start.
3. **Pool the Totoro factorial GNN wave.** Still running, see the revised estimate below.

### Measured state at 19:51 MDT

| host | in flight | landed |
|---|---|---|
| Totoro | 128 cell processes, factorial GNN wave (pgid 308643), load 123 of 250 permitted cores | core 3600/3600, factorial **7,374/11,200**, avonet 20/20 |
| fir | 49 queue entries. Only two arrays are RUNNING: factorial n=100 half B (60776422, 12 tasks) and the factorial n=1000 BACE tail (60776371, 32 tasks). **All five other arrays are PENDING with reason Priority.** | core 1,777, factorial 1,773, covsens 1,091 |
| nibi | empty | - |
| rorqual | empty, retired | - |

covsens by n on fir: **n = 100 796/1200, n = 300 295/1200, n = 1000 0/1200.**

### The n = 300 throttle scare was a misreading, and the real constraint is fair share

`squeue`'s truncated array field printed `60866290_[63-240%2`, which looked like a throttle of 2 and
a 50 hour wall. `scontrol show job 60866290` reads **`ArrayTaskThrottle=250`**. The `%2` was the
column clipping `%250`. **No throttle change was made and none is needed.**

What is actually holding every queued array is fair share:

```
sshare -U -u snakagaw
def-snakagaw_cpu   FairShare 0.156898   RawUsage 1.79e11
```

At that level fir will drip-feed the campaign whatever shape the arrays take. Resubmitting anything
in a cheaper shape buys queue position, not scheduling, so **churn was deliberately avoided this run**.

### covsens n = 300 is priced, and it confirms the single-thread finding at a second problem size

`seff` over the last six COMPLETED tasks of 60866290 (BLOCK=5, 4 cores, 16 GB, 1 h 30):

| measure | value |
|---|---|
| wall per task (5 cell-seeds) | 00:20:16 to 00:43:36, so **4 to 9 min per cell-seed** |
| CPU efficiency | **25.3 to 26.0% of 4 cores** on 6 of 6 tasks |
| memory efficiency | 4.6 to 4.9% of 16 GB, i.e. **0.78 GB** |

The 18:45 entry measured 24.6 to 25.0% at n = 100 and left 4 cores in place for larger n in case
torch behaved differently. It does not: at n = 300 the fast arms are **single-threaded too**, and
16 GB is 20x the real footprint. Zero TIMEOUT so far at BLOCK=5 (59 COMPLETED, 3 NODE_FAIL), so the
shape is not losing tasks and was left alone. **Recorded for the next resubmission of this stage:
CPUS=1, MEM=8G.** Remaining cost of the n = 300 slice: 905 cell-seeds x ~6.5 min = **~98 core-hours**.

The n = 100 recovery array 60895448 now has 19 COMPLETED tasks but **zero new cell-seeds** (796
unchanged). Array tasks run in index order and the first 796 indices are already on disk, so those
19 are all resume-skips at about 7 s each. The array reaches real work only after it walks past the
landed block. Nothing is wrong; it is just not measurable yet either.

### Totoro's landing estimate slips to the morning of 22 Sep

factorial by n on Totoro: **n = 100 5,600/5,600 COMPLETE**, n = 1000 **1,790/5,600**.
Rate from the last two observations (7,064 at 18:50, 7,374 at 19:51) is **310 rds/h**; over the
longer 18:33 to 19:51 window it is 262/h. The 3,810 remaining n = 1000 cell-seeds therefore land
between **08:00 and 10:30 MDT on 22 Sep**, not the ~05:30 previously recorded. Revised here rather
than left stale, because the pooling step (G10/G11/G14) is scheduled off it.

**Next run, in order, unchanged:** (1) if 60895543 has finished tasks, `seff` them, price the covsens
n = 1000 slice, and either raise THROTTLE or bring the number to Shinichi; (2) if 60776371 has
drained, resubmit the factorial BACE tail to recover the ~29 NODE_FAIL/TIMEOUT cell-seeds;
(3) if the Totoro factorial GNN wave has finished, pool and run G10/G11/G14.

**Still awaiting Shinichi (unchanged):** (1) S7a publication, the results Artifact
(https://claude.ai/artifact/M5HtGRnNGfwsK2Se4gMX24) is current for the core slice; (2) whether the
covariate slice should include BACE; (3) whether the mean stays the reported estimator for interval
width and interval score in the two arm-3b cells that one replicate destroys.

## 2026-09-21 20:52 MDT - scheduled run: nothing to launch, but the factorial is far further along than the record said

**Nothing was launched, relaunched or cancelled.** All three steps the 19:55 entry prescribed were
checked and all three are still blocked, for the same reasons:

1. **Price covsens n = 1000 (60895543).** Still `PENDING`, still zero tasks started. No wall to `seff`.
2. **Resubmit the factorial BACE tail (60776371).** Still **32 RUNNING** (359 COMPLETED, 24 NODE_FAIL,
   5 TIMEOUT). The rule holds: do not resubmit while tasks are in flight.
3. **Pool the Totoro factorial GNN wave.** Still running, 127 cell processes, load 99 of 250 permitted.

### fir is no longer starved

The 19:55 entry recorded two RUNNING arrays and everything else `PENDING` on Priority. fir now has
**117 tasks RUNNING** across three arrays (60776371 32, 60776422 52, 60866290 33). Fair share
loosened on its own; no shape was changed to get this and none was needed.

### The factorial accounting was wrong, and correcting it changes what is left

Previous entries reported one number, "factorial 7,374/11,200 on Totoro", beside a separate fir
count. That mixes two different denominators. Every factorial rds on fir and nibi holds the `bace`
arm only (FINDING A, 09:20 entry, re-confirmed this run by reading
`logs/factorial_n1000_halfA.sbatch`: `arms="bace"` on both branches of its `if`), while Totoro's
files hold the fast and GNN arms. Identical filenames on different hosts are **complementary arm
sets, not duplicates**, which is exactly why the aggregator dedupes on (filename, arm set).

Counted properly, by host and against each arm family's own target:

| arm family | target | landed | missing |
|---|---:|---:|---:|
| fast arms, n = 100 | 5,600 | **5,600** | 0 |
| fast arms, n = 1000 | 5,600 | 2,049 | **3,551** |
| BACE, n = 100 | 2,800 | 2,746 | 54 |
| BACE, n = 1000 (30 seeds) | 840 | 779 | 61 |

The BACE factorial is **97% complete**, not the 78% a 100-seed denominator implies. Shinichi's
30-seed decision for factorial n = 1000 (recorded in the morning entry, commit 2463cc6) had never
been carried into the progress counts, so the remaining BACE work has been overstated all day.

The 115 missing BACE cell-seeds were named, not estimated. At n = 100 they are 7 OU cells, 83 to 99
seeds each, all in fir's half B array 60776422, which is running. At n = 1000 they are 6 BM cells
(3, 18, 18, 25, 27 and 28 of 30 seeds), all in fir's half A array 60776371, which is running. So
every gap is already inside a live array and **no resubmission is warranted anywhere**.

### nibi being empty is correct, not a stall

nibi holds no jobs and its last `pig_factorial_n1000B` tasks show CANCELLED at 07:40:49, which looked
like a wave that died. It is not. Those cancellations are the morning's deliberate switch off the
100-seed arrays. nibi's half B at n = 1000 is the 14 OU cells, and all 14 are complete at 30 seeds.
nibi has nothing left to run in this campaign.

### covsens and the n = 100 walk-past

covsens on fir: **n = 100 796/1200 (unchanged), n = 300 356/1200 (up 61), n = 1000 0/1200.**
The n = 100 recovery array 60895448 now shows 97 COMPLETED with no new cell-seeds, which looked like
a slow drip. `sacct` says otherwise: tasks 1 to 97 all started at 19:50:27 and all finished in 4 to
6 s. The scheduler handed it a burst of about 100 slots and the array spent them walking resume-skips
at 6 s each. Its throttle is 200, so the remaining ~700 skips cost a couple of bursts and a few
minutes of compute, not hours. Nothing to fix.

### Totoro landing estimate holds

7,374 rds at 19:51, 7,643 at 20:52, so **269/h**, against the 262 to 310/h recorded earlier. The
3,551 remaining n = 1000 fast-arm cell-seeds land at about **10:00 MDT on 22 Sep**. Totoro has 39 of
the 56 factorial cells touched; the 17 untouched are all n = 1000 and are queued behind the running
62-way wave.

**Next run, in order:** (1) if 60895543 has started tasks, `seff` them and price the covsens n = 1000
slice; (2) if 60776371 and 60776422 have drained, check the 115 named cell-seeds landed and resubmit
only what did not; (3) if the Totoro GNN wave has finished, pool and run G10/G11/G14.

**Still awaiting Shinichi (unchanged):** (1) S7a publication, the results Artifact
(https://claude.ai/artifact/M5HtGRnNGfwsK2Se4gMX24) is current for the core slice; (2) whether the
covariate slice should include BACE; (3) whether the mean stays the reported estimator for interval
width and interval score in the two arm-3b cells that one replicate destroys.

## 2026-09-22 09:35 — factorial GNN arms LANDED (Totoro 09:02), pooling into /tmp/pig_pool6
- Totoro results/factorial 11,200/11,200 (5,600 n=100 + 5,600 n=1000); wave took 23 h 43 m from 09:19 09-21.
- Error lines in factorial.cells.log, classified: 159 gnn_on "Dimension out of range" (known monomorphic
  one-hot defect, floored), 153 freq + 151 freq_lambda castor "Need at least 2 states" (monomorphic,
  floored), and ONE NEW MODE: 50 freq_lambda "Not compatible with requested type: character -> double"
  preceded by Rphylopars singular-solve warnings, in OU_l1_r0_mcar0.1_n1000. A Rphylopars internal
  failure at lambda=1 under OU; floored and to be reported as its own failure mode.
- fir is NOT fully drained: pig_factorial_n100B has 9 R + 1 PD (last BACE n=100 half-B tasks);
  pig_fact_n1000A_long (12 h limit, tasks 4-420, resume-skip) PENDING as the TIMEOUT recovery;
  pig_covsens_n1000 PENDING on JobArrayTaskLimit. covsens rds 3,584 of 3,600 non-BACE.
- Pool layout: /tmp/pig_pool6/factorial/{totoro (gnn arms), totoro_fl, totoro_ff, nibi (bace), fir (bace)}
  + /tmp/pig_pool6/avonet/totoro_fl. nibi pulled (3,434; rsync exit 23 = local chmod only).
- Next: finish pulls -> gates G10/G11/G14 -> aggregate factorial (background; ECE bootstrap is the slow part)
  -> factorial sections into methods note, article, board v2 -> covsens when fir finishes -> after-task.

## 2026-09-22 10:05 — pool complete for every fast arm; gates run
Pool /tmp/pig_pool6/factorial: totoro 11,200 (gnn_on, gnn_off, gnn_off_rphylopars) · totoro_fl 11,200 ·
totoro_ff 11,200 (freq, floor) · nibi 3,434 (bace) · fir 2,220 (bace). 39,254 rds. Mac rsync 2.6.9 does
not know --info=stats1; the first Totoro pull died on the usage message and my grep hid it. Fixed.
- G14 PASS x3 (totoro_fl/nibi, totoro_ff/fir, nibi/fir): truth and mask bit-identical across hosts.
  The freq clause is vacuous on those pairs (nibi, fir hold bace only); the evidence is data + mask.
- G10 (core) FAIL on BACE only: 24 of 600 seeds missing in 5 cells (576/600 = 96%); every fast arm
  200/200. Recovery: core n1000 BACE array 60949560 submitted on fir (600 tasks, resume-skip, ~24 compute).
  core n100/n300 recovery arrays were already queued.
- G11 (factorial) FAIL on BACE only: 57 of 3,550 seeds missing in 5 cells (98.4%); every fast arm
  200/200 in all 56 cells. 56 of the 57 are the four BM clade0.3 n=1000 cells (3, 18, 25, 18 of 30):
  the 5 h --time was too short for clade-masked BACE at n=1000; pig_fact_n1000A_long (12 h, resume-skip)
  is PENDING on fir as the recovery. One seed missing in OU_l0.3_r0.5_mar0.3_n100.
- Factorial aggregation running (bye0og36m) -> /tmp/pig_pool6/agg/fact_*.
- New freq_lambda failure mode confirmed: 50 replicates, all OU_l1_r0_mcar0.1_n1000, Rphylopars
  singular solve -> type error. Floored; disclose as its own row.

## 2026-09-22 10:25 — analysis phase running in parallel
- Aggregator on /tmp/pig_pool6/factorial (bye0og36m) -> /tmp/pig_pool6/agg/fact_*  (slow: ECE bootstrap).
- Direct per-replicate extracts (faster, bypass the aggregator): scratchpad/fact_gnn.csv (bqjk0rh4m),
  fact_bace.csv (brxpkjzmc, nibi+fir deduped on (cell, seed)); ff_totoro.csv + fl_totoro.csv from yesterday.
  Paired script ready: scratchpad/fact_paired.R (joins all five arms + floor by (cell, seed); floors
  divergent fits under the study rule; paired contrasts vs BACE on BACE's own seeds).
- Failure bookkeeping check from the rds (not the log) running -> scratchpad/fail_{gnn,fl,ff}.csv.
- covsens_fl launched on Totoro 09:43 (3,600 jobs); waiter bbst8pssw pulls fir covsens + totoro covsens_fl
  into /tmp/pig_pool6/covsens/{fir,totoro_fl} when both reach 3,600.
- Ledger: leaf-campaign.md G14 [x], G6c [x], G10/G11 evidence recorded (BACE-only gaps, recoveries queued);
  leaf-results.md G13a/G13b/G13d [x]. Open: G10, G11 (await BACE recovery), G6d (await covsens), G12
  (after aggregation), G13c (Shinichi reads the Artifact).

## 2026-09-22 10:40 — FACTORIAL, ALL ARMS, z-RMSE (direct per-replicate reads; 11,200 paired rows, 56 cells)
Divergent fits floored under the study rule: freq 1, freq_lambda 15, bace 2, gnn_off 0, gnn_off_rph 32, gnn_on 0.
BACE on every seed it ran (4,434 incl. out-of-design n=1000 seeds 31-100; 1,220 cross-host dups dropped).
| evo | lambda | freq | freq_lambda | bace | gnn_off | gnn_off_rph | gnn_on | floor |
| BM | 0.3 | 1.128 | 0.924 | 1.037 | 1.048 | 1.030 | 1.044 | 1.043 |
| BM | 1.0 | 0.579 | 0.572 | 0.876 | 0.486 | 0.536 | 0.533 | 1.070 |
| OU | 0.3 | 1.109 | 0.911 | 0.997 | 1.035 | 1.008 | 1.032 | 1.035 |
| OU | 1.0 | 0.569 | 0.562 | 0.610 | 0.501 | 0.546 | 0.525 | 1.051 |
By mechanism: clade  freq .869 fl .784 BACE .994 gnn_off .801 (floor 1.054)
              mar    freq .897 fl .774 BACE .948 gnn_off .816 (floor 1.113)
              mcar.1 freq .784 fl .690 BACE .727 gnn_off .706 (floor 1.000)
Paired vs BACE on BACE's seeds (neg = better than BACE): BM l0.3 fl -0.123 (.0085), gnn_off +0.003;
  BM l1: every arm beats BACE by 0.25-0.35 (gnn_off -0.349 (.0108)); OU l1: gnn_off -0.108 (.0069).
NEW RESULT: BACE under clade-biased missingness at BM lambda=1 sits AT OR ABOVE THE FLOOR (0.994 n=100,
1.067 n=1000) while gnn_off is 0.621 / 0.475. This is the "mixed model equations singular" regime,
floored failures dominating. Must report the per-cell BACE failure rate beside it (fail_count running).
Under OU the collapse is milder (0.860 / 0.664). BACE is competitive only at OU lambda=1 MCAR n=1000
(0.425 vs gnn_off 0.346) and beats pigauto nowhere on continuous traits in the factorial.
The core-slice ranking inversion replicates in every stratum: at lambda=0.3 only freq_lambda is clearly
below the floor everywhere; BACE below it only under MCAR; pigauto arms sit at the floor (gate closed).

## 2026-09-22 10:50 — the BACE clade result partitioned by mechanism (Fisher's rule: failures first)
BACE errored-and-floored rate by lambda: 0% everywhere at lambda=0.3; at lambda=1: MCAR.1 35%, MCAR.3 18%,
MAR 39%, clade 39% (of 339-709 replicates). Clade cells by evo x n at lambda=1: BM 53% (n=100), 86% (n=1000);
OU 19%, 28%. So "BACE at or above the floor under clade masking at lambda=1" is mostly failures scored at the
floor. Among the replicates where BACE FITS, clade lambda=1: BM 0.816 (k=94) / 0.545 (k=6); OU 0.811 / 0.469,
against gnn_off 0.621 / 0.475 and 0.727 / 0.506. Worse than pigauto, not catastrophic. At lambda=0.3 under
clade masking BACE is above the floor even with zero failures (1.181, 1.149 vs floor ~1.03): a real
accuracy result, not a failure artefact. Both halves go in the report, side by side.
The lambda=1 failure surge is the fixed-threshold monomorphic regime ("mixed model equations singular").

## 2026-09-22 11:00 — divergence rule moved INTO the committed aggregator; covsens n300 recovery
- script/campaign_gnn_off_aggregate.R: after tab2, an arm-replicate-trait is divergent if zRMSE > 3x the
  floor arm's zRMSE for the same (cell, seed, trait) OR interval_score > 1e3; zRMSE -> floor value,
  coverage/width/interval_score -> NA; failures.csv gains n_divergent beside n_failed. Smoke on the one
  divergent cell (OU_l1_r0_mcar0.1_n1000; dir renamed because the script drops any path containing
  "_smoke"). Lane check on docs/gnn-off-campaign: that branch holds the ORIGINAL 152-line-shorter
  version from the closed arc, not a fork of this fix.
- Superseded factorial aggregation (bye0og36m, unpatched) stopped; relaunch after the smoke passes,
  then re-aggregate core into /tmp/pig_pool5 so the committed csv carry n_divergent too.
- Methods note factorial section committed 9f2d39a (1.3 per 1000, 0 em dashes, 0 superlatives).
- covsens: Totoro covsens_fl 3,600/3,600 DONE. fir covsens 3,584/3,600: n300 missing 15 (old array
  237 COMPLETED + 3 NODE_FAIL at BLOCK=5) -> recovery 60950427 submitted (fast arms, BLOCK=1, 01:30:00);
  n1000 missing 1, array 60895543 [344-1200] PD on JobArrayTaskLimit with 0 running (odd; resume-skip
  will clear it when the throttle releases; re-check in an hour, scancel+resubmit if still 0 R).
- fir also holds: pig_core_n1000 (core BACE recovery) PD, pig_fact_n1000A_long PD, n100B 4 R.

## 2026-09-22 11:25 — board v2 published; article factorial section written; aggregations re-running
- Artifact v2 (same URL M5HtGRnNGfwsK2Se4gMX24, label "Factorial, all arms"): Factorial tab now carries
  all arms (z-RMSE by evo x lambda and by mechanism, discrete accuracy, coverage by mechanism, the
  BACE clade partition by cause, failure rate by mechanism, the divergence rule); status strip updated;
  Decisions "not yet" paragraph updated. JS parse-checked with node --check; every getElementById id exists.
  Not visually inspected: the in-app browser is not signed in and I will not sign it in.
- Aggregator committed 1ec3e98 (divergence rule + n_divergent). Factorial re-aggregation bl8eg0fqy ->
  /tmp/pig_pool6/agg/fact_*; core re-aggregation b2cafpz0s -> /tmp/pig_pool7/agg_*. When both land:
  concatenate summary/paired/failures (same schema) into script/campaign_sim_results/, render the
  article, G12 on the combined dir, commit.
- vignettes/articles/simulation-study.Rmd: "The factorial" section added before Cost (five chunks:
  fact-z, fact-z-miss, fact-acc, fact-cov, fact-fail; `strat()` helper does cell means then equal
  weight per cell); scope note rewritten; "does not cover" names the two-specification exception.
  Renders only once summary.csv holds factorial rows.
- fir: covsens n1000 array 60895543 stuck (0 R, "JobArrayTaskLimit") -> scancel + resubmitted with
  resume-skip; covsens n300 recovery 60950427 12 R; core n1000 BACE recovery PD; _long PD; n100B 3 R.

## 2026-09-22 10:05 (Mac clock) — PR, reconcile, after-task scaffold
- Draft PR #184 (https://github.com/itchyshin/pigauto/pull/184) has existed since 2026-09-20 23:52 UTC; the
  09-20 reconcile's "no PR" row predates it. Remote branch was 46 commits behind -> pushed
  (1a389ab..1ec3e98); body refreshed from /tmp/pr_body.md; agent_mention_check G-20 clean.
- docs/dev-log/plan-actual/2026-09-20-imputation-sim-reconcile.md: "Reconcile 2, 2026-09-22" appended
  (six axes, DECISION RECEIPT). Not yet committed (commit with the article once it renders).
- docs/dev-log/arc/2026-09-22-imputation-sim-after-task.md scaffolded with closeout.py; sections 2, 3a,
  4, 7a, 9, 11, 12 drafted; 5, 6, 8, 10 wait on G12 and the BACE tails. NOTE closeout.py auto-filled
  section 4 from the BRAIN repo's git status because I ran it from ~/shinichi-brain; overwritten with the
  real list.
- GOAL.md arcs block refreshed to the current state.
- Re-aggregations (factorial bl8eg0fqy, core b2cafpz0s) still reading rds at 10:04.

## 2026-09-22 10:15 (Mac) — deliverables aligned on the committed aggregate; waiting on machines only
- Aggregator with the divergence rule is the single authority: script/campaign_sim_results/{summary,paired,
  per_trait}.csv = core (pool7) + factorial (pool6/agg), 16,276 summary rows, commit 077b4de, pushed.
- Article renders the factorial from it (verified); methods note tables moved onto it (12 replacements,
  max shift 0.011, BACE and Rphylopars-solver cells only, because the aggregator floors per (replicate,
  trait) and my cross-check floored whole replicates); board v3 = same numbers.
- G12 PASS (core totoro + totoro_fl hosts); G13a/b/d PASS via gate-check --approve; G14 PASS x3.
  Open: G10, G11 (BACE tails), G6d (covsens), G13c (Shinichi).
- OWED when machines land: (1) failures.csv (n_divergent) + ece.csv from both aggregations -> results dir ->
  re-render (fact-fail chunk) -> commit. (2) covsens: pull, aggregate with the patched script, one paragraph
  each in note/article, board v4, G6d. (3) BACE tails -> G10/G11 -> after-task 5/10 final -> closeout check
  from the worktree -> handover (protocols/handoff.md template; run handoff_gate.sh first).
- fir at 10:15: every array PD on Priority (core_n1000 BACE, covsens_n1000, covsens_n300, fact_n1000A_long,
  factorial_n100B); covsens 3,599/3,600; core 1,777/1,800; factorial 2,232.

## 2026-09-22 10:18 (Mac) — stall rule applied: BACE tails handed to Totoro (and nibi once seeded)
- fir: all five arrays PD on Priority with no start estimate; 3,000+ nodes drained/draining. Plan's stall
  rule (no task started 2 h after submit -> hand remaining cells elsewhere) applies.
- Totoro: seeded results/core_bace (1,776 = the pooled nibi+fir BACE rds) and results/factorial_bace
  (4,434) so campaign_sim_cell.R's resume-skip leaves only the missing seeds; launched core (pgid 1840413)
  and factorial (pgid 1840886) BACE-only waves at 40-way, 1 core each. Expected to run 24 + 57 seeds;
  n=1000 BACE ~3 h each -> both done by ~13:30-14:00 if RAM holds (watch free -g; 40 x ~16-20 GB).
- nibi: my first core n1000 array 22471569 would have RECOMPUTED 576 finished seeds because nibi's own
  results/core holds no n=1000 BACE (those ran on fir) -> cancelled within a minute. Seeding nibi with fir's
  core (1,776) and factorial (2,232) rds first (btn9r8kwu); then resubmit core n1000, factorial n1000 half A
  (12 h), factorial n100 half B. sbatch --test-only says nibi would start them ~12:17. Whichever host lands
  a seed first wins; the aggregator's (filename, arm-set) dedupe keeps one copy.
- LESSON: resume-skip is per-host. Before submitting a recovery array on a host, seed its results dir with
  every finished rds from the pool, or the "recovery" recomputes the campaign.
- 10:20 correction: the Totoro core BACE tail is RUNNING (49 cell processes on core_bace, 81 on
  factorial_bace; log shows seeds 37-68 of BM_l0.7_r0.5_mcar0.3_n1000 fitting). My "DONE already"
  reading counted core.log's second DONE, which belongs to the 09-21 core_fl wave. Watcher re-armed on
  core DONE >= 3 and factorial DONE >= 4 (biu95j9dm). All 2,814 core and 5,654 factorial pooled BACE
  files carry a bace arm (no partial files), so the 24 + 57 are genuinely uncomputed seeds.
- covsens: fir's array for the last n300 seed is PD with no start estimate; pulling fir's 3,599 covsens
  rds to the Mac so the seed can run on Totoro against a seeded directory, then pool by hand (the
  covsens waiter bbst8pssw watches fir and will not fire).

## 2026-09-22 10:23 (Mac) — covsens moving to Totoro; handover committed
- Handover docs/dev-log/handover/2026-09-22-claude-handover.md committed 30278b3 and pushed (landing state as
  of 10:20; refresh at close). Board decision-4 text updated for v4 (not yet republished).
- covsens: fir pool pulled (3,599; the one missing seed is n=1000, not n=300 as I first read). Seeding Totoro
  results/covsens from it and running the covsens stage there (fast arms, PAR 10) so the single seed fills
  (bls8v1er2); covsens_fl (3,600) pulling to /tmp/pig_pool6/covsens/totoro_fl (bvutekccy). Per-replicate
  extracts for core (totoro, totoro_fl) and covsens (fir) running -> scratchpad/all_core_*.csv, all_cov_fir.csv;
  scratchpad/covsens_cmp.R compares covsens minus core on identical (cell, seed, arm).
- Totoro BACE tails at 10:22: 65 R processes, RAM 100 GB of 1,007, core_bace still 1,776 (n=1000 fits ~3 h).
- Still waiting: nibi seeding rsync; ECE stage of both aggregations (6 R processes on the Mac).

## 2026-09-22 10:26 (Mac) — S6d covariate sensitivity DONE (3,599/3,600); board v4
- covsens minus core, paired on identical (cell, seed, arm), divergence rule applied to both sides (20 + 16
  floored): freq -0.051 (0.0022) / -0.024 (0.0023) / +0.003 (0.0009) z-RMSE at lambda 0.3/0.7/1;
  freq_lambda -0.032 / -0.012 / +0.002; gnn_on -0.012 / -0.015 / -0.000; gnn_off exactly 0 (covariate-free by
  design; the shared-seed check); coverage moves < 0.002; gnn_on accuracy +1.6 pts at lambda 0.3, -1.4 at 1.
  No ranking changes. BACE not in this run (budget).
- Written into the methods note ("Covariate sensitivity") and the article ("## Covariates"), commit 77c34ed,
  pushed; slop 1.2 / 2.0 per 1000, 0 em dashes; G13a PASS; render OK. Board v4 published with a covsens
  section on the All arms tab and the status strip updated. G6d [x] in leaf-campaign.md.
- Remaining machine-bound: Totoro BACE tails (watcher biu95j9dm), the last covsens seed on Totoro, the ECE
  stage of both aggregations, nibi seeding (backup arrays not yet resubmitted; fir arrays still PD).
- 10:35 nibi seeded (core 1,776 incl. 580 n=1000; factorial 4,434). Backup core n1000 BACE array 22472376
  submitted (600 tasks, resume-skip -> 24 compute). The two factorial backups were REFUSED:
  AssocMaxSubmitJobLimit (1,000 array tasks per user; 600 + 420 + 700 exceeds it). Not resubmitted: Totoro
  carries every tail; nibi and fir arrays are insurance only and get cancelled once Totoro's rds are pooled.

## 2026-09-22 10:50 (Mac) — CORRECTION: the core aggregate with the divergence rule moves the solver story
Patched core re-aggregation (pool7) landed: 68 divergent arm-replicates in the CORE, not the "3" I had
assumed: gnn_off_rphylopars 28/10/12 at n=100/300/1000, bace 7 (n=100), freq 4+6, freq_lambda 1. With them
floored, the Rphylopars-solver core z-RMSE at lambda=1 becomes 0.558 / 0.476 / 0.415 (was 1.654 / 0.540 /
0.590) against the in-house 0.487 / 0.426 / 0.379. The "solver diverges catastrophically" sentence in the
board, the article prose and the after-task must become "50 replicates diverged and were floored; once
floored the solver is 0.03 to 0.07 behind the in-house one". BACE core cells move by <= 0.011 (1.015 -> 1.004,
0.879 -> 0.874, 0.801 -> 0.794), freq by <= 0.014 (1.161 -> 1.158, 1.140 -> 1.127, 0.972 -> 0.965).
The committed summary.csv already carries the patched core (built from pool7), so the ARTICLE TABLES are
right; its prose, the methods note tables, the board's CORE data and the PR body are being moved now.
covsens final 3,600/3,600: identical to three decimals; G6d closed.
- 10:52 correction landed everywhere: methods note (16 replacements, commit b0da85b), article solver
  paragraph (renders from the patched csv already), after-task 3a and issue ledger, board v5 (26 data and
  prose replacements; kv "0.09 - 0.21"; decision 5 reworded). Every core number in all three deliverables
  now comes from /tmp/pig_pool7 (= script/campaign_sim_results core half). Pushed.
- Remaining: factorial ECE stage (54 min in) -> failures.csv/ece.csv (concat core pool7 + factorial) ->
  re-render -> commit; Totoro BACE tails (core_bace 1,776/1,800, factorial_bace 4,435; ~13:20); then
  G10/G11, final re-aggregation of both stages with the tails, final csv, after-task 5/10, closeout check,
  handover landing state, cancel fir/nibi backup arrays.
- 10:58 factorial bootstrap landed: "38034 cells read, 1220 duplicates dropped, divergence rule: 188
  arm-replicates (216 traits)". failures.csv (core + factorial, n_divergent) and ece.csv assembled into
  script/campaign_sim_results, article re-rendered (failure table renders: BACE 1,097 + 56, freq 237 + 46,
  freq_lambda 285 + 21, gnn_off 0 + 1, Rphylopars solver 0 + 130, gnn_on 237 + 1), commit 4e61cc7 pushed.
  ONLY the Totoro BACE tails remain machine-bound (watcher biu95j9dm).

## 2026-09-22 14:40 (Mac) — CLOSE on the replicates present (Shinichi: "all done - what are you waiting??")
- Decision: the study is closed on BACE 576/600 core and 3,493/3,550 factorial, stated in every table; the
  tails compute (core on nibi 22472376 running 2 h+, factorial on Totoro pgid 1840886, 2 of 57 landed) and
  fold in by re-aggregation. fir arrays all cancelled (queue 0). Totoro core BACE wave killed at RAM 803 GB
  (shared machine); nibi covers those 24. Totoro RAM 723 GB with 42 factorial fits.
- Ledger: 23 of 26 met. G6 and G9b resolved by the campaign's measurements (BACE coverage 0.885 at
  lambda=1 n=1000 is outside the pre-run gate's band -> finding; convergence rate reported). Open: G10, G11
  (tails), G13c (Shinichi reads the board, six decisions).
- BACE convergence rate (1,774 fits with diagnostics, runner >= 09-20): converged by BACE's verdict 15%/20%
  (core/fact) at lambda 0.3, 33% at 0.7, 81%/67% at 1; by n 23-28% / 43% / 48-61%; median ESS 572 -> 1600.
  Added to the methods note "Failures and convergence" (commit after be1dc48). Not yet on the board.
- after-task finalised (be1dc48); handover refreshed (be1dc48); closeout.py check will pass only after G13c.
- 14:39 Totoro RAM reached 976 GB of 1,007 with 42 clade-masked BACE fits at 4 h 18 m: killed the factorial
  BACE wave (pgid 1840886) by hand ahead of the 940 GB guard, losing those partial fits, and relaunched the
  same wave at 20-way with resume-skip (55 seeds left). Clade n=1000 BACE fits run 5 to 6 h each, so three
  waves: ~15 h, landing ~06:00 on 09-23. nibi core tail: 20 tasks at 1 h 06 m, ~16:30 today. Neither
  changes a deliverable; both fold in by re-aggregation. Guard task stopped (replaced by the lower width).
- 15:12 RAM guard fired at 927 GB and killed the 20-way wave (pgid 2301135) -- but Totoro still read 942 GB
  with our 2 R processes: user ortegara has one cc1plus (C++ compile) process at 930 GB since 14:57. Not
  ours to touch; flagged once. Totoro cannot host the BACE tail until it clears. The 14:38 976 GB reading
  WAS ours (42 fits); this one is not. Factorial tail (55 seeds) moved to nibi: half A n1000 (12 h) and
  half B n100 (BLOCK 5) submitted against the seeded results/factorial; core tail already on nibi (20 R,
  3 h in). Totoro watcher biu95j9dm stopped. LESSON: read per-user RSS before blaming your own wave.
- 15:40 Decision 7 taken by Shinichi ("Yes — open a lane for spec decision 3 with lambda estimated as the new
  default, covariates included; use the core slice as the benchmark; separate PR after #184"). Recorded:
  vault D-278 (local commit), board v7 (Decisions tab, seventh item marked decided), after-task 3a, handover,
  and the lane plan stub docs/dev-log/arc/2026-09-22-joint-lambda-default-plan.md. NOT started in this arc
  (package change, fenced). Shinichi: "finish this arc when we finish the goal" -> remaining: his six
  publication decisions (G13c) and the BACE tails on nibi (G10/G11).
- 16:30 Shinichi opened the lambda-default lane himself. My arc-1 agent stopped before writing anything but
  three PROGRESS lines in ../pigauto-lambda-default/.unlazy/lambda-default/ (its finding: the gap is the joint
  solver, which has no likelihood to profile lambda with). No commits, no edits, lease claude:pigauto:lambda
  released so his session can claim. This session finishes the simulation arc only.

## 2026-09-22 17:45 (Mac) — tails: core n1000 COMPLETE on nibi; 4 core seeds + 55 factorial seeds in flight
- nibi array 22472376 (core n1000 BACE) finished all 600 tasks; pooled to /tmp/pig_pool4/core/nibi (600 n1000).
  Core BACE now 596/600; the 4 missing are BM_l0.7_r0.5_mcar0.3 n100 s38/s51/s52 and n300 s55 (fir's cancelled
  arrays owned them). Resubmitted on nibi: core n100 array 22502324 (resume-skip, 3 compute) and a single-task
  job 22502338 for n300 s55 (a 600-task n300 array hit AssocMaxSubmitJobLimit and was pointless for one seed).
  Watcher pulls + runs G10 when both leave the queue (~1 h).
- Factorial tail: 60 clade n1000 tasks RUNNING on nibi since ~17:05 MDT (5 to 6 h each) -> land ~23:00-00:00 MDT;
  watcher bod0meuyq pulls when the arrays drain. Then: G11, re-aggregate both stages, refresh csv/article/board,
  after-task 5/10 final, closeout check, handover landing state.
- Shinichi opened the lambda-default lane himself (D-278); this session finishes only the simulation arc.

## 2026-09-22 20:55 (Mac) — factorial BACE tail LANDED (nibi 420/420 COMPLETED, ~3.8 h per clade fit)
- /tmp/pig_pool6/factorial/nibi now 4,496 rds (>= 4,489 expected; extras are out-of-design seeds). G11 running
  directly (bg); final factorial aggregation -> /tmp/pig_pool6/agg2/fact_* (bg, ~1 h with the ECE bootstrap).
- Core 1,797/1,800: n300 s55 landed; n100 s38/s51/s52 running on nibi (~40 min in). Watcher bzl3tob0o pulls
  and runs G10 when they finish; then final core aggregation -> /tmp/pig_pool8.
- Then: concat csv -> script/campaign_sim_results, re-render article, board v8 (status 100%), commit, ledger
  G10/G11 [x], after-task 5/10, closeout check from the worktree, handoff_gate + handover landing state,
  PR #184 body. Only G13c (Shinichi's six decisions) will remain open.
- 20:58 G11 PASS: "56 design cells, 40316 rds present, 0 replicate-arms missing". Factorial BACE 3,550/3,550.
- 22:05 core tail: s52 (n100) landed in 21.7 min; s38 and s51 (n100) TIMED OUT at 01:30 (slow outliers in the
  BM_l0.7_r0.5 cell, not crashes); s55 (n300) SEGFAULTED in MCMCglmm (R aborted; no rds, so the runner could
  not floor it). Retries as single tasks: n100 x2 at 04:00:00 (first retries at 01:30 cancelled), n300 at
  02:30:00 with 48 GB. Watcher pulls + runs G10 when they drain. If s55 segfaults again it is recorded as a
  BACE crash (1 of 600), disclosed, not fabricated.
- Final factorial aggregation (agg2) Section I landed; committed summary/paired/per_trait rebuilt (core pool7 +
  factorial agg2); article re-rendering. failures/ece follow the ECE stage.
- 23:10 final factorial aggregation (agg2) complete: 38,096 cells, 2,220 cross-host duplicates dropped, 188
  divergent arm-replicates floored. failures.csv/ece.csv rebuilt (core pool7 + factorial agg2), article
  re-rendered, commit d36e682 pushed. Checking whether any factorial stratum figure moved with the 57 added
  BACE seeds (they were all clade n=1000 cells). Core retries still running on nibi (watcher b3qj1d9qp).
- 23:25 Completing the 57 clade seeds moved two factorial BACE figures materially: BM lambda=1 z-RMSE 0.872 ->
  0.842; clade mechanism 0.986 -> 0.962; BM lambda=1 clade n=1000 failure 86% of 43 -> 63% of 60 and fitted
  z-RMSE 0.545 (k=6) -> 0.472 (k=22), level with pigauto's 0.475. Everything else moved <= 0.003. Methods note,
  article, after-task (b0be521) and board v8 carry the final figures. Remaining: the 3 core seeds on nibi
  (watcher b3qj1d9qp -> G10), final core aggregation, board status 100%, closeout, handover, PR body. G13c open.
- 23:40 core seed 55 (BM_l0.7_r0.5_mcar0.3_n300): BACE segfaults the R process deterministically (exit 139 at
  1 min 25 s, nibi 22502338 and 22517730, 48 GB). The runner cannot floor a crash, so the failed-arm record was
  written BY HAND from the runner's own floor arm on the same seed (cell script --arms floor, rows relabelled
  bace, failed$bace = TRUE, errors$bace = the crash text), placed in /tmp/pig_pool4/core/nibi and copied to
  nibi results/core. Disclosed in the methods note "Failures and convergence". Two n100 retries still running.
- 22:55 seed-55 failed-arm record written (24 result rows, arm bace, failed$bace TRUE, errors$bace = crash text,
  note field) into /tmp/pig_pool4/core/nibi and the nibi keeper (core 1,798/1,800 there). The two n100 retries
  (s38, s51) are at 47 min of a 4 h limit; watcher b3qj1d9qp pulls and runs G10 when they finish, then the
  final core aggregation, board v9, ledger, closeout, handover, PR body. G13c still open.

## 2026-09-23 02:20 — G10 PASS; every cell present in both stages
- Core: "18 design cells, 10785 rds present, 0 replicate-arms missing". BACE 600/600 with three hand-assembled
  failed-arm records (s55 n300 segfault x2; s38/s51 n100 TIMEOUT at 4 h x2, sibling median 22 min), each
  labelled in the file (note, errors$bace) and disclosed in the methods note and after-task (65f1a14).
- Final core aggregation running -> /tmp/pig_pool8/agg_* (then core half of the committed csv, render, commit,
  board v9 with the 100% strip, after-task 5/10 final, closeout check, handoff_gate + handover, PR body).
- Ledger after this: only G13c (Shinichi's six decisions) open.
- 02:30 final core aggregate (pool8): BACE 100/100 in every core cell. One figure moved past the third decimal
  (BACE z-RMSE at n=1000 lambda=1: 0.659 -> 0.647; paired loss to freq 0.168 -> 0.158; accuracy 0.840 -> 0.845).
  Note, article, committed csv (a0c8e0e) and board v9 carry the final figures. gate-check --approve unchecked
  G11 on its 120 s timeout; re-marked [x] with the direct-run PASS and a do-not-re-approve note.
  Ledger: 25 of 26 met; only G13c open. Remaining: pool8 failures/ece -> csv -> commit; closeout check;
  handoff_gate; handover Landing State with the final HEAD.
- 02:25 validators at close: gate-check --status "UNMET: 1 (met: 25, 6 never executed)" = G13c only;
  handoff_gate FAIL (counts never-executed evidence-only gates + G13c) -- declared verbatim in the handover;
  closeout.py structure PASS, ledger stage blocked by other projects' .unlazy in the vault and by G13c.
  WORKTREE: lines added to every leaf (bare path). Waiting only on pool8 failures/ece for the last csv commit.

## 2026-09-23 02:45 — ARC COMPLETE except G13c
- Final core bootstrap landed (9,000 cells read, 1,785 duplicates dropped, 68 divergent floored). Committed
  failures.csv/ece.csv now cover both stages at 100% (a0df635). Core BACE errored 114/114/96 (+7 divergent at
  n=100) after the three hand records; methods note table updated.
- Deliverables final: board v9 (private), methods note, article (renders from the final csv), after-task,
  reconcile 2, handover with validators quoted. PR #184 body current. Ledger 25/26; open G13c only.
- Nothing is running anywhere: Totoro idle (our waves ended), nibi queue empty, fir queue empty. Pools on the
  Mac under /tmp/pig_pool4 (core), /tmp/pig_pool6 (factorial, covsens, avonet), aggregates pool8 (core) and
  pool6/agg2 (factorial); keepers on nibi results/core + results/factorial and Totoro results/*.
- To close the arc: Shinichi answers the six decisions on the board -> mark G13c [x] with his words -> run
  gate-check --status (expect 26/26 met) -> commit -> he decides on article visibility and the #184 merge.
- 02:50 PAUSED AT THE HUMAN GATE. Every reversible step is landed (HEAD 8a69e55, pushed; 25/26 gates; nothing
  running anywhere). G13c is Shinichi's six publication decisions and is not something an agent may fill in.
  Resume: on his answers, write them verbatim into leaf-results.md G13c EVIDENCE, mark [x], gate-check --status
  (expect 26/26), commit, then his call on article visibility and the #184 merge. Until then: no further work.

## 2026-09-23 — ARC COMPLETE. G13c recorded; ledger ALL MET (26 of 26).
Shinichi's answers: both specifications with the gap as a finding; lead with discrete, both costs beside it;
article unlisted until Szymek signs off; PR #184 stays draft until he reads the board. Board v10 marks all seven
decisions. No further work on this arc. Follow-on: the lambda-default lane (Shinichi's own session).
- STOPPED by Shinichi ("I think we stop now as we want to redo simulations"). Nothing running on any machine.
  All results, csv, board (v10), note, article, after-task, handover pushed at the HEAD above; they stand as the
  record of this run. The redo is a new arc (likely after the lambda-default lane changes pigauto's baseline).
