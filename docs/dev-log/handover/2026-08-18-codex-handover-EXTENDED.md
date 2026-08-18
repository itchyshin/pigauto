# EXTENDED HANDOVER: pigauto, 15–18 Aug 2026 (Claude → Codex)

Meta: 2026-08-18 · `AUTHOR = claude` · `TARGET = codex` · supersedes/extends
`2026-08-18-codex-handover.md` (that one is the 5-minute version; **this one is lossless**).

**You inherit no chat.** Four days, three arcs, six merged PRs, ~2,600 simulation jobs, and
one unmerged branch. Everything below is grounded in committed artefacts you can open.

---

# PART 1 — THE ONE-PARAGRAPH STORY

pigauto lost to BACE and to raw Rphylopars on continuous traits. Over three arcs we found out
why, and it was **not** the neural network. The joint baseline **estimated a cross-trait
covariance Σ and then discarded it at prediction time** (`max_iter = 0` short-circuits to
per-column BM), so "joint" was joint in estimation only. Fixing the prediction step — by
computing the *exact* matrix-normal conditional in the precision form, using
Hadfield & Nakagawa (2010) sparse precision — recovers **50–60% of the AVONET gap on 3 of 4
continuous traits with no new dependency**. That work is on `arc/exact-conditional`, green,
pushed, and **not yet PR'd**. Separately, the GNN's headline claim was narrowed to something
defensible, the conformal-interval story was diagnosed and partially repaired, and the first
working external comparisons in the repo's history were produced.

---

# PART 2 — REPO STATE (exact, verified at handover)

| Item | State |
|---|---|
| `origin/main` | **`f925e17`**, version `0.10.0.9000` |
| Open PRs | **0** |
| Merged this arc | **#167, #168, #169, #170, #171, #172** |
| Docs branch | `handover/2026-08-09-cursor`, fully pushed |
| Unmerged branch (ready) | **`arc/exact-conditional`** — 6 commits, pushed, remote == local (`ba255ee`), suite **FAIL 0 / PASS 2180** |
| Unmerged branch (benches) | `arc/bace-comparators` (pushed), `arc/pertype-benches` (`4586b62`, **unpushed**) |
| Worktrees | `.worktrees/arc-exact`, `arc-bace`, `arc-pertype` |
| Dirty protected files | **exactly 15** — must stay unstaged |
| Running jobs | **none** (local and Totoro both idle) |

**The 15 protected files** (never stage; another lane owns them): `.gitignore`, `AGENTS.md`,
`CLAUDE.md`, `README.md`, `_pkgdown.yml`, `dev/gnn_attribution_low_lambda_smoke.{R,log,md}`,
`script/bootstrap_gnn_attribution_fir.sh`, `script/collect_gnn_attribution_low_lambda.R`,
`script/gnn_attribution_low_lambda_ladder.R`, `script/returned_gnn_attr/`,
`script/rsync_gnn_attribution_low_lambda.sh`,
`script/run_gnn_attribution_low_lambda_totoro.sh`,
`script/submit_gnn_attribution_low_lambda_fir.sh`.

**`docs/` is gitignored** → commit docs with `git add -f` on explicit paths. **Never
`git add -A`.**

---

# PART 3 — WHAT WAS MERGED, AND WHY IT MATTERS

| PR | What | Why it matters |
|---|---|---|
| **#167** | Silent-fallback honesty | λ silently dropped under covariates; covariates silently ignored by joint baselines (P1-8); the `n_val < 19` warning understated an *arithmetic impossibility*; `conformal_split_val` exposed in `impute()` |
| **#168** | `conformal_method = "mondrian"` | Locality-stratified conformal intervals; repairs the confirmed MAR/MNAR exchangeability failure at n≥1000. Falls back below 19 residuals/stratum (that floor is load-bearing — see §6) |
| **#169** | `joint_solver = "rphylopars"` | Opt-in converged-REML joint solver; the **measured yardstick** that quantified the gap |
| **#170** | Per-type λ dispatch | `lambda_mode != "fixed_1"` was disabling joint/OVR for *all* traits; cost 19 pp of Trophic.Level accuracy (0.789 → 0.600). Now discrete keeps its joint/OVR baseline |
| **#171** | UQ tail | Tree-MI conformal draws (P1-11); zi_count **conditional-on-non-zero** intervals (deliberately *not* an E[X] interval); one definitive `$se`-by-type doc section |
| **#172** | Refinement variance + guard | `.mvn_estep_refine` pooled two **dependent** posteriors with the independent-sources precision formula → coverage 0.925 → 0.618. Fixed; plus a divergence guard making `max_iter > 0` safe after its 2026-05-17 disablement |

---

# PART 4 — THE EVIDENCE BASE (open these, don't re-derive)

**Campaign results** (all pre-registered *before* results existed):
- `2026-08-16-regime-map-results.md` — 5,400 jobs. The GNN recovers **nonlinear** cross-trait
  structure the baseline structurally cannot (F2, 58% of cells at λ=1); the linear control F1
  is null. **The large low-λ gains are baseline misspecification, not GNN skill.**
- `2026-08-16-lambda-attribution-results.md` — 1,920 jobs. Confirms it: fitting λ in the
  **baseline** reproduces **100–117%** of those low-λ gains with no neural network. The F2
  lift **survives** a λ-fitted baseline (5.9–8.5 MCSE) — that is the claim that stands.
  Also: `lambda_mode="bayes"` eliminates the ML boundary collapse (P(λ̂<0.05) 0.47 → 0.00);
  NEWS's planned "v0.11 CV-λ" remedy is measurably the wrong one.
- `2026-08-16-mechanism-coverage-results.md` — 200 jobs. **Exchangeability failure CONFIRMED**:
  conformal coverage drops under MAR_phylo (−3.4 pp) and MNAR, and does **not** wash out with n.
  MCAR control healthy, so the pipeline is sound; the distribution shift is real.
- `2026-08-16-coverage-results.md` — **contains a CORRECTION of my own overclaim.** Two
  distinct defects: the n=100 arithmetic ceiling (ordinary) and large-n real-data
  undercoverage (0.87–0.91 at n=4k–10.6k) which the ceiling *cannot* explain.
- `2026-08-16-external-comparison-results.md` — **first working BACE head-to-head in repo
  history** (every prior one said "BACE skipped" — an API-signature bug, not a missing
  install). Plus a 5-seed Rphylopars/phylolm bench: **raw Rphylopars beat pigauto on all four
  AVONET continuous traits.** Mixed-type remains pigauto's unambiguous win (+30 pp vs BACE).
- `2026-08-16-continuous-gap-diagnosis.md` — the four-arm decomposition that **exonerated the
  gate and the mixed-type path** and convicted the solver.
- `2026-08-17-sigma-recovery-results.md` — the identical-predictions finding (**two different
  Σ estimators → numerically identical predictions in 12/12 cells**, because Σ was discarded).
  Also **retracts my own ill-posed Frobenius-vs-phylopars metric** — do not cite that column.
- `2026-08-17-refinement-results.md` — refinement helps at λ=0.2 but **fails on AVONET**
  (λ 0.993–0.998). Honest negative.
- `2026-08-17-exact-conditional-{design,results}.md` — the branch you're inheriting.

**After-task reports**: `after-task/2026-08-15-p0-land.md`,
`.../2026-08-16-earn-the-claims-arc.md`, `.../2026-08-17-arc3-solver-and-tail.md` (§9 = ranked
open items).

**Lesson draft (for the brain, not yet written there)**:
`2026-08-16-lesson-compute-lane-preflight.md`.

---

# PART 5 — THE UNMERGED BRANCH (your main job)

`arc/exact-conditional` — new file `R/exact_conditional.R` + `predict_method` threaded through
`fit_baseline` / `fit_pigauto` / `impute` and every joint/OVR path.

**The method.** For a Gaussian, `precision(L_unk | L_obs) = P[unk,unk]` — **no inversion of the
observed block**. With `P = Σ⁻¹ ⊗ S⁻¹` and `S⁻¹` the H&N sparse precision (`nnz = 4n − 6`),
a sparse Cholesky is the whole computation. It is the Kronecker lift of
`henderson_bm_predict()`, which has done exactly this for K=1 in production.

**Correctness, three independent ways:**
1. **7e-09** relative error vs a dense conditional built from `Σ ⊗ R`, mean *and* variance, on
   **both** the `A` and `R` scales.
2. **Given the true Σ it beats phylopars** (λ=0.2, n=300, 20 reps): 0.639 vs 0.645 (exch07);
   0.870 vs 0.879 (hetero). The conditional is optimal; residuals are Σ *estimation*.
3. Recovery sim at λ=1.0: 0.419 vs phylopars 0.413 (shipped path 0.498); coverage 0.963/0.934.

**Speed** (K=4, 30% MCAR): 9× dense at n=300, **86× at n=600**, dense infeasible by n=1200
while this stays at 0.10 s. Without H&N this feature could not exist in-house.

**G4 real-data gate — passes 3 of 4** (AVONET300, 5 masks, **paired**):

| trait | per_column | exact | paired Δ | \|Δ\|/SE | gap closed |
|---|---:|---:|---:|---:|---:|
| Beak.Length_Culmen | 0.892 | 0.717 | −0.175 | 3.67 | 60% |
| Tarsus.Length | 1.170 | 1.023 | −0.148 | 2.95 | 50% |
| Wing.Length | 0.629 | 0.523 | −0.106 | 2.64 | 59% |
| Mass | 1.468 | 1.465 | −0.004 | 0.10 | 2% |

Reference on the same masks: `joint_solver="rphylopars"` = 0.602 / 0.873 / 0.449 / 1.295.
**It does not reach parity.** Mass is unmoved and unexplained — flagged, not glossed.

**Your steps** (suite is already green — do NOT re-run unless you change code):
```sh
cd "/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/arc-exact"
export NOT_CRAN=true
Rscript -e 'torch::torch_is_installed()'   # libtorch must be present
Rscript -e 'r <- rcmdcheck::rcmdcheck(args="--as-cran", error_on="never", quiet=TRUE);
            cat("E:",length(r$errors),"W:",length(r$warnings),"N:",length(r$notes),"\n")'
# expect 0 / 0 / 1  (known NOTE: version string + BACE not in mainstream repos)
```
Then open the PR against `main`, framed as **opt-in, not parity, default unchanged**.
**Do not auto-merge. Do not propose flipping the default.**

---

# PART 6 — FAILURE PATTERNS THAT WILL BITE YOU (read before coding)

**1. Silent parameter drops — FIVE instances in `fit_baseline.R`'s call chain.**
`lambda_mode` (×2, PRs #126/#127), baseline config, covariates-into-joint (#167), and
`predict_method` in **two** places this week. The threading pattern in that file is
structurally unsafe: a parameter added to a signature does not reliably reach its consumer.
**Every new argument needs an end-to-end test that exercises every dispatch branch it can
reach** — a continuous-only fixture takes the *continuous-joint* path and will pass while the
*threshold-joint* path is broken.

**2. Broad `tryCatch` turns crashes into silent quality loss — twice.**
`fit_ovr_categorical_fits()` is called inside `tryCatch(error = function(e) NULL)`. When it
received an argument it didn't accept, the error was swallowed, `probs` became `NULL`, and
**every categorical trait silently fell back to label propagation** — no error, no warning.
Same shape as the "BACE skipped" string that hid an API-signature bug for months.
**Recommendation: those handlers should log what they caught.**

**3. Verify the treatment applied before reporting a negative.** The first G4 run reported
all-null effects and would have been reported as "the exact conditional fails on real data."
It was void — the option was a silent no-op. Two hours of compute, a believable false
negative, caught only by checking the treatment actually ran. Log kept:
`dev/bench_exact_avonet_VOID.log`.

**4. Bench scripts resume from a committed `.rds`.** In a fresh worktree they silently reuse
**old** cells and regenerate a summary that looks like a genuine re-run
(`Resuming: 35/35 cells already done`). Move the `.rds` aside to force a real run.

**5. Comparisons to the April per-type references are confounded.** The *baseline algorithm*
changed since (Level-C joint MVN landed after), and zi_count also carries the P1-9 split fix.
Only `multi_proportion`'s baseline reproduces April exactly.

**6. Operational** (cost real time this week): `pkill -f <pat>` over SSH **kills your own
session** when the pattern appears in your own remote command line — use `pkill -f 'foo[2]'`.
GPU/MPS contention between concurrent torch jobs produces `Internal Error (00000206)` storms —
serialise torch work or force CPU. Concurrent runs in the *same worktree* overwrite each
other's outputs.

---

# PART 7 — FENCES (do NOT do these without a fresh decision)

- **No default change** to `predict_method` or `joint_solver` on this evidence (one dataset).
- **No new EM/refinement loop.** Tried in #172: helped a low-λ sim, failed AVONET.
- **No GNN architecture work / transformer pretraining.** The GNN is not the bottleneck.
- **Do not replace the BM kernel with phylolm** — it is monotype continuous and would
  dismantle unified mixed-type imputation, pigauto's actual selling point.
- **No Σ-estimation work at low λ** without a decision — separate project, unclear payoff.
- **No CRAN submission** while `Suggests: BACE` (BACE still 404 on CRAN).
- **No BACE-tree edits** (either tree). **Never `git add -A`.** 15 dirty files stay unstaged.
- **Do not cite** the regime map's low-λ GNN gains as accuracy, the retracted
  Frobenius-vs-phylopars column, or any pre-P0 benchmark number.

---

# PART 8 — OPEN ITEMS, RANKED

1. **`--as-cran` → PR** for `arc/exact-conditional` (mechanical; suite already green).
2. **Shinichi's decisions** (not yours): paper framing; and A/B/C from his 2026-08-17 phone
   Doc — **C** document the gap and move on, **B** recommend rphylopars for continuous,
   **A** keep building in-house. **A has now partially succeeded.** Claude's recommendation
   on record: **A-as-shipped (opt-in) + C for the paper**; do not chase full parity.
3. **Multi-seed external comparison** (BACE + Rphylopars + phylolm + missForest) — run
   **after** the PR lands so it benches the shipped solver. Totoro-sized.
4. **Real-data Mondrian coverage** re-run (fishbase / PanTHERIA) — the second coverage defect.
5. **Push `arc/pertype-benches`** (`4586b62`, unpushed) and land its completion doc.
6. **OVR categorical head-to-head** vs label-propagation and BACE — the 72% Trophic.Level
   figure in the docs is **BACE's own number**, not a pigauto measurement.
7. **A recovery sim for `suggest_next_observation()`** — the feature no other phylogenetic
   imputer exposes, and it has no ADEMP campaign. If the paper leads on mixed-type + design,
   this is the missing figure.
8. **Liability-decode calibration** (binary/ordinal): are `Φ(μ/√(1+σ²))` probabilities
   calibrated once per-tip SEs feed them? Unmeasured.
9. **Packaging**: capability-surface rows for the joint baseline and conformal are **stale**
   (they still imply an optional Rphylopars dependency for the joint path — untrue since
   PR #100); every manuscript number needs a committed regenerating script; version strings
   drift across DESCRIPTION / AGENTS.md / README.

---

# PART 9 — DEPENDENCY STATUS (asked directly, answered precisely)

**pigauto has NO Rphylopars dependency.** It is in `Suggests:`, there is exactly **one**
`Rphylopars::` call in the package (`R/joint_mvn_solver.R`, in `.fit_mvn_bm_rphylopars()`),
and it is unreachable unless a user explicitly passes `joint_solver = "rphylopars"`.
`joint_mvn_available()` returns `TRUE` unconditionally — the joint path has needed no external
solver since PR #100. If the exact conditional ever reaches parity, that branch and the
`Suggests` entry can both be deleted.

---

# PART 10 — COMPUTE STATE

**Totoro** (`snakagaw@totoro.biology.ualberta.ca`, via the standing `~/.ssh/cm-*` socket —
no Duo needed). `~/pigauto_regime_map/` holds: `results/regime_map` (5,400),
`results/lambda_attr` (1,920), `results/mech_cov` (200), `results/mech_cov_mondrian`
(200 valid + `_bad300` superseded), `lib_mondrian/`, plus runners `regime_cell.R`,
`lambda_cell.R`, `mech_cell.R`, `run_campaign*.sh` and their summarisers. **Nothing running.**

**D-143 is an ACCOUNT-level cap (≤150 cores), not per-lane.** On 2026-08-16 three lanes — two
of them other projects — each individually compliant, summed to **192 cores**. Measure what
the account is already using *before* launching:
```sh
ssh totoro 'ps -eo user,pcpu --no-headers | awk "{a[\$1]+=\$2} END {for(u in a) if(a[u]>50) print u, int(a[u]/100)\" cores\"}"'
```

---

# PART 11 — RESUME PROMPT

```text
Read AGENTS.md and docs/dev-log/handover/2026-08-18-codex-handover-EXTENDED.md. Run the handover rehydration steps, reconcile them with the current git state, then continue only the OWED Next Immediate Steps.
```
