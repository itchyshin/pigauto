# Session Handoff: pigauto — exact Σ⊗R conditional, unmerged and gated (Claude → Codex)

Meta: 2026-08-18 · from Claude (`AUTHOR = claude`) · `TARGET = codex`
**This file is the authoritative input. You inherit no chat.** Everything needed is here.

---

## 0. READ FIRST — the one thing that is unfinished

Branch **`arc/exact-conditional`** (worktree `.worktrees/arc-exact`, 6 commits ahead of
`origin/main`, **pushed**, **no PR opened**) implements the exact matrix-normal conditional.
It is **evidence-complete and behaviour-verified**, but its **final full-suite run was still
executing when this handover was written**. Re-run it before opening a PR.

```sh
cd "/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/arc-exact"
NOT_CRAN=true nice -n 15 Rscript -e 'devtools::test(stop_on_failure = FALSE)'
```

Expected: **FAIL 0**. The last observed count was FAIL 4, but **all four were traced to a
regression that has since been fixed** (see §3). The three directly-affected files
(`test-lambda-per-type.R`, `test-exact-conditional.R`, `test-ovr-categorical.R`) each pass
individually post-fix, and default output is bit-identical to `origin/main`. If FAIL > 0
persists, **stop and diagnose — do not open the PR.**

---

## 1. Repo state

- `origin/main` = **`f925e17`**, version `0.10.0.9000`, **zero open PRs**. PRs #167–#172
  merged 16–17 Aug (honesty warnings; `conformal_method="mondrian"`; `joint_solver` switch;
  per-type λ dispatch; UQ tail; refinement variance fix).
- Docs branch **`handover/2026-08-09-cursor`** carries every results doc. Push with
  `git add -f` on explicit paths (`docs/` is gitignored).
- **15 carried-over dirty files must stay unstaged** (`.gitignore`, `AGENTS.md`, `CLAUDE.md`,
  `README.md`, `_pkgdown.yml`, `dev/gnn_attribution_low_lambda_*`,
  `script/*gnn_attribution*`, `script/returned_gnn_attr/`). **Never `git add -A`.**
- Worktrees: `arc-exact` (this work), `arc-bace` (comparator benches, pushed),
  `arc-pertype` (per-type benches, branch `arc/pertype-benches`, commit `4586b62`, unpushed).
- **No BACE-tree edits. No CRAN submission** (`Suggests: BACE`, BACE still 404 on CRAN).

## 2. What `arc/exact-conditional` does

`predict_method = c("per_column", "exact")` — **opt-in, default unchanged.** Computes the
exact conditional under `vec(L) ~ MVN(0, Σ ⊗ R)` instead of predicting traits independently.

Conditioning is done in the **precision form**: `precision(L_unk|L_obs) = P[unk,unk]` with
`P = Σ⁻¹ ⊗ S⁻¹`, so the `|O|×|O|` observed block is never inverted. `S⁻¹` is the
**Hadfield & Nakagawa (2010) eq. 29** sparse precision over the extended tree
(`build_henderson_S_inv()`, `nnz = 4n − 6`). New file `R/exact_conditional.R`; it is the
Kronecker lift of `henderson_bm_predict()`, which has done this for K=1 in production.

**Motivation** (`docs/dev-log/2026-08-16-continuous-gap-diagnosis.md`): the in-house joint
solver estimated Σ and then discarded it (`max_iter = 0` short-circuits to per-column BM),
which is why raw Rphylopars beat pigauto on AVONET continuous traits.

## 3. Evidence — and the two false results that nearly shipped

**Correctness, three independent ways:**
1. Matches a dense conditional built from `Σ ⊗ R` to **7e-09** relative error, mean and
   variance, on **both** the `A` and `R` scales (testing one scale hides errors: means are
   scale-invariant, variances differ by tree height).
2. **Given the true Σ it beats phylopars** (λ=0.2, n=300, 20 reps): 0.639 vs 0.645 (exch07),
   0.870 vs 0.879 (hetero) — so the conditional is optimal; residuals are Σ *estimation*.
3. Recovery sim λ=1.0: 0.419 vs phylopars 0.413 (shipped path 0.498). Coverage 0.963/0.934.

**G4, the real-data gate — PASSES on 3 of 4 traits** (AVONET300, 5 masks, **paired**):

| trait | per_column | exact | paired Δ | \|Δ\|/SE | gap closed |
|---|---:|---:|---:|---:|---:|
| Beak.Length_Culmen | 0.892 | 0.717 | −0.175 | 3.67 | 60% |
| Tarsus.Length | 1.170 | 1.023 | −0.148 | 2.95 | 50% |
| Wing.Length | 0.629 | 0.523 | −0.106 | 2.64 | 59% |
| Mass | 1.468 | 1.465 | −0.004 | 0.10 | 2% |

**Does NOT reach parity** with `joint_solver="rphylopars"` (0.602/0.873/0.449/1.295). Mass is
unmoved and unexplained — flagged, not glossed.

**Two false results caught, both recorded in the results docs:**
- The **first G4 run was VOID**: `predict_method` reached `fit_joint_threshold_baseline()` but
  was never forwarded to `fit_joint_solver()`, so on any dataset with binary/ordinal traits it
  was a silent no-op. Two hours of compute produced a believable *false negative*. Caught only
  by verifying the treatment applied (tracing showed `exact_conditional_mvn` never called).
  Void log preserved: `dev/bench_exact_avonet_VOID.log`.
- The **4 suite failures were a real regression**: `fit_baseline` passed `predict_method` to
  `fit_ovr_categorical_fits()`, which did not accept it — and that call sits inside
  `tryCatch(error = function(e) NULL)`. The error was **swallowed**, `probs` became `NULL`, and
  **every categorical trait silently fell back to label propagation**. Fixed in `ba255ee`;
  default output re-verified bit-identical to main (max|diff| 1.509 → 0).

**This is the fifth silent-parameter failure in this call chain** (`lambda_mode` ×2, baseline
config, covariates-into-joint, `predict_method` ×2 places) and the **second** where a broad
`tryCatch` converted a crash into a silent quality loss. The Mission Control board carries a
standing high-priority recruit item for exactly this.

## 4. Your job (Codex runs the LIVE toolchain)

Live env for this repo:
```sh
cd "/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/arc-exact"
export NOT_CRAN=true            # torch/testthat gates
# torch + libtorch must be present: skip_if_not_installed("torch") does NOT check the backend
Rscript -e 'torch::torch_is_installed()'
```

1. **Full suite** (above). Expect FAIL 0. If not, diagnose before anything else.
2. **`--as-cran`** — expect the baseline **0 errors / 0 warnings / 1 known NOTE**
   (version string + `BACE` not in mainstream repos):
   ```sh
   Rscript -e 'r <- rcmdcheck::rcmdcheck(args = "--as-cran", error_on = "never", quiet = TRUE);
               cat("E:", length(r$errors), "W:", length(r$warnings), "N:", length(r$notes), "\n")'
   ```
3. **Open the PR** (base `main`) only if 1–2 pass. Frame it as an **opt-in improvement that
   does NOT reach parity and does NOT change the default**. Do **not** auto-merge; do not
   propose flipping the default — that needs multi-dataset evidence and is Shinichi's call.
4. **`devtools::document()`** if you touch roxygen; `man/` is already regenerated.

## 5. Rehydration recipe (Codex-native)

- `AGENTS.md` at repo root is native to you — read it first; it holds the architecture,
  the eight trait types, the UQ design (three distinct mechanisms — do not conflate), and
  the gotchas (`.Rbuildignore` is wide; `k_eigen="auto"` scales with tree size).
- `.codex/agents/*.toml` — the team. **Rose is mandatory** for any claim-gate/after-task work.
- Then, in order: this file → `docs/dev-log/2026-08-17-exact-conditional-results.md` →
  `...-exact-conditional-design.md` (pre-registered gates) →
  `...-2026-08-16-continuous-gap-diagnosis.md` (why any of this exists).
- Repo state is truth; this doc is a pointer.

## 6. Fences — do NOT do these

- **No default change** to `predict_method` / `joint_solver` on this evidence.
- **No new EM/refinement loop.** It was tried (PR #172): helped a low-λ sim, failed AVONET.
- **No GNN architecture work.** The GNN is not the bottleneck; the measured claim that
  survives is nonlinear cross-trait structure at correctly-specified λ (F2, 5.9–8.5 MCSE).
- No BACE-tree edits, no CRAN submission, no `git add -A`, 15 dirty files stay unstaged.
- **Do not chase Σ estimation at low λ** without a fresh decision — it is a separate project.

## 7. Open decisions (Shinichi's, not yours)

Emailed himself a Doc 2026-08-17 with three options; **A has since partially succeeded**:
- **C** = document the covariance gap and move on (mixed-type is the paper).
- **B** = recommend rphylopars for continuous traits.
- **A** = keep building the in-house prediction step ← **this branch; 50–60% of the gap**.

Claude's recommendation on record: **A-as-shipped (opt-in) + C for the paper.** Do not chase
full parity. **Also open:** paper framing; and whether `joint_solver="rphylopars"` stays.

## 8. Still owed (unstarted)

- Real-data Mondrian coverage re-run (fishbase / PanTHERIA).
- Multi-seed external comparison (BACE + Rphylopars + phylolm + missForest) — run it
  **after** the PR lands so it benches the shipped solver.
- Per-type benches: `arc/pertype-benches` (`4586b62`) has proportion/zi_count/multi_proportion
  results, **unpushed**. Caveat from that run: comparisons to the April references are
  confounded — the *baseline algorithm* changed since (Level-C joint MVN), and zi_count also
  carries the P1-9 split fix. Only multi_proportion's baseline reproduces April exactly.
- **Bench caching gotcha:** the per-type scripts resume from a committed `.rds`. In a fresh
  worktree they silently reuse OLD cells and regenerate a summary that looks like a re-run.
  Move the `.rds` aside to force a genuine run.

## 9. Resume prompt

```text
Read AGENTS.md and docs/dev-log/handover/2026-08-18-codex-handover.md. Run the handover rehydration steps, reconcile them with the current git state, then continue only the OWED Next Immediate Steps.
```
