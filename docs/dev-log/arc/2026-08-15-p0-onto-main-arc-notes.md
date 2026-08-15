# Arc notes — P0 → main merge (scratch; becomes the after-task receipt)

Branch: `arc/p0-onto-main`, cut from `origin/main` `416561b`.
Merged: `origin/fix/ci-install-libtorch` `21d2ea6` (carries PR #155).
Merge commit `8df0b0e`; doc regen `f5e8416`.
**Not landed. `main` untouched. Rung 3 gated on Shinichi's G0.**

## Conflict resolutions (7, exactly as the pre-merge probe predicted)

| File | Resolution | Judgement required? |
|---|---|---|
| `.Rbuildignore` | union — both sides added different ignore lines | no |
| `tests/testthat/test-gbif-centroids.R` | took P0 side; difference was **two blank lines only** (verified by printing both sides) | no |
| `tests/testthat/test-joint-baseline.R` | took P0 side — strictly more informative test name + the `3c31f1d` contract comments | no |
| `.github/workflows/R-CMD-check.yaml` | took main side — main added a macOS MPS regression canary step; P0 side was empty there. P0's libtorch-install change auto-merged elsewhere in the file. | no |
| `AGENTS.md` | kept P0's LOAD-FIRST manifest **and** main's cleaned prose; **dropped** the `~/.Codex/projects/…` host note and the "Codex.ai/code" framing | **yes — see F3** |
| `NEWS.md` | union under `# pigauto 0.10.0.9000 (development)` | **yes — see F2** |
| `R/fit_helpers.R` | took main's `calibrate_gates()` signature (`seed = NULL`) | **yes — see F1** |

## Findings

### F1 — `seed` default: main's value is the CRAN-requested one

Merge base had `calibrate_gates(seed = 1L)`. **main changed it to `NULL`** and added
`is.null(seed)` guards through the body. P0 branched earlier, kept `1L`, and only
reformatted the signature — so P0 never touched the body, and **main's NULL-handling
code auto-merged silently**. Taking P0's signature would have left code written to
service a `NULL` default paired with a default that is never `NULL`.

`cran-comments.md` records this as a CRAN reviewer correction for the 0.10.0
resubmission: *"changes public random-seed defaults to `NULL`, retaining
reproducibility only when the user explicitly supplies a seed."*

Scope, stated precisely — the public halves were never at risk:

| Function | base | main | merged |
|---|---|---|---|
| `impute()` | `1L` | `NULL` | `NULL` (auto-merged) |
| `multi_impute()` | `1L` | `NULL` | `NULL` (auto-merged) |
| `calibrate_gates()` internal | `1L` | `NULL` | `NULL` (**conflicted; resolved**) |

Consequence of the wrong choice: an *inconsistent* package (public API honouring
`NULL`, internal layer reverted to `1L`), and force-seeding on direct
`pigauto:::calibrate_gates()` calls — i.e. `test-safety-floor.R`. **Not** a production
behaviour change: `fit_pigauto.R:956` passes `seed = seed` explicitly.

### F2 — NEWS.md placement (an error I made and caught)

First union appended P0's four Fix sections *after* main's block, which ends by opening
the `# pigauto 0.10.0` heading — filing four never-released fixes under the **released
CRAN tarball**. Corrected: they now sit under `# pigauto 0.10.0.9000 (development)`,
matching their placement on P0's own branch. Verified absent from `origin/main:NEWS.md`,
confirming they never shipped in 0.10.0.

### F3 — AGENTS.md: P0 predates PR #136

P0's side would have reintroduced a `~/.Codex/projects/-Users-…/memory/` path that
**never existed** (journal 2026-07-08 item #9: FIXED in PR #136) plus a "Codex.ai/code"
framing main had replaced. Kept P0's intentional addition (the LOAD-FIRST compute
manifest, commit `18ff66c`); dropped both defects.

### F4 — pre-existing, NOT introduced here

`R/fit_pigauto.R` roxygen still shows `seed = 1)` in an example block while the actual
default is `NULL`. Present on `origin/main` already. Flagged, not fixed — out of scope
for a merge resolution.

## Verification status

- [x] no conflict markers anywhere; `R/fit_helpers.R` parses
- [x] `devtools::document()` — regenerated `man/fit_pigauto.Rd` only (P0's `refine_steps`
      roxygen; the `.Rd` was stale relative to merged source)
- [x] `devtools::test()` — **1 failure**, see F5. Everything else green.
- [ ] `rcmdcheck --as-cran` — **NOT RUN.** Deliberately skipped: a 45-min check on a branch
      with a known unexplained failure spends the time and adds nothing.
- [x] claim-gate — 1 BLOCKING, 5 SHOULD-FIX, 3 NOTE → `CLAIM-GATE-FINDINGS.md`

## F5 — BLOCKER: P0 breaks the strict val-floor invariant on main's data

`test-safety-floor.R:712` — *"strict val-floor: pigauto val-loss <= baseline val-loss per
trait, all types"*, trait `x_cont` (continuous).

```
Expected loss_blend <= loss_bm * 1.05 + 1e-10
run 1:  0.3434 > 0.3410
run 2:  0.3442 > 0.3410
```

**Controlled comparison — same test, same machine, same data shape:**

| Branch | Result |
|---|---|
| pristine `main` `416561b` (detached worktree `.worktrees/p0-control`) | **PASSES** — no Failed section, no Skipped section |
| `arc/p0-onto-main` (`main` + P0) | **FAILS**, reproducibly, twice |

The control provably *ran* the test rather than skipping it: it emitted the identical
`Small validation set for 4 trait(s): x_cont (n=10), x_bin (n=10), x_cat (n=13),
x_cnt (n=15)` line as the failing run.

`loss_bm * 1.05 + 1e-10` is byte-identical (`0.3410`) on every run → the BM baseline is
deterministic. Only `loss_blend` moves (`0.3434`, `0.3442`) → the nondeterminism is in the
GNN/torch path, but the *overshoot* is systematic, not a coin-flip flake.

**Leading hypothesis (NOT yet confirmed — this is the open question):** P0 fix #4 stops the
DAE seeing held-out val/test truth as encoder *context* during training. Pre-fix the GNN
trained with val cells visible, so its validation loss was flattered. Post-fix it is
genuinely worse on validation — ~5.7% vs a 5% tolerance. On this reading the failure is the
fix working, and the threshold was calibrated against the leaky pipeline.

Independent corroboration: the claim-gate agent, working only from the diff with no
knowledge of this test, wrote that *"any benchmark number that depended on a full
`fit_pigauto()` GNN training run pre-fix was generated on a training loop with this
held-out-context leak."* This test is such a number.

**Why that is not sufficient to land on.** The competing explanation — P0 genuinely degrades
the blend past its own safety invariant — produces the same signature. The strict val-floor
is a *safety property* (`r_cal = 0` must always remain a valid fallback), so a real breach is
serious, not cosmetic. One test cannot separate the two readings.

**Ruled out:**
- *My `seed` resolution (F1).* The test seeds itself: `set.seed(20260429L)` and
  `seed = 20260429L` passed explicitly. The `seed = NULL` default never fires here.
- *A merge-resolution error.* `test-safety-floor.R` auto-merged; main's only changes to it
  are `skip_if_no_libtorch()` refactors, no logic.
- *Pure MPS flakiness.* Plausible a priori — main pins a *neighbouring* test to CPU citing
  exactly this ("MPS kernels can vary enough ... to move a single masked cell across the
  deliberately tight smoke threshold") — but ruled out as sufficient, because the control on
  `main` passes and the merged branch fails twice.

**Next step for whoever picks this up:** determine whether the degradation is confined to
fix #4 by fitting with `model_config$train_mask_heldout` forced FALSE on the merged branch.
If the test passes with the leak restored, the hypothesis is confirmed and the threshold
needs re-deriving on honest data (with the reason recorded), not loosening.

### F5-RESULT — experiment run 2026-08-15. Hypothesis WRONG; the real cause is worse.

Script: `docs/dev-log/arc/2026-08-15-train-mask-heldout-experiment.R`.
`train_mask_heldout` in `model_config` is a **record, not a switch** (`fit_pigauto.R:1099`);
the behaviour lives at the single call site `fit_pigauto.R:439`
`X_fill_heldout <- mask_heldout_with_baseline(X_fill, MU, hold_idx)`, which feeds training,
val eval, **and gate calibration**. The leak is restored by mocking that function to identity.
Same fixture, same seed, same device; only the leak differs.

**Losses (`x_cont`, 4 reps each; `bm = 0.3248` in every run; threshold `0.3410`):**

| Condition | blend across reps | Result |
|---|---|---|
| LEAKED (pre-fix) | 0.3248, 0.3248, 0.3299, 0.3248 | 4/4 PASS |
| FIXED (P0) | 0.3422, 0.3403, 0.3330, 0.3415 | **2/4 FAIL** |

**Calibrated gate (`x_cont`, 2 reps each) — the decisive measurement:**

| Condition | `r_cal_bm` | `r_cal_gnn` | `r_cal_mean` |
|---|---|---|---|
| LEAKED | **1.0** | **0.0** | 0.0 |
| FIXED | **0.9** | **0.1** | 0.0 |

**What this means.** Under LEAK the gate closes completely, so `blend == bm` exactly and the
floor passes trivially. Under P0 the gate **opens 10% to the GNN, and that contribution makes
the prediction worse than pure BM.** The original hypothesis — "the GNN is honestly worse now,
so the threshold is stale" — is **wrong**: on the pre-fix path the GNN was not used at all here.

This is a **safety-invariant failure, not a stale threshold.** `calibrate_gates()` is supposed
to open the gate only when it improves validation loss (relative-gain floor + half-A/half-B
cross-check). It opened anyway and degraded the result. `r_cal = 0` stopped being the
protective fallback the architecture documents.

**Inferred mechanism — NOT confirmed, do not cite as fact.** Calibration evaluates with
held-out cells baseline-filled (P0's masking), while `predict()` later sees those same cells as
*observed* context. The gate would then be chosen on one surface and applied to another — an
asymmetry, in the fix named "train/cal/predict symmetry". Confirming this requires comparing
the delta surface at calibration against the one at predict, which was not done.

**Regime (do not over-generalise).** n = 80, one synthetic 4-trait BM fixture, `x_cont` only,
30 epochs, single-obs, 4 reps, one machine. This does **not** show P0 breaks the gate in
general. It shows P0 breaks it reproducibly on the fixture pigauto's own safety test uses.

**Also learned:** the failure is *marginal in isolation* (2/4 reps) but was 2/2 inside the full
test file, consistent with accumulated torch/MPS state pushing a borderline case over — the
same effect `main` already documents when it pins a neighbouring test to CPU. The borderline-ness
is caused by P0 opening the gate; the run-to-run noise only decides which side of the line a
given run lands on.

**Verdict: do not land P0 as-is.** The next question is why calibration opens a harmful gate —
not how to make the test green.

### F5-MECHANISM — surface comparison run 2026-08-15. The asymmetry is MEASURED.

Script: `docs/dev-log/arc/2026-08-15-delta-surface-compare.R`. One fit, no refit, so model and
calibrated gate are held fixed and **only the input context varies**. `.mask_observed_idx`
un-pins chosen cells at predict, reproducing the calibration-time context.

Gate on this fit: `r_cal_bm = 0.85`, `r_cal_gnn = 0.15`.

| Surface | val MSE (`x_cont`, 10 cells) | vs floor `0.3410` |
|---|---|---|
| pure BM baseline | 0.3248 | — |
| **B** calibration-like (val truth HIDDEN) | **0.3315** | within floor |
| **A** predict surface (val truth PINNED) | **0.3429** | **BREACHES** |

`corr(A, B) = 0.9996`; `max|A − B| = 0.0723`; `mean(A − B) = −0.0179`. Nearly the same shape,
a few cells moving enough to cross the line.

**Gate calibration scores surface B; the package delivers surface A.** The code asymmetry:

- **calibration** `fit_pigauto.R:819` — `X_init = t_X_eval` (val cells start at BM),
  `pin_mask = t_M_obs`, and `M_obs_mat[c(val_idx, test_idx)] <- FALSE` (`:429`) → val truth
  never enters as DAE context.
- **predict** `predict_pigauto.R:339` — `observed_mask <- !is.na(object$X_scaled)`. Val cells
  are genuinely non-NA (held out only for *fitting*), so predict **pins them to their own
  truth** for every refine step.

**Counter-intuitive and probably the real story:** pinning a cell's own truth makes the
prediction *worse* (0.3429 vs 0.3315). **Inferred, NOT measured:** P0 trains the DAE with
held-out cells at baseline, so real data in those positions at predict time is
**out-of-distribution**; pre-fix `main` trained *with* that truth visible, so pinning it is
in-distribution and free.

If so, P0 fixed train↔cal symmetry and thereby **exposed a pre-existing cal↔predict
asymmetry** — invisible on `main` only because the gate closes to 0, and at `r_cal_gnn = 0`
the input context cannot change the output at all. **Cheap test of that claim (NOT run):** on
pristine `main`, force the gate open and check whether A ≠ B there too. If it does, the
asymmetry predates P0 and P0 is the messenger.

**Which surface is correct? B.** A genuinely missing cell has no truth to pin, so B is what
real users get. That indicts the test: it evaluates val cells through a plain `predict(fit)`,
leaking their truth into context. Passing `.mask_observed_idx = splits$val_idx` would measure
the surface that corresponds to use.

**Do not treat that as the fix without a decision.** Two restraints: (1) even on the correct
surface B, the blend (0.3315) is still **worse** than pure BM (0.3248) — the gate opened for a
GNN that does not help, merely not badly enough to breach; (2) whether the val-floor invariant
should be asserted on the pinned or the unpinned surface is a design question about what the
safety property *means*. That is Shinichi's call, not a merge-resolution decision.

## F6 — pre-existing oddity, not introduced here

`tests/testthat/test-safety-floor.R` ends with a bare `skip_if_no_libtorch()` at top level,
outside any `test_that()` block (last line). Present on `origin/main`; looks like a stray
from main's `skip_if_*` refactor. Harmless but almost certainly unintended. Not fixed here.

## Residual risk the tests must clear

Three files auto-merged **clean** despite heavy independent movement on both sides.
Textual cleanliness is not semantic correctness. Prime suspect first:

| File | main churn | P0 churn |
|---|---:|---:|
| `R/multi_impute.R` | +63 / −95 | 38 lines |
| `R/fit_pigauto.R` | +25 / −7 | 134 lines |
| `R/impute.R` | +19 / −49 | 6 lines |

**Risk branch (declared before running):** if failures are not confined to the seven
conflicted files, STOP and report the diagnosis rather than pushing through to Rung 1.

### F5-CONTROL — gate-open check on pristine `main`, 2026-08-15. **Asymmetry PREDATES P0.**

Script: `docs/dev-log/arc/2026-08-15-gate-open-check-on-main.R`, run on detached `416561b`.

```
AS-CALIBRATED gate on main: r_cal_bm=0.90  r_cal_gnn=0.10
  gate as-is : A=0.3301  B=0.3266  max|A-B|=0.0289   (BM=0.3248, floor=0.3410)
FORCED gate 0.85/0.15 on the SAME model:
  A predict (truth PINNED) 0.3349 within floor
  B calib   (truth HIDDEN) 0.3297 within floor
  max|A-B|=0.0434  corr=0.9999  mean(A-B)=-0.0059
```

**A != B on `main`. The cal<->predict asymmetry is PRE-EXISTING, not introduced by P0.**

#### Two corrections to earlier entries in this file

1. **F5-RESULT's "on main the gate closes (r_cal_gnn = 0)" is WRONG.** That value came from the
   *mocked* LEAKED condition, which I wrongly treated as standing in for `main`. Pristine `main`
   calibrates to **`r_cal_gnn = 0.10`** — the GNN does contribute there.
   Why the mock misleads: mocking `mask_heldout_with_baseline()` reverts only ONE of fix #4's
   three parts (the masking), leaving the phylo-signal-gate fallback and the calibration-time
   refine loop active. LEAKED was a P0/main chimera, not pre-fix behaviour. Any inference drawn
   from LEAKED alone is unsafe; this control supersedes it.
2. **P0 AMPLIFIES the asymmetry, it does not create it.** Like-for-like at the same forced gate:

| | mean(A-B) | max abs(A-B) | A at gate 0.15 | vs floor |
|---|---|---|---|---|
| `main` | -0.0059 | 0.0434 | 0.3349 | passes |
| `main` + P0 | -0.0179 | 0.0723 | 0.3429 | **breaches** |

~3x larger under P0. `main` has enough headroom to stay inside the 5% floor; **P0 consumes it.**

Consistent with the OOD reading (now the third independent measurement pointing the same way):
P0's model never sees real data in held-out positions during training, so pinning truth there at
predict hurts it more than it hurts `main`'s model, which trained with that truth visible.

**Reframing for the decision.** The question is not "is P0 broken". `main` already ships this
asymmetry and already runs an open gate; it merely has margin. P0 spends the margin. The
asymmetry itself is the more interesting defect and is **independent of P0** — worth its own
issue whether or not P0 ever lands.
