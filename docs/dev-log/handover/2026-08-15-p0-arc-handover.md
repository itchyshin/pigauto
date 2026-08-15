# Session Handoff: pigauto P0 merged and verified on a branch, NOT landed — one blocker

Meta: 2026-08-15 · from Claude (`AUTHOR = claude`) · `TARGET = any`
Authoritative copy: this file. Chat is disposable. You inherit **no** authoring chat.

Shinichi is **stepping away from pigauto for a while.** Nothing here is urgent. Do not start
this lane unless he asks.

## Critical Context

1. **`origin/main` is unchanged at `416561b`** (Version `0.10.0.9000`, Suggests `BACE`).
   The 2026-08-14 wrap merge (#156) is still the head. Nothing this session touched it.
2. **P0 is merged and verified on `origin/arc/p0-onto-main` `f5e8416`** — pushed, **no PR**.
   Seven conflicts resolved, `document()` clean, full suite run.
3. **ONE BLOCKER stops the land.** `test-safety-floor.R:712` (strict val-floor invariant)
   **passes on `main` and fails on `main` + P0**, reproducibly, under a clean control.
   This is the whole reason P0 is not on `main`. See "The blocker" below.
4. **Shinichi locked the G0** ("land it, and fix the false doc claim too") **conditional on a
   clean check.** The check is not clean, so the land did not happen. The lock does **not**
   carry forward to a future session — re-confirm before landing.
5. **CRAN pigauto 0.10.0 live; CRAN BACE still 404.** Do not CRAN-submit from this `main`
   while Suggests includes `BACE`.
6. **PROTECTED / do not touch:** standalone BACE `@ce8bc87`; in-tree `pigauto/BACE` `@de87d8c`;
   the 15 dirty uinit/GNN items on the `handover/2026-08-09-cursor` checkout; the 44 foreign
   unpushed commits; EM (`max_iter > 0`); DRM.jl. Never `git add -A`. Never archive this project.

## The blocker (read this before doing anything)

`test-safety-floor.R:712` — *"strict val-floor: pigauto val-loss <= baseline val-loss per trait"*,
trait `x_cont`:

```
Expected loss_blend <= loss_bm * 1.05 + 1e-10
run 1:  0.3434 > 0.3410
run 2:  0.3442 > 0.3410
```

| Branch | Result |
|---|---|
| pristine `main` `416561b` | **PASSES** (test provably ran, not skipped — same `x_cont (n=10), x_bin (n=10), x_cat (n=13), x_cnt (n=15)` line) |
| `arc/p0-onto-main` | **FAILS** twice |

`loss_bm * 1.05 + 1e-10` is byte-identical `0.3410` every run → BM baseline deterministic.
Only the GNN blend moves, and the overshoot is systematic (~5.7% vs a 5% tolerance).

**Leading hypothesis — the fix working, not a bug.** P0 fix #4 stops the DAE seeing held-out
val/test truth as encoder *context* during training. Pre-fix the GNN trained with val cells
visible, so its validation loss was flattered. Post-fix it is honestly worse on validation, and
this test's 5% threshold was calibrated against the leaky pipeline. Independently corroborated
by the claim-gate agent, which (working only from the diff, with no knowledge of this test)
noted that any benchmark from a pre-fix `fit_pigauto()` training run carries that leak.

**Why that is not enough to land on.** The competing reading — P0 genuinely degrades the blend
past its own safety invariant — has the identical signature. The strict val-floor is a *safety
property* (`r_cal = 0` must always remain a valid fallback), so a real breach matters. One test
cannot separate the two.

**Ruled out:** the `seed` resolution (the test seeds itself explicitly: `set.seed(20260429L)`,
`seed = 20260429L`); a merge-resolution error (`test-safety-floor.R` auto-merged; main's only
changes to it are `skip_if_no_libtorch()` refactors); pure MPS flakiness (the control on `main`
passes).

**The next experiment, and it is cheap:** on `arc/p0-onto-main`, fit with
`model_config$train_mask_heldout` forced `FALSE` (restoring the leak) and re-run that one test.
- **Passes with the leak restored** → hypothesis confirmed. The threshold was calibrated on
  leaked evidence and must be **re-derived on honest data with the reason recorded in the test**
  — not quietly loosened until it goes green.
- **Still fails** → P0 degrades the blend for some other reason. Do not land; investigate.

## What was accomplished

- 7/7 conflicts resolved (merge `8df0b0e`), doc regen (`f5e8416`), branch pushed.
- Full test suite: 1 failure (above), everything else green — `covariate-alignment` 13/13,
  `henderson-s-inv` 18/18, `bm-internal` 58/58.
- Controlled comparison against pristine `main` proving the failure is P0-attributable.
- Claim-gate audit of `main`'s public surface: 1 BLOCKING, 5 SHOULD-FIX, 3 NOTE.
- Earlier the same session: rehydration receipt for the 2026-08-14 Cursor handover
  (`cb285ad`), all claims reconciled, zero discrepancies.

## The BLOCKING claim-gate finding (independent of the blocker above)

`vignettes/gnn-architecture.Rmd` §5 states `pred$se` is **"Exact under BM, model-dependent."**
On `main`, the default multi-trait path computes:

```r
se[miss_idx_local] <- sqrt(max(sigma2, eps))   # sigma2 <- var(y[obs_idx_local])
```

— one scalar broadcast to every missing tip, no phylogenetic position. Verified in code, not
taken from the agent. P0 replaces it with per-tip `σ√v_i`. Because this SE feeds
`decode_binary_liability()`'s `Φ(μ/√(1+σ²))`, pre-fix it degrades **predicted binary/ordinal
probabilities**, not just interval widths.

**`main` ships a documented guarantee its code does not satisfy.** That is an argument for
landing P0 sooner, not for landing it carelessly. Shinichi authorised fixing the vignette row in
the same PR; do that when the PR happens.

## Landing state

| Artifact | Committed | Pushed | PR | State |
|---|---|---|---|---|
| `origin/main` `416561b` | y | y | — | **UNCHANGED.** Do not land without re-confirming the G0. |
| `origin/arc/p0-onto-main` `f5e8416` | y | y | **none** | **CARRIED-OVER.** Merge + doc regen, clean tree, PR-ready once the blocker resolves. |
| `handover/2026-08-09-cursor` docs | y | y | none | Rehydration receipt `cb285ad` + this handover. |
| 15 dirty uinit/GNN items | n | n | — | **CARRIED-OVER**, untouched. Do not stage. |
| 44 foreign unpushed commits | mixed | n | — | **CARRIED-OVER.** Not this lane. |
| `.worktrees/p0-onto-main`, `.worktrees/p0-control` | — | — | — | Local worktrees; `p0-control` removed at close. |

## Resume

```bash
cd "/Users/z3437171/Dropbox/Github Local/pigauto"
git fetch origin
git worktree add .worktrees/p0-onto-main arc/p0-onto-main   # if not present
cd .worktrees/p0-onto-main
# Re-run the blocker:
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-safety-floor.R")'
# Then the train_mask_heldout experiment described above.
```

Still owed if and when it lands: `rcmdcheck --as-cran` (never run — expect the pre-existing
2 WARN / 3 NOTE and no new ones), the vignette §5 doc fix, and the 5 SHOULD-FIX claim-gate items.

## Gotchas

- `docs/` is `.gitignore`d — force-add specific paths only.
- `git merge-tree --write-tree <a> <b>` predicts the conflict set read-only, in seconds. It was
  exactly right here (7/7) and is what made the arc estimate honest. Use it before sizing a merge.
- Do not loosen the val-floor threshold to make the suite green. That is the one move that would
  convert an honest finding into a hidden regression.
- Do not run this on DRAC/Totoro. It is a ~45-min local check, and the toolchain is provisioned here.

> Related: `docs/dev-log/after-task/2026-08-15-p0-arc.md` ·
> `docs/dev-log/arc/2026-08-15-p0-onto-main-arc-notes.md` ·
> `docs/dev-log/arc/2026-08-15-p0-claim-gate-findings.md` ·
> `docs/dev-log/handover/2026-08-14-claude-handover.md` · `docs/dev-log/coordination-board.md`
