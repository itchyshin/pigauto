# Tier-1 bench re-run — the leak had no detectable effect on the per-type numbers

2026-08-15/16. Re-ran the per-type bench suite on the corrected pipeline
(`main` @ `c655d75`: leak-free GNN training, two-layer #157 gate floor) and compared
against the pre-fix outputs (May 30 / Jun 11 / Apr 28) that the publication roadmap had
flagged as leak-tainted and unusable.

**Headline: they are usable.** The 2026-08 GNN train/cal-symmetry fix — which removed
held-out truth from the training-time input context — changed the aggregate benchmark
metrics by less than one Monte Carlo standard error in every trait type measured.

## Paired per-rep comparison

Old and new runs share seeds (confirmed: the `baseline` and `mean` rows are
byte-identical across old/new, since those paths are deterministic and untouched). The
comparison is therefore **paired per (scenario, rep, trait)**, and the MCSE is
`sd(Δ)/√n_pairs` on the paired difference — the right error bar for this question.

| Trait type | Metric | n pairs | mean Δ (new − old) | MCSE | \|Δ\| / MCSE |
|---|---|---:|---:|---:|---:|
| continuous | RMSE | 2720 | −0.00168 | 0.00243 | **0.69** |
| binary | accuracy | 640 | +0.00382 | 0.00363 | **1.05** |
| ordinal | RMSE | 360 | +0.00048 | 0.00631 | **0.08** |
| count | RMSE | 360 | +0.00206 | 0.00527 | **0.39** |
| categorical | accuracy | 280 | +0.00253 | 0.00563 | **0.45** |
| proportion | RMSE | 640 | +0.00629 | 0.00553 | **1.14** |

Every type sits below 1.15 MCSE. Nothing here is distinguishable from zero.

Per-scenario means tell the same story: continuous deltas within ±0.004 RMSE; binary
accuracy mostly *improved* (largest +0.017 at signal_1.0); count within +0.006 RMSE.

**Categorical is the sharpest evidence:** 6 of its 7 scenarios reproduce to
**exactly 0.0000**. The calibrated gate closes completely for categorical traits, so the
GNN contributes nothing there — and a leak in GNN *training* cannot move a number the GNN
does not influence. The architecture's own safety mechanism made those cells immune.

## What this does and does not license

**Does:** the existing per-type benchmark numbers can be cited, because they have now been
reproduced on the corrected pipeline within noise. The roadmap's blanket rule ("no GNN
number measured before the P0 land may appear in the manuscript") can be relaxed to: cite
the **re-run** values, and note that the pre-fix values agree within MCSE.

**Does not:**

- This is *not* "the leak was harmless." It is "the leak did not move these aggregate
  metrics at this resolution." The `test-safety-floor.R:712` invariant — a tighter,
  per-cell measurement on a single fixture — *did* detect a difference. Aggregate means
  over hundreds of cells wash out effects that a strict per-trait floor catches.
- **Power:** with MCSE ≈ 0.002–0.006, effects smaller than roughly 0.01–0.02 in the metric
  are invisible here. A real-but-small leak effect cannot be excluded, only bounded.
- Regime: the bench suite's own DGPs, n, missingness and rep counts (5 reps per cell,
  3 for multi_proportion). Not a statement about other regimes — that is what the
  regime-map campaign is for.

## Status

Complete: continuous, binary, ordinal, count, categorical, **proportion**.
Pending: zi_count, multi_proportion (running).

**Read `zi_count`'s eventual result with care.** P1-9 landed tonight (observed zeros could
never enter val/test), so its *current* `main` behaviour differs from both the pre-fix and
the re-run numbers. The re-run is executed on `c655d75`, which predates P1-9, so the
comparison below stays a clean leak-only contrast — but neither side reflects today's
`main`. `script/bench_zi_count.*` needs a third run once P1-9 is in the bench worktree.

## Operational note (cost this lane paid)

The bench drivers use PSOCK clusters because torch is not fork-safe. Two failure modes
cost ~90 minutes:

1. **`nice`-at-launch kills PSOCK.** Wrapping the launch in `nice -19` made workers miss
   `makePSOCKcluster()`'s connect window ("4 of 16 workers failed to connect"). Renice
   *after* startup instead — a 30 s watcher loop achieves the same politeness.
1b. **`parLapply` logs nothing until every cell finishes**, so a master at 0% CPU with a
   log frozen at "Dispatching cells" is *normal*, not hung. I killed at least two healthy
   drivers on that misread. `bench_proportion` needed **103 min** where `categorical` needed
   46 — my 60-minute timeout was simply too short for it.
2. **torch hangs the R process AFTER the work is done.** `bench_categorical` wrote its
   `.rds`/`.md` and logged `done`, then sat at 0% CPU forever, blocking a chain that
   waited on process exit. Drive such chains off the driver's own completion marker in the
   log, not off `wait`.

Both are recorded as candidate brain lessons in
`docs/dev-log/2026-08-15-brain-lesson-drafts.md`.
