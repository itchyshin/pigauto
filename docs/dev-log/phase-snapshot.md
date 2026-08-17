# Phase Snapshot — pigauto

> **ONE entry. Replace it; never prepend.** Before adding a new entry, move the current one —
> **verbatim** — to `docs/dev-log/phase-snapshot-archive.md`. Repo state is truth; this block
> is a pointer, not a log.
>
> **Multi-lane:** do not replace this with a single-lane START HERE. Rehydrate from the
> coordination board’s Active Lane Split (every row).

- **2026-08-17 evening (Arc 3 CLOSED, slice 1b landed)** · **START HERE: `docs/dev-log/handover/2026-08-16-claude-handover.md` — UPDATE 2** · `main` `f925e17`, **PRs #167–#172 all merged, zero open** · **Arc 3 headline: the joint solver estimated Σ and DISCARDED it** (`max_iter=0` short-circuits to per-column BM; two estimators → identical predictions in 12/12 sim cells) · **slice 1b (G0 granted) fixed a real variance bug** — `.mvn_estep_refine` added precisions of two DEPENDENT posteriors, coverage 0.925→0.618; now `refine_variance="conservative"` + divergence guard + opt-in `joint_refine_iter` (default 0L). Sim gate PASSED at λ=0.2 (RMSE −12%, coverage 0.943→0.964); **real-data gate FAILED on AVONET** (λ 0.993–0.998, refinement adds nothing there; all deltas within MCSE) · **NO Rphylopars dependency** — Suggests-only, one call site, unreachable by default · **open for Shinichi: (1) exact Σ⊗R conditional in-house vs keep the rphylopars yardstick vs document-and-move-on; (2) paper framing** · results `2026-08-17-refinement-results.md` + `...-sigma-recovery-results.md`, after-task `2026-08-17-arc3-solver-and-tail.md`
