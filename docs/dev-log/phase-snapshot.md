# Phase Snapshot — pigauto

> **ONE entry. Replace it; never prepend.** Before adding a new entry, move the current one —
> **verbatim** — to `docs/dev-log/phase-snapshot-archive.md`. Repo state is truth; this block
> is a pointer, not a log.
>
> **Multi-lane:** do not replace this with a single-lane START HERE. Rehydrate from the
> coordination board’s Active Lane Split (every row).

- **2026-08-17 (Arc 3 closed; Arc 2 before it)** · **START HERE: `docs/dev-log/handover/2026-08-16-claude-handover.md` — read UPDATE 2 at the top** · `main` `584ea81`, PRs #167–#170 merged, **#171 open at the gate** · **Arc 3 headline: the in-house joint solver ESTIMATES Σ AND DISCARDS IT** (two estimators → numerically identical predictions in 12/12 sim cells; `max_iter=0` short-circuits to per-column BM). That is the AVONET continuous gap; Fisher-ML port does NOT ship · **FENCED: enabling refinement recovers ~85% of the gap but breaks SE coverage 0.925→0.618 — that is the EM restore needing its own G0 (slice 1b, scoped, NOT started)** · per-type λ dispatch merged (#170) removes the 19pp categorical regression · after-task `2026-08-17-arc3-solver-and-tail.md` (§9 ranked) · **open for Shinichi:** G0 for slice 1b, paper framing
