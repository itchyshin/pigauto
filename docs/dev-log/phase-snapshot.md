# Phase Snapshot — pigauto

> **ONE entry. Replace it; never prepend.** Before adding a new entry, move the current one —
> **verbatim** — to `docs/dev-log/phase-snapshot-archive.md`. Repo state is truth; this block
> is a pointer, not a log.
>
> **Multi-lane:** do not replace this with a single-lane START HERE. Rehydrate from the
> coordination board’s Active Lane Split (every row).

- **2026-08-15** · `origin/main` **UNCHANGED** `416561b` Version `0.10.0.9000` · rehydration of the 2026-08-14 handover reconciled with **zero discrepancies** (`cb285ad`) · P0 **merged + verified on `origin/arc/p0-onto-main` `f5e8416`, NOT landed, no PR** — 7/7 conflicts resolved, one **BLOCKER**: `test-safety-floor.R:712` strict val-floor passes on pristine `main`, fails on `main`+P0 (controlled) · claim-gate found **1 BLOCKING doc defect on `main`** (`gnn-architecture.Rmd` §5 promises `pred$se` is "Exact under BM"; default path broadcasts one `sd(observed)` per tip) · **defect filed: GitHub [#157](https://github.com/itchyshin/pigauto/issues/157)** — gate calibration scores a different delta surface than `predict()` delivers; **on `main` today, independent of P0** · Shinichi **stepping away — do not resume unasked** · **publication roadmap:** `docs/dev-log/2026-08-15-publication-readiness-roadmap.md` (tiered gaps: #157 decision → P0 land → bench re-runs → regime map) · **START HERE:** `docs/dev-log/coordination-board.md` · handover `docs/dev-log/handover/2026-08-15-p0-arc-handover.md` · P0 handover `docs/dev-log/handover/2026-08-08-p0-rose.md`
