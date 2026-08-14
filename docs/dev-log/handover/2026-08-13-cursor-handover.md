# Session Handoff: pigauto wrap LANDED on main (Cursor)

Meta: 2026-08-13 ~18:00 MDT · from Cursor · `TARGET = cursor`  
Authoritative copy: this file. Chat is disposable.

You are **Cursor**. Remit is pigauto only. Wrap merge G0 is **DONE**. Do not re-merge #156. Do not CRAN-submit from `main`. Do not start DRM.jl from this file.

## Critical Context

1. **Wrap is on `origin/main`.** Merge commit `416561b` (PR [#156](https://github.com/itchyshin/pigauto/pull/156) MERGED). Version **0.10.0.9000**. Suggests `BACE`. `R/fit_baseline_bace.R` present. This is GitHub-dev, **not** a CRAN tarball.
2. **CRAN pigauto 0.10.0** published 2026-07-30. **CRAN BACE** still 404. Next CRAN cut must drop Suggests BACE or wait until BACE is on CRAN.
3. **P0 sibling** #155 remains merged only to `fix/ci-install-libtorch`, **not** `main`. Not this G0.
4. **PROTECTED trees.** Standalone BACE `@ce8bc87`. In-tree `pigauto/BACE` — no vendor-sync. No EM. Never `git add -A`. Never archive this project.

## Goals / mission

Nothing wrap-OWED. Option B-minus already on main: `final_imp = FALSE` default, `n_final = 15L`.

## Current Working State

- **Working:** `origin/main` `416561b` Version 0.10.0.9000 with wrap.
- **In progress:** nothing wrap-OWED.
- **Not working / blocked:** CRAN-submit from this `main` (Suggests BACE). P0 parent still off `main` until a separate G0.

## Next Immediate Steps (OWED only)

1. **Rehydrate.** Lane preflight. `gh pr view 156` should be MERGED. Read the Active Lane Split.
2. **Do not re-merge, re-restore, or CRAN-submit.**
3. **If Shinichi names a new lane** (P0 parent → main, next CRAN cut, DRM.jl), that is a **new G0 that day**. Until then stop.

## How to Resume

```text
Read AGENTS.md and docs/dev-log/handover/2026-08-13-cursor-handover.md. Run lane preflight. Continue only OWED steps. Wrap merge is DONE.
```
