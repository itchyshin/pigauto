# S0 leftover inventory — pkgdown-011-repair

Date: 2026-08-28. Worktree `pigauto-011-ship` @ `6fddd79` plus uncommitted repair edits.

## Named defects (G0) and source

| ID | Defect | File | Status after S1–S4 |
| --- | --- | --- | --- |
| 1 | Step 7 GNN-beats-BM table | `vignettes/getting-started.Rmd` | Replaced with gate-closed interpretation + `evaluate_imputation` / `r_cal` |
| 2 | Live-looking validation suite | `pkgdown/assets/validation_suite.html` | Tombstoned in place (`noindex` + withdrawn banner) |
| 3 | Body link contrast | `pkgdown/extra.css` | `main` / `.contents` / `.page-header` links `#2c5e4f` |
| 4 | Raw `[!WARNING]` | `README.md` (kept GFM tag; CSS + after_body paint cream/amber); `_pkgdown.yml` sidebar now uses the same GFM warning block | Source tagged; render check is S6 |
| 5 | BACE 42%→72%; stale `phylopars()`; dense §5 table | `vignettes/gnn-architecture.Rmd` | Attribution removed; threshold-joint text uses in-house solver; §5 is three prose blocks |
| 6 | AGENTS.md / CLAUDE.md | — | Out of scope (leftover 42%/72% in AGENTS.md **not** touched) |

## Leftover grep (in-scope only)

- `vignettes/`: no remaining `42%` / `72%` / `GNN reduces test RMSE`.
- `README.md` / `_pkgdown.yml`: `[!WARNING]` remains as the GFM marker the after_body script strips.
- `gnn-architecture.Rmd` still mentions `Rphylopars::phylopars()` once as **historical** (“earlier versions delegated…”). That is the correction, not the stale claim.
- `NEWS.md` / `AGENTS.md` still mention phylopars and 42%/72% — **out of lease**.

## Out of scope leftovers (do not “fix”)

- `script/validation_suite.html` still has the old suite (not in lease).
- Parked branches’ LOOP files (`evidence/gnn-sentinel-prerun`, `handover/2026-08-09-cursor`).
