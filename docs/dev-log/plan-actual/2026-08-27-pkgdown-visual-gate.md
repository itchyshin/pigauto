# Plan vs actual — pkgdown-011-repair (2026-08-28)

## Planned

S0 inventory → parallel S1 CSS, S2 gnn-architecture, S3 getting-started, S4 tombstone → S5 `build_site` + `check_pkgdown` (+ R CMD check if feasible) → S6 see-pass → S7 close. #174 stays draft.

## Actual

- S0–S4 done in the `pigauto-011-ship` worktree under the expanded visual-audit lease (vignettes + README included by this G0).
- S5: `pkgdown::build_site(devel = TRUE)` + `check_pkgdown()` clean. `R CMD check` skipped (no `R/` / `tests/` edits).
- S6: Read rendered HTML for the named pages (not screenshots of a live browser).
- S7: this file + `docs/dev-log/2026-08-27-pkgdown-visual-gate.md`.

## Drift

- LOOP was still the older visual-gate-only contract; rewritten to this repair contract before edits.
- `CreateGoal` was not available in this runtime; LOOP on disk is the durable goal.
- Parallel `pigauto_finishing` visual-audit CSS (overflow / mobile) was already in `pkgdown/extra.css`; kept, plus link-contrast rules.
- Full R CMD check not run.

## Still true

#174 draft. #175 parked. No push. No new public capability claim.
