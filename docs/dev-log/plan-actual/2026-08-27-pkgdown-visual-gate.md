# Plan vs actual — pkgdown visual gate (2026-08-28)

## 1. Goal

Written Florence/Tufte visual verdict for pigauto 0.11 user-facing pkgdown pages; fail-to-ship defects fixed in CSS/yml only; #174 stays draft.

## 2. Implemented

LOOP/ kit rewritten for this lane (replaced a stale P0 kit). Local `build_site(devel = TRUE)`. See-pass on home + 5 articles + reference + methodology layout/tombstones. Rose fence. CSS/yml fixes for GFM warning token, logo overflow, flex `min-width`. Verdict + inventory written.

## 3a. Decisions and Rejected Alternatives

- `devel = TRUE` rather than the letter of `devel = FALSE`, so the site is this worktree’s 0.11 pages, not the installed library.
- Did not knit `eval=FALSE` plot chunks (G0 default).
- Did not add `_brand.yml`.
- Did not edit README / vignettes (Codex lease). JS strip in `_pkgdown.yml` instead of a `pkgdown/index.md` fork of README.
- Did not treat 390×844 headless crops as overflow after DOM widths showed no element past the layout viewport.

## 4. Files Touched

- `LOOP/GOAL.md`, `LOOP/arcs.md`, `LOOP/checkpoint.md`, `LOOP/ultra-plan.md`
- `pkgdown/extra.css`, `_pkgdown.yml`
- `docs/dev-log/2026-08-27-pkgdown-visual-gate.md`
- `docs/dev-log/2026-08-27-pkgdown-visual-inventory.md`
- `docs/dev-log/plan-actual/2026-08-27-pkgdown-visual-gate.md`
- `.unlazy/pkgdown-visual-gate/` (gitignored)

Not touched: `R/`, `tests/`, `NEWS.md`, `DESCRIPTION`, `vignettes/`, `man/`, PR #175.

## 5. Checks Run

| Check | Result |
| --- | --- |
| `pkgdown::build_site(devel = TRUE, preview = FALSE)` | exit 0; 2 interactive() example warnings |
| `pkgdown::check_pkgdown()` | No problems found |
| Headless screenshots + DOM `scrollWidth` | see-pass; 390 crop is a Chrome artifact |
| `gh pr view 174` | draft=true state=OPEN |

## 6. Tests of the Tests

`check_pkgdown()` is structural, not visual. Visual verify was Read-the-PNG plus measured widths. A false mobile-FAIL from cropped 390 PNGs was rejected by `innerWidth` / `scrollWidth` dump.

## 7a. Issue Ledger

None opened. Codex ticket suggested (not filed): README `[!WARNING]` → normal markdown so the after_body strip can go.

## 8. Consistency Audit

Navbar, reference groups, and tombstone wording match `_pkgdown.yml` and the G9 claim tests. No #175 strings on in-gate pages.

## 9. What Did Not Go Smoothly

- First `build_site` aborted because `docs/` already held tracked dev-log files.
- `working_directory` / sandbox mapped to the Dropbox checkout; all writes used the ship worktree with `all` permissions.
- Headless Chrome `--window-size=390` still laid out at 500px and cropped — wasted a CSS refine loop.

## 10. Known Residuals

Orphan agent markdown pages; empty logo alts; search placeholder; live pkgdown badge may show failing; Visualise sections are code listings; README still contains `[!WARNING]` source.

## 11. Team Learning

pkgdown does not render GitHub GFM alerts. Bootstrap `.row` flex + long `<pre>` is the overflow class to fix with `min-width: 0`, not a bigger logo hack.

## 12. Cross-Product Coverage

Reused Florence/Tufte method from sister pkgdown audits. Did not re-run May 2026 claim audits.

## Melissa reconcile

| Plan slice | Actual | Drift |
| --- | --- | --- |
| S0 inventory | written | none |
| S1 build | done, devel=TRUE | documented |
| S2 Emmy | in verdict | none |
| S3 Florence see | screenshots + widths | none |
| S4 Rose | in verdict | none |
| S5 CSS/yml | extra.css + after_body | none |
| S6 reverify | check_pkgdown + home rebuild | none |
| S7 close | this file | none |

Adaptive changes: 1 (`devel=TRUE`). Drift requiring G0 reopen: 0.
