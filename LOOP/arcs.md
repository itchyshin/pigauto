# Arcs — pigauto BACE wrap (Option B-minus)

Branch `handover/2026-08-09-cursor`. Plan:
`docs/dev-log/handover/2026-08-09-bace-wrap-g0-ultra-plan.md`.

| Slice | Name | Bar | Status |
| --- | --- | --- | --- |
| SCAFFOLD | LOOP kit for this lane (P0 kit preserved to `docs/dev-log/loop-archive/2026-08-08-p0-lane/`) | — | DONE |
| S0 | Install standalone BACE `@ce8bc87` | hand off / local | DONE — SHA gate passed |
| R1 | Recon `bace_final_imp()` return shape | Cursor Models | DONE |
| S1 | Implement `final_imp` + `n_final` in `fit_baseline_bace()` | Cursor Models | DONE |
| S2 | OVR footnote (docs-only) | Other Models | DONE |
| M2 | New focused test file | Cursor Models | DONE — 18 pass / 0 fail / 1 expected skip |
| M1 | Mechanical verify: paths-scoped drift, focused tests, `devtools::check()`, bit-identity | hand off / local | **PARTIAL** — drift + tests + bit-identity DONE; `check()` in flight (~40 min) |
| S3 | Re-bench wrap (both paths) vs BACE `@ce8bc87` | hand off / local | DONE — coverage 0.672 → 0.940 (nominal 0.95); 1/5 seeds fail in `bace_final_imp` |
| S3b | Robustness follow-up: contextual error instead of bare MCMCglmm message | Cursor Models | DONE (uncommitted at checkpoint time) |
| N1 | NEWS entry | Other Models | DONE (+ robustness caveat) |
| X1 | Melissa reconcile (plan vs git reality) | Other Models | DONE — `docs/dev-log/plan-actual/2026-08-09-bace-wrap-reconcile.md` |
| X2 | After-task report + brain-write **proposal** (not a write) | Other Models | DONE — `docs/dev-log/2026-08-09-bace-wrap-brain-write-proposal.md`; vault untouched |

## Not green, stated plainly

`devtools::check()` died at `checking Rd \usage sections`; the full testthat suite was still
running at hand-off. Neither is claimed as passing. Both are Codex hand-off items.

## Blocking discovery (raised during M1)

`origin/main` **deleted** `R/fit_baseline_bace.R` in `b615579 docs: prepare v0.10.0 CRAN
release surface` (Shinichi, 10 Jul 2026) — removed from `R/`, `NAMESPACE`, `man/`, and
`_pkgdown.yml`. This branch is 43 behind and still carries the pre-deletion file. The wrap
exists on this branch only. See the OPEN GATE in `LOOP/checkpoint.md`.
