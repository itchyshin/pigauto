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
| M1 | Mechanical verify: paths-scoped drift, focused tests, `devtools::check()`, bit-identity | hand off / local | **PARTIAL** — drift + tests + bit-identity DONE; `check()` in flight |
| S3 | Re-bench wrap (both paths) vs BACE `@ce8bc87` | hand off / local | TODO — smoke-first |
| N1 | NEWS entry | Other Models | DONE |
| X1 | Melissa reconcile (plan vs git reality) | Other Models | TODO |
| X2 | After-task report + brain-write **proposal** (not a write) | Other Models | TODO |

## Blocking discovery (raised during M1)

`origin/main` **deleted** `R/fit_baseline_bace.R` in `b615579 docs: prepare v0.10.0 CRAN
release surface` (Shinichi, 10 Jul 2026) — removed from `R/`, `NAMESPACE`, `man/`, and
`_pkgdown.yml`. This branch is 43 behind and still carries the pre-deletion file. The wrap
exists on this branch only. See the OPEN GATE in `LOOP/checkpoint.md`.
