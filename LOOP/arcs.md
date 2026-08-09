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
| M1 | Mechanical verify: paths-scoped drift, focused tests, `devtools::check()`, bit-identity | hand off / local | **PARTIAL** — wrap Rd `\usage` OK; focused tests 18/0/1; full check still 2 WARN / 3 NOTE (pre-existing). Log `/tmp/pigauto_m1c3/pigauto.Rcheck/00check.log` |
| S3 | Re-bench wrap (both paths) vs BACE `@ce8bc87` | hand off / local | DONE — coverage 0.672 → 0.940 (nominal 0.95); 1/5 seeds fail in `bace_final_imp` |
| S3b | Robustness follow-up: contextual error instead of bare MCMCglmm message | Cursor Models | DONE (`3bfd740`) |
| KEEP | Shinichi 2026-08-09: keep wrap; landing path TBD; no merge | this chat | **LOCKED** — landing executed on `feat/bace-wrap-restore`; still do not merge |
| N1 | NEWS entry | Other Models | DONE (+ robustness caveat) |
| X1 | Melissa reconcile (plan vs git reality) | Other Models | DONE — `docs/dev-log/plan-actual/2026-08-09-bace-wrap-reconcile.md` |
| X2 | After-task report + brain-write **proposal** (not a write) | Other Models | DONE — `docs/dev-log/2026-08-09-bace-wrap-brain-write-proposal.md`; vault untouched |
| RESTORE | New branch from `origin/main`: restore wrap + folded API; GitHub-dev only | this chat | DONE — `feat/bace-wrap-restore` @ `b180555` (do not merge) |
| CLOSEOUT | After-task + Melissa LIGHT + draft DO-NOT-MERGE PR; stop until CRAN | this chat | **DONE** — PR [#156](https://github.com/itchyshin/pigauto/pull/156) draft; STOP until CRAN |

## Not green, stated plainly

Wrap Rd `\usage` **OK**. Full `--as-cran` check is **not** 0/0/0: 2 WARNINGs + 3 NOTEs,
all pre-existing (bench rds in `data/`, undeclared `jsonlite`, `.uinit`, `results.tsv`,
BACE not on CRAN). Do not claim check-green. Full testthat suite still not claimed.

## KEEP wrap / landing (Shinichi 2026-08-09)

Keep the wrap. Landing lives on `feat/bace-wrap-restore` @ `b180555` with draft
DO-NOT-MERGE PR [#156](https://github.com/itchyshin/pigauto/pull/156). Do **not**
merge to `main` until BACE is on CRAN / the v0.10.0 CRAN cut ships. Do not rebase
the handover branch. Slack / public pigauto-vs-BACE / vendor-sync / EM remain
parked. See `LOOP/checkpoint.md`.
