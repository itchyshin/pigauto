# Plan: merge pigauto wrap #156 now that CRAN 0.10.0 is live

Meta: 2026-08-13 · Cursor · **G0 LOCKED and executed.**

## 🎯 GOAL (locked)

Land GitHub-dev BACE wrap (Option B-minus) on `origin/main` via PR [#156](https://github.com/itchyshin/pigauto/pull/156): `final_imp = FALSE` default, `n_final = 15L`, Version `0.10.0.9000`, Suggests `BACE`. After merge, `main` is **not** a CRAN tarball. Do not vendor-sync BACE. Do not restore EM. Do not Slack Dan. Do not merge the P0 parent (`fix/ci-install-libtorch`) in the same G0.

## Executed

- NEWS lead-in on `feat/bace-wrap-restore`: `a54e6a4`
- #156 undrafted, retitled, merge commit `416561b` on `origin/main` (2026-08-13 23:59 UTC)
- Verified: Version `0.10.0.9000`, Suggests `BACE`, `R/fit_baseline_bace.R` present
- CI on the NEWS commit was UNSTABLE (checks re-running); merge used the 2026-08-09 green matrix plus NEWS-only prose. Not a CRAN submit.

## Still out of scope

P0 parent → `main` · vendor-sync · EM · Slack · DRM.jl · CRAN-submit from this `main`
