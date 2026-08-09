# Checkpoint — pigauto BACE wrap (Option B-minus)

**Date:** 2026-08-09 · **Branch:** `handover/2026-08-09-cursor`
**Repo:** `/Users/z3437171/Dropbox/Github Local/pigauto`

> Previous LOOP kit belonged to the **closed P0 lane** and is preserved verbatim at
> `docs/dev-log/loop-archive/2026-08-08-p0-lane/`. Nothing was destroyed.

## STATE

Wrap implemented and verified locally. Default path proven **bit-identical**. Blocked at an
OPEN GATE that predates this lane: `origin/main` no longer ships the file we just extended.

## ARCS DONE (verified, evidence in-repo)

- **S0** BACE installed from standalone `/Users/z3437171/Dropbox/Github Local/BACE`.
  `git rev-parse HEAD` = `ce8bc87d5b1dd4059bacfba23ab42c3f4dfe6080`, tree clean, SHA gate
  passed. In-tree `pigauto/BACE` untouched. `packageVersion("BACE")` 0.0.0.9000, built
  2026-08-09 14:35:46 UTC.
- **R1** `bace_imp()` returns class `"bace"` (`bace_imp.R:407`), so `bace_final_imp()`
  consumes it directly. `all_datasets` = list of length `n_final`, each a full data.frame of
  the same shape the wrapper already re-encodes — existing re-encode loop works unchanged.
  Source comment confirms each final run starts independently from the converged chain
  ("truly independent posterior draws suitable for Rubin's rules pooling").
- **S1** `final_imp = FALSE` + `n_final = 15L` added; branch on dataset selection only;
  validation errors for bad flags and for a BACE too old to export `bace_final_imp()`.
- **S2** `R/ovr_categorical.R` footnote: pigauto's OVR is no longer described as "the OVR
  strategy BACE uses"; records BACE's 2026-08 multinomial default and that pigauto's
  rank/conditioning motivation is unchanged. Docs-only.
- **M2** `tests/testthat/test-fit-baseline-bace-final-imp.R`:
  `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 18 ]` (the one skip is the inapplicable
  "BACE lacks bace_final_imp" branch).
- **N1** NEWS entry written. Explicitly states this is **not** an imputed-as-observed fix.
- **M1 (partial)**
  - Bit-identity vs `git show HEAD:R/fit_baseline_bace.R`, same seed:
    `mu identical TRUE`, `se identical TRUE`, `max|dmu| = 0`, `max|dse| = 0`.
  - `test-shipping-coverage.R`: `FAIL 0 | WARN 10 | SKIP 0 | PASS 54` — the T4
    `fit_baseline_bace` smoke now **runs** instead of skipping (it skipped before S0).
  - `test-bace-compat-eval.R`: `FAIL 0 | WARN 0 | SKIP 0 | PASS 20`.
  - Paths-scoped drift diff run (no rebase) — surfaced the OPEN GATE below.
  - `devtools::check()` still in flight.

## MEASUREMENT (toy regime — not a claim)

n=20 species, 2 traits (1 continuous + 1 binary), `runs=3`, `nitt=300`, one seed:
mean `se` default **0.0175** vs `final_imp=TRUE` (`n_final=5`) **0.2118** — about 12x larger.
Direction matches the G0 argument that chain-sweep SD understates imputation uncertainty.
**Single toy seed. Not a benchmark. Do not quote as a result.**

## ARC IN PROGRESS

M1 `devtools::check()` (background). Then S3 re-bench.

## OPEN GATES

1. **`origin/main` deleted `R/fit_baseline_bace.R`.** Commit `b615579 docs: prepare v0.10.0
   CRAN release surface` (Shinichi, 10 Jul 2026) removed the file plus its `NAMESPACE`
   export, `man/fit_baseline_bace.Rd`, and its `_pkgdown.yml` entry. This branch is 43
   behind main and still carries the pre-deletion file, so the wrap exists **only here**.
   Neither handover mentioned this. Needs a Shinichi decision before any landing path.
2. **S3 is an open gate for public claims.** Numbers may be produced; no pigauto-vs-BACE
   capability sentence ships without Shinichi.

## TRUTH

Working tree scoped to: `R/fit_baseline_bace.R`, `R/ovr_categorical.R`,
`man/fit_baseline_bace.Rd`, `NEWS.md`, `tests/testthat/test-fit-baseline-bace-final-imp.R`,
`LOOP/*`. `man/fit_pigauto.Rd` regeneration drift belonged to the P0 lane and was reverted.
uinit files and `gnn_attribution` artefacts remain unstaged. PR #155 untouched.

## RESUME

```
Wrap lane: S0–S2, M2, N1 done; M1 partial (check in flight); S3, X1, X2 open.
Default path proven bit-identical. OPEN GATE: main deleted R/fit_baseline_bace.R in b615579
— ask Shinichi whether the wrap should exist at all before planning any landing.
```
