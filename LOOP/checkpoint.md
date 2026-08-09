# Checkpoint — pigauto BACE wrap (Option B-minus)

**Date:** 2026-08-09 · **Branch:** `handover/2026-08-09-cursor` · **HEAD:** `b57da54` (+ uncommitted S3 follow-up)
**Repo:** `/Users/z3437171/Dropbox/Github Local/pigauto`

> Previous LOOP kit belonged to the **closed P0 lane** and is preserved verbatim at
> `docs/dev-log/loop-archive/2026-08-08-p0-lane/`. Nothing was destroyed.

## STATE

Wrap implemented, measured, and verified. Default path proven **bit-identical** twice (before
and after the robustness follow-up). S3 measurement done and it supports the G0 decision.
One OPEN GATE that predates this lane blocks any landing conversation.

## ARCS DONE (verified, evidence pasted below)

- **S0** BACE installed from standalone `/Users/z3437171/Dropbox/Github Local/BACE`.
  SHA gate: `ce8bc87d5b1dd4059bacfba23ab42c3f4dfe6080`, tree clean. In-tree `pigauto/BACE`
  untouched. Built 2026-08-09 14:35:46 UTC.
- **R1** `bace_imp()` returns class `"bace"`; `bace_final_imp()` consumes it directly;
  `all_datasets` matches the shape the wrapper already re-encodes. BACE source confirms each
  final run starts independently from the converged chain.
- **S1** `final_imp = FALSE` + `n_final = 15L`; dataset-selection branch only.
- **S2** `R/ovr_categorical.R` footnote (docs-only).
- **M2** `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 18 ]`.
- **N1** NEWS written; explicitly **not** an imputed-as-observed fix claim.
- **S3** measured (see below).
- **M1** drift + focused tests + bit-identity DONE; `devtools::check()` still in flight
  (~40 min, buffered — no partial output available).

## S3 RESULT — wrap default vs `final_imp` (the measurement that matters)

Simulated BM, n=100 tips, 3 traits (2 correlated continuous + 1 binary), 5 seeds,
`runs=5 nitt=2000 burnin=500 thin=10`, held-out val+test cells, latent scale.
Harness: `docs/dev-log/2026-08-09-bace-wrap-s3-harness.R`.

| path | RMSE | cover95 | mean se | secs | seeds ok |
|---|---|---|---|---|---|
| default (chain averages) | 0.2274 | **0.672** | 0.0646 | 2.0 | 5/5 |
| `final_imp`, `n_final=15` | 0.2162 | **0.940** | 0.1766 | 7.9 | 4/5 |

Nominal target 0.95. The default path undercovers by ~28 points; the proper-MI path lands
essentially on nominal. RMSE also improves slightly (~5%). Cost ~3.9x runtime.

**Caveat, and it is a real one:** 1 of 5 seeds failed inside `bace_final_imp()` with
*"Mixed model equations singular: use a (stronger) prior"*. The chain path succeeded on all
5. The final phase refits MCMCglmm per draw and is genuinely less robust.

**Regime fence:** simulated BM only, one n, one trait mix, 5 seeds. Not AVONET, not
PanTHERIA. This is wrap-config vs wrap-config — it is **not** a pigauto-vs-BACE comparison.

## FOLLOW-UP APPLIED (uncommitted at time of writing)

`bace_final_imp()` is now wrapped in `tryCatch` and re-thrown with context pointing at
`final_imp = FALSE`, instead of surfacing a bare MCMCglmm message. Deliberately **not** a
silent fallback: a caller who asked for proper MI must not receive chain averages unknowingly.
Roxygen + NEWS document the failure mode. Bit-identity re-verified after this change
(`identical()` TRUE, max diff 0).

## OPEN GATES

1. **`origin/main` deleted `R/fit_baseline_bace.R`.** Commit `b615579 docs: prepare v0.10.0
   CRAN release surface` (Shinichi, 10 Jul 2026) removed the file plus its `NAMESPACE`
   export, `man/` page, and `_pkgdown.yml` entry. This branch is 43 behind and still carries
   the pre-deletion file, so the wrap exists **only here**. Neither handover mentioned it.
   Needs Shinichi before any landing path is designed.
2. **S3 gate for public claims.** Numbers above are internal wrap-config comparison only. No
   pigauto-vs-BACE capability sentence ships without Shinichi.
3. **Scope fence held:** the S3 harness lives in `docs/dev-log/` (gitignored). Promoting it
   to `script/bench_bace_wrap_final_imp.R` so it is reproducible for others needs a one-line
   scope OK, since the approved scope was `R/` + `man/` + `tests/` + `NEWS.md`.

## TRUTH

Committed `b57da54` (9 files, explicit paths). Uncommitted: the robustness follow-up in
`R/fit_baseline_bace.R` + `man/fit_baseline_bace.Rd` + `NEWS.md`, plus `LOOP/` updates and
`docs/dev-log/coordination-board.md`. uinit files and `gnn_attribution` artefacts remain
unstaged and untouched. PR #155 untouched. Neither BACE tree touched.

## RESUME

```
Wrap lane: S0-S3, M2, N1 done; M1 check in flight; X1/X2 open.
Default path bit-identical. S3: coverage 0.672 -> 0.940 (nominal 0.95), 1/5 seeds fail in
bace_final_imp. OPEN GATE: main deleted R/fit_baseline_bace.R in b615579 — ask Shinichi
whether the wrap should exist before designing any landing.
```
