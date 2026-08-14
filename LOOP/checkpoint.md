# Checkpoint — pigauto BACE wrap (Option B-minus)

**Date:** 2026-08-13 · **Landed:** `origin/main` `416561b` (PR #156 MERGED)
**Wrap worktree:** `.worktrees/bace-wrap-restore` (`feat/bace-wrap-restore` @ `a54e6a4`)
**Handover checkout:** `handover/2026-08-09-cursor` (docs; dirty uinit + gnn unstaged)
**PR:** https://github.com/itchyshin/pigauto/pull/156 (**MERGED**)
**Repo:** `/Users/z3437171/Dropbox/Github Local/pigauto`

> Previous LOOP kit belonged to the **closed P0 lane** and is preserved verbatim at
> `docs/dev-log/loop-archive/2026-08-08-p0-lane/`. Nothing was destroyed.

## STATE

**WRAP MERGED TO MAIN (2026-08-13).** `origin/main` `416561b`, Version `0.10.0.9000`,
Suggests `BACE`. PR [#156](https://github.com/itchyshin/pigauto/pull/156) MERGED.
This is GitHub-dev, **not** a CRAN tarball. Do not CRAN-submit. Do not re-restore.
Do not merge P0 parent in this lane. Slack / vendor-sync / EM / public
pigauto-vs-BACE still parked. D-43: M1c full check WARN/NOTE remains pre-existing.

## KEEP-WRAP DECISION (Shinichi, 2026-08-09)

- **Keep** `fit_baseline_bace()` + opt-in `final_imp` / `n_final = 15L`.
- Landing: new branch from `origin/main` (executed). **Do not merge to `main`.**
- `n_final = 15L` is Shinichi's explicit AskQuestion pick (also in the `/goal` paste), not
  merely Ada's IF-YOU-DO-NOT-MIND default.
- Slack to Dan remains **parked** per the `/goal` paste ("No Slack to Dan").

## ARCS DONE (verified, evidence pasted below)

- **S0** BACE installed from standalone `/Users/z3437171/Dropbox/Github Local/BACE`.
  SHA gate: `ce8bc87d5b1dd4059bacfba23ab42c3f4dfe6080`, tree clean. In-tree `pigauto/BACE`
  untouched. Built 2026-08-09 14:35:46 UTC.
- **R1** `bace_imp()` returns class `"bace"`; `bace_final_imp()` consumes it directly;
  `all_datasets` matches the shape the wrapper already re-encodes.
- **S1** `final_imp = FALSE` + `n_final = 15L`; dataset-selection branch only.
- **S2** `R/ovr_categorical.R` footnote (docs-only).
- **M2 / CLOSEOUT re-verify** `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 18 ]` (fresh 2026-08-09 closeout run on restore worktree).
- **N1** NEWS written; explicitly **not** an imputed-as-observed fix claim.
- **S3** measured (see below).
- **M1** drift + focused tests + bit-identity DONE. M1c wrap Rd `\usage` **OK**;
  full check still 2 WARN / 3 NOTE (pre-existing). Log
  `/tmp/pigauto_m1c3/pigauto.Rcheck/00check.log`.
- **X1** `docs/dev-log/plan-actual/2026-08-09-bace-wrap-reconcile.md`.
- **X2** brain-write **proposal** (not a write) at
  `docs/dev-log/2026-08-09-bace-wrap-brain-write-proposal.md`.
- **RESTORE** `feat/bace-wrap-restore` @ `b180555` from `origin/main` `bf46991`.
- **CLOSEOUT** after-task + Melissa LIGHT + draft PR #156. Then stop.

## M1c — wrap Rd green; full check still red (pre-existing)

Wrap Rd is **not** what killed the earlier driver. `man/fit_baseline_bace.Rd` `\usage`
matches formals exactly (`final_imp`, `n_final` included).
`tools::checkRd("man/fit_baseline_bace.Rd")` is silent. No wrap-only Rd/test edit was needed.

Finished check (existing tarball, `_R_CHECK_FORCE_SUGGESTS_=false`,
`--as-cran --no-tests --no-vignettes --no-manual`), log
`/tmp/pigauto_m1c3/pigauto.Rcheck/00check.log`:

```
* checking for missing documentation entries ... OK
* checking for code/documentation mismatches ... OK
* checking Rd \usage sections ... OK
* checking Rd contents ... OK
* checking examples ... OK
* DONE
Status: 2 WARNINGs, 3 NOTEs
```

WARNINGs (pre-existing, not wrap): `data/bench_*.rds` not allowed in `data/`;
`jsonlite` `::` import undeclared in tests. NOTEs: CRAN incoming (version
`0.10.0.9000` + **BACE not in mainstream repos** + README URL 301); `.uinit`
hidden dir; top-level `results.tsv` / `run.log`.

Focused test (pasted this closeout, restore worktree `/tmp/pigauto-bace-wrap-restore`):

```
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-fit-baseline-bace-final-imp.R")'

ℹ Loading pigauto
SKIP: 'test-fit-baseline-bace-final-imp.R:127:3' ----------
Reason: installed BACE does export bace_final_imp(); nothing to assert

[ FAIL 0 | WARN 0 | SKIP 1 | PASS 18 ]
```

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

**Caveat:** 1 of 5 seeds failed inside `bace_final_imp()` with *"Mixed model equations
singular"*. The chain path succeeded on all 5. Attributed in `3bfd740` (fail loud, no silent
fallback). **Regime fence:** simulated BM only, one n, one trait mix, 5 seeds. Wrap-config vs
wrap-config — **not** a pigauto-vs-BACE comparison.

## LANDING PATH (executed + draft PR — not merged to main)

Why main deleted the wrap: `b615579` (Shinichi, 10 Jul 2026) *docs: prepare v0.10.0 CRAN
release surface* removed `R/fit_baseline_bace.R`, `man/fit_baseline_bace.Rd`, the NAMESPACE
export, the `_pkgdown.yml` entry, the `[T4]` smoke, and **BACE from Suggests**.
Current `origin/main` is `0.10.0` (`bf46991`). BACE is not on CRAN, so Suggests BACE is a
CRAN blocker.

**LOCKED:** option **A after #155**. Executed: `feat/bace-wrap-restore` @ `b180555`.
Draft PR [#156](https://github.com/itchyshin/pigauto/pull/156) is the wait-for-CRAN
pointer. **Do not merge. Do not `gh pr merge`.**

## OPEN GATES

1. **Wait for CRAN (human).** Do not merge #156 / `feat/bace-wrap-restore` to `main`
   until BACE is on CRAN or the v0.10.0 CRAN cut has shipped.
2. **S3 gate for public claims.** Numbers above are internal wrap-config comparison only. No
   pigauto-vs-BACE capability sentence ships without Shinichi.
3. **Scope fence held:** S3 harness lives in `docs/dev-log/` (gitignored). Promoting it to
   `script/bench_bace_wrap_final_imp.R` needs a one-line scope OK.

## PARKED (unchanged)

- No public pigauto-vs-BACE sentence.
- No Slack to Dan.
- No vendor-sync of in-tree `pigauto/BACE`.
- No PR #155 work.
- No edit to either BACE tree.
- No EM restore.
- No rebase / merge to `main`.
- No DRM.jl. No DESCRIPTION 95% claim-gate.
- uinit files and `gnn_attribution` artefacts stay unstaged.

## TRUTH

Restore: `feat/bace-wrap-restore` @ `b180555` (pushed). Wrap source: `b57da54`, `3bfd740`.
Draft PR https://github.com/itchyshin/pigauto/pull/156 (DO NOT MERGE). After-task
`docs/dev-log/after-task/2026-08-09-bace-wrap-closeout.md`. Melissa LIGHT
`docs/dev-log/plan-actual/2026-08-09-bace-wrap-closeout-reconcile.md`. Handover dirty
uinit / gnn unstaged. Neither BACE tree touched. No merge to `main`.

## RESUME

```
STOP until CRAN. Wrap closeout DONE. feat/bace-wrap-restore @ b180555.
Draft DO-NOT-MERGE PR https://github.com/itchyshin/pigauto/pull/156 → main.
Suggests BACE would block CRAN if merged now. Do not merge. Do not re-restore.
Do not start DRM.jl / DESCRIPTION claim-gate / EM / Slack / vendor-sync /
public pigauto-vs-BACE. Focused tests this closeout: FAIL 0 WARN 0 SKIP 1 PASS 18.
D-43: ultra-plan not wholly done (M1c WARN/NOTE + wrap→main deferred).
```
