# Checkpoint — pigauto BACE wrap (Option B-minus)

**Date:** 2026-08-09 · **Landing branch:** `feat/bace-wrap-restore` @ `b180555`
**Worktree:** `/tmp/pigauto-bace-wrap-restore` (from `origin/main` `bf46991`)
**Handover checkout (untouched dirty tree):** `handover/2026-08-09-cursor` @ `a5976b0`
**Repo:** `/Users/z3437171/Dropbox/Github Local/pigauto`

> Previous LOOP kit belonged to the **closed P0 lane** and is preserved verbatim at
> `docs/dev-log/loop-archive/2026-08-08-p0-lane/`. Nothing was destroyed.

## STATE

Shinichi 2026-08-09 chose **KEEP the wrap**. Landing timing: **A after #155**. #155 is now
**MERGED** to parent `fix/ci-install-libtorch` (not to `main`). Wrap restore **executed** on
isolated worktree `/tmp/pigauto-bace-wrap-restore`, branch `feat/bace-wrap-restore` @
`b180555`, pushed to `origin/feat/bace-wrap-restore`. GitHub-dev-only: **do not merge to
`main`** until BACE is on CRAN / v0.10.0 ships. No PR this turn. Handover checkout left
dirty (uinit + gnn artefacts unstaged). PARKED stays parked.

## KEEP-WRAP DECISION (Shinichi, 2026-08-09)

- **Keep** `fit_baseline_bace()` + opt-in `final_imp` / `n_final = 15L`.
- Design a landing path onto current `main` (`b615579` deleted the wrap for the v0.10.0 CRAN
  surface). **Do not merge this turn.**
- `n_final = 15L` is Shinichi's explicit AskQuestion pick (also in the `/goal` paste), not
  merely Ada's IF-YOU-DO-NOT-MIND default.
- Slack to Dan remains **parked** per the `/goal` paste ("No Slack to Dan").

## ARCS DONE (verified, evidence pasted below)

- **S0** BACE installed from standalone `/Users/z3437171/Dropbox/Github Local/BACE`.
  SHA gate: `ce8bc87d5b1dd4059bacfba23ab42c3f4dfe6080`, tree clean. In-tree `pigauto/BACE`
  untouched. Built 2026-08-09 14:35:46 UTC.
- **R1** `bace_imp()` returns class `"bace"`; `bace_final_imp()` consumes it directly;
  `all_datasets` matches the shape the wrapper already re-encodes. BACE source confirms each
  final run starts independently from the converged chain.
- **S1** `final_imp = FALSE` + `n_final = 15L`; dataset-selection branch only.
- **S2** `R/ovr_categorical.R` footnote (docs-only).
- **M2** `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 18 ]` (re-verified 2026-08-09 this turn).
- **N1** NEWS written; explicitly **not** an imputed-as-observed fix claim.
- **S3** measured (see below).
- **M1** drift + focused tests + bit-identity DONE. M1c wrap Rd `\usage` **OK**;
  full check still 2 WARN / 3 NOTE (pre-existing). Log
  `/tmp/pigauto_m1c3/pigauto.Rcheck/00check.log`.
- **X1** reconcile at `docs/dev-log/plan-actual/2026-08-09-bace-wrap-reconcile.md`.
- **X2** brain-write **proposal** (not a write) at
  `docs/dev-log/2026-08-09-bace-wrap-brain-write-proposal.md`.

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
hidden dir; top-level `results.tsv` / `run.log`. First `--as-cran` run without
`FORCE_SUGGESTS=false` ERRORed earlier on missing Suggests `RSpectra` /
`Rphylopars` / `rgbif` and never reached Rd. PDF-manual run hit pre-existing
LaTeX `Σ` (U+03A3) in non-wrap Rd and was killed mid "without index".

Focused test (pasted):

```
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-fit-baseline-bace-final-imp.R")'

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

## LANDING PATH (executed 2026-08-09 — not merged to main)

Why main deleted the wrap: `b615579` (Shinichi, 10 Jul 2026) *docs: prepare v0.10.0 CRAN
release surface* removed `R/fit_baseline_bace.R`, `man/fit_baseline_bace.Rd`, the NAMESPACE
export, the `_pkgdown.yml` entry, the `[T4]` smoke, and **BACE from Suggests**.
`cran-comments.md` in that commit states user-facing pages contain **no live BACE
integration** and BACE pages are historical benches only. Current `origin/main` is `0.10.0`
(`bf46991`), Suggests has `jomo` not `BACE`, and `0cc11c5` records a CRAN resubmission
candidate. BACE is not on CRAN, so Suggests BACE is a CRAN blocker.

Options (do **not** rebase this 43-behind branch onto main):

1. **New branch from current `origin/main`** that restores the wrap file + folded wrap API
   (`b57da54` + `3bfd740`), tests, NEWS, NAMESPACE, pkgdown. Cherry-pick of wrap commits
   alone will not apply — the file is absent on main; restore first, then apply the delta.
2. **Wait until after #155** is off the `R/` lock, then do (1). #155 does not restore the
   wrap; waiting only sequences traffic.
3. **Rebase / merge `handover/2026-08-09-cursor` onto main.** Forbidden; 12 ahead / 43
   behind; mixed P0 history + uinit dirt.

**LOCKED (Shinichi AskQuestion 2026-08-09):** option **A after #155**. Executed: new branch
`feat/bace-wrap-restore` from `origin/main`, wrap + folded `b57da54`/`3bfd740` API, T4
smoke, OVR footnote, NEWS 0.10.0.9000, BACE in Suggests (no Remotes). Keep wrap
**GitHub-dev-only** until BACE is on CRAN / the v0.10.0 cut ships. Do not merge wrap to
`main` as part of the CRAN surface. Do not rebase the handover branch. No PR into `main`
this turn (branch on origin is enough).

Focused restore-branch test (pasted):

```
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-fit-baseline-bace-final-imp.R")'

SKIP: 'test-fit-baseline-bace-final-imp.R:127:3' ----------
Reason: installed BACE does export bace_final_imp(); nothing to assert

[ FAIL 0 | WARN 0 | SKIP 1 | PASS 18 ]
```

## OPEN GATES

1. **Landing (in progress; wait for CRAN cut).** Wrap is on `origin/feat/bace-wrap-restore`
   @ `b180555`. Do **not** merge to `main` until BACE is on CRAN / v0.10.0 ships. Draft PR
   only if marked DO NOT MERGE; prefer wait.
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
- uinit files and `gnn_attribution` artefacts stay unstaged.

## TRUTH

Restore commit: `feat/bace-wrap-restore` @ `b180555` (pushed). Source wrap commits on
handover: `b57da54`, `3bfd740`. Handover checkout `a5976b0` left dirty; uinit /
`gnn_attribution` unstaged. Neither BACE tree touched. No merge to `main`. No PR.

## RESUME

```
KEEP wrap. Landing EXECUTED on feat/bace-wrap-restore @ b180555 (worktree
/tmp/pigauto-bace-wrap-restore, pushed). GitHub-dev-only until BACE on CRAN / v0.10.0
ships. Do not merge to main. Do not rebase handover. Focused tests 18 pass / 1 skip.
Slack / public BACE claim / vendor-sync / EM parked.
```
