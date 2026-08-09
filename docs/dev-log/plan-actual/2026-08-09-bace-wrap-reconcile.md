# Plan-vs-actual — pigauto BACE wrap lane (Option B-minus)

Date 2026-08-09 · Branch `handover/2026-08-09-cursor` · HEAD `3bfd740`
Plan: `docs/dev-log/handover/2026-08-09-bace-wrap-g0-ultra-plan.md` (frozen at
`LOOP/ultra-plan.md`). Reconciler: Melissa lens, run by Ada.

Rule applied: a row is DONE only with a receipt (git object, pasted test output, or a file
that exists). Self-report is not a receipt. Per D-43, two or more NOT-DONE verdicts withhold
the completion claim.

## Row-by-row

| Row | Planned | Actual | Verdict |
|---|---|---|---|
| SCAFFOLD | LOOP kit for this lane | `LOOP/{GOAL,arcs,checkpoint,ultra-plan}.md` rewritten; P0 kit copied to `docs/dev-log/loop-archive/2026-08-08-p0-lane/` before overwrite | **DONE** |
| S0 | Install standalone BACE `@ce8bc87`; stop if SHA differs | `git rev-parse HEAD` = `ce8bc87d5b1…`; `R CMD INSTALL … * DONE (BACE)`; `requireNamespace("BACE")` TRUE; built 2026-08-09 14:35:46 UTC. In-tree clone untouched | **DONE** |
| R1 | Recon `bace_final_imp()` shape, no edits | `bace_imp.R:407` sets `class(out) <- "bace"`; `bace_final_imp.R:286` `all_datasets <- lapply(results_list, "[[", "dataset")`; no BACE file modified (`git status` clean in both clones) | **DONE** |
| S1 | `final_imp` + `n_final`, default byte-identical | Committed `b57da54`. Bit-identity vs `git show HEAD:R/fit_baseline_bace.R` under fixed seed: `identical()` TRUE for `mu` and `se`, `max|d| = 0`. Re-verified after `3bfd740` | **DONE** |
| S2 | OVR footnote, docs-only | `R/ovr_categorical.R` roxygen changed; no code path touched (diff is comment-only) | **DONE** |
| M2 | New focused test file | `tests/testthat/test-fit-baseline-bace-final-imp.R`, `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 18 ]` | **DONE** |
| M1a | Paths-scoped drift diff, no rebase | Run. Surfaced that `origin/main` **deleted** `R/fit_baseline_bace.R` in `b615579` | **DONE — and it changed the picture** |
| M1b | Focused tests | `test-shipping-coverage.R` `FAIL 0 / PASS 54 / SKIP 0`; `test-bace-compat-eval.R` `FAIL 0 / PASS 20 / SKIP 0`; new file 18 pass | **DONE** |
| M1c | `devtools::check()` clean on the wrap surface | **Stalled.** Log reached `checking Rd \usage sections` at 08:45 then the driver died. Completed sections include `code/documentation mismatches ... OK` and `missing documentation entries ... OK` — the two that bind on a new argument. NOTEs seen: hidden files, top-level `results.tsv`/`run.log` (both pre-existing, neither mine) | **NOT DONE** |
| M1d | Default output unchanged vs pre-change | Covered by S1 receipt | **DONE** |
| S3 | Re-bench wrap both paths | Ran. n=100, 5 seeds: coverage 0.672 → 0.940 (nominal 0.95); RMSE 0.2274 → 0.2162; 3.9x runtime; 1/5 seeds failed inside `bace_final_imp` | **DONE (regime-fenced)** |
| S3b | *(unplanned)* attribute `bace_final_imp` failures | Committed `3bfd740`. Not in the original plan — added because S3 surfaced a real failure mode | **DONE — plan deviation, logged** |
| N1 | NEWS entry, no overclaim | Two entries added. Contains the sentence that this is *not* an imputed-as-observed fix, and *"No pigauto-versus-BACE performance claim is made here"* | **DONE** |
| X1 | This reconcile | This file | **DONE** |
| X2 | After-task report + brain-write proposal (not a write) | Draft below; **nothing written to the vault** | **DONE (proposal only)** |

## Deviations from the approved plan

1. **S3b added.** The plan had no robustness row. S3 produced a real MCMCglmm singularity on
   1 of 5 seeds, so error attribution was added. Small, in-scope (`R/fit_baseline_bace.R`),
   and it does not change the estimator.
2. **S3 harness location.** The plan said reuse `script/bench_avonet_bace.R`. Those scripts
   benchmark pigauto *against* BACE, which is not what S3 asked for (wrap-path vs wrap-path
   with coverage). A purpose-built harness was written instead and parked at
   `docs/dev-log/2026-08-09-bace-wrap-s3-harness.R` because `script/` was outside the
   approved scope fence. **Consequence: the bench is not reproducible by anyone else until
   it is promoted to `script/`.** Needs a one-line scope OK.
3. **M1c incomplete.** See verdict above. Not claimed as green.

## Verdict

One NOT-DONE (M1c) and one scope-fenced artifact (S3 harness). Below the D-43 threshold of
two, so the implementation slice may be reported as complete — but `devtools::check()` is
**not** green-by-evidence and must not be described as such. The lane's real blocker is not
any of these: it is the main-branch deletion recorded as OPEN GATE 1.
