## 1. Goal

Reconcile the pigauto CRAN 0.11 audit after PR #234 merged and the current-source GitHub Actions jobs remained queued; preserve the NOT_READY verdict until current-source site and exact-artifact gates pass.

## 2. Implemented

Updated the gate ledger and machine-readable release ledger with PR #234's merged source commit, current queued workflow evidence, live-site mismatch, a fresh full-site build, corrected candidate-site verifier, visual review of three key pages, and the fact that PR #228 is historical and merged. G5, G5b, and G6 now have current local evidence; G7, G8, and G9 remain open. Corrected the ledger's visual-review field to reflect the full local review.

## 3a. Decisions and Rejected Alternatives

Do not freeze another tarball while the deployed site is stale. Do not infer the cause of queued Actions from the public status page or other repositories' active jobs. Do not treat PR #228 as the required final unmerged evidence PR because it is now merged. Correct the candidate-site verifier's 568 expected nonempty paths to 567, because 46 of its 613 records contain `path: []`.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-15-queue-reconciliation.md`
- `docs/dev-log/cran-0.11-audit/provenance/` site build, crawl, verifier, negative-control, receipt, and Chrome review records
- `script/cran-0.11-site/verify_candidate_site.py`
- `script/cran-0.11-site/test_candidate_site.py`

## 5. Checks Run

- Parsed `release-ledger.json` with Python's JSON tool: passed.
- Fresh full-site build from `8082ef8273177ff664d49ac4bd643f7def671f20`: 62 HTML pages; `pkgdown::check_pkgdown()` reported no problems.
- Six crawler tests passed; crawler checked 3,584 local references and 34 retired pages with zero errors; build retirement assertions passed.
- Old candidate-site verifier reproduced its defect (`expected 568 paths; found 567`). The corrected verifier passed all 16 tests and returned `CANDIDATE_SITE_OK` against the fresh output.
- Chrome reviewed the fresh local Getting Started, tree-uncertainty, and multiple-imputation pages. A cache-busted live-page check still showed predecessor content.
- GitHub Actions browser check confirmed the three R-CMD-check jobs and pkgdown job remain queued.
- `cran_release_gate.py` returned `READY FOR CLAIMED RUNG` for `source-clean`; this does not claim tarball or release readiness.

## 6. Tests of the Tests

The JSON syntax check, crawler, candidate-site verifier, and verifier tests are directly reproducible. `check-after-task.R` passed the required 12-section shape check. The scoped Unlazy status check for `.unlazy/cran-0.11-audit` reports G3, G4, and G5 unmet: checksum-bound release/devel results, final independent artifact/site review, and the final unmerged evidence PR. Those are genuine later release gates. The after-task helper passed its structure check. Its repository-wide acceptance recheck reports the scoped audit ledger G3-G5 unmet and also flags the separate `imputation-sim` arc; that separate arc was not changed. A direct `gate-check --status --root <audit worktree> --scope cran-0.11-audit` confirmed only G3-G5 are unmet here, without executing or modifying gates. The release validator records only the `source-clean` rung.

## 7a. Issue Ledger

- Fixed in this checkpoint: stale current-state notes, an outdated browser-review status field, and the candidate verifier's incorrect expected search-path count.
- Deferred: obtain runner results, verify deployment, freeze a new exact artifact, obtain checksum-bound platform checks and independent verdicts, then open a new unmerged final evidence PR.
- No unrelated source issues were changed.

## 8. Consistency Audit

The gate file and JSON ledger agree that PR #234 is merged at `8082ef8273177ff664d49ac4bd643f7def671f20`, current local G5/G5b/G6 passed, the live deployment is stale, the previous tarball predates the correction, and the verdict is NOT_READY. Historical site and artifact passes remain labeled as historical rather than promoted to current-source evidence.

## 9. What Did Not Go Smoothly

GitHub Actions remained queued although GitHub Status reported Actions operational. The full build was initially thought blocked because R was absent from the default PATH; locating the installed R runtime and supplying the integrity-checked cache enabled the fresh build. The after-task helper is anchored to the brain checkout and does not support this detached worktree; its failed invocation made no report write. The report was then written directly in the owned checkout.

## 10. Known Residuals

G7 is unmet because the corrected page is not deployed. G8 and G9 are also unmet; no new tarball was frozen. The exact cause of the runner queue is unknown. The scoped audit acceptance ledger still has G3-G5 unmet: checksum-bound platform results, exact-artifact/site review, and the final unmerged evidence PR. The separate imputation-sim recheck is outside this workstream. This checkpoint does not establish any CRAN readiness rung beyond `source-clean`; `tarball-clean` is next.

## 11. Team Learning

Keep the deployment check separate from source merge: a merged reader correction can coexist with stale live content. Mark historical evidence explicitly when new source changes invalidate its scope. The queue status page alone cannot diagnose runner allocation.

Memory receipt: applied the pigauto release scope and exact-artifact gate context from the project memory pass.

Golden Set: not checked in this checkpoint.

## 12. Cross-Product Coverage

Covers the audit ledger, current full local website build and review, the merged reader correction's deployment state, the verifier correction, and current Actions queue evidence. Does NOT cover deployed-site verification, a new frozen tarball, platform checks on that tarball, broad statistical recovery campaigns, or CRAN submission.
