## 1. Goal

Record the completed candidate-source R CMD check matrix for source/docs PR #231 and preserve the boundary between candidate-source checks and the exact-artifact gate.

## 2. Implemented

Added the latest source-pinned matrix result to the CRAN 0.11 gate ledger and this dated after-task record. The receipt identifies source commit `c5a1cb971f57646f7ad9873be5d337bb2b940dee` and workflow run #37872160178.

## 3a. Decisions and Rejected Alternatives

Counted the successful Ubuntu release, Ubuntu devel, and macOS jobs as candidate-source evidence only. Kept G0, G7, G8, and G9 open because the rights/provenance basis, deployed site, final frozen artifact, and its independent review remain unverified. Did not treat this workflow as the required force-Suggests check on the final tarball.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-candidate-matrix.md`

## 5. Checks Run

- Chrome, GitHub Actions run #37872160178: Success; all three matrix jobs completed. Ubuntu R release: 11m12s. Ubuntu R-devel: 12m09s. macOS R release: 19m06s, including the focused MPS prediction test (4m05s) and `R CMD check` (10m55s). Total workflow runtime: 19m16s.
- `node ~/.codex/skills/unlazy/scripts/gate-check.mjs --status docs/dev-log/cran-0.11-audit/GATES.md`: 11 gates; 7 met and 4 open (G0, G7, G8, G9).
- `git diff --check -- docs/dev-log/cran-0.11-audit/GATES.md docs/dev-log/cran-0.11-audit/after-task-2026-10-08-candidate-matrix.md`: passed (exit 0; no whitespace errors).
- `python3 ~/shinichi-brain/tools/slop_check.py docs/dev-log/cran-0.11-audit/after-task-2026-10-08-candidate-matrix.md`: 0 hits in 601 words.
- `Rscript --vanilla ~/shinichi-brain/tools/check-after-task.R /Users/z3437171/.codex/worktrees/cran-011-pr231-verified/pigauto/docs/dev-log/cran-0.11-audit/after-task-2026-10-08-candidate-matrix.md`: required sections and negative-space check passed; the subsequent whole-ledger re-verification was interrupted because G0/G7/G8/G9 are still open. It produced no whole-ledger acceptance verdict.

## 6. Tests of the Tests

No test code or CI workflow was changed. The three independently reported matrix jobs ran the package's configured checks on the candidate source commit. No negative-control test was needed for this evidence-only update.

## 7a. Issue Ledger

- Resolved for this candidate-source snapshot: Ubuntu R release, Ubuntu R-devel, and macOS R release CI jobs all succeeded.
- Still open: BirdTree redistribution basis and remaining G0 provenance details; deployed-site verification (G7); the final post-merge tarball and exact checks (G8); independent review of that artifact and site evidence (G9).

## 8. Consistency Audit

Compared the PR head and run identity shown in Chrome with the source commit recorded in the ledger. Confirmed the run completed all three named matrix jobs. The ledger explicitly limits the claim to candidate source and retains the later release gates as open. The source commit tested here is unchanged by this evidence-only report.

## 9. What Did Not Go Smoothly

The closeout helper initially resolved a relative output path under the brain's docs tree. I removed that misplaced draft and recreated the report at the explicit path in the attached pigauto worktree. Chrome then stopped responding while I tried to inspect an already completed job's detailed log; the run summary and saved job result still provided the needed completion evidence.

## 10. Known Residuals

This candidate-source matrix does not prove the exact frozen tarball passes local or platform checks. The live website has not been checked after deployment of the current source changes. The rights basis for the bundled BirdTree-derived files is still not independently documented. No merge, deployment, or CRAN submission occurred.

## 11. Team Learning

Memory receipt: loaded the pigauto LOAD-FIRST manifest with `route.py pigauto`; the current lane preflight confirmed one active Codex lane and inspected the ledger path before editing. The main release lesson applied here is to bind each result to its exact source commit and keep candidate-source evidence separate from frozen-artifact evidence.

Golden Set: not in scope for this evidence-only update.

## 12. Cross-Product Coverage

Covers: the three configured GitHub Actions R CMD check matrix jobs for the candidate source commit, including the focused macOS MPS prediction test.

Does NOT cover: exact frozen-tarball checks, force-Suggests validation, optional-model interoperability beyond the existing independently recorded gates, deployed-site state, CRAN acceptance, or inferential validity beyond the bounded recovery evidence already recorded.
