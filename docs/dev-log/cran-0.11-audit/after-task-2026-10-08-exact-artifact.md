## 1. Goal

Record and independently review one exact locally checked pigauto 0.11.0 candidate tarball in the unmerged CRAN evidence branch. It records one candidate; final release approval remains open.

## 2. Implemented

Added the exact-artifact receipt, ordered archive inventory, full `R CMD check` log, and raw testthat output to the release-evidence branch. Bound them in `release-ledger.json` to the read-only tarball with SHA-256 `1bfed5ad7be61437f4fdb09ece053d6a40211b5a5b7da4b2c947c3343493b719`, built from a clean detached worktree at merged `main` commit `0b0f71fee838c6ed51ef832ed819270eeafaf29b`.

The exact tarball is 5,128,791 bytes with 248 entries. Its full inventory hash is `52192b529b009a47b9bc3e7f1687f543194994de42e9e508ee536ebbacc0deb5`. Its macOS/R 4.6.0 `R CMD check --as-cran --no-manual` ended with exit 0 and `Status: OK`. Testthat reports 3,150 passes, 0 failures, 161 warnings, and 86 skips.

Updated G8 to distinguish this new local candidate from the older archives and left G8 open. Updated G6 with Chrome's direct `file:` rejection and left visual review open. Marked earlier panel READY entries as applying only to their prior candidate. The release ledger remains `NOT_READY`.

## 3a. Decisions and Rejected Alternatives

Kept the raw `testthat.Rout` byte-for-byte. Git's generic whitespace check flags four trailing-whitespace lines in that captured R output; editing the log would invalidate its evidence hash. Ran the diff check on all other staged evidence files and retained the original log hash.

Recorded the tarball as `candidate_not_release_final`. It does not contain the pending source and site corrections, has no Windows result bound to this hash, and has no exact-hash panel verdict. Did not upload it to another service in this slice.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `docs/dev-log/cran-0.11-audit/provenance/exact-artifact-0b0f71f-2026-10-08.md`
- `docs/dev-log/cran-0.11-audit/provenance/exact-artifact-0b0f71f-2026-10-08.inventory.txt`
- `docs/dev-log/cran-0.11-audit/provenance/exact-artifact-0b0f71f-2026-10-08-R-CMD-check.log`
- `docs/dev-log/cran-0.11-audit/provenance/exact-artifact-0b0f71f-2026-10-08-testthat.Rout`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-exact-artifact.md`

The pre-existing unstaged `retirement-manifest.md` change was not modified.

## 5. Checks Run

- Recomputed archive SHA-256: `1bfed5ad7be61437f4fdb09ece053d6a40211b5a5b7da4b2c947c3343493b719`; it matches the read-only copy and the checked build output.
- Recomputed inventory SHA-256: `52192b529b009a47b9bc3e7f1687f543194994de42e9e508ee536ebbacc0deb5`; the independent reviewer confirmed all 248 ordered members match the tarball.
- Exact command `env -u NOT_CRAN R CMD check --as-cran --no-manual pigauto_0.11.0.tar.gz`: exit 0, incoming feasibility passed, `Status: OK`, on macOS Tahoe 26.7, arm64, R 4.6.0. The check log SHA-256 is `4798992b15bcc845b7d99ca7b331d17cffcf72342882da9f4cda5a7b685568d1`.
- `python3 -m json.tool docs/dev-log/cran-0.11-audit/release-ledger.json`: parsed successfully.
- `git diff --cached --check -- . ':(exclude)docs/dev-log/cran-0.11-audit/provenance/exact-artifact-0b0f71f-2026-10-08-testthat.Rout'`: passed. The unfiltered check reports whitespace in the immutable raw test output as described above.
- `node .../unlazy/scripts/gate-check.mjs --status docs/dev-log/cran-0.11-audit/GATES.md`: parsed 11 gates, 5 marked met and 6 unmet. Two marked-met gates have runnable CHECK commands with no approval record; status mode runs nothing, so those commands remain unverified by Unlazy re-execution.
- `python3 ~/shinichi-brain/tools/closeout.py check <absolute report path>`: failed after the structural report check passed. The acceptance-ledger audit found five unmet gates under `.unlazy/imputation-sim/` whose worktree is the separate `pigauto-imputation-sim` project. Those campaign gates are outside this lane; they were left untouched.
- `python3 ~/shinichi-brain/tools/slop_check.py <exact-artifact receipt>`: findings 0.
- Chrome loaded the official CRAN pigauto index, which lists version 0.10.0 published 2026-07-30. Chrome rejected direct opening of the local candidate site's `file:` URL under browser security policy.

## 6. Tests of the Tests

No test code changed in this evidence slice. The exact package check ran the tests shipped in the tarball. The independent artifact reviewer separately recomputed the archive hash, member inventory, retained log hashes, and test counts. No deliberately corrupted tarball or planted artifact-control test was run.

## 7a. Issue Ledger

- Fixed: GATES.md said the release ledger had empty artifact/evidence fields after this slice populated them. The two statements now identify the earlier snapshot as historical.
- Fixed: prior panel READY values could be mistaken as verdicts on this tarball. They now say `READY_FOR_PRIOR_CANDIDATE_ONLY`; review for the 1bfed5ad hash remains pending.
- Open: Windows results are not bound to this hash.
- Open: the basis for CRAN redistribution of bundled BirdTree-derived data remains unresolved.
- Open: the local visual-review attempt is blocked by Chrome's explicit protocol policy.
- Open: release-site corrections have not been merged or deployed, and exact-hash panel verdicts are outstanding.

## 8. Consistency Audit

Compared the receipt, ledger, full local check log, raw test output, frozen archive, and archive-order inventory. The independent `artifact_gate_status` reviewer confirmed the archive, inventory, log and test output match byte-for-byte; confirmed the exact check result and test counts; and found no remaining factual inconsistency after the stale-snapshot and panel-scope corrections.

The receipt distinguishes the tarball's verified bytes from its source-commit provenance, which rests on the recorded clean detached worktree and build command. The tarball itself does not encode a Git commit. Its DESCRIPTION omits `drmTMB` and `gllvmTMB` dependency fields, but optional-backend workflow checks are separate from this artifact receipt. Chrome showed the public CRAN page still lists 0.10.0. The release candidate version is 0.11.0. The live-site and redistribution-rights gates remain open.

## 9. What Did Not Go Smoothly

The browser rejected a direct Chrome `file:` URL and prohibited local-server, alternate-browser, or indirect workarounds. The first sandboxed write attempt to the separate worktree was denied; retrying the leased paths with scoped approval succeeded. The closeout helper resolved a relative report path inside the Shinichi brain repository and created a scaffold there. I removed that accidental scaffold, then reran the helper with the absolute pigauto worktree path. No brain change remains.

## 10. Known Residuals

The tarball is a locally checked candidate; a post-documentation tarball remains to be built. No Windows result is tied to this hash. Earlier Win-builder logs cannot be attributed to it because they do not identify an archive checksum. BirdTree redistribution rights are unresolved. G6 visual review, G7 final deployment and route checks, and G9 exact-hash independent panel review remain open. The release ledger remains `NOT_READY`; no merge, deployment, or CRAN submission occurred.

## 11. Team Learning

Memory receipt: loaded the repo's LOAD-FIRST manifest with `route.py pigauto`. Its guidance to preserve scoped evidence and separate tested evidence from broader claims shaped this receipt. A vault-first Shinichi-brain search and a follow-up all-project search returned no result; repo records and the current worktrees supplied the evidence. No brain memory files were changed.

Golden Set: not run; this slice changed release evidence, not package behavior or the known-mistake classes covered by the Golden Set.

For the closeout helper, pass an absolute report path when the target is a worktree outside the brain repository. Candidate and release-final artifact records must remain visibly distinct, and a previous panel READY verdict does not transfer to a new archive hash.

## 12. Cross-Product Coverage

Covers: the exact locally checked tarball identity, archive inventory, macOS/R 4.6.0 check, testthat summary, and evidence-ledger linkage for this candidate.

Does NOT cover: pending package source or documentation corrections, Windows results for this hash, rights to redistribute bundled tree data, visual inspection of the local site, post-correction deployment and full live route/sitemap checks, exact-hash Grace/Rose/Pat verdicts, CRAN submission, or scientific validity beyond the existing bounded recovery evidence.

**Style:** 2/10, moderate confidence; technical closeout with direct evidence labels and bounded claims. **Genre/coverage:** entire report. **Evidence/repair:** the candidate-versus-final distinction and open-gate list are specific; no change needed after checking the retained archive and logs. **Gates:** science = pass for this evidence report's scope; facts = pass against the current files and independent artifact review; references = pass for the official CRAN page linked in the gate ledger, with no scientific literature claims here. **Provenance:** self-review after independent artifact review, draft `exact-artifact-0b0f71f`, 2026-10-08; prior examples were seen and were not treated as current artifact evidence.
