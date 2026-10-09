# After-task: local website gate re-verification

## 1. Goal

Reverify the current candidate source's local website build and retirement controls, then retain the evidence in the CRAN 0.11 release-evidence checkout. This is one pigauto audit lane with multiple checkouts.

## 2. Implemented

Updated G5's output handling so it prints only the decisive completion marker while retaining the full build output. Re-ran G5 and G5b through Unlazy against source commit `cf88d78e37d08130bd33823e1fd4594f3afa2d4f`. Both passed. Retained the wrapper, pkgdown render, and crawler logs with SHA-256 values in GATES.md.

## 3a. Decisions and Rejected Alternatives

Kept the release ledger's gate states controlled by Unlazy. G5 and G5b are met; G7 deployment, G8 final post-merge artifact, and G9 independent final review remain open. Followed Shinichi's direction that these worktrees belong to one pigauto lane. No merge, deployment, artifact freeze, or CRAN submission occurred.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/provenance/site-build-2026-10-09-cf88d78.log`
- `docs/dev-log/cran-0.11-audit/provenance/pkgdown-render-2026-10-09-cf88d78.log`
- `docs/dev-log/cran-0.11-audit/provenance/site-crawl-2026-10-09-cf88d78.log`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-gate-reverify.md`

## 5. Checks Run

- `python3 /Users/z3437171/shinichi-brain/tools/route.py pigauto` loaded the LOAD-FIRST manifest. The project routing directs release checks to source-bound evidence and preserves separate deployment and artifact gates.
- `node /Users/z3437171/.codex/skills/unlazy/scripts/gate-check.mjs --approve --timeout 900 /Users/z3437171/.codex/worktrees/cran-011-gate-reconcile/pigauto/docs/dev-log/cran-0.11-audit/GATES.md` ran G5 and G5b sequentially. G5 passed with 62 HTML pages and 613 search entries. G5b passed with 62 HTML pages, 3,584 local references, 34 retired pages, zero errors, and `SITE_CRAWL_OK`. The command returned nonzero only because G7, G8, and G9 are still open.
- `pkgdown::check_pkgdown()` reported no problems. Two non-interactive `@examplesIf interactive()` conditions evaluated false as expected.
- Raw log hashes: wrapper `af0193c9eb2f5c970553df996965ab800cb7d7adff5d73d34cb0719268e06660`; pkgdown `b6c576cad42c65fd7b7cfea06aa8b4f958d6c62ff665fabf399d9125fb0fc056`; crawler `666b188ae633c31388a96b76c6b704ae575a69281e6606c220539e00f8f159d4`.

## 6. Tests of the Tests

The G5 expectation requires the builder's final source and retired-route marker. G5b reads the exact fresh G5 output, confirms its embedded source commit, and checks routes, links, anchors, search, sitemap, and retirement controls. Unlazy reported both checks as passed. The prior wrapper failure was reproduced as an output-volume failure and corrected by printing only the marker; the full log remains available for inspection.

## 7a. Issue Ledger

Resolved: G5/G5b evidence had been pinned to an older website source and the G5 wrapper exceeded the gate output limit. Updated G5/G5b checks to current source `cf88d78` and bounded the wrapper output.

Open: G7 deployed-site verification, G8 final post-merge tarball and platform results, and G9 independent review of the exact final artifact.

## 8. Consistency Audit

Confirmed the source checkout remained clean at `cf88d78e37d08130bd33823e1fd4594f3afa2d4f`. The gate edits and evidence live in the release-evidence checkout; no source checkout files changed. GATES.md now records the passing local results and leaves deployed-site, exact-artifact, and independent-review gates open. The crawler's 34 retired-page controls and sitemap/search checks agree with the recorded local site output.

## 9. What Did Not Go Smoothly

The first Unlazy attempt printed the entire build log, exceeded its output pipe, and failed before writing a site-root receipt. G5b then failed because that prerequisite was absent. The corrected wrapper emitted only the success marker. The sandbox also required a scoped write escalation for the managed worktree and Unlazy approval receipt. `git diff --cached --check` flags two trailing-space lines repeated in the wrapper and pkgdown logs because the raw output is preserved byte-for-byte for its SHA-256; the edited Markdown files pass the whitespace check.

## 10. Known Residuals

This verifies the local candidate site only. It does NOT establish the deployed site, a post-merge frozen tarball, Windows or independent platform checks, final independent review, or CRAN acceptance. The three release gates G7-G9 remain unmet.

## 11. Team Learning

Memory receipt: loaded the pigauto LOAD-FIRST manifest through `route.py` and checked the pigauto CRAN 0.11 record in `MEMORY.md`; source-bound gates and explicit remaining release rungs shaped this work.

Golden Set: no package-code mistake class changed in this evidence-only slice, so the Golden Set was not in scope.

Team lesson: keep large raw logs as file evidence and make an automated gate print only a concise marker. Treat checkout count as an operational observation; the maintainer's clarification defines this audit as one lane.

## 12. Cross-Product Coverage

Covers: local pkgdown build, local route/link/anchor checks, sitemap/search retirement checks, and the local website source commit.

Does NOT cover: live deployment, exact post-merge source artifact, cross-platform package checks, independent final-artifact verdict, or submission.
