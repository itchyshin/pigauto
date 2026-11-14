# Submission preparation checkpoint

This receipt records preparation for the remaining gates of PR #236. The original c2 archive remains platform-clean for its recorded source. The updated installation guidance changes shipped README and vignette bytes, so the next candidate requires a clean post-merge freeze and fresh artifact-dependent checks. Overall NOT_READY. No merge, deployment or CRAN submission has occurred. A separate source-only followup branch 259d7e98c86dd5aa7680367d9e621ad76037d323 contains exactly 11 reviewed files from main `780c084`; its approval and CI are separate from the unmerged evidence PR.

## Primary records refreshed

On 2026-10-17, the primary CRAN PACKAGES.rds record lists pigauto 0.10.0, published 2026-07-30. No direct reverse dependencies occur in Depends, Imports, LinkingTo, Suggests or Enhances among 25,421 indexed packages. The package page also returned HTTP 200. The archive directory returned 404; this does not establish absence of any pending submission. Remote GitHub matching refs for tags/v0.11 returned an empty list. Keep version 0.11.0; do not create a public tag as part of this preparation.

The refreshed [CRAN policy](https://cran.r-project.org/web/packages/policies.html), revision 6875, asks for short checks, examples lasting only a few seconds, and at most two simultaneous cores. It does not state a universal ten-minute total-check limit. The separate local release protocol calls for timing adjudication. The submission checklist and URL guidance were also refreshed. `primary-cran-refresh.json` records URLs, retrieval results, sizes and hashes; `cran-publication-and-revdeps.json` retains the relevant index result.

Shinichi is the sole Authors@R author, maintainer and copyright holder. The existing tree-use warranty with citation/credit and the distinction in NOTICE remain the rights basis. No newly discovered collaborator-consent requirement or independent BirdTree redistribution licence is asserted.

## Local warning capture

The exact c2 archive was installed with --install-tests into a task library. Only its pristine tests tree was extracted into a separate .Rcheck-style directory. NOT_CRAN remained unset, force-Suggests was enabled, testthat 3.3.2 CheckReporter was retained, and an actual CPU torch tensor succeeded. All process thread settings were capped at one.

The corrected replay took 183.448 seconds for tests and 205.568 seconds including installation. It returned 3,195 passing expectations, zero failures, 161 warnings and 86 skips, matching the original CRAN-mode totals. All 161 warning objects and their row-level messages, test names, source files, source lines and final calls are retained. The classifier rejects an injected unknown warning.

| Local warning category | Count |
|---|---:|
| Small validation sets |119|
| Additional conformal sample-size limits |2|
| Undefined correlation on constant values |9|
| Constant matched traits |4|
| Few tree-sensitivity draws |5|
| Covariates ignored by the joint baseline |5|
| Covariates ignored with the GNN off |3|
| Ordinal imbalance |12|
| Opt-in signal-gate fallback |1|
| Deliberately unmatched species dropped |1|

Independent review matched every RDS message to its classified CSV row, checked complete unique IDs and callers, and found no source blocker established by these warnings. They retain real consequences: tiny calibration sets do not establish coverage, constant correlations are undefined, warned covariates contribute nothing on those routes, few-tree summaries can be unstable, and dropped species are not modeled. See `independent-warning-review.md`.

This is a local replay only. Equal Windows counts do not establish equal warning identities. The final candidate requires renewed capture or demonstrated relevant equivalence; the final panel must adjudicate the retained scope. Do not describe the package as warning-free or broadly inferentially validated.

The first diagnostic used the installed-tests directory, whose DESCRIPTION file made source-only tests run against an installed layout. It returned 3,217 pass, 4 fail, 165 warnings and 80 skips. Those failures are retained in `attempt1-installed-layout/`; they are not a native CRAN check or a package-regression claim. The corrected harness did not change package tests or suppress warnings.

## Timing and package size

The final c2 Windows checks took 536.636 seconds on release and 629.689 seconds on devel. Tests took 333/407 seconds, examples 17/16 seconds, and all six vignette rebuilds 19/21 seconds. Earlier successful devel timing 656.047 seconds and the earlier 5.92-second example observation remain in the platform receipt. These update-check costs require an explicit maintainer decision after final-candidate timings exist.

The c2 archive is 5,136,713 bytes. Its uncompressed regular files total 7,148,405 bytes, including 2,313,444 bytes of data, 1,497,652 bytes under inst, and 1,026,612 bytes under man. No src files or forbidden docs/script/.github paths are present. These are predecessor measurements; recompute for the post-merge archive.

## Remaining sequence and approvals

1. Review and pass CI on the separate source-only followup, then obtain explicit approval to merge that source PR and permit its automatic main-push Pages deployment. PR #236 remains the unmerged evidence PR required by the approved plan. The plan reserves merges to Shinichi.
2. Verify deployment commit, all 44 retired routes, 62 retained routes, sitemap/search and the changed pages in the browser. Search counts must be measured from the new build rather than pinned to the old 613-entry site.
3. Rebuild from clean merged source and freeze the final archive. R 4.6.0 writes a fresh Packaged timestamp during build; do not assume byte identity or silently normalize it. Rerun all affected artifact gates if bytes differ.
4. Fill the four submission-ready conditional evidence rows from final proof: reverse dependencies, compiled diagnostics or justified non-applicability, external-service behaviour, and size/vignette budget. Obtain fresh Grace/Rose/Pat verdicts on that artifact, then run the submission-ready validator.
5. Present the exact final archive path, hash, size, cran-comments and timing decision for separate CRAN-upload and maintainer-confirmation approval. Nothing is submitted by this checkpoint.

The scope does NOT cover current candidate deployment, final post-merge artifact identity, Windows warning-message classification, final release-panel READY votes, CRAN submission, confirmation or acceptance.

Raw CSV/session bytes are preserved in matching `.raw.gz` files. The readable copies use LF endings and trim display-only trailing whitespace; the classifier writes LF. The exact installation patch is retained as `installation-guidance.patch.gz`. No warning message or classification was changed.
