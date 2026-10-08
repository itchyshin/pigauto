# Pre-PR source check, 2026-10-06

The audit-branch build is a pre-PR rehearsal. The exact release candidate is to be frozen from the merged source after the maintainer's merge gate.

The worktree was copied to `/private/tmp/pigauto-cran-011-src/pigauto` with `.git`, `graft`, `_site`, and `.unlazy` excluded. The copy had no symlinks. `R CMD build pigauto` produced `pigauto_0.11.0.tar.gz` (SHA-256 `a9c7ce7b65aec751edd0c8246697ecb51eda9ffe19c874566ad246b4b67c4852`). Its DESCRIPTION says version 0.11.0 and omits drmTMB and gllvmTMB dependencies. The tarball contains no `graft/`, `.unlazy/`, `_site/`, `script/`, `dev/`, `docs/`, `pkgdown/`, or `BACE/` subtree.

On macOS Tahoe with R 4.6.0, `R CMD check --as-cran --no-manual pigauto_0.11.0.tar.gz` finished with **Status: OK**. The retained `pre-pr-check.log` is the check's `00check.log` (SHA-256 `1b3d54d54d7c4cc7341ead4feb3b2ef320b4bc1ea0c1d5cc423bdc6968610346`). It includes successful incoming-feasibility, dependency, installation, documentation, example, test, and vignette checks. The source defaults audit also passed 126 assertions. The four isolated installed-library cells passed their fixture and real adapter checks; those cells installed from the worktree rather than this tarball.

The local pkgdown build passed a crawl of 65 rendered pages and 3,621 local references, with all 34 retired direct pages absent. Desktop and mobile views were inspected. Public deployment will be checked after merge.

The remaining gates are the nine gllvmTMB recovery seeds, which require approval for the measured run above three hours; the maintainer's audit-PR merge; a frozen post-merge tarball and platform checks; the open BirdTree redistribution-rights point; and the separate fresh release panel. Nothing has been submitted to CRAN.



## Superseding post-merge status (2026-10-06)

This is a pre-PR rehearsal record; its remaining-gates paragraph was true only at that
point in the audit and is superseded by
`../post-merge-verification-2026-10-06.md`. PR #226 is merged, the live deployment is
verified, and the exact merged-source tarball plus its full macOS/R 4.6.0 check and
merged-commit Ubuntu/macOS check jobs are recorded there. The earlier pre-PR artifact
remains a distinct historical artifact. The BirdTree rights issue and fresh release
panel remain open. Nothing has been submitted to CRAN.

The 2026-10-07 source/site checks and their raw logs are retained in this file's prior
merge history and summarized in `../GATES.md`. Their statements about the then-current
live site and absence of a post-edit tarball are historical; the 2026-10-08 ledger
records subsequent main, candidate, and site checks. The nine-seed recovery count above
is likewise a pre-PR snapshot; `../recovery/RESULTS.md` records the later ten-seed result.
