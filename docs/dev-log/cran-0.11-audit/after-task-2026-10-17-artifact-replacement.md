> Superseded for current status by `after-task-2026-10-17-winbuilder-and-ledger.md`. The completed exact-hash as-cran check and upload acknowledgements below are newer; this report remains as the incident record.

# After-task: replacement frozen artifact and exact-hash check, 2026-10-17

## 1. Goal

Continue the approved pigauto 0.11 audit after the prior frozen tarball was overwritten during a local reproducibility build. Preserve the release boundary: no merge, public upload, deployment, or CRAN submission.

## 2. Implemented

- Rebuilt and froze a replacement tarball from clean source commit `780c08470b263cc26b4139c711a5c144fd29389c` at `/private/tmp/pigauto-release-candidate/sha256-fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a/pigauto_0.11.0.tar.gz`.
- Bound the archive hash, size, inventory, source-file comparison, local check outcomes, rights-related file hashes, and exact-hash MI integration matrix in `provenance/exact-artifact-fe1600d7-source-binding.md` and the updated release ledger.
- Marked the old `0fccf4bb0a179434fe15fa3a7ff514112bc5220e60c4b3e025fa00536f5ca29d` archive superseded and its local check, MI, Windows upload, and review receipts historical only.
- Requested and received an independent artifact identity/check review from Grace; requested an independent platform-binding recheck from Rose. The replacement-hash MI review passed independently within the bounded scope recorded in `provenance/independent-review-fe1600d7-2026-10-17.md`.

## 3a. Decisions and Rejected Alternatives

- Highest proven rung is `source-clean`; next is `tarball-clean`. The local ordinary check used `--no-manual`; the as-cran incoming stage failed on DNS, so the former hash's `tarball-clean` claim cannot transfer.
- Kept the overall verdict `NOT_READY`; no current-hash Windows result exists and no server digest attests to the prior upload.
- Did not upload to R-hub, because that publishes the package and results. No publication approval was given.
- No source behavior or scientific defaults were changed.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-16-blocker-refresh.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-17-artifact-replacement.md`
- `docs/dev-log/cran-0.11-audit/provenance/exact-artifact-fe1600d7-source-binding.md`
- `docs/dev-log/cran-0.11-audit/provenance/exact-artifact-fe1600d7-macos-check.log`
- `docs/dev-log/cran-0.11-audit/provenance/exact-artifact-fe1600d7-testthat.Rout`
- `docs/dev-log/cran-0.11-audit/provenance/exact-artifact-fe1600d7-as-cran-network-block.log`
- `docs/dev-log/cran-0.11-audit/provenance/mi-matrix-fe1600d7-provenance-complete/` (manifest and 12 logs)

## 5. Checks Run

- Replacement tarball SHA-256: `fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a`; size 5,136,083 bytes; 250 archive entries; sorted path inventory SHA-256 `d7b9b0e59ea051b3a22d7f347d8d925e6265b076b238ac893410e79d5b61329d`.
- Source checkout HEAD is `780c08470b263cc26b4139c711a5c144fd29389c` and `git status --short` is empty. Archive comparison: 215 source files byte-identical, one build normalization (`DESCRIPTION`), and 20 build-generated regular files.
- Exact-file ordinary `R CMD check --no-manual --run-donttest --timings`, `NOT_CRAN` unset and force-Suggests true: `Status: OK`; testthat 3,194 pass, 0 fail, 161 warnings, 86 skips. R 4.6.0, macOS arm64.
- Exact-file `R CMD check --as-cran --run-donttest --timings`: stopped at CRAN incoming feasibility because DNS could not resolve CRAN/Bioconductor. It did not complete package checks.
- Earlier no-`--no-manual` run ended with a manual-PDF warning/error because `pdflatex` is absent; the successful rerun explicitly used `--no-manual`.
- Exact-hash MI matrix: four backend states, 12/12 install/fixture/real steps passed; archive hash was unchanged before/after each step. Manifest SHA-256 `e81c14f45d3d3837a25999f21d4f83f60276a918f532adfd5c9f2ddba0bb1991`. Summed step execution time: 217.7 seconds.
- Direct archive streaming confirmed `inst/NOTICE` and six data-object hashes match the earlier receipt.
- Grace independently confirmed new archive hash, size, inventory, source comparison, local check, and retained log hashes. Rose independently confirmed Windows evidence is absent for the replacement hash; prior response size 5,136,073 differs from replacement size 5,136,083.

## 6. Tests of the Tests

The exact-MI runner independently checks the required tarball SHA before starting and before/after each of its 12 steps, stores per-step exit codes and hashes, and verifies package availability in each cell. All 12 steps report exit code 0 and unchanged archive hash. The independent MI reviewer rehashed the archive, runner, harness, five executed sources, manifest, and all 12 logs and returned PASS within scope; no material discrepancy was found.

## 7a. Issue Ledger

- Resolved for this candidate: frozen archive identity and clean-source binding are recorded; ordinary exact-archive R CMD check and exact-hash MI integration matrix passed within their stated scopes.
- Open: `--as-cran` incoming feasibility, because this host cannot resolve CRAN/Bioconductor.
- Open: manual PDF check, because `pdflatex` is unavailable.
- Open: Windows/platform evidence; previous upload acknowledgement reported a different size and no server-side hash.
- Open: final release closure remains unmet; independent reviews confirm the artifact, MI matrix, and bounded site receipt while platform checks remain open.
- Open: merge, deployment, and CRAN submission remain separate approval gates.

## 8. Consistency Audit

The current release ledger, GATES opening checkpoint, source-binding receipt, no-manual check log, test output, blocked as-cran log, and replacement MI manifest now name `fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a`. The previous blocker-refresh report is visibly superseded. Site evidence remains tied to the unchanged source commit `780c08470b263cc26b4139c711a5c144fd29389c`. The old artifact's Win-builder, as-cran, and MI evidence is retained only as historical and is not counted for this hash.

## 9. What Did Not Go Smoothly

The local reproducibility build reused the conventional filename of the old frozen tarball and overwrote it. I immediately disclosed the collision, searched the available temp, worktree, and repository locations, and found no copy matching the old digest. I have not represented the old checks as applying to the replacement. The host also lacks `pdflatex` and DNS access to CRAN/Bioconductor.

## 10. Known Residuals

The replacement passes ordinary package check with manual PDF generation excluded and passes the bounded exact-hash MI integration matrix. It is not `tarball-clean`: CRAN incoming feasibility and manual PDF checks remain unmet, and independent MI review is pending. It is not `platform-clean`: Windows evidence is absent for the new digest. The evidence PR remains unmerged; nothing has been published or submitted.

## 11. Team Learning

A hash-named frozen artifact must never be rebuilt in place or share a conventional output path. Build to a fresh staging directory, hash first, and copy only into a new digest-named freeze after checking the destination is absent. Exact-artifact evidence must become stale immediately when the archive digest changes.

## 12. Cross-Product Coverage

This checkpoint covers the replacement archive, local package check, exact-hash optional-backend MI workflow, and reviewer status. It does NOT cover as-cran incoming success, manual PDF generation, Windows/platform readiness, statistical recovery/coverage, CRAN acceptance, merge, deployment, or publication.
