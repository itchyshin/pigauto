# Win-builder results for `fe1600d7`, 2026-10-17

## Result

Both official Win-builder result logs completed with `Status: OK` for pigauto 0.11.0:

| Lane | Result | Environment | Start (UTC) | Runtime |
|---|---|---|---|---:|
| R-release | [00check.log](https://win-builder.r-project.org/ys47lC1Ld715/00check.log) | R 4.6.1 (2026-06-24 ucrt), Windows Server 2022 x64 | 2026-10-17 00:55:24 | 271 s |
| R-devel | [00check.log](https://win-builder.r-project.org/m53luBDZqnM7/00check.log) | R-devel r90686 ucrt, Windows Server 2022 x64 | 2026-10-17 01:04:23 | 271 s |

Both logs report successful incoming-feasibility, dependency/install/load, R code, examples, tests, rebuilt vignettes, PDF/HTML manual, and final checks. The R-release test stage took 116 seconds and vignette rebuild 12 seconds; R-devel took 114 and 11 seconds, respectively. No warning, error, or note is indicated in either final log.

## Artifact association and limit

The matching upload acknowledgement is `provenance/winbuilder-upload-fe1600d7-2026-10-17.md`. It records the exact frozen local input, SHA-256 `fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a`, unchanged local before/after hash, 5,136,083-byte size, and successful R-release and R-devel form responses echoing filename and size. Both result runs began after that upload window. This makes the results strongly operationally attributable to the frozen candidate and distinguishes them from the earlier same-name, 10-byte-smaller upload.

Neither the upload acknowledgement nor result logs expose a server-side SHA-256 or checked-archive byte size. Therefore this is a bounded operational attribution, not cryptographic proof that Win-builder checked the exact frozen bytes. Record both Windows outcomes as `Status OK; operationally attributed to fe1600d7; digest unverified`. The strict checksum-bound `platform-clean` gate remains open. The earlier form's `Currently unavailable` notice for R-devel was superseded by the completed R-devel result log.

## Independent review

Grace independently opened both official logs and confirmed their final statuses, environments, stages, and lack of a server digest. Rose independently rehashed the frozen local archive, compared the upload receipt and result timestamps, confirmed the size distinction from the earlier candidate, and reached the same bounded attribution. Both reviewers retain NOT READY for a checksum-bound platform-clean verdict.
