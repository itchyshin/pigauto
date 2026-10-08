# Exact-tarball Windows builder uploads, 2026-10-06

At 2026-10-06 21:26 UTC, the exact frozen audit tarball was uploaded through the
R-release and R-devel forms on [win-builder](https://win-builder.r-project.org/upload.aspx).
Both postbacks displayed the accepted filename `pigauto_0.11.0.tar.gz`, size
5,122,013 bytes, and content type `application/x-gzip`.

| Field | Value |
|---|---|
| Source commit | `bb5835d1b214d783da7b0b414df99aa6ba926bc7` |
| Tarball SHA-256 | `d34d469981386546b4ad03b9c7aab81674556708eeaf524c1c7365b6c68dec26` |
| R-release builder | Windows Server 2022, R-release 4.6.0 |
| R-devel builder | Windows Server 2022, R-devel 4.7.0 |
| Builder result | Pending; no check logs or result email observed yet |

The builder's main page says it emails the package maintainer when checks finish
and keeps result directories for about 72 hours. The source package is public;
the builder warns that uploaded files and logs have no confidentiality guarantee.
The upload requests pre-release package checks. No CRAN submission has occurred.
Exact-artifact Windows evidence remains an open G7 gate until both results are
retrieved and reviewed. The browser showed the file metadata after each form
submission; it did not provide a build ID or result URL at submission time.

## Targeted result search

On 2026-10-06, after Shinichi authorized a narrowly scoped search, Gmail was
searched for the exact filename `pigauto_0.11.0.tar.gz` and for the package and
R-release/R-devel builder terms. No matching message IDs were returned, and no
email messages were opened. This records only that the two result messages were
not found in the connected account at search time; it does not establish that
either builder failed or completed. G7 remains partial pending the two logs.

A second exact-filename search at 2026-10-06 19:52 MDT again returned no message
IDs. No other messages or mailbox surfaces were inspected.

## Recheck of older Windows results, 2026-10-08

The current PR #228 conversation records two Win-builder result notices dated
2026-10-07 at 17:40 and 17:48 UTC for older R-release and R-devel uploads. The
R-release result at
https://win-builder.r-project.org/98A2ziQk5Ob4/00check.log identifies pigauto
0.11.0 on Windows Server 2022 with R 4.6.1 and ends with `Status: 1 ERROR`.
Its test output reports two failures in `test-check-pigauto.R` lines 47-48:
the default `check_pigauto()` result was not ready. The paired R-devel log
directory is recorded in the PR conversation as
https://win-builder.r-project.org/DHFjzzWXlqui; its result notice also reports
one error. I opened and inspected the R-release check log in the browser. The
R-devel result is recorded from the PR conversation; direct navigation to its
log was blocked by the browser.

These results cannot be tied to a tarball checksum or source commit. The PR
conversation says both predate the later `5fbbec8c...` uploads; they therefore
also predate the `10573d4f...` candidate. The current source on `main` has the
`check_pigauto(gnn = FALSE)` default and skips the torch runtime probe unless
`gnn = TRUE`, recorded in the merged PR #229. That makes the earlier failure
consistent with pre-fix source, but does not prove which source was checked.

Do not use these older results as either a pass or a failure verdict for the
`10573d4f...` candidate. The exact-artifact Windows-results gate remains open:
rerun R-release and R-devel on the
final post-merge tarball, and retain the upload and result receipts with its
SHA-256. The existing 10573 artifact is itself not final while documentation
and provenance edits remain outstanding.
