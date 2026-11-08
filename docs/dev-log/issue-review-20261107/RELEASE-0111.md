# pigauto 0.11.1 release preparation

Status: candidate source only. Nothing has been submitted to CRAN.

The source is based on merge commit 0731a3bf4dc2a2136aedba7a6cb2fa5e7eb2be59.
Its release metadata changes set 0.11.1, credit Daniel Noble as a contributor,
update the installation guidance and prepare CRAN comments. The development
source diff at preflight had SHA-256
cc17d43b7de10fe5d438dba7b9d78ad013e8cb8e300c75da1576d6efab19bd92.

A local development preflight finished with zero errors, warnings or notes;
it preceded the final contributor and heading edits. A separate no-tests
check of the updated source also finished with zero errors, warnings or notes.
The candidate archive it built had SHA-256
97df5aea6c3749ab7edaea5f37e23505c4bbc31aaeb0ce2d063f4ee82ad0d967 and size
5,419,087 bytes. That predecessor has not been frozen as a release artifact.
Focused documentation tests passed. The current CRAN index returned no
reverse dependencies on 7 November 2026. The release-gate selftest passed.

Still required: merge reviewed release source; build one fresh artifact from
clean merged source; capture its inventory, hash and size; run exact-artifact
local and external platform checks; reconcile rights, maintainer consent,
timings and submission text; obtain the required independent release review;
and distinguish submission, confirmation, acceptance and availability.
Do not infer any of these from the development preflight.

The merge gate for PR 262 reported success after seeing only the skipped
pkgdown job, before platform checks registered. GitHub records the merge at
23:36:34 UTC on 7 November 2026; the platform jobs started seconds later.
Main-branch platform checks are being observed separately. This gate failure
must be corrected separately before it is trusted for another merge.
