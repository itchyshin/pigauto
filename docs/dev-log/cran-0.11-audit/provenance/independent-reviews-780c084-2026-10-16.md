# Independent review receipts for commit 780c084

These are bounded reviews. They do not establish whole-package CRAN
submission readiness.

## Grace: exact artifact and local check

Verdict: local exact-artifact evidence verified; broader release gates remain
open. Grace verified the tarball's 5,136,073-byte size, SHA-256
`0fccf4bb0a179434fe15fa3a7ff514112bc5220e60c4b3e025fa00536f5ca29d`,
250-entry inventory and inventory digest. She independently reconstructed the
Git tree from the retained GitHub archive and matched it to commit
`780c08470b263cc26b4139c711a5c144fd29389c`. All 236 regular package files
matched the retained `R CMD check` input byte-for-byte. The retained macOS/R
4.6.0 check log says `Status: OK`; Grace verified its hash but did not rerun
the check. She did not pass the pending main-branch job or the pending Windows
results.

## Rose: Windows result binding

Verdict: Windows results are not bound to the frozen archive digest. The local
upload helper checked the expected SHA-256 before sending; Win-builder's HTTP
response acknowledged only `pigauto_0.11.0.tar.gz` and 5,136,073 bytes. No
server-side checksum or result log has arrived. The checksum-bound platform
gate is therefore not met by these receipts.

## Pat: current deployed site

Verdict: pass within the supplied current-site receipt. The fresh audit reports
44/44 retired routes as 404, 62/62 retained routes as 200, a 62-entry sitemap
with no retired target, and 613 search entries with no retired target. Pat
opened the live NEWS page and confirmed the flagged historical development
entry is beneath “Earlier development changes included in 0.11.0” and states
that no CRAN submission or public release is implied. Pat did not replay the
106 URL requests or recheck all retained-page content.

## Combined limit

The current deployment and exact local tarball check are verified. Windows
R-release/R-devel results remain pending and are not checksum-bound. Ubuntu
R-release source CI has since completed successfully, but is not a substitute
for exact-tarball platform evidence. No reviewer recommends claiming
platform-clean or submission-ready from this evidence set.
