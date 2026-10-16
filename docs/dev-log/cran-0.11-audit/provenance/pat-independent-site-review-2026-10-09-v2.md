# Independent site re-review: Pat, 2026-10-09

Verdict: **G7 PASS for the checked deployment**, superseding the earlier NOT ASSESSED verdict about whether durable receipts were available.

Pat independently rechecked the committed v2 receipts at evidence commit `7752956` and verified the Pages deployment run `37975894397` succeeded on `main` commit `d76804e768bf77f430f43dcf591633f9cd900dab`. The receipt hashes match the committed files. The site checks cover 44 retired routes returning 404, all 62 sitemap-derived URLs, 613 search entries with no retired target, and three passing verifier tests.

In Chrome, Pat confirmed the homepage shows version 0.11.0 and the literal `[!WARNING]` token; the Getting Started content and defaults are current. Direct live sitemap access was blocked in this re-review, so the sitemap claim relies on the durable verifier capture and its deployment-bound receipts.

This is a point-in-time deployment result. It does not establish future permanence, fix the homepage warning rendering, or close the exact-artifact Windows gate. Grace and Rose remain NOT READY because the R-release/R-devel diagnostics report failures and are not checksum-bound to the frozen tarball. Overall audit status remains NOT READY; PR #228 remains Draft and unmerged; no CRAN submission has occurred.
