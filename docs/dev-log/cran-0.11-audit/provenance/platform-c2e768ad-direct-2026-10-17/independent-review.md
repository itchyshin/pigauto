# Independent review of the original-archive platform matrix

Reviewed 2026-10-17 by a separate read-only platform results reviewer. Scope is run 39070922656 on frozen c2e768ad, with separate static checker review before each retry.

Verdict: PASS for the three-cell exact-original-archive matrix. The separate candidate-source CI result is recorded in source-ci-359.json and verified by the driver; all three jobs succeeded. This additional readback closes the source-CI prerequisite for platform-clean. This is not the final release panel or a submission approval.

The reviewer independently reran verify_results.py and obtained DIRECT_ARCHIVE_MATRIX_VERIFIED_3_OF_3. Raw Windows release, Windows devel and Linux check logs each report Status OK. Incoming checks, examples, vignette rebuilding, indexed PDF manuals and HTML manuals pass. Actual manual logs read Rd2.ind; the Windows manuals have 84 pages and the Linux manual has 90.

All 29 dependencies match the frozen DESCRIPTION, are available and installed. Download/pre/post identities match the pinned original gzip SHA-256 and size. The run and workflow merge SHA agree across all records. The fetched synthetic merge has parents 780c084 and reviewed head a550071; executable checker/workflow files have an empty diff. The executed commit's seven checker file hashes match the binding receipt.

All 96 retained compressed files were checked against their compressed hashes, decompressed hashes and downloaded originals. The matrix summary agrees with raw records. Source comparison found 215 byte-identical regular files, 20 generated files, and only expected DESCRIPTION normalization and generated metadata fields.

Timing review: Windows release 599.377 seconds, devel 656.047 seconds, Linux 324.048 seconds. Tests took 379, 417 and 192 seconds respectively. Retain Windows devel's 10.934-minute total for the separate submission timing decision. The earlier 5.92-second compare_methods example observation remains historical evidence, even though the necessary indexed-manual repair run ended Status OK.

Every runner reports 3,195 passing expectations, zero failures, 161 unclassified warnings and 86 identical skips. The CRAN-mode reporter suppresses warning messages. The skipped features and dependency states remain outside this matrix's claim. Raw Windows OS labels differ between Python and R; preserve the build 26100 and both labels without inventing distinct machines.

This review does NOT cover deployed-site parity, a new post-merge archive, submission timing approval, final release closure, CRAN submission or acceptance. Source CI is verified separately. Historical failed attempts retain their original verdicts.

A separate reader/governance pass found five receipt issues before publication: one stale nested rung, a pending validator label, a duplicate tree hash typo, stale panel wording and a pending assessment label. All were corrected and independently rechecked. The reviewer rehashed all six data streams plus NOTICE from the frozen archive, verified the actual validator log and its digest, and checked the report assessment's exact draft hash. Bounded verdict: PASS for publishing the evidence commit; final published-head readback remains required.
