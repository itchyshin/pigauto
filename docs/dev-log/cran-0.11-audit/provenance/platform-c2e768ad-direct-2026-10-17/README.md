# Original gzip platform evidence

Candidate: `pigauto_0.11.0.tar.gz`, source `9a457eb6b3363d413dbc502de88ff5bc602ed98a`, SHA-256 `c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2`, 5,136,713 bytes, 250 entries. Sorted inventory SHA-256: `d7b9b0e59ea051b3a22d7f347d8d925e6265b076b238ac893410e79d5b61329d`.

The original gzip passed the complete direct-archive matrix in [run 39070922656](https://github.com/itchyshin/pigauto/actions/runs/39070922656), with Status OK on all three runners. Candidate-source [CI 359](https://github.com/itchyshin/pigauto/actions/runs/39070922883) also passed all three jobs at head a5500712b771cb4e34401fb5d9b7681a13219cf4, including macOS prediction tests. Independent results review returned PASS for this bounded exact-archive matrix. Highest proven rung is platform-clean; next is submission-ready. Overall release verdict remains NOT_READY.

The checker downloads the public upstream R-hub object, rejects any mismatch in pinned identity, and passes the original absolute compressed archive to `R CMD check --as-cran --run-donttest --timings`. It does not invoke R CMD build. It reads a separate DESCRIPTION copy to install dependencies. Before and after identities include the run ID, workflow SHA, original path, SHA-256, size and inventory. Checks keep NOT_CRAN unset, force-Suggests enabled, manuals enabled, and a functioning CPU libtorch backend.

## Attempts retained

| Attempt | Workflow run | Verdict and cause |
|---|---|---|
| 1 | [39068078482](https://github.com/itchyshin/pigauto/actions/runs/39068078482) | Windows stopped before package checking because a trailing DESCRIPTION comma produced an empty dependency token. Remaining Linux work was cancelled. Checker parser repaired and tested against R's parser. |
| 2 | [39068501072](https://github.com/itchyshin/pigauto/actions/runs/39068501072) | Windows release/devel ended with 1 ERROR, 1 WARNING, 2 NOTEs; first manual errors name missing Courier metrics pcrr8t/pcrr8c. Linux ended with 1 NOTE for missing HTML Tidy and V8. All three passed tests, examples and vignette rebuilding. These are failed checks. |
| 3 | [39069727324](https://github.com/itchyshin/pigauto/actions/runs/39069727324) | Linux Status OK. Windows PDF without index and HTML manual passed, but indexed PDF still failed. Windows release: 1 WARNING, 1 NOTE; devel: 1 WARNING, 2 NOTEs, including compare_methods 5.92 elapsed seconds (2.82 CPU seconds). These are not platform passes. |

Attempt 4, [39070922656](https://github.com/itchyshin/pigauto/actions/runs/39070922656), completed successfully at checker commit a5500712b771cb4e34401fb5d9b7681a13219cf4. All three runners first recorded missing makeindex (empty executable path, command exit 127); targeted installation produced a working nonempty index and an indexed PDF through R. The actual package manuals also read their generated Rd2.ind files. No fallback-only manual result is counted.

| Platform | R | Check status | Total check seconds |
|---|---|---|---:|
| Windows x64, build 26100 | 4.6.1 release | OK | 599.377 |
| Windows x64, build 26100 | devel r90686 | OK | 656.047 |
| Ubuntu 24.04.5 | 4.6.1 release | OK | 324.048 |

Windows OS labels differ between reporting tools: Python identifies Windows Server 2025, while R reports Windows Server 2022 x64 (build 26100). The raw records retain both; these are the same named GitHub runner legs. Windows total time is close to or above ten minutes, so this receipt does not establish submission timing approval. Attempt 3's 5.92-second example observation remains in the record despite the later Status OK results; the retry was required to repair manual indexing, not to obtain a faster timing observation.

The final matrix retains 3,195 passing expectations, zero failures, 161 unclassified test warnings and 86 skips on each runner. All 29 archive-declared dependencies are available and installed; NOT_CRAN is unset and force-Suggests is enabled. Raw check, test, manual, vignette, session, dependency and probe records are compressed in attempt4/. The retained manifest supplies original paths and raw/compressed hashes. Full job logs retain the identity-control outputs and setup commands. Server artifact metadata retains the three artifact archive digests.

The recorded workflow SHA f78a42c4329e3249ed91ed6bb54ab4c39a4da8b8 is GitHub's synthetic PR merge. Its parents are base 780c08470b263cc26b4139c711a5c144fd29389c and reviewed head a5500712b771cb4e34401fb5d9b7681a13219cf4. Fetching that object and comparing the workflow/checker files gives an empty diff; workflow-source-binding.json records the parents and file hashes.

Attempt 2's compressed raw logs, dependency/session manifests, commands and identities are in `attempt2/`; `attempt2-retained-manifest.json` identifies raw and compressed hashes. `attempt2-artifact-metadata.json` records the server artifact IDs and archive digests. The first Windows failure log is retained separately. No failed attempt is counted as a pass.

The current-source comparison finds 215 common regular package files byte-identical to the frozen archive. DESCRIPTION has only expected build normalization; the recorded comparison preserves the list of build-generated files. The checker, scripts and development logs are build-excluded. The archive was not rebuilt in this slice.

## Test and scope limits

The CRAN-mode test reporter reports 3,195 passes, zero failures, 161 warnings and 86 skips in attempt 2, matching the exact local macOS receipt. It suppresses warning text in CRAN mode; the 161 messages remain unclassified. The warnings are not described as benign or absent.

Independent review classified the 86 skips: 54 explicit CRAN exclusions, 2 unavailable BIEN-cache checks, 3 slow/live integration checks, 24 source-only documentation/helper checks, and 3 missing-dependency branches skipped because those dependencies are installed. The test counts describe executed expectations, not complete feature validation. The separate exact-artifact optional drmTMB/gllvmTMB matrix retains its own bounded interoperability claim.

This receipt does NOT cover deployment, final post-merge archive freezing, final release closure, CRAN submission or acceptance. Historical R-hub Windows results still checked rebuilt source, and historical Win-builder results remain associated with the predecessor archive. Neither is relabelled by this new direct-archive work.


## Reverification

Download the three named artifacts from run 39070922656 into separate named directories, then run `python script/cran-0.11-platform/verify_results.py RESULT_DIRECTORY`. `matrix_controls.py RESULT_DIRECTORY` exercises the matrix oracle on isolated evidence copies: the positive matrix passes; wrong hash, failed status, wrong CRAN mode and wrong check command each fail. Original downloaded evidence is unchanged.

The portable release-validator regression tests and selftest were rerun successfully. Its claimed-rung result and the prepublication unlazy reverification are retained with this receipt. The publication/readback gate is closed separately in the final PR description. The release plan still has six met gates, three unmet gates (deployed site, final post-merge freeze, final independent closure), and zero abandoned gates. Final independent closure includes the remaining submission timing decision. No merge, deployment or CRAN submission is performed.

Reader reconciliation corrected a stale nested status and a one-character duplicate trees300 hash typo. All six current data hashes and NOTICE were recomputed from the frozen gzip and match the current authoritative hash map. Older embedded snapshots remain verbatim and are explicitly corrected by current-data-hash-audit.json.
