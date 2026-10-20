# Pat independent readback: frozen successor archive

Verdict: **PASS** for the frozen source/archive identity and rights/inventory scope.

I independently checked the merged source checkout at `b2e6054685e6bdedb073ba46e058af2c0a619a18` and confirmed the worktree is clean. The retained main CI record for run `39090594763` has conclusion SUCCESS at the same source commit; all three recorded jobs succeeded (Ubuntu R release, macOS R release, Ubuntu R devel).

The frozen archive is `/private/tmp/pigauto-release-candidate/sha256-a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb/pigauto_0.11.0.tar.gz`. Its actual SHA-256 is `a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb`; size is 5,137,837 bytes; mode is 0o444; it contains 250 entries. I independently verified the sibling freeze receipt equals `FINAL.json`, the content-addressed directory name, and that the inventory file hash is `d7b9b0e59ea051b3a22d7f347d8d925e6265b076b238ac893410e79d5b61329d`. Recomputed sorted tar member names match the inventory byte-for-byte. All members are regular files or directories; no symlink or special member was found.

Of 216 packaged tracked files compared against the exact source commit, 215 are byte-identical. DESCRIPTION is the sole difference. I independently parsed both DCF files: all 18 source fields match after whitespace normalization, and the built DESCRIPTION adds only `Author`, `Maintainer`, `NeedsCompilation`, and `Packaged`. The 7 data/NOTICE streams have hashes matching both the frozen receipt and predecessor archive receipt; the rights-preserved flag is true and the forbidden-entry list is empty.

The concrete identity/inventory/rights checks all passed. This review does not cover results of the already-running native check or MI fit runner, the separate human assessment of the long Mac timing, or authority/submission/acceptance at CRAN. The prior 14ee verdict and artifact were not modified.

Machine-readable checks: [actual-verdict.json](actual-verdict.json).
