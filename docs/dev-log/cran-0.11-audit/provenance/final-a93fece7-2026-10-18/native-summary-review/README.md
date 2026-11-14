# Pat independent native-summary readback

**Verdict: PASS for summary-to-raw-evidence consistency.** No native checks or warning replays were rerun. The readback covers four retained native `00check.log` files, their timing receipts, dependency inventories and test output counts. It does not waive or adjudicate the recorded warnings.

All four native logs report `Status: OK`, each declared-dependency inventory lists 29 available packages, and each test output reports 0 failures, 161 warnings, 86 skips and 3,195 passes. Each raw `00check.log` yields the exact five reported stage lines in the summary. Each per-platform elapsed value matches its `check-exit.txt` receipt. The test summary appears twice in each `testthat.Rout`; both occurrences agree on every platform. The parser accepts identical duplicate counts and refuses absent or conflicting counts.

| Platform | Check elapsed s | Examples | Tests | Vignettes | PDF | HTML |
|---|---:|---|---|---|---|---|
| macos-arm64 | 452.674 | 9s CPU / 10s elapsed | 204s CPU / 256s elapsed | 12s CPU / 15s elapsed | 6s CPU / 7s elapsed | 3s CPU / 4s elapsed |
| ubuntu-latest-release | 534.689 | 17s CPU / 17s elapsed | 333s CPU / 333s elapsed | 19s CPU / 19s elapsed | 3s CPU / 3s elapsed | 4s CPU / 4s elapsed |
| windows-latest-devel | 516.487 | 16s reported | 324s reported | 17s reported | 3s reported | 3s reported |
| windows-latest-release | 627.141 | 16s reported | 406s reported | 20s reported | 4s reported | 12s reported |

The Windows stage numbers are preserved as single reported seconds; the readback does not infer CPU versus elapsed time. The macOS whole-runner receipt reports 459.703 seconds; the native check itself reports 452.674 seconds in `check-exit.txt`. Those measure different scopes.

Identity layers remain distinct:

- Package source commit: `b2e6054685e6bdedb073ba46e058af2c0a619a18` (local source checkout clean at this commit).
- Frozen package archive SHA-256: `a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb`.
- Local native `check.R` SHA-256: `1e641c87e699fba2846aaa9c8958fa77348ecd22221b5ed08e1ccf31c2de9fb0`.
- CI workflow source SHA: `999645fae22107d838e73e90e33996c025f2575a`; platform artifact head SHA: `5e1014afadbae91dd2883d3a074e3e9d253b33cd`; direct-platform run ID: `39091958419`.
- Three downloaded CI artifact ZIP digests are recorded separately from the package tarball hash in `actual-readback.json`.
- Thirteen known-answer parser controls are synthetic input-only checks, separately hashed; they are not platform or package results.

All four actual archive identity pre/post receipts bind the package digest, size, inventory hash and 250-entry count. All platform rows tie back to the same archive digest. The parser AST for `parse_stage_timing` and `parse_native_stages` is unchanged from the prior reviewed copy.

Raw log hashes, `testthat.Rout` hashes, dependency-name-set hashes, artifact ZIP hashes, the parser controls and machine-readable comparisons are retained in [actual-readback.json](actual-readback.json).
