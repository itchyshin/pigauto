# Final archive local macOS evidence

This package preserves the regular-file evidence from the native macOS `R CMD check` and its separately run warning capture for frozen pigauto 0.11.0 archive `14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb` (5,137,648 bytes), built from merged source commit `8e275ca738129c8d608f06ffc407b8141c1a2636`.

The native check finished with `Status: 1 NOTE`; the sole NOTE is the examples timing record for `compare_methods` (user 2.532 s, system 0.992 s, elapsed 18.406 s). The R check exit was 0. Its wrapper exited 1 because the wrapper explicitly refuses check findings. The timing NOTE remains **UNADJUDICATED**. This package does not claim native `Status: OK`, release readiness, or CRAN submission readiness.

The warning capture was a separate diagnostic on the package installed by the completed native check, with `NOT_CRAN=UNSET` and `_R_CHECK_FORCE_SUGGESTS_=true`. It exited 0 after 687.196 s at the process driver; its R-side reporter records 669.225 s. Those diagnostic times are separate from the native check time. The replay counts match the native test summary (0 fail, 161 warnings, 86 skips, 3195 pass). The counts and raw warning objects are preserved for review; they do not adjudicate the warnings.

The four-cell warning audit did not pass. It returned `NOT_DONE` because its environment guard required `error` to be null while all replay receipts encode an empty object (`{}`). The independent self-test passed. No checker or contract change is included here.

`local-macos/platform-evidence/` contains the complete evidence directory. `local-macos/frozen-check/pigauto.Rcheck/` contains each regular file directly at the Rcheck root, including the manual PDF and check logs, plus only `tests/testthat.Rout`. The installed package tree, source tree, fixtures, symlinks, and frozen archive gzip are excluded. All included source files are copied as original bytes; line endings are not normalized.

The internal member manifest records each included file's byte length, SHA-256, and source path where applicable. The outer `publication-index.json` binds the deterministic tarball, manifest, producer identities, validation receipts, and native/replay receipts. Run `python3 readback.py <archive> --selftest` for the positive roundtrip and planted hash/path rejection checks.
