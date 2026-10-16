# Exact tarball check receipt: pigauto 0.11.0

**State:** locally checked, hash-bound candidate. The CRAN release ledger remains `NOT_READY`.

## Source and artifact identity

- Package: `pigauto` 0.11.0.
- Source commit: `0b0f71fee838c6ed51ef832ed819270eeafaf29b` on merged `main`.
- A fresh detached worktree at that commit had empty `git status --porcelain=v1` before building.
- Build command: `R CMD build /Users/z3437171/.codex/worktrees/cran-011-frozen-source/pigauto`.
- Hash-qualified local copy: `/private/tmp/pigauto-release/sha256-1bfed5ad7be61437f4fdb09ece053d6a40211b5a5b7da4b2c947c3343493b719/pigauto_0.11.0.tar.gz`.
- SHA-256: `1bfed5ad7be61437f4fdb09ece053d6a40211b5a5b7da4b2c947c3343493b719`.
- Size: 5,128,791 bytes; tar inventory: 248 entries.
- The frozen copy and checked build output have the same SHA-256. The frozen copy is read-only.
- Complete archive-order inventory is in `exact-artifact-0b0f71f-2026-10-08.inventory.txt`, SHA-256 `52192b529b009a47b9bc3e7f1687f543194994de42e9e508ee536ebbacc0deb5`.
- The scan found no `.git`, `BACE`, `script`, `dev`, `data-raw`, checkpoint, `docs`, `.unlazy`, `cran-comments.md`, `.Rhistory`, `.RData`, log, output, or nested tarball paths. The `.rds` files present are shipped data, test fixtures, or generated vignette metadata and were not treated as forbidden.

## Exact local check

Command: `env -u NOT_CRAN R CMD check --as-cran --no-manual pigauto_0.11.0.tar.gz` from the directory containing this exact hash. `NOT_CRAN` was unset; network access was enabled for incoming feasibility. R 4.6.0, macOS Tahoe 26.7, `aarch64-apple-darwin23`.

Result: exit 0, `Status: OK`. CRAN incoming feasibility passed. The check log is retained as `exact-artifact-0b0f71f-2026-10-08-R-CMD-check.log`, SHA-256 `4798992b15bcc845b7d99ca7b331d17cffcf72342882da9f4cda5a7b685568d1`. It reports tests at `[300s/391s]` and vignette rebuilding at `[15s/21s]`. `--as-cran` ran the `--run-donttest` examples successfully.

The exact `testthat.Rout` is retained as `exact-artifact-0b0f71f-2026-10-08-testthat.Rout`, SHA-256 `db050f787770be1588b532d7bd53285f80d00cd2ea02ae047c385c5df0a72c20`. It reports 3,150 passes, zero failures, 161 test warnings, and 86 skips. These test-level warnings and skips coexist with the check's `Status: OK` and remain visible here.

The tarball's `DESCRIPTION` omits `drmTMB` and `gllvmTMB` from both `Imports` and `Suggests`. This verifies package dependency independence for this artifact; optional-backend integration remains a separate bounded claim.

## Gates that remain open

- No Windows R-release or R-devel result is bound to this SHA-256 yet. Results for earlier hashes remain predecessor evidence.
- The package contains BirdTree-derived data. Citation guidance is documented, but the evidence for redistribution rights is unresolved; this receipt does not resolve that question.
- No Grace, Rose, or Pat verdict has been issued on this exact hash. Earlier panel verdicts do not transfer.
- The live Pages workflow succeeded on merged `main` at `0b0f71fee838c6ed51ef832ed819270eeafaf29b`, and Chrome confirmed the deployed getting-started page reflects the merged defaults. Full retained/retired route and sitemap verification, local candidate visual review, and all other release-site checks remain open.
- The release ledger remains `NOT_READY`. No CRAN upload or submission has occurred.

This receipt binds one locally checked tarball. It does not establish platform-clean or submission-ready status.
