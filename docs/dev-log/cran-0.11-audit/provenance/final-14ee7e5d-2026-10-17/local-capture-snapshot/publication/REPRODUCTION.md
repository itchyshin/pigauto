# Reproducing or reading this local evidence package

The source evidence is read-only and remains at:

- Frozen archive: `/private/tmp/pigauto-release-candidate/sha256-14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb/pigauto_0.11.0.tar.gz`
- Native Mac check tree: `/private/tmp/pigauto-release-candidate/sha256-14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb/local-macos/`
- Frozen archive identity record: `/private/tmp/pigauto-final-closure-20261017/FINAL.json`
- Published warning-capture script: `/private/tmp/pigauto-pr236-validator-20261017/script/cran-0.11-platform/capture-warnings.R`

The packaged raw evidence is under `local-macos/`: the full `platform-evidence/`, each regular file directly under `frozen-check/pigauto.Rcheck/` (including manuals and check logs), `frozen-check/pigauto.Rcheck/tests/testthat.Rout`, plus the four named native/replay receipts and logs. Source and installed-package trees, fixture trees, symlinks, and the frozen archive bytes themselves are intentionally excluded; the frozen archive identity is SHA-bound in the records.

All copied inputs retain their original bytes. No CRLF normalization was done, so an `original.raw.gz` duplicate was unnecessary. `publication/member-manifest.json` records each archive member's size and SHA-256. `publication-index.json` binds that manifest and archive; `publication/producer-identity.json` binds the source archive, local receipts/logs, and producer tools.

Verify the existing archive and run both planted-failure controls with:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -B publication/readback.py pigauto-final-local-macos-evidence.tar.gz --selftest
```

To reproduce the archive from original evidence, copy only this `publication/` helper directory to a new empty parent directory; do not copy an existing tarball or `.staging-final-local` directory (the builder refuses to overwrite either). Keep the three absolute scratch inputs above available, then run from the copied `publication/` directory:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -B publication/build.py
```

The resulting tarball is written beside `publication/`. The builder copies regular files only, sorts member paths, fixes tar/gzip metadata, and checks every copied source hash. It does not run R, tests, warning capture, or the audit.

Native result stays `Status: 1 NOTE`, R CMD check exit 0, wrapper exit 1; the one examples timing NOTE remains UNADJUDICATED. The separate warning replay exited 0 (driver 687.196 s; R reporter 669.225 s). The warning audit remains NOT_DONE because the current reader guard rejects the replay receipt's empty `error` object; no warning adjudication or release readiness is claimed.
