# Exact-artifact multiple-imputation backend matrix, 2026-10-09

**Result: PASS for the bounded optional-backend integration checks below.**

The tested archive was `pigauto_0.11.0.tar.gz`, SHA-256 `f8fee9f631460a1a49fa7ef78b4256e0861ea1929f25cba5093d0ea034f5a9c4`, size 5,135,324 bytes, 250 archive entries, built from source commit `d76804e768bf77f430f43dcf591633f9cd900dab`. The archive hash was checked before the matrix and again after it. The package was installed from that archive. Each isolated test library loads pigauto 0.11.0 from its own `pigauto` directory; a fresh `Rscript --vanilla` process per cell confirmed the optional package visibility shown below.

The matrix runs `script/cran-0.11-integration/check-adapters.R` in installed mode, with the real workflow in `test-real-backends.R` and the fixture/provenance tests in `tests/testthat/`. The test script runs the `multi_impute_analysis()` → `with_imputations()` → model fitting → `pool_mi()` workflow, automatic fixed-effect adapters, independent arithmetic checks, and saved-fit reloads.

| Installed cell | Verified backend visibility | Result |
|---|---|---|
| Neither | `drmTMB=FALSE`, `gllvmTMB=FALSE` | 50 adapter-fixture expectations and 55 provenance expectations passed. Saved optional objects report the missing-package condition; the cell completed with `NEITHER_BACKEND_OK`. |
| drmTMB only | `drmTMB=TRUE`, `gllvmTMB=FALSE` | Three fits had convergence 0 and positive-definite Hessians. Gaussian and Rubin arithmetic oracles passed; 28 expectations passed. |
| gllvmTMB only | `drmTMB=FALSE`, `gllvmTMB=TRUE` | Three fits had convergence 0 and positive-definite Hessians. Rubin arithmetic passed for eight terms; 22 expectations passed. |
| Both | `drmTMB=TRUE`, `gllvmTMB=TRUE` | Three fits per backend had convergence 0 and positive-definite Hessians. Gaussian and Rubin arithmetic passed for drmTMB and Rubin arithmetic for eight gllvmTMB terms; saved-fit reloads passed; 50 expectations passed. |

The run-level manifest and original cell logs are also retained. Manifest SHA-256: `fb9af9cbd7337358497ac8c72ecfcb8256f9504d01f8f696fac8ccdd39fb062e`. The manifest records archive hashes before and after every cell, isolated library paths, the installed-mode test script path and hash, and each original run-log hash. The test script hash is `d2b0b5f01e3470ce1f0f0d2116cadbd1a9a2beffcd864a5644e71c0b7d53fb7d` and matches `script/cran-0.11-integration/check-adapters.R` at evidence commit `2c212ab`.

Original run logs and hashes:

- `exact-mi-original-neither-2026-10-09.log`: `4bb5735cd580de0410b16ec3268b1f52cb756393f858cafe66aee8c264985d0f`.
- `exact-mi-original-drmTMB-only-2026-10-09.log`: `fc0385d9aa199457fe021079acfb37d582db7762de2a71e69e9f9d86738c8e1c`.
- `exact-mi-original-gllvmTMB-only-2026-10-09.log`: `67c1e5126ba261fa75e682da9ab7f3b7e8e5f7c2b00dd202550bb81951a982a7`.
- `exact-mi-original-both-2026-10-09.log`: `45dc91f0acb44796fbb92054ad797cf24d191d8a4fb35b914b336f34b4ac58b9`.

The retained logs and SHA-256 values are:

- [`exact-mi-both-backends-2026-10-09.log`](exact-mi-both-backends-2026-10-09.log): `dcba11ad9f11fe5740af2f1301bc894fb7196ac2082e2483cb189cb8612cbac9`.
- [`exact-mi-neither-fixtures-2026-10-09.log`](exact-mi-neither-fixtures-2026-10-09.log): `cb398da193bd8c6eb8c8ecf46cdf7b115e328991a7f5aad39c66efcc529ecb4e`.
- [`exact-mi-neither-reload-2026-10-09.log`](exact-mi-neither-reload-2026-10-09.log): `334e8ca94f3b17898c1b414ca275907ced1e7b85330b272464073a5646634615`.
- [`exact-mi-drmTMB-only-2026-10-09.log`](exact-mi-drmTMB-only-2026-10-09.log): `4d4cca3f2a713b1fb8ddc631347aab04484f9dc655d6335cb65e4ccab0b4ba6b`.
- [`exact-mi-gllvmTMB-only-2026-10-09.log`](exact-mi-gllvmTMB-only-2026-10-09.log): `e32344192ae01b7d8a0a6085b8271b180421cde77e6906023e2b4077f1c3729e`.

This establishes installation independence and bounded interoperability of the optional fixed-effect adapters for this archive. It does not validate every downstream model family or pooled estimand, establish Windows compatibility, or replace a full platform check. Windows logs currently report failures but are not checksum-bound to this archive. G8 and G9 therefore remain open, and the release verdict remains NOT READY.
