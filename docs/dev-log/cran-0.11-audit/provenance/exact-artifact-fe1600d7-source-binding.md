# Exact artifact source and content binding: `fe1600d7d562`

Date: 2026-10-16 MDT / 2026-10-17 UTC
Source commit: `780c08470b263cc26b4139c711a5c144fd29389c`
Archive: `/private/tmp/pigauto-release-candidate/sha256-fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a/pigauto_0.11.0.tar.gz`
SHA-256: `fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a`
Size: 5,136,083 bytes
Tar entries: 250; regular files: 236
Sorted archive path inventory SHA-256: `d7b9b0e59ea051b3a22d7f347d8d925e6265b076b238ac893410e79d5b61329d`.

The tarball was rebuilt from the clean detached checkout `/private/tmp/pigauto-clean-780c084`; Command Line Tools Git reports HEAD `780c08470b263cc26b4139c711a5c144fd29389c` and an empty `git status --short`. A raw tar-member comparison against that checkout found 215 byte-identical regular package files, one content difference (`DESCRIPTION`, generated/normalized during `R CMD build`), and 20 build-generated files absent from the source checkout. This binds the new archive to the same source commit. It does **not** prove byte identity with the earlier archive.

The prior recorded archive SHA-256 was `0fccf4bb0a179434fe15fa3a7ff514112bc5220e60c4b3e025fa00536f5ca29d`. The local rebuild accidentally replaced its conventional `/private/tmp/pigauto_0.11.0.tar.gz` path. A search found no surviving copy of that old hash. The old hash is therefore superseded and unavailable; its checks and reviews do not validate this replacement archive.

The new archive has the same sorted 250-entry path inventory as the old receipt. Direct streaming from the new tarball confirms the same SHA-256 values for `inst/NOTICE` and all six shipped `.rda` data objects as the old receipt. This is a narrow content comparison only; it does not transfer MI, Windows, or independent-review evidence.

## Exact-file package checks

- `R CMD check --as-cran --run-donttest --timings` with `NOT_CRAN` unset and `_R_CHECK_FORCE_SUGGESTS_=true` completed on this exact hash on macOS arm64 / R 4.6.0. CRAN incoming feasibility, installation, examples, tests, vignette rebuilds, and PDF/HTML manual checks completed. Final status was `1 NOTE`, solely because `compare_methods` took 6.142 seconds (the timing-note threshold is 5 seconds). See `exact-artifact-fe1600d7-as-cran-check.log` and `exact-artifact-fe1600d7-as-cran-testthat.Rout`.
- Testthat: 3,194 pass, 0 fail, 161 warnings, 86 skips; the as-cran run output is `exact-artifact-fe1600d7-as-cran-testthat.Rout`.
- Reproducible invocation: `env -u NOT_CRAN _R_CHECK_FORCE_SUGGESTS_=true R CMD check --as-cran --run-donttest --timings --output=/private/tmp/pigauto-ascran-fe1600d7-network /private/tmp/pigauto-release-candidate/sha256-fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a/pigauto_0.11.0.tar.gz`. The frozen file still hashes to `fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a`; an independent review confirmed all 236 regular archive files byte-match the check's `00_pkg_src` tree.
- Evidence hashes: `exact-artifact-fe1600d7-as-cran-check.log` SHA-256 `35524a843ebe50d4b85f1653767597960a7c5ea81ff238de651582a20459d2c0`; `exact-artifact-fe1600d7-as-cran-testthat.Rout` SHA-256 `f868df437cb7cc69d474baec4321878541c59fed1fb0433d45d279374fca93af`; generated 83-page manual PDF SHA-256 `7f35acb0001c21f22ea05c78d69f2b301b2715152a74cdcc70a464960f6f50e9`.
- An earlier ordinary check with `--no-manual` returned `Status: OK`; a still earlier attempt without `--no-manual` failed while generating the manual because `pdflatex` was unavailable. The later complete as-cran run demonstrates that the PDF manual check now completed successfully in the active environment.
- Testthat: 3,194 pass, 0 fail, 161 warnings, 86 skips; output is `exact-artifact-fe1600d7-testthat.Rout`.
- A preceding ordinary check without `--no-manual` completed package checks but ended with one manual-PDF warning and one error because `pdflatex` was not installed. It did not identify an Rd/package defect.
- `R CMD check --as-cran --run-donttest --timings` on this same hash stopped during CRAN incoming feasibility because this host could not resolve CRAN/Bioconductor hosts. It did not complete an as-cran check. Captured output: `exact-artifact-fe1600d7-as-cran-network-block.log`.

The complete as-cran run supports the exact-tarball check rung, with one timing NOTE to explain in release comments if this artifact is submitted. It does not close platform checks. The optional-backend MI matrix on this hash has now passed independently within its bounded scope; checksum-bound Windows results remain open. Grace independently rehashed the archive and check outputs, matched the full archive contents to the check input tree, and confirmed incoming feasibility, manual checks, and the sole timing NOTE.

## Optional-backend MI integration matrix on the replacement hash

The existing exact-artifact runner was rerun against `fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a` in a new output directory. All four availability cells (`neither`, `drm-only`, `gllvm-only`, `both`) completed install, fixture, and real-workflow steps: 12/12 exit codes were 0, and the archive hash was checked before and after every step. The manifest is `provenance/mi-matrix-fe1600d7-provenance-complete/manifest.json` (SHA-256 `e81c14f45d3d3837a25999f21d4f83f60276a918f532adfd5c9f2ddba0bb1991`). Summed step time was 217.7 seconds; all 12 logs are retained beside the manifest. The harness, runner, and five executed R sources are hashed in the manifest.

This is bounded installed-workflow interoperability evidence. Shared dependency symlinks mean it is not clean-room dependency installation. It covers termwise fixed-effect pooling and the recorded Gaussian/Rubin arithmetic checks; it does not establish broad recovery/coverage, joint covariance/loading/variance-component pooling, Windows portability, or release readiness. The independent MI review passed within this scope; see `independent-review-fe1600d7-2026-10-17.md`.


## Win-builder submission acknowledgement (2026-10-17 UTC)

The frozen archive was submitted without rebuilding to the official Win-builder R-release and R-devel upload forms. Both HTTP responses returned the exact filename `pigauto_0.11.0.tar.gz`, size 5,136,083 bytes, and octet-stream type. The R-devel page also showed “Currently unavailable.” These are upload acknowledgements only; they do not constitute Windows check results, and the service provides no server-side SHA-256. The local archive remained SHA-256 `fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a` after submission. Receipt: `winbuilder-upload-fe1600d7-2026-10-17.md`.

The seven checked rights/data paths were rehashed directly from this exact archive; individual SHA-256 values are recorded in `exact-artifact-fe1600d7-rights-data-hashes.txt`. This checks shipped bytes only and does not establish third-party redistribution terms.
