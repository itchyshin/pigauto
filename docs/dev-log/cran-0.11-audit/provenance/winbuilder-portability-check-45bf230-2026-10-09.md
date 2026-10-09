# Source-package check of the Windows portability patch

This check used a temporary archive of source commit `45bf230c7a10f56d0357b7fd5f5403a263eb72d4` from the separate portability branch. It did not alter that branch or the frozen release tarball.

- Environment: macOS arm64, R 4.6.0, 2026-10-09.
- Build: `R CMD build --no-build-vignettes .` from the extracted source snapshot. Build log: `winbuilder-portability-build-45bf230-2026-10-09.log` (SHA-256 `ac4a4636459c529341daf67858aae891358e667377b552014d479a15d197994b`).
- Built source archive: `pigauto_0.11.0.tar.gz`, SHA-256 `58e4b5628f0f58339fd49e50f5aa8dabed9ed8e635b7b0564d310a5fd185d405`, 5,052,462 bytes. This is a temporary branch-source archive, not the frozen release artifact.
- `R CMD check --as-cran --no-manual --run-donttest` stopped at incoming feasibility because CRAN and Bioconductor package indexes could not resolve in this environment. It did not install or test the package. Raw log: `winbuilder-portability-as-cran-network-45bf230-2026-10-09.log` (SHA-256 `4767d92b9b45d2339b2a9f7eab62c1aeb0bbb00d5f4254be9d9fce253b04edba`).
- Offline-capable package check: `_R_CHECK_FORCE_SUGGESTS_=true R CMD check --no-manual --run-donttest pigauto_0.11.0.tar.gz`. Installation and all runnable tests passed. Testthat: 3,190 passes, 0 failures, 161 warnings, 86 skips. The overall check ended with 2 WARNINGs because `--no-build-vignettes` left `inst/doc` absent; this is not a clean CRAN check. Raw log SHA-256: `ff78131befefea08b35112e85cc3d094ef8d5204c58cb09e620534b3b82d940c`. Test output SHA-256: `cecb5f6af6b5f2f4e780064027a4f7da60b00fa6c825f7bf2b6247d94538c811`.

This broadens the earlier three-file macOS smoke to the full runnable package test suite for the patch source snapshot. It does not establish Windows or libtorch-enabled behavior, does not bind the Win-builder failures to a source hash, does not validate the frozen release tarball, and does not close G8.

## Preserved logs

- `winbuilder-portability-build-45bf230-2026-10-09.log.gz` stores exact raw bytes; uncompressed SHA-256 `ac4a4636459c529341daf67858aae891358e667377b552014d479a15d197994b`, compressed SHA-256 `87847ab7d7062f0bc4572c62c6f28244d35c74bfeeb2834571d92d88bcf259f4`.
- `winbuilder-portability-as-cran-network-45bf230-2026-10-09.log.gz` stores exact raw bytes; uncompressed SHA-256 `4767d92b9b45d2339b2a9f7eab62c1aeb0bbb00d5f4254be9d9fce253b04edba`, compressed SHA-256 `c172d4d0025a6a649c8d3efddb787e2c76fabb9848aa20bda726b68ec151961c`.
- `winbuilder-portability-package-check-45bf230-2026-10-09.log.gz` stores exact raw bytes; uncompressed SHA-256 `ff78131befefea08b35112e85cc3d094ef8d5204c58cb09e620534b3b82d940c`, compressed SHA-256 `4d9a1c14cd98234077d8303c096b6320da750c8f2371f2872155643c1cb2c8dd`.
- `winbuilder-portability-testthat-45bf230-2026-10-09.Rout.gz` stores exact raw bytes; uncompressed SHA-256 `cecb5f6af6b5f2f4e780064027a4f7da60b00fa6c825f7bf2b6247d94538c811`, compressed SHA-256 `4d2372879cf42971c526063dd066c99ad745e17a63e2362aedf65d2b37dc3226`.
