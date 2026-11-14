# Current artifact correction and public receipt map

Checked 2026-10-17 against frozen `pigauto_0.11.0.tar.gz` and PR #236 head `e2fd7f72a09a90a7f0c85315042e94b53c25b0e7`. This file supersedes inconsistent current-state fields in the embedded evidence bundle, `release-ledger.json`, and the opening checkpoint in `GATES.md`. Older dated entries remain historical evidence for their stated hashes.

## Authoritative artifact identity

- Source commit: `9a457eb6b3363d413dbc502de88ff5bc602ed98a`.
- Tarball SHA-256: `c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2`.
- Size: 5,136,713 bytes; 250 entries; sorted inventory SHA-256 `d7b9b0e59ea051b3a22d7f347d8d925e6265b076b238ac893410e79d5b61329d`.
- Exact macOS arm64/R 4.6.0 `R CMD check --as-cran --run-donttest --timings`: Status OK; 3,195 pass, 0 fail, 161 warnings, 86 skips. Check-log SHA-256: `93180d66bf6ee3ddeed26806ef1a998849d9759bc6eea5e2abda6292e4d48478`.

The bundle at `provenance/current-artifact-c2e768ad-evidence-bundle.md` contains the receipts embedded as text. Earlier ledger fields pointing to `/private/tmp` or to standalone files absent from the PR are provenance references from the producer checkout, not public links. Use the bundle and the hashes above as the public receipt source.

## Windows result scope

The R-hub uploaded object was read back byte-identical to the frozen tarball (object SHA-256 `c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2`; 5,136,713 bytes). R-hub extracted that object, committed the extracted source as `9386c07717bc7a2dc63f072fcfbea2345b6e1e15`, then rebuilt and checked that source on Windows R-release and R-devel. Both checks reported Status OK, 2,442 pass, 0 fail, 28 warnings, and 178 skips. This is a byte-verified input-object-to-source chain followed by fresh Windows source builds. The Windows jobs did not execute the original compressed tarball in place, and R-hub did not attest a server-side tarball digest. The previous Boolean `windows_results_bound_to_hash: true` overstated direct archive execution; strict archive-in-place binding is false.

The older Win-builder URLs under `win_builder_logs` are the separate 5,136,083-byte `fe1600d7` attempt and are not Windows results for `c2e768ad`. Keep them as predecessor evidence only. The current artifact's Windows evidence is the R-hub chain above.

## Optional MI integration

The current installed optional-backend matrix is the `c2e768ad` matrix: four availability cells (neither, drmTMB only, gllvmTMB only, both), 12 successful steps, independent arithmetic and serialization/provenance checks. Its manifest hash is `bc29ddc1ba91dd7619ac44e2553c19588bd3cacceed1477fce68ffb84c94c525`. This supports bounded optional model-object interoperability and termwise fixed-effect Rubin pooling. The `fe1600d7` matrix fields in the older ledger are superseded. This does not show clean-room dependency installation, broad statistical recovery or interval coverage, or joint pooling of covariance, loadings, or variance components.

## Shipped tree data

The NOTICE and six shipped R data streams were rehashed from `c2e768ad`. Their bytes match the hashes listed in `rights_exact_artifact` and `current_shipped_data_hashes.shipped_rda_sha256`; the archive identity for those hashes must be `c2e768ad`, not `fe1600d7`. This confirms package contents only. The rights basis remains the maintainer's explicit authorization and warranty to use the trees with citation and credit; no separate BirdTree redistribution licence or independent legal opinion is claimed.

## Site and gate status

The local site build is from source `9a457eb6b3363d413dbc502de88ff5bc602ed98a`: 62 HTML pages, 614 search entries, pkgdown check passed, and the structural crawl found 3,588 references and 34 retired routes with zero errors. Browser review checked rendered content and accessibility-tree text for the regenerated example page. Pixel appearance was not assessed. This candidate has not been deployed, and the older live deployment does not prove candidate parity.

Using the approved plan's gate definitions, G0-G4 are met within their recorded scopes (5 met). G5 is not fully met because pixel appearance review is outstanding; G6 deployed-site verification, G7 post-merge artifact freeze, and G8 independent closure/release-validator reverification are also unmet (4 unmet). None are abandoned. Overall verdict: NOT READY. Highest proven rung: tarball-clean. Next unproven rung: platform-clean under the exact artifact gate, then the local-site visual, deployed-site, post-merge freeze, and closure gates. Current PR checks on `e2fd7f72` were still running at the time of this correction; the pkgdown PR workflow was skipped by its configured guard.

No merge, deployment, or CRAN submission occurred.
