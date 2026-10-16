# Independent review of frozen pigauto 0.11.0 evidence

Date: 2026-10-16
Artifact: `pigauto_0.11.0.tar.gz`
SHA-256: `4bdf4de71650716b158824f545ad92c4d125fe1601adefab97eeb5f301c5debe`
Size: 5,135,816 bytes
Source commit: `8082ef8273177ff664d49ac4bd643f7def671f20`

## Verdicts

- **Grace — PASS, bounded exact-hash review.** Independently rehashed the archive and checked its size, version, absence of drmTMB/gllvmTMB dependency fields, source binding, local check and testthat log hashes, optional-backend MI manifest and individual logs, and Windows R-release receipt. Confirmed the MI claim is termwise fixed-effect interoperability only. Limits: Windows R-release log does not echo the archive digest; Windows R-devel and GitHub R-devel remain pending. Site evidence is a separate gate.
- **Pat — PASS, bounded exact-artifact review.** Independently hashed the archive and checked source identity, inventory, local check, testthat, Windows R-release, and twelve MI log hashes. Checked fresh site HTTP receipt/output hashes separately. Limits: Windows R-release association is traceable through unique result URL, matching filename/size, and local pre/post-upload digest, but is not a server-side digest assertion. R-devel remains pending; this does not establish platform-clean.
- **Rose — PASS, bounded adversarial reconciliation.** Verified archive identity, local as-cran check, MI integration, attribution reconciliation, receipt/output hashes, and the current site evidence after requested corrections. Limits: R-devel/platform completion and the final evidence PR remain open.

Preservation note: the fetched R-release log bytes had SHA-256 `056ea23d31efc6fcaebd06530ffbdf51c92c8ca6c6f4adb3450fb1c86c1d0340` and 3,128 bytes. They are retained as deterministic gzip at `provenance/winbuilder-release-4bdf4de7-2026-10-16.log.gz` (compressed SHA-256 `a770d0f1e26748fd7951ea8aff1004431b3157cde93376fb215f90bb66c358d6`, 1,034 bytes); decompression reproduces the fetched bytes exactly.

## Scope and release status

The three reviews support the exact local artifact, bounded MI integration, current site receipt, and R-release upload association within the stated limits. They do not establish Windows R-devel, GitHub R-devel, full platform-clean status, or CRAN acceptance. The release ledger therefore remains `NOT_READY` at `tarball-clean`; `platform-clean` is the next unproven rung. The final unmerged evidence PR is deferred until the remaining platform results are available. Nothing has been merged or submitted to CRAN.

## Draft assessment

**Style:** 4/10, medium confidence; the repeated “Limits” phrasing is formulaic, while the short parallel verdicts make scope easy to compare. No change needed for this evidence-record genre.

**Genre/coverage:** concise technical review receipt; the complete file.

**Evidence/repair:** artifact identity, check, MI, site, and Windows statements were checked against the cited ledger fields, receipts, log bytes, and independent reviewer messages. The three separate platform limitations are retained because they differ materially.

**Gates:** science = not applicable (no new method claim); facts = pass for the bounded statements above; references = not applicable (no external literature citations).

**Provenance:** self-review by Codex, this file version, 2026-10-16; no style-model examples or controls were used.
