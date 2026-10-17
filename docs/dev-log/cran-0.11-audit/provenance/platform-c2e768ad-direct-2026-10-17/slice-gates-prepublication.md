# Exact frozen-artifact platform closure

OWNS: .github/workflows/frozen-tarball-check.yaml, script/cran-0.11-platform/, docs/dev-log/cran-0.11-audit/

- [x] P0: Same frozen artifact is locally intact and public object readback matches.
  EVIDENCE: c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2; 5136713 bytes; inventory 250; independent identity and downloaded object both match.
- [x] P1: Checker has positive and mutation controls and independent pre-launch review.
  EVIDENCE: Positive and five negative identity controls passed; Python/R/YAML parsed; stdout/stderr streaming canary passed; independent platform reviewer approved bounded design.
- [x] P2: Latest applicable source CI passes.
  EVIDENCE: SourceCI35939070922883 completed success at a5500712,all3jobs. Final evidence tip CI is checked as part of P6 readback.
- [x] P3: Direct archive checks complete on Windows release/devel and Linux release.
  CHECK: PYTHONDONTWRITEBYTECODE=1 python3 /private/tmp/pigauto-pr236-validator-20261017/script/cran-0.11-platform/verify_results.py /private/tmp/pigauto-platform-c2e768ad-20261017/results
  EXPECT: DIRECT_ARCHIVE_MATRIX_VERIFIED_3_OF_3
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/private/tmp/pigauto-platform-c2e768ad-20261017; path=e0fb875d37d1/39 entries; output=] | DIRECT_ARCHIVE_MATRIX_VERIFIED_3_OF_3
- [x] P4: Independent completed-log review approves the platform claim with timings and skips adjudicated.
  EVIDENCE: Fresh read-only reviewer PASS on all3 rawStatusOK logs,all29dependencies,96 retained file hashes,actualindexedmanuals,provenance;161warningsunclassified/86skips retained;599/656secWindows time keepssubmissiontiminggateopen.
- [x] P5: Release ledger at platform-clean passes the landed validator.
  CHECK: python3 /private/tmp/pigauto-pr236-validator-20261017/docs/dev-log/cran-0.11-audit/provenance/validator-2026-10-17/cran_release_gate.py /private/tmp/pigauto-pr236-validator-20261017/docs/dev-log/cran-0.11-audit/release-ledger.json --root /private/tmp/pigauto-pr236-validator-20261017
  EXPECT: READY FOR CLAIMED RUNG
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/private/tmp/pigauto-platform-c2e768ad-20261017; path=e0fb875d37d1/39 entries; output=READY FOR CLAIMED RUNG
- [ ] P6: Reviewed receipts are published in the existing unmerged PR and externally read back.
