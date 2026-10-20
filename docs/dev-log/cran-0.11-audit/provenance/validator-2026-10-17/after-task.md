> Historical snapshot of the brain-fix report at `docs/dev-log/cran-0.11-audit/validator-fix-2026-10-17.md` in local brain commit `1b5598eb05971ccde3020dd07295891b9b3caee8`. Paths and staging statements below describe that fix at its original checkpoint, not the current pigauto checkout.

# Portable CRAN ledger validator fix, 2026-10-17

## 1. Goal

Land the validator path fix, obtain independent review, and rerun pigauto PR #236's public ledger at its existing tarball-clean claim.

## 2. Implemented

Relative evidence and artifact paths resolve from --root, or from the nearest .git/DESCRIPTION ancestor of the ledger. Without a marker, the CLI uses the ledger's directory. Direct API calls require a root for relative paths. Absolute paths retain their previous behavior.

## 3a. Decisions and Rejected Alternatives

Keep checksum, byte-size, evidence-presence, and release-rung checks intact. Do not convert ledger evidence to machine-specific absolute paths. Assumptions made without asking: the existing ancestor markers define the repository root; this is tested for ordinary repositories, worktrees, and R-package snapshots.

## 4. Files Touched

- tools/cran_release_gate.py
- tools/test_cran_release_gate.py
- docs/dev-log/cran-0.11-audit/validator-fix-2026-10-17.md (this report)

## 5. Checks Run

The corrected validator's built-in selftest passed. The new offline unittest suite passed all eight tests. The original HEAD validator rejected the unchanged public ledger's repository-relative evidence paths. The corrected candidate accepted that same ledger and evidence bundle with an explicit repository root, returning READY FOR CLAIMED RUNG for tarball-clean.

Validator SHA-256: 81a897cc6c3351e36d082099663aa3255f780d698a377c35f201c1b41768572f.
Test SHA-256: 12542f73787a9f25d6eee0dd22145977cd22d004eb5dd996a77dc9894736acab.

## 6. Tests of the Tests

Against the initial buggy patch, the new suite failed exactly one regression: a relative artifact was accepted without a root when the working directory held the file. The corrected patch rejects it from both tested working directories. Missing relative evidence, wrong roots, bad hashes and sizes, and directories also reject.

## 7a. Issue Ledger

Fixed: public relative evidence rejected solely for lacking absolute paths.
Fixed during independent review: rootless relative artifact acceptance depended on the working directory.
Independent Codex reviewer returned PASS after reproducing the regression, reviewing the full corrected delta, rerunning eight tests and the selftest, and checking negative controls. The review was read-only and bound to the hashes above.

## 8. Consistency Audit

Checked every verify_evidence caller, inventory evidence, artifact path, and later submission evidence. CLI discovery covers .git directories, worktree .git files, DESCRIPTION, fallback, and explicit roots. Absolute paths remain supported. No package source or scientific default changed.

## 9. What Did Not Go Smoothly

Default Git prompted for an Xcode licence; the command-line tools Git worked. The initial sandboxed lease attempt could not write the registry. A scoped elevated attempt succeeded and its lease file was read back. Only the authorized validator paths are staged.

## 10. Known Residuals

This closes the validator path defect. It does not establish platform-clean, deployed-site parity, post-merge artifact identity, or submission readiness. The validator checks file presence and the selected artifact's hash and size; it does not independently adjudicate scientific claims inside evidence files. Evidence paths are not confined to the repository, consistent with existing absolute-path support.

## 11. Team Learning

Test the direct API as well as the CLI, and vary the working directory. Routing receipt: pigauto LOAD-FIRST, the CRAN release protocol, unlazy, and verification-before-completion shaped the bounded claim and unchanged-artifact check. Golden-Set statement: the built-in CRAN partial-green negative controls ran; the broader memory regression suite was not run. No brain-memory files were edited.

## 12. Cross-Product Coverage

Covers relative/absolute file evidence, artifact hash/size checks, CLI root discovery, and direct API calls. It does NOT cover independent validation of evidence contents, external platform execution, publication, merging, deployment, or CRAN submission.

Assessment: self-review of this engineering report, 2/10 perceived AI-like style, medium confidence. Sections 5 and 10 give measured outcomes and bounded claims; no prose repair was needed. Science gate: not applicable, no new statistical claim. Facts gate: pass within the recorded tests and review. References gate: not applicable, no literature citations. Examples and test outputs were already seen during implementation.
