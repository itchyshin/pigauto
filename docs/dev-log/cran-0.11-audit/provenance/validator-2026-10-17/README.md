# Validator fix and public-ledger rerun, 2026-10-17

The portable path fix is landed in the local-only brain repository at `1b5598eb05971ccde3020dd07295891b9b3caee8`. These source snapshots are byte-identical to that commit. The fix was independently reviewed and the unchanged public ledger from PR head `3bba584c53dc2ff93f2dbed2812faa8fe64b4515` now returns `READY FOR CLAIMED RUNG` at its existing `tarball-clean` claim. The metadata-only ledger update in this evidence commit was also rerun successfully.

This resolves the relative-path validator defect. The overall audit remains **NOT READY**. Platform-clean, deployed-site verification, final post-merge artifact checks, and independent release closure remain open.

## Source and execution identity

- Validator SHA-256: `81a897cc6c3351e36d082099663aa3255f780d698a377c35f201c1b41768572f`.
- Regression suite SHA-256: `12542f73787a9f25d6eee0dd22145977cd22d004eb5dd996a77dc9894736acab`.
- Unchanged input ledger SHA-256: `7b0f512905b86d087fab09b422377256e1464be63a93e5816668beffa856811b`.
- Frozen archive SHA-256: `c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2`; 5,136,713 bytes.
- Verification log SHA-256: `fcfbde9ecfd80fc085e5ee5be9cf65b9827c5ca2105ad12715a6c969e34c950b`.

## Review and regression evidence

An independent Codex reviewer first returned NOT READY for a regression: direct `validate()` calls could accept relative artifact paths from the working directory without a root. The corrected code rejects that case. The reviewer reran all eight regression tests and the built-in selftest, reproduced the negative control against the buggy version, and returned PASS with no remaining actionable findings, bound to the two source hashes above. Review was read-only.

The original committed validator rejected the public ledger's relative file evidence. The buggy intermediate patch failed exactly one of the new regression tests. The landed fix passes all eight tests, the built-in selftest, and the public-ledger rerun. Coverage includes explicit roots, automatic discovery through .git directories, worktree .git files and DESCRIPTION, ledger-directory fallback, absolute-path compatibility, missing files, wrong roots, directories, and bad archive hashes/sizes.

## Reproduction

From a checkout of this evidence branch, run:

```sh
python3 docs/dev-log/cran-0.11-audit/provenance/validator-2026-10-17/test_cran_release_gate.py -v
python3 docs/dev-log/cran-0.11-audit/provenance/validator-2026-10-17/cran_release_gate.py --selftest
python3 docs/dev-log/cran-0.11-audit/provenance/validator-2026-10-17/cran_release_gate.py docs/dev-log/cran-0.11-audit/release-ledger.json --root "$PWD"
```

The ledger command requires the frozen archive at `artifact.path`. On another machine, provide that exact archive and change only its local path in a copy of the ledger. The validator independently checks its SHA-256 and byte size. Missing archives fail closed; the archive is not included in this receipt. Relative evidence paths resolve from the selected root. Absolute paths and evidence outside the repository remain supported; this fix does not impose repository containment.

The validator tests structural requirements and evidence-file presence. It does not independently validate every claim inside those files. Its success at tarball-clean does not promote the claim to platform-clean or submission-ready. No package source, installed bytes, scientific default, deployment, or submission changed in this continuation.

See [verification.log](verification.log), [updated-ledger.log](updated-ledger.log), and [after-task.md](after-task.md).
