# After-task: defaults and MI source-equivalence check

## 1. Goal

Check whether the package source revision used for the defaults and focused multiple-imputation audits matches the merged source used to build the frozen 0.11.0 tarball.

## 2. Implemented

Compared `R`, `tests`, `DESCRIPTION`, `NAMESPACE`, and `man` between tested source `cf88d78e37d08130bd33823e1fd4594f3afa2d4f` and merged source `d76804e768bf77f430f43dcf591633f9cd900dab`. Git reported no differences in those paths. Recorded the exact comparison and bounded result in `provenance/source-equivalence-d768-cf88-2026-10-09.log` and added the interpretation to the release ledger.

## 3a. Decisions and Rejected Alternatives

Use the source-equality result to connect the defaults and MI source-mode evidence to merged package code. Do not treat it as a substitute for installation tests on the frozen archive or platform checks.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/provenance/source-equivalence-d768-cf88-2026-10-09.log`
- This after-task report.
- No package code, tests, dependency fields, namespace, or generated help changed.

## 5. Checks Run

- `git diff --quiet d76804e768bf77f430f43dcf591633f9cd900dab cf88d78e37d08130bd33823e1fd4594f3afa2d4f -- R tests DESCRIPTION NAMESPACE man`: exit 0, no differences.
- The defaults receipt remains 198 passed, 0 failed, 0 warnings, 0 skipped at the tested source revision. The focused MI receipt remains 234 passed, 0 failed, 9 expected small-sample warnings, 0 skipped. Their raw logs and hashes are listed in G1.
- This was a read-only revision comparison. No test suite was rerun.

## 6. Tests of the Tests

No test code was authored or changed. The check compares all package source and test paths relevant to the prior receipts, plus dependency metadata, namespace, and generated help.

## 7a. Issue Ledger

- Resolved: source-level defaults and focused MI receipts are applicable to the package code and tests in the merged source used for the frozen artifact.
- Open: installed exact-tarball checks and platform results remain artifact-specific and are tracked under G8.

## 8. Consistency Audit

The compared commits are not in a parent-child relationship, so the result is stated as a scoped path comparison rather than ancestry. `R`, `tests`, `DESCRIPTION`, `NAMESPACE`, and `man` are identical. The audit scripts and evidence files under `script/` and `docs/` were outside this comparison and are not used to assert package-code identity.

## 9. What Did Not Go Smoothly

No comparison errors occurred. The commits have different ancestry, which makes the explicit path scope important.

## 10. Known Residuals

This does not establish that every declared default combination was executed, that defaults are scientifically optimal, or that the frozen archive passes Windows checks. G8 and G9 remain open.

## 11. Team Learning

When test receipts and release artifacts name different revisions, compare the exact shipped package and test paths directly before either rerunning a suite or treating the receipts as equivalent.

## 12. Cross-Product Coverage

This comparison covers package source, tests, dependency metadata, namespace, and generated help. It does NOT cover the audit scripts, a fresh package installation from the tarball, Windows or other platform behavior, rendered-site output, statistical recovery beyond the existing regimes, or CRAN acceptance.
