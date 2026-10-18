# Successor actual warning-object audit

Verdict: **PASS** for warning identity and consequence classification across all four actual cells.

The audit is bound to source `b2e6054685e6bdedb073ba46e058af2c0a619a18` and archive `a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb`, direct run 39091958419. The frozen warning contract was copied byte-for-byte from the prior adjudicated audit. Strict auditor self-tests passed; local macOS, Windows release, Windows devel, and Ubuntu release all passed.

Native and replay counts matched: 0 failures, 161 warnings, 86 skips, and 3,195 passes. Warning category counts are in `actual-successor-a93fece7/report.md`. Native check status was `Status: OK`, exit 0.

Two successful warning captures are retained, with coherent supervisor/diagnostic pairs of 171.199674/168.258 seconds and 167.394193/164.807 seconds. The prior derived receipt incorrectly combined the second supervisor time with the first diagnostic time. The audit's local replay JSON is hash-bound to capture A's diagnostic record; the raw RDS bytes are identical in the surviving output of capture B, so unique RDS process attribution cannot be distinguished. The corrected provenance is documented in `timing-provenance-amendment/amendment.md`; original derived files are retained there before amendment.

The PASS covers warning identities and bounded consequences only. It does not establish general inference quality, coverage, latent-parameter recovery, GPU capability, timing approval, CRAN upload or CRAN acceptance.
