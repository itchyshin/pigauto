# Nine-repair integration checkpoint

Frozen source `3015dbe18f242a1f8bd7011f6d40d40ed82961ee` combines the coverage tests with eight runtime repairs and the separate Bayesian cleanup in draft PR 294. It does not contain PR 295 or the later backend and CPU diagnostic tests. This is evidence for the integration branch; no coverage percentage is attributed to main.

The full instrumented run covers **12,825 of 13,025 executable locations (98.464491%)**, leaving **200 uncovered locations**. All 253 tracked locations in the four neural methods are covered. Collection took 17 minutes 14 seconds on one CPU thread, within its 30-minute cap. The runtime receipt records torch 0.17.0, CPU execution, and unavailable CUDA and MPS.

| Prior catalog disposition | Locations |
|---|---:|
| Unchanged mapped locations covered by named test sets | 12,777 |
| Unchanged mapped locations uncovered | 200 |
| Changed or removed locations pending separate repair approval | 25 |
| Total prior catalog | 13,002 |
| New current locations, all covered | 48 |
| Current denominator: 13,002 minus 25 plus 48 | 13,025 |

The original 1,904 misses reconcile as **1,696 covered by named test sets, 200 uncovered, and eight changed locations pending approval**. All original entries remain in the accounting. No previously covered unchanged location regressed. The two report-label regressions in the preceding f780 snapshot are now covered by the deterministic report test. The reachable forest-panel guard is also covered; its real-input test uses missing trait labels with na.action set to "na.pass".

Of the current covered locations, 12,813 have single named-test witnesses and 12 require collective witnesses. Their distinction is preserved in the compressed evidence. The raw counter RDS remains on Totoro and its byte hash is bound by the receipt, alongside the CSV, manifest, extractor and source identities.

The instrumented suite passed 5,123 assertions with 192 warnings and 33 skips. The exact-source package check passed with zero errors, warnings or notes, including rebuilt vignettes and excluding the PDF manual. Its installed suite passed 5,120 assertions with 192 warnings and 34 skips. Every condition is inventoried and classified. The installed suite's additional skip is the NEWS source test; the three-assertion difference matches that test, which has no skip in the instrumented run. The separate source-only gate passed 192 assertions without conditions. Missing optional fixtures and HTML Tidy remain disclosed.

Independent Sol review recomputed all artifact hashes and checked every prior catalog key, mapped Git line, counter value and new-location complement. It independently parsed both raw test logs and inspected the check, source and CPU receipts. The raw RDS bytes were retained remotely rather than downloaded for that review. G5 and G6 support this exact partial snapshot. The strict verifier still rejects full source identity for the nine modified runtime files; the separate conservative reconciliation never reports completion and leaves changed entries unresolved.

The remaining current locations are prediction/baselines 118, inputs/evaluation 28, MI 30 and GNN/fitting 24. The residual review CSV preserves provisional source-based classifications and their evidence limits. A failed search for a fixture does not establish that a branch is unreachable. Named tests for newer changes await measurement on candidate `0aad929c263419ab238df326951d1f2ab186651a`, which also includes the separately reviewed PR 295 cleanup. That run is active with the same time and CPU limits.

G7 remains unmet because 200 current locations and eight changed original misses remain unresolved. Final G8 acceptance is pending. Coverage does not establish calibration, GPU correctness, singleton resolution or inferential validity of optional-backend analysis fixtures. No PR merge, public issue update or CRAN release action has been taken.
