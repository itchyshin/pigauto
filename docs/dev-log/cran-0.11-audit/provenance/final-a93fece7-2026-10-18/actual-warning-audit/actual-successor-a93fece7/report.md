# Final-archive warning audit

PASS warning identities and consequences only

Archive `a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb`; direct run 39091958419; 4/4 warning audits passed.

## local-macos
PASS

{"status": "OK", "adjudication": "native_Status_OK", "exit": 0}
{"fail": 0, "warnings": 161, "skips": 86, "pass": 3195}
{"small_validation": 119, "undefined_constant_correlation": 9, "constant_matched_trait": 4, "few_tree_draws": 5, "covariates_joint_baseline": 5, "covariates_gnn_off": 3, "imbalanced_ordinal": 12, "signal_gate_fallback": 1, "conformal_sample_limit": 2, "species_dropped": 1}

## windows-latest-release
PASS

{"status": "OK", "adjudication": "native_Status_OK", "exit": 0}
{"fail": 0, "warnings": 161, "skips": 86, "pass": 3195}
{"small_validation": 119, "undefined_constant_correlation": 9, "constant_matched_trait": 4, "few_tree_draws": 5, "covariates_joint_baseline": 5, "covariates_gnn_off": 3, "imbalanced_ordinal": 12, "signal_gate_fallback": 1, "conformal_sample_limit": 2, "species_dropped": 1}

## windows-latest-devel
PASS

{"status": "OK", "adjudication": "native_Status_OK", "exit": 0}
{"fail": 0, "warnings": 161, "skips": 86, "pass": 3195}
{"small_validation": 119, "undefined_constant_correlation": 9, "constant_matched_trait": 4, "few_tree_draws": 5, "covariates_joint_baseline": 5, "covariates_gnn_off": 3, "imbalanced_ordinal": 12, "signal_gate_fallback": 1, "conformal_sample_limit": 2, "species_dropped": 1}

## ubuntu-latest-release
PASS

{"status": "OK", "adjudication": "native_Status_OK", "exit": 0}
{"fail": 0, "warnings": 161, "skips": 86, "pass": 3195}
{"small_validation": 119, "undefined_constant_correlation": 9, "constant_matched_trait": 4, "few_tree_draws": 5, "covariates_joint_baseline": 5, "covariates_gnn_off": 3, "imbalanced_ordinal": 12, "signal_gate_fallback": 1, "conformal_sample_limit": 2, "species_dropped": 1}

## Consequences and limits
- small_validation: Small calibration samples leave gates/scores noisy; no coverage certification.
- conformal_sample_limit: The stated finite-sample size cannot attain 95% split-conformal coverage.
- constant_matched_trait: The deliberately constant matched trait does not identify phylogenetic signal.
- few_tree_draws: Few tree draws can make tree-sensitivity summaries unstable.
- covariates_joint_baseline: Covariates do not contribute through the selected joint-baseline route.
- covariates_gnn_off: Covariates do not contribute when the GNN is disabled.
- imbalanced_ordinal: Ordinal class imbalance limits the tested prediction/calibration regime.
- signal_gate_fallback: The opt-in signal gate selected its documented BM fallback.
- species_dropped: The three deliberately unmatched species are dropped and not modeled.
- undefined_constant_correlation: Correlation is undefined for the constant values in this fixture.

Separate actual warning-object review. Not a final release-panel READY vote, coverage certification, timing approval or submission authorization.

## Capture provenance amendment

Two successful diagnostic captures are retained: A supervisor/diagnostic 171.199674/168.258 s; B 167.394193/164.807 s. The earlier derived verdict combined B supervisor time with A diagnostic time; that combination is invalid. Exact replay JSON used by this audit is capture A's derived record. The warning RDS used by the audit is SHA-identical to both successful output copies, so unique run attribution for that byte-identical RDS is unresolved. No replay was run for this amendment. See `/private/tmp/pigauto-final-closure-20261017/successor-b2e6054/actual-warning-audit/timing-provenance-amendment/amendment.md`.
