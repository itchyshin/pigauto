# Final-archive warning audit

NOT_DONE

Archive `14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb`; direct run 39082452558; 0/4 warning audits passed.

## local-macos
NOT_DONE
Replay error or changed CRAN environment
{"status": "NOTE", "adjudication": "UNADJUDICATED", "exit": 0, "typed_note": {"kind": "exact_local_example_timing_note", "method": "compare_methods", "user_seconds": 2.532, "system_seconds": 0.992, "elapsed_seconds": 18.406, "raw_section": "* checking examples ... [13s/63s] NOTE\nExamples with CPU (user + system) or elapsed time > 5s\n                 user system elapsed\ncompare_methods 2.532  0.992  18.406"}, "scope": "Allowed only to classify warning objects. Timing/native release verdict remains UNADJUDICATED."}

## windows-latest-release
NOT_DONE
Replay error or changed CRAN environment
{"status": "OK", "adjudication": "native_Status_OK", "exit": 0}

## windows-latest-devel
NOT_DONE
Replay error or changed CRAN environment
{"status": "OK", "adjudication": "native_Status_OK", "exit": 0}

## ubuntu-latest-release
NOT_DONE
Replay error or changed CRAN environment
{"status": "OK", "adjudication": "native_Status_OK", "exit": 0}

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
