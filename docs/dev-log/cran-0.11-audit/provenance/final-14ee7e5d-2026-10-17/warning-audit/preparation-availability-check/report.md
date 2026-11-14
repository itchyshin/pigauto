# Final-archive warning audit

NOT_DONE

Archive `14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb`; direct run 39082452558; 0/4 cells passed.

## local-macos
NOT_AVAILABLE


## windows-latest-release
NOT_AVAILABLE


## windows-latest-devel
NOT_AVAILABLE


## ubuntu-latest-release
NOT_AVAILABLE


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
