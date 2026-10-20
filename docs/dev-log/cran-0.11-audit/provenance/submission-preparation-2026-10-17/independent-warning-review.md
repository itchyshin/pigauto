# Independent local-warning review

The platform reviewer inspected the corrected local c2 replay, warning RDS, details CSV, classifier and unknown-warning control. Verdict: PASS for local classification, with no source blocker established by these 161 warnings. All 161 messages matched the CSV exactly, IDs were complete and unique, and callers supported every category.

The 121 calibration/sample-size warnings arise in calibrate_gates and retain their stated uncertainty limits. Nine stats::cor warnings identify undefined constant correlations; four other warnings concern deliberately constant matched traits. Eight warnings disclose ignored covariates under the selected joint-baseline or GNN-disabled route. Twelve disclose ordinal imbalance, five disclose too few tree-sensitivity draws, one identifies the tested opt-in signal-gate BM fallback, and one identifies three deliberately unmatched species.

This is not a blanket benign-warning verdict. The review does NOT cover Windows warning identities, a later archive, broader calibration/coverage, final timing acceptance, or final release-panel approval. Matching platform totals do not prove matching messages. Preserve the failed installed-layout diagnostic separately and repeat relevant capture on the final artifact.

Reviewer: platform_gate_review, native read-only independent review, 2026-10-17. This was a preparation review, not the fresh final Grace/Rose/Pat release panel.
