# Independent pilot evidence review before production

Reviewer was separate from the runner implementation. Read-only raw-record inspection; no fits run by reviewer.

Original pilot: eighteen exact job IDs; thirteen ok, two failed, three unavailable. Rows: 203 ok, two failed, 45 unavailable. Follow-up: seventeen exact IDs, all envelopes ok; rows: 86 ok, fifteen unavailable, zero failed. Unavailable metrics remain visible within successful jobs.

All twelve analysis cells used m=50 with 50/50 successful downstream fits and smoke=FALSE. Posterior replay used four chains, 5,000 iterations, 1,000 burn-in, thinning20, 1,000 retained draws and m=50. Sampler convergence was TRUE and downstream fits 50/50. Posterior x intervals retained the complete 83/83 masked-cell denominator; y had zero requested cells and is unavailable.

The MAR case retained sixty reliability-bin rows. ECE is descriptive fixed-bin calibration, not a coverage guarantee. AVONET and BIEN mask counts match originally observed-cell denominators. Raw source identities remain separate, with evaluator provenance recorded independently.

Scalar summaries are descriptive or indeterminate. Paired contrasts have at most one independent pair. No replicated calibration or scientific superiority conclusion is supported by this pilot.

This review covers only the pilot evidence and production estimate. Production approval, execution, and scientific outcomes are recorded in `review-production.md` and `scientific-report.md`.
