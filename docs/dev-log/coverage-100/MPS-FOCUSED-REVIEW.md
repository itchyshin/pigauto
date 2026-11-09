# Real MPS fixed-effect prediction review

A real Apple M1 Ultra MPS runtime passes the targeted stored-state prediction test with 12 assertions. The coordinator reran the final RNG-isolated test on the actual device. The test uses a 12-tip, one-epoch fit with two numeric covariates, a one-column response and a length-one direct-effect bias. A known weight matrix and bias isolate the baseline-plus-fixed-effect contract while the fitted model and tensor operations remain real.

The expected latent matrix is computed independently in R from the fit's baseline and preprocessed covariates. Assertions check the full matrix, missing cells, finite output and observed-value preservation when completing data. Increasing the stored bias by 0.4 must increase missing-cell predictions by 0.4. No numerical or device-selection mock is used. R random state and torch thread count are restored.

Final test SHA-256: `50eec2bf6c3230beeb66ff222e724af6e5509d356f75423db481c84cd490ce91`. Sol independently reviewed the oracles and final isolation change. Focused runs on the primary source and a fresh f96a export both passed; these have different complete source identities. The frozen f96a predictor SHA-256 is `78dd630a9bb5c55ccabd16f9120f61d644a65a75a29fcfff7a1912a5355d2762`.

A CPU-only Totoro run explicitly skips with `real MPS device unavailable`. Missing torch runtime uses the separate exact reason `real MPS runtime unavailable`; unknown MPS reasons remain unclassified by the audit parser.

This test checks a narrow stored-parameter prediction and scalar-bias contract. Learned coefficient recovery, general GPU support and statistical calibration remain outside its evidence. The test was added after the f96a freeze. A source-bound instrumented MPS artifact and independently verified multi-runtime accounting remain pending, so no coverage gain is credited to the frozen CPU run.
