# PR 283 review and verification

Reviewed on 9 November 2026. The changes add development scripts and evidence only; package runtime and CRAN release operations remain separate.

One input-validation gap was reproduced: the post-hoc subgroup script accepted records with the wrong source or runner identity, and accepted altered truth values shared across methods. It now requires the same complete-roster validation as the primary summarizer before computing subgroup metrics. The primary validator also rejects missing, non-finite or reversed interval endpoints, matching the runner's existing contract. A synthetic regression test failed before the fix and passed afterwards for source, runner, truth and interval corruption. No model fits or reported scientific results changed.

Verification: the seed harness, tie diagnostic and stored-draw subset self-tests all passed. The review checked target masking before preprocessing, input-order mapping, fixed seed groups, complete attempt accounting and the distinction between row holdouts and singleton species. No additional blocking findings were identified. The subgroup limitation remains explicit; these results do not establish a general fallback, interval calibration or downstream-inference validity.

The latest successful main-branch coverage artifact, run 42537328573 at source 033649b9, reports 10,845 covered lines out of 12,749 (85.0655%), leaving 1,904 uncovered lines. This PR's development scripts are outside the package coverage denominator. A package-wide coverage effort must use tests of behavior without hiding uncovered lines or treating line execution as scientific validation.
