# Cache-busted public page recheck, 2026-10-06

After the independent reader review surfaced older search-indexed content, the four current pages were loaded directly in the primary app browser with the unique query `?audit=2026100622`. The pages rendered successfully. Their navigation identifies version 0.11.0.

| Page | Observed current content |
|---|---|
| Home | Navigation is Get started, Articles, Reference, Changelog. The page says the baseline is the default, `gnn = FALSE`, and the GNN is optional. No retired per-trait benchmark links appeared. |
| Getting started | `torch::install_torch()` is conditional on opting into `gnn = TRUE`; the GNN-off default and narrowed uncertainty claims match the merged source. No blanket 95% coverage guarantee appeared. |
| Multiple-imputation article | The article distinguishes posterior imputations from diagnostic draws, limits coverage statements to named simulations, and retains the downstream model examples with scope caveats. |
| `multi_impute()` reference | The usage lists `draws_method = "auto"` and `gnn = FALSE`; the description explains the posterior route and the diagnostic-only fallback. |

This direct browser observation resolves the earlier conflicting search-index fetch for the reader review. The browser did not expose response headers or a raw-response hash, so this receipt records rendered content checks rather than byte-level HTTP evidence. A shell `curl` attempt could not resolve `itchyshin.github.io` in the local environment; the browser fetch itself succeeded.

## Re-verification

On the same audit checkout, the CRAN release-gate selftest passed, including every planted negative control. The post-merge R CMD check log SHA-256 reverified as `657a19dc6584f1475f78f35d96824304be45b1279459f72bfa5f1c72113d41a0`; the test output as `2f97450cb9b7eebbf4c89d0ece75b59f3447f778382901294d916f2e9db29613`; and the artifact identity receipt as `d668c8cd985ceaa0af56cb4766f83e3a5496f9eab3c65de108bfed7d1ee7853f`. The checkout HEAD is the recorded merge commit `bb5835d1b214d783da7b0b414df99aa6ba926bc7`.

The reader reviewer revised its reader-surface verdict to READY after this recheck. This does not resolve the BirdTree redistribution-rights question or grant CRAN submission authorization.
