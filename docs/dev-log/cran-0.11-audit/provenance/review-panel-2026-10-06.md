# Fresh review panel, 2026-10-06

Three independent, read-only reviews were run against the merged-source audit
packet. The panel assessed bounded evidence and user-facing surfaces. It did not
authorize a CRAN release.

| Reviewer | Lens | Final verdict | Basis and scope |
|---|---|---|---|
| Grace | Release and CI evidence | READY | Re-verified merge, tarball and check-log identities; confirmed the exact-artifact check receipt. |
| Rose | Claims, evidence and provenance | READY | Re-verified the artifact/check receipts and negative controls; treated missing raw web response hashes as a provenance limitation after the direct cache-busted check. |
| Pat | Reader workflow | READY | Revised the initial NOT_READY after current pages were fetched directly with a cache-busting query and matched the merged navigation, defaults and scoped MI claims. |

The first Pat review surfaced older search-indexed content that conflicted with
the retained deployment receipt. A direct browser recheck with
`?audit=2026100622` confirmed the current homepage, getting-started article,
multiple-imputation article and `multi_impute()` reference. The detailed markers
and limitation are in `post-review-live-page-recheck-2026-10-06.md`. Pat revised
the reader-surface verdict to READY. The other two reviewers independently
re-verified the changed evidence packet and also returned READY.

The panel's READY applies only to G8 audit evidence. BirdTree underlying-data
redistribution rights remain unresolved. The release ledger remains
`NOT_READY`, and this review does not permit upstream contact or CRAN submission.
