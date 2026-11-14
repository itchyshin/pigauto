# Installation-guide local site review

Snapshot base: d0b3b2516b4ab8f60d60a07c2ffe2b83cd21f313. The README/Rmd patch hash is 9d37f7d54558672eaa9687aa6afe45728f3286e0eac0f79679886ef34c8196e2. The later generated vignette R script reproduces the Rmd installation chunks; it is not a different article input.

The final fresh build completed with install=FALSE and one-thread settings. Its recorded font overrides use system fonts, and CRAN NEWS-date lookup is disabled. Those local build choices do not validate production typography or deployment. R/ matched all 49 files in the checked c2 source; the installed package was 0.11.0.

Cleaner and six crawler controls passed. The structural crawler measured 62 HTML pages, 3587 local references and 34 archived retired pages, with zero errors. The retirement audit passed all 44 retired routes; sitemap 62, search 614 entries and 568 path values were measured. The verifier's old 613/567 expectations were replaced by the observed build counts for that audit; this verifies route/discovery constraints, not a prediction of those counts.

In Codex's own browser over localhost, the driver inspected the changed homepage and Getting Started layouts and read the rendered version guidance. Both visibly ask readers to check packageVersion and explain GitHub installation if CRAN is older. The novice code and warning callout remained legible. The earlier candidate review also sampled the three changed reference pages and NEWS. The UI screenshots are in this task's tool transcript; no screenshot-file hash is claimed.

The first scratch copy put the Rmd at the wrong directory level. The rendered-text check detected the missing article note; that output was retained separately and the corrected snapshot was rebuilt. It was not counted as a pass.

This review does NOT cover the actual main deployment or a visual inspection of every page. After merge, verify the deployment SHA, all live route/discovery checks and the changed-page delta sample, including Getting Started using the production website.
