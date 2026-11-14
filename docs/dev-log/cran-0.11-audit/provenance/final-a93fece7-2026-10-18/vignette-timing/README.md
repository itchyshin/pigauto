# Successor archive vignette timing diagnostic

**Verdict: PASS for this bounded diagnostic.** All six packaged vignettes rendered from a fresh extraction of the exact successor archive using the package installed in the successor `mi-matrix/neither` library. This is separate vignette evidence; it is not a native check result or a release-wide readiness verdict.

Archive: `/private/tmp/pigauto-release-candidate/sha256-a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb/pigauto_0.11.0.tar.gz`  
SHA-256: `a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb`; size: 5,137,837 bytes  
Source commit: `b2e6054685e6bdedb073ba46e058af2c0a619a18`  
Installed pigauto: `/private/tmp/pigauto-final-closure-20261017/successor-b2e6054/mi-matrix/neither/pigauto`, version 0.11.0  
R: 4.6.0 on aarch64 macOS; torch installed; Pandoc 3.12; `NOT_CRAN` unset; all numerical thread caps set to 1.  
The install receipt records a successful install from this exact archive digest. `drmTMB` and `gllvmTMB` were absent from the library.

Each render used the default rmarkdown/knitr behavior. The option trace hook recorded chunk `eval` settings and returned options unchanged. The 20 teaching chunks in `getting-started.Rmd` and seven analysis chunks in `multiple-imputation.Rmd` remained disabled. No global eval override, retry, or warning replay was used. Render errors and nonzero exits remained failures.

| Vignette | Render seconds | Whole R process seconds | Enabled chunks | Disabled chunks | Exit |
|---|---:|---:|---:|---:|---:|
| common-pitfalls.Rmd | 0.771 | 4.065 | 1 | 8 | 0 |
| getting-started.Rmd | 2.719 | 5.502 | 7 | 20 | 0 |
| gnn-architecture.Rmd | 0.505 | 3.353 | 1 | 0 | 0 |
| mixed-types.Rmd | 1.656 | 4.992 | 11 | 1 | 0 |
| multiple-imputation.Rmd | 0.809 | 3.966 | 1 | 8 | 0 |
| tree-uncertainty.Rmd | 0.668 | 3.434 | 1 | 3 | 0 |

Total process time was 25.340 seconds, below the 600-second total cap and 120-second per-vignette cap. The prior estimate was 5–8 minutes; this run finished much faster (about 25 seconds). The render timer excludes R startup and preflight; process time includes those costs.

Archive SHA-256 was unchanged before and after the six renders. All six extracted Rmd files retained their recorded SHA-256 values before and after. Raw process logs, render-level JSON, HTML outputs, chunk traces, commands, per-script hashes, and the installation receipt reference are retained in this directory.

Exact launcher: `python3 /private/tmp/pigauto-final-closure-20261017/successor-b2e6054/vignette-timing/run-timings.py`. Per-vignette commands and timings are recorded in `execution.json` and `timing-summary.json`.

Limits: this measured documentation rendering only. It does not establish native or Windows `R CMD check` timing, comprehensive teaching-code execution, or release/upload/CRAN acceptance.
