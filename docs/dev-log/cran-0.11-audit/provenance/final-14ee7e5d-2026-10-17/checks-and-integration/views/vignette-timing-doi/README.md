# Final archive per-vignette timing diagnostic

All six vignettes rendered successfully from a fresh exact-archive copy using final installed pigauto, verified CPU torch, explicit library paths, NOT_CRAN unset and one numerical thread. Original source Rmd bytes and archive SHA-256 agree before/after. Teaching chunks retain their disabled settings; chunk-option traces are retained.

| Vignette | Render elapsed s | Whole R process s | Exit |
|---|---:|---:|---:|
| common-pitfalls.Rmd | 5.812 | 39.061 | 0 |
| getting-started.Rmd | 23.682 | 64.075 | 0 |
| gnn-architecture.Rmd | 9.878 | 43.365 | 0 |
| mixed-types.Rmd | 16.645 | 64.054 | 0 |
| multiple-imputation.Rmd | 10.701 | 52.372 | 0 |
| tree-uncertainty.Rmd | 8.168 | 45.354 | 0 |

Render time begins after package/runtime preflight; whole-process time also includes startup/preflight. These are separate local diagnostic measurements under current host load. They do not replace native R CMD check timing, Windows timing, or manual/platform proof. No global eval override or retry was used. Setup extraction was still running when an initial launch attempt found no script; no R render ran on that failed prelaunch command. The ordinary five-minute estimate was exceeded and reported; authorized upper limits remained in force.

The known DOI returned403 to the exact-source URL checker. One Crossref metadata request returned200 and matches DOI/title/Hackett/Science320(5884)/1763–1768/2008-06-27. See doi-adjudication.md and raw response/metadata. Bibliographic identity is supported; automated access remains403, and URL-clean/CRAN acceptance is not claimed.

All scripts, commands, versions, raw logs, per-document start/elapsed/exit, hashes, outputs and chunk traces are retained. No source/archive/repository changes, native-check/MI directory writes, campaign, commit/push, merge or submission. Final package/platform/timing approval belongs to the driver. Formal scratch-draft completion assessment not performed: NOTASSESSED.
