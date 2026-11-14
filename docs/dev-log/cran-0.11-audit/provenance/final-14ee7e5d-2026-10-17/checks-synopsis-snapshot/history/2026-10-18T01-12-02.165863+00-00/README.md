# Final checks evidence synopsis

This is an evidence snapshot for pigauto 0.11.0, not a release verdict. The exact archive is 5,137,648 bytes, SHA-256 `14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb`, generated from clean merged source `8e275ca738129c8d608f06ffc407b8141c1a2636`. It contains 250 entries/236 regular files and 7,150,500 regular uncompressed bytes. `synopsis.json` retains exact values; `inputs.json` binds them to retained receipts.

All three remote native checks report Status:OK. Each has 29 available declared dependencies and 3,195 passes, 161 warnings, 86 skips and no failures. The 161 warnings and 86 skips remain observations pending final consequence review; their counts are not a blanket benign classification. Separate warning replays matched these counts and are excluded from native totals.

| Remote platform | Native total s | Examples s | Tests s | Vignette rebuild s | PDF s | HTML s |
|---|---:|---:|---:|---:|---:|---:|
| ubuntu-latest-release | 336.330 | 14.0/14.0 | 206.0/206.0 | 11.0/11.0 | 2.0/2.0 | 2.0/2.0 |
| windows-latest-devel | 622.432 | 16.0 | 404.0 | 22.0 | 4.0 | 4.0 |
| windows-latest-release | 645.813 | 17.0 | 415.0 | 22.0 | 4.0 | 13.0 |

Linux stage pairs mean rounded CPU/elapsed seconds. Windows single values remain labelled reported seconds; unrounded Rout process timings are separate in JSON. Windows totals 645.813 and 622.432 seconds exceed 600 seconds, while vignette rebuilding is 22 seconds on each. This does not establish that CRAN will accept the timings. The protocol calls 600 seconds an observed incoming signal, not an immutable policy limit. Its strict first-submission rule is distinct: the retained CRAN snapshot shows 0.10.0, so this is an update. Maintainer/final-panel timing review stays open.

The local Mac native check is `COMPLETED_UNADJUDICATED`; warning replay is `PENDING_NOT_YET_PROVEN`. The known compare_methods NOTE reports user 2.532/system 0.992/elapsed 18.406 seconds (CPU 3.524 seconds). It is unadjudicated and must not be rerun until faster. Observed local test counts are retained separately from final check completion. The guarded refresh helper preserves native exit versus strict checker-driver rejection and leaves any completed warning replay unadjudicated.

Native installed-size INFO reports 5.9Mb on all three remote platforms, with data 2.2Mb, extdata 1.1Mb and help 1.1Mb; the local snapshot reports 7.0Mb, with data 3.2Mb, extdata 1.1Mb and help 1.1Mb. These are raw platform observations, not identical byte-size measurements or native warnings. Archive bytes are unchanged. No package-owned src files ship; compiled dependencies remain separate requirements.

All six individual vignette diagnostics exited 0, with render-time sum 74.8865 seconds and whole-R-process sum 308.2825 seconds. The enclosing sequential loop took 310.1696 seconds, including between-call overhead. Package/runtime startup is outside the render interval. These separate local measurements do not replace native or Windows timing. Thirty-two explicit evalFALSE chunks and the multiple-imputation global evalFALSE setting remained disabled.

| Vignette | Render elapsed s | Whole R process s |
|---|---:|---:|
| common-pitfalls.Rmd | 5.8116 | 39.061 |
| getting-started.Rmd | 23.6822 | 64.075 |
| gnn-architecture.Rmd | 9.8780 | 43.365 |
| mixed-types.Rmd | 16.6452 | 64.054 |
| multiple-imputation.Rmd | 10.7013 | 52.372 |
| tree-uncertainty.Rmd | 8.1682 | 45.354 |

The isolated missing-rgbif/missing-terra proof records 60 passes, 0 failures, 0 warnings and 8 expected skips: 2 rgbif-dependent tests, 3 terra-dependent tests and 3 live/slow integrations. It exercises helpful absence branches and bounded cache/helper/failure tests. Dependency symlinks plus base R exclude user/site-library fallthrough; this is not a clean-room dependency build. Live GBIF/WorldClim availability, downloads and skipped raster/stored-point paths are not proven.

The exact-source URL check retains one distinct DOI 403 in two locations. The retained [Crossref metadata response](https://api.crossref.org/works/10.1126/science.1157704) returned 200 and identifies Hackett's 2008 Science paper, volume 320(5884), pages 1763–1768, with the exact DOI/title/date. Identifier identity is supported; automated DOI access remains 403. No URL-clean verdict is made. Initial sandbox DNS/export failures remain in original receipts.

Spelling review covers 404 flagged tokens, 1,886 contexts and 43 archive-bound source files. It identified no ordinary-word typo or uncertain spelling token among these flags. This does not establish unflagged prose, grammar, scientific claims or reference quality. The category/scope counts and original contexts remain in the spelling receipt.

Installed MI verification passed all 12 steps across neither, drm-only, gllvm-only and both optional-backend cells. Both measured backend versions are 0.7.1; their source commits are not embedded in the installed packages. Dependency symlinks and three downstream fits per real fixture bound this proof; it is not a clean-room build or a recovery campaign. Artifact source and checker commits are recorded separately.

Earlier Linux “timings absent” wording was a derived parser error and is superseded by paired timings. The original service receipt's pending broad spelling review is superseded only for its 404 flags. STATE's predecessor warning-capture and old running-build/checker wording must not be promoted into current final-artifact evidence. `contradictions_and_repairs` records these distinctions.

No panel verdict or release-rung promotion is made. Local completion/replay, timing and warning-consequence review remain open. PR 237 is recorded merged; PR 236 remains unmerged under retained authority. CRAN upload/confirmation are separately gated.

After native completion, update with:

```
python3 /private/tmp/pigauto-final-closure-20261017/final-checks-synopsis/build-synopsis.py --refresh-local
```

The helper reads proofs only and preserves prior outputs in history. It refuses a changed or extra native NOTE, error/warning, hash mismatch or invalid completed replay. A later refresh changes the exact draft: repeat its assessment before describing it as assessed. Assessment completion is recorded separately, never implied by this synopsis.
