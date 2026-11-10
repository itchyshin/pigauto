# Prepared pull requests

All listed code PRs are drafts and require separate approval before merging. The final evidence closeout is a docs-only draft stacked on #288; its PR body links this report. Code heads below are independent of the documentation publication commit.

## Coverage and test infrastructure

| PR | Exact code head | CI snapshot | Scope |
|---|---|---|---|
| [#287](https://github.com/itchyshin/pigauto/pull/287) | `bd6a643a9be222a23af54b311df374c98e2f5cdd` | Green; pkgdown intentionally skipped | Audit coverage gaps and test presentation and evaluation contracts |
| [#288](https://github.com/itchyshin/pigauto/pull/288) | `638445a002ce5118f78561c81cf35c2f313a514c` | Pending runner/checks | test: extend prediction fitting and imputation coverage |

## Separate runtime repairs

| PR | Exact code head | CI snapshot | Scope |
|---|---|---|---|
| [#284](https://github.com/itchyshin/pigauto/pull/284) | `08b7cdc6df9fbacea3dbf305328f28b57d67f02d` | Green; pkgdown intentionally skipped | Fix uncertainty plots to use sorted prediction ranks |
| [#285](https://github.com/itchyshin/pigauto/pull/285) | `81e5da2c636190c086ed31f079a15dd5ead30537` | Green; pkgdown intentionally skipped | Fix phylogenetic masking with one eligible clade |
| [#286](https://github.com/itchyshin/pigauto/pull/286) | `ae2f1c70fd0dec3eb4c032a503e18c692dfd5490` | Green; pkgdown intentionally skipped | Refuse Gaussian interval displays for discrete uncertainty scores |
| [#289](https://github.com/itchyshin/pigauto/pull/289) | `8bc94aaae17d9a0eed4af934ed5b61e9a55be5f9` | Green; pkgdown intentionally skipped | Fix Blomberg K signal gate |
| [#290](https://github.com/itchyshin/pigauto/pull/290) | `e43655e1d2cbd493db69a9bb2dcefb4cca2cde99` | Green; pkgdown intentionally skipped | fix: preserve covariance dimensions for one missing cell |
| [#291](https://github.com/itchyshin/pigauto/pull/291) | `34dd16617b1f1aa91aa8480fabee208a981490af` | Green; pkgdown intentionally skipped | fix: label comparisons with the supplied method |
| [#292](https://github.com/itchyshin/pigauto/pull/292) | `9ed71bb8c014e707f38c0170f44ac950a2931066` | Green; pkgdown intentionally skipped | Fix validation-best GNN checkpoint restoration |
| [#293](https://github.com/itchyshin/pigauto/pull/293) | `5e8deabbec07c4999097a8fe01eda75a549b49b3` | Green; pkgdown intentionally skipped | Preserve closed gates and reject non-finite GNN blends |
| [#294](https://github.com/itchyshin/pigauto/pull/294) | `22f8a5fa20ecbfed906944a3c9435e9667006879` | Green; pkgdown intentionally skipped | Remove unreachable Bayesian weight-sum fallback |
| [#295](https://github.com/itchyshin/pigauto/pull/295) | `5f2562c15d4c650aa57ce7893ef20378cfcccdab` | Green; pkgdown intentionally skipped | Remove fallbacks ruled out by earlier validation |
| [#296](https://github.com/itchyshin/pigauto/pull/296) | `c50f80f41bc4cf4b5a966acd2c7a94406a5a9fc5` | Green; pkgdown intentionally skipped | Remove duplicate analysis-MI completeness check |
| [#297](https://github.com/itchyshin/pigauto/pull/297) | `592df7340edc018f10ad929fc231d4dcf27fa99b` | Green; pkgdown intentionally skipped | Preserve causes when benchmark pipeline attempts all fail |
| [#298](https://github.com/itchyshin/pigauto/pull/298) | `946c33cb5014d276d2a501622665d7472eb03359` | Green; pkgdown intentionally skipped | Simplify redundant imputation guards |
| [#300](https://github.com/itchyshin/pigauto/pull/300) | `16555adc0abb7179e8f9d5407653411b686c4e34` | Green; pkgdown intentionally skipped | Check private conditional-draw observations before estimation |
| [#301](https://github.com/itchyshin/pigauto/pull/301) | `96cd4dbeca19ce47d80f1ee4da556fe866768afe` | Green; pkgdown intentionally skipped | Remove unreachable empty GBIF cache-key fallback |
| [#302](https://github.com/itchyshin/pigauto/pull/302) | `f75d2c994f90c67f3dee49dad597cc7a07c4058e` | Green; pkgdown intentionally skipped | Retain last accepted joint refinement after divergence |
| [#303](https://github.com/itchyshin/pigauto/pull/303) | `aa2ccddb51a46d11022c5a4f61b63f60ae767e3a` | Green; pkgdown intentionally skipped | Refuse non-finite exact conditional variances |
| [#304](https://github.com/itchyshin/pigauto/pull/304) | `ffe108cdf72f5abde733e870ade4b0d800e3e728` | Green; pkgdown intentionally skipped | Reject overflowed numeric covariate encodings with an early diagnostic |
| [#305](https://github.com/itchyshin/pigauto/pull/305) | `0ba4ec0df598a43c29515b9e84a990b645411b8f` | Green; pkgdown intentionally skipped | Remove unreachable empty calibration fold fallbacks |
| [#306](https://github.com/itchyshin/pigauto/pull/306) | `b913036942ad0aa73ee48c62c08343bf66c82ab7` | Green; pkgdown intentionally skipped | Validate prediction counts before device selection |

Snapshot UTC: 2026-11-10T00:30:11.586019+00:00.

PR #288 stacks on #287. The combined verification branch includes all twenty repairs but is separate from main and the coverage-only PR head. The docs-only closeout does not change runtime files, tests, coverage configuration or workflows. Scientific calibration, general GPU correctness, singleton limitations and CRAN release remain separate.
