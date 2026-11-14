# Prepared pull requests

All listed PRs are drafts and require separate approval before merging. This table records exact heads; the final CI snapshot will be refreshed after final verification. Coverage-only PRs and runtime repairs have different source identities.

## Coverage and evidence

| PR | Exact head | Scope |
|---|---|---|
| [#287](https://github.com/itchyshin/pigauto/pull/287) | `bd6a643a9be222a23af54b311df374c98e2f5cdd` | Audit coverage gaps and test presentation and evaluation contracts |
| [#288](https://github.com/itchyshin/pigauto/pull/288) | `b142267903b2b8ab64e5cecb9adc3ec06d00d67f` | test: extend prediction fitting and imputation coverage |

## Separate repairs

| PR | Exact head | Scope |
|---|---|---|
| [#284](https://github.com/itchyshin/pigauto/pull/284) | `08b7cdc6df9fbacea3dbf305328f28b57d67f02d` | Fix uncertainty plots to use sorted prediction ranks |
| [#285](https://github.com/itchyshin/pigauto/pull/285) | `81e5da2c636190c086ed31f079a15dd5ead30537` | Fix phylogenetic masking with one eligible clade |
| [#286](https://github.com/itchyshin/pigauto/pull/286) | `ae2f1c70fd0dec3eb4c032a503e18c692dfd5490` | Refuse Gaussian interval displays for discrete uncertainty scores |
| [#289](https://github.com/itchyshin/pigauto/pull/289) | `8bc94aaae17d9a0eed4af934ed5b61e9a55be5f9` | Fix Blomberg K signal gate |
| [#290](https://github.com/itchyshin/pigauto/pull/290) | `e43655e1d2cbd493db69a9bb2dcefb4cca2cde99` | fix: preserve covariance dimensions for one missing cell |
| [#291](https://github.com/itchyshin/pigauto/pull/291) | `34dd16617b1f1aa91aa8480fabee208a981490af` | fix: label comparisons with the supplied method |
| [#292](https://github.com/itchyshin/pigauto/pull/292) | `9ed71bb8c014e707f38c0170f44ac950a2931066` | Fix validation-best GNN checkpoint restoration |
| [#293](https://github.com/itchyshin/pigauto/pull/293) | `5e8deabbec07c4999097a8fe01eda75a549b49b3` | Preserve closed gates and reject non-finite GNN blends |
| [#294](https://github.com/itchyshin/pigauto/pull/294) | `22f8a5fa20ecbfed906944a3c9435e9667006879` | Remove unreachable Bayesian weight-sum fallback |
| [#295](https://github.com/itchyshin/pigauto/pull/295) | `5f2562c15d4c650aa57ce7893ef20378cfcccdab` | Remove fallbacks ruled out by earlier validation |
| [#296](https://github.com/itchyshin/pigauto/pull/296) | `c50f80f41bc4cf4b5a966acd2c7a94406a5a9fc5` | Remove duplicate analysis-MI completeness check |
| [#297](https://github.com/itchyshin/pigauto/pull/297) | `592df7340edc018f10ad929fc231d4dcf27fa99b` | Preserve causes when benchmark pipeline attempts all fail |
| [#298](https://github.com/itchyshin/pigauto/pull/298) | `946c33cb5014d276d2a501622665d7472eb03359` | Simplify redundant imputation guards |
| [#300](https://github.com/itchyshin/pigauto/pull/300) | `16555adc0abb7179e8f9d5407653411b686c4e34` | Check private conditional-draw observations before estimation |
| [#301](https://github.com/itchyshin/pigauto/pull/301) | `96cd4dbeca19ce47d80f1ee4da556fe866768afe` | Remove unreachable empty GBIF cache-key fallback |
| [#302](https://github.com/itchyshin/pigauto/pull/302) | `f75d2c994f90c67f3dee49dad597cc7a07c4058e` | Retain last accepted joint refinement after divergence |
| [#303](https://github.com/itchyshin/pigauto/pull/303) | `aa2ccddb51a46d11022c5a4f61b63f60ae767e3a` | Refuse non-finite exact conditional variances |
| [#304](https://github.com/itchyshin/pigauto/pull/304) | `ffe108cdf72f5abde733e870ade4b0d800e3e728` | Reject overflowed numeric covariate encodings with an early diagnostic |
| [#305](https://github.com/itchyshin/pigauto/pull/305) | `0ba4ec0df598a43c29515b9e84a990b645411b8f` | Remove unreachable empty calibration fold fallbacks |
| [#306](https://github.com/itchyshin/pigauto/pull/306) | `b913036942ad0aa73ee48c62c08343bf66c82ab7` | Validate prediction counts before device selection |

PR #288 stacks on #287. Its head will receive the final evidence report, so refresh that head before approval. The final combined source includes the prepared repairs for verification only; no merge or CRAN release is authorized. Scientific calibration, general GPU correctness and the singleton limitation remain separate from these engineering checks.
