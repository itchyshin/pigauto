# pigauto: Mondrian conformal on real data, and conformal MI under phylogenetic GLS

Draft for the brain vault (not written there; needs Shinichi's approval). Links:
[[pigauto-split-conformal-explainer]], [[PROJECTS]].

- Mondrian conformal, tested on masked observed cells in PanTHERIA, AVONET and FishBase
  under a pre-registered two-arm design (2026-09-23): pooled over traits it lifts
  far-stratum coverage from 0.92-0.93 to 0.94-0.96 and trims near-stratum
  over-coverage (0.97-0.98 to 0.95-0.97). The pre-registered rule kept split as the
  default: its near non-inferiority condition penalised removal of over-coverage. A new
  rule would need a new confirmation. Evidence: pigauto PR #188,
  `docs/dev-log/mondrian-realdata/results.md`.
- Conformal multiple-imputation draws (split or Mondrian) halve a phylogenetic GLS slope
  in simulation while OLS stays unbiased; an oracle proper imputation is unbiased. Causes
  on one tree: draw spread about 1.3 times the proper conditional spread, and the default
  prediction route ignores co-observed traits. Pre-existing. Lane:
  `arc/mi-gls-attenuation` (fir array 61137481). Lesson: validate MI with a pooled GLS
  downstream, not only marginal draw calibration.
