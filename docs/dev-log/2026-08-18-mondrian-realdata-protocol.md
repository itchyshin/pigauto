# Stage C: real-data Mondrian confirmation protocol

Independently of Stage B, FishBase and PanTHERIA will mask only cells that were
originally observed.  For each trait, compare `conformal_method = "split"` and
`"mondrian"`, retaining masked-cell count, calibration/validation count,
activation or fallback, coverage, interval width, MCSE, runtime, and errors.

The only promotable result is empirical calibration improvement with transparent
fallback and a stated width cost.  This study cannot establish an unconditional
95% guarantee or justify changing the default.
