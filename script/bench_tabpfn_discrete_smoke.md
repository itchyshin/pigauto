# TabPFN AVONET Discrete Benchmark

- Generated: 2026-06-10 17:27:15 MDT
- Git commit: `9f9ddf2cf1a835d3eb05494f56e8fd04eb5d4c92`
- Scales: `50`
- Variants: `lappe_nfa`
- Replicates: `1`
- Dry run: `FALSE`

## Design

AVONET categorical/ordinal targets are held out with pigauto's `make_missing_splits()` and scored as classification tasks on the same test rows used for baseline and pigauto. TabPFN feature sets use same-row latent traits, optional Laplacian phylogenetic eigenvectors, and optional nearest-training-target class-proportion features. This benchmark reports accuracy and balanced accuracy only; it does not attempt classification prediction sets.

## Status Counts

           method status Freq
         baseline     ok    3
         majority     ok    3
 tabpfn_lappe_nfa     ok    3

## Test Summary

           method scale_n             trait        type accuracy
         majority      50         Migration     ordinal   1.0000
 tabpfn_lappe_nfa      50         Migration     ordinal   1.0000
         baseline      50         Migration     ordinal   0.7143
         baseline      50 Primary.Lifestyle categorical   0.7500
         majority      50 Primary.Lifestyle categorical   0.7500
 tabpfn_lappe_nfa      50 Primary.Lifestyle categorical   0.7500
         baseline      50     Trophic.Level categorical   0.6667
 tabpfn_lappe_nfa      50     Trophic.Level categorical   0.6667
         majority      50     Trophic.Level categorical   0.5000
 balanced_accuracy
            1.0000
            1.0000
            0.7143
            0.3333
            0.3333
            0.3333
            0.5556
            0.5000
            0.3333
