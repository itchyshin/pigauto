# pigauto 0.11.0 candidate component ledger

This ledger defines the local candidate's data and asset boundary. It does not
claim platform-clean or submission-ready status. The frozen tarball inventory
must prove every `SHIP` and `EXCLUDE` decision below.

| Component | Source and holder | Permission | Transformation / consumer | Decision |
|---|---|---|---|---|
| `data/avonet300.rda`, `data/avonet_full.rda` | AVONET, Tobias et al. (2022) | CC BY; attribution in `inst/NOTICE` | Subsets and transformations used by package examples and documentation | SHIP |
| `data/tree300.rda`, `data/trees300.rda`, `data/tree_full.rda` | `megatrees` 1.0.0 release; Li; underlying Jetz et al. posterior | MIT upstream distribution; provenance and upstream digest in `inst/NOTICE` | Pruned example phylogenies | SHIP |
| `data/ctmax_sim.rda` | pigauto-generated simulation | Original package material under MIT | Reproducible generator in `data-raw/make_ctmax_sim.R`; repeated-observation example | SHIP |
| `inst/extdata/legacy_fit_v091.rds` | pigauto-generated compatibility fixture | Original package material under MIT | Reproducible generator in `data-raw/make_legacy_fit_v091.R`; saved-fit compatibility tests | SHIP |
| `inst/extdata/novice_traits.csv`, `inst/extdata/novice_tree.nwk`, `inst/examples/novice-workflow.R` | pigauto-generated teaching fixture and script | Original package material under MIT | Eight-tip, three-trait novice workflow with two missing cells | SHIP |
| `man/figures/logo.png` | pigauto repository artwork | Original package documentation under MIT | README and package identity | SHIP |
| `BACE/**` | External comparator package retained in repository | Separate upstream terms; not adjudicated for pigauto distribution | Comparator development only | EXCLUDE |
| `avonet/**` | Raw AVONET and phylogeny inputs | Third-party sources described in `inst/NOTICE`; raw redistribution not part of the installed product | Source-only data preparation inputs | EXCLUDE |
| `data-raw/**` | First-party generators plus references to upstream inputs | Development material; generated SHIP objects adjudicated above | Data construction scripts | EXCLUDE |
| `useful/**` | Historical PDFs, data snapshots and development results from multiple holders | Not adjudicated for redistribution in the installed package | Research/development archive | EXCLUDE |
| `pkgdown/**` | Generated and historical website assets | Reader surface reviewed separately by G9; not installed-package content | Local website inputs | EXCLUDE from tarball |
| `script/**`, `dev/**`, `docs/**`, `checkpoints/**`, `.github/**`, `misc/**` | Development and evidence material | Not part of installed product | Build exclusions | EXCLUDE |

`Authors@R` contains one `aut`, `cre` and `cph` person, Shinichi Nakagawa. This
local candidate records the repository metadata as supplied by that author; it
does not substitute for a later submission-time consent and policy refresh.
