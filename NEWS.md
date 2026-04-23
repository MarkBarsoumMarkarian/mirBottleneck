# mirBottleneck 0.99.0

## New features

* Added modular exported APIs: `compute_vss()`, `compute_cis()`, and `classify_archetypes()`.
* Preserved backward compatibility: `score_vss()`, `score_coherence()`, and `classify_bottleneck()` remain available and now delegate to modular APIs.
* Added `SummarizedExperiment` support in modular compute functions with assay-name overrides (`mirna_assay`, `mrna_assay`) and input validation helpers.
* Added `plot_archetype_landscape()` returning a customizable `ggplot2` object.
* Updated `run_mirBottleneck_project()` orchestration to use modular APIs while keeping existing workflow and interface behavior.
* Added initial test coverage for modular outputs, wrapper smoke parity, `SummarizedExperiment` path, and validation error messaging.
* Added vignette scaffold `vignettes/modular-workflow.Rmd` and pkgdown scaffold `_pkgdown.yml`.
