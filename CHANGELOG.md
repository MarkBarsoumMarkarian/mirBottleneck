# Changelog

## [v0.2-tooling] - in progress

### Added
- `v0.2-tooling` branch for package hardening and tooling infrastructure
- Planned: `testthat` unit tests for core scoring functions
- Planned: `pkgdown` documentation site
- Planned: GitHub Actions CI workflow running tests against demo data

## [v0.1.0] - 2026-03-20

### Initial release
- `build_network()` — miRNA-target network via miRTarBase/multiMiR
- `score_vss()` — Variance Suppression Score
- `score_coherence()` — Coherence Induction Score  
- `classify_bottleneck()` — miRNA classification into Silencer/Conductor/Dual/Weak
- `composite_score()` — patient-level weighted bottleneck index
- `survival_model()` — Cox models + log-rank test
- Harmonized TCGA-PAAD data objects (26 `.rds` files)
- Kaplan-Meier figure (High vs Low bottleneck score)
