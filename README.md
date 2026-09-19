# FiO₂ de-escalation: landmark and measurement analysis

Analysis code accompanying **Documented FiO₂ De-escalation and In-Hospital Mortality Across Three Critical-Care Databases: A Landmark and Measurement Analysis**.

This repository contains extraction SQL, cohort construction, statistical estimation, diagnostic and figure-generation code for MIMIC-IV, eICU-CRD and SICdb. It is a code-only release of the analysis underlying the September 2026 manuscript revision. Repository creation is not prospective study registration.

## Scientific scope

The reference comparison is the first classifiable opportunity per ICU stay/case: documented FiO₂ reduction of at least 10 percentage points versus stability within ±5 points. Outcome follow-up begins after the classification interval. The reference endpoint is in-hospital mortality; source-specific recorded 28-day death is secondary. Treatment comparability, recording selection, and SICdb hourly/minute measurement definitions are examined separately.

The study is observational and exploratory. This release does not claim a sustained-strategy target trial, formal longitudinal clone-censor-weight analysis, externally validated individualized treatment rule, or clinically adjudicated minute-level gold standard. Post hoc support and measurement analyses do not replace the reference estimate based on significance.

## Start here

1. Read [data access and governance](docs/DATA_ACCESS.md).
2. Follow [environment configuration and execution order](docs/REPRODUCING.md).
3. Consult [manuscript-to-code mapping](docs/ANALYSIS_MAP.md).
4. Run the data-free checks: `python tools/validate_release.py` and `python analysis/audited_reanalysis_v2/test_estimand.py`.

The release contains no patient records, source archives, fitted individual predictions, credentials, private signal-review packets, manuscript files or clinical data outputs. Credentialed source access is required for statistical reproduction. The original computation environment is recorded in `analysis/software_versions.json`. Full analyses are computationally intensive and write controlled outputs next to their scripts; execute a copy in a private, unshared working directory.

## Layout

The `analysis/` stage folders must remain siblings because later stages import earlier estimator code. Stage names preserve analysis provenance; they are not competing versions from which to choose a favorable result. The v2 estimator remains the reference. Later stages provide the specified post hoc analyses and figure refinements.

Publication preparation scripts that modify author Word documents or Desktop delivery folders are intentionally excluded. The conceptual Figure S5 renderer is retained because it is independent of those author files.

## Verification and limitations

Release preparation checks Python syntax, permitted file types, hard-coded credential patterns, source hashes, and the existing synthetic AIPW formula fixtures. No patient-data reanalysis was run solely to publish this repository. Passing release checks does not independently validate the source databases, clinical phenotypes, or causal assumptions. Figure scripts require outputs from the documented preceding stages; no automatic end-to-end runner or bundled data is claimed.

## Citation and reuse

Until the manuscript has a verified publication record, cite the repository URL together with the exact Git commit used. No publication DOI, invented author list or registration identifier is supplied. Public visibility does not confer access to the governed databases. No open-source license has been selected by the authors for this initial release; public inspection should not be mistaken for an unrestricted reuse license.
