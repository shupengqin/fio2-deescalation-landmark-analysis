# Recording availability and time aggregation in documented FiO₂ decreases

Analysis code accompanying **Recording Availability and Time Aggregation Shape Documented FiO₂ Decreases Across Three Critical-Care Databases: A Measurement and Selection Analysis**.

This repository contains extraction SQL, cohort construction, statistical estimation, diagnostic and figure-generation code for MIMIC-IV, eICU-CRD and SICdb. It is a code-only release. Repository creation and later code commits are not prospective study registration.

## Scientific scope

The primary analyses examine recording availability, cohort selection and exposure reclassification. The comparison uses the first classifiable opportunity per ICU stay/case: a documented between-bin FiO₂ decrease of at least 10 percentage points versus a change within ±5 points. These labels describe selected recorded changes, not verified bedside treatment strategies.

Minute-level comparisons come only from SICdb, a single-centre source. Comparisons of hourly means, last-valid settings and minute-event rules quantify the size of differences between definitions. They do not validate one definition as clinically correct or establish a clinical gold standard.

MIMIC-IV mortality is retained as an exploratory specification example within the classifiable landmark population. Outcome follow-up begins after the classification interval. eICU-CRD and SICdb mortality estimates are support and model diagnostics because small effective samples and residual imbalance limit interpretation. Source-specific recorded 28-day death is secondary. Treatment comparability, recording selection and measurement definitions are examined separately.

The study is observational and exploratory. Its initial mortality-oriented question was narrowed after recording, support and measurement audits. Historical code and folder names are retained to make that development traceable; the repository slug `fio2-deescalation-landmark-analysis` is retained for stable links. The current manuscript's contribution is the measurement and selection audit. Neither repository naming nor a later revision should be read as prospective registration. Post hoc specifications are not selected according to statistical significance.

## Start here

1. Read [data access and governance](docs/DATA_ACCESS.md).
2. Follow [environment configuration and execution order](docs/REPRODUCING.md).
3. Consult [manuscript-to-code mapping](docs/ANALYSIS_MAP.md).
4. Run the data-free checks: `python tools/validate_release.py` and `python analysis/audited_reanalysis_v2/test_estimand.py`.

The release contains no patient records, source archives, fitted individual predictions, credentials, private signal-review packets, manuscript files or clinical data outputs. Credentialed source access is required for statistical reproduction. The original computation environment is recorded in `analysis/software_versions.json`. Full analyses are computationally intensive and write controlled outputs next to their scripts; execute a copy in a private, unshared working directory.

## Layout

The `analysis/` stage folders must remain siblings because later stages import earlier estimator code. Stage names preserve analysis provenance; they are not competing versions from which to choose a favorable result. The v2 estimator remains the reference. Later stages provide the specified post hoc analyses and figure refinements.

The October 2026 calendar-group extension is in `analysis/posthoc_calendar_groups_20261002/`. It refits the reference spline models within the five existing MIMIC-IV de-identified admission-year groups, with patient-cluster bootstrap intervals and support diagnostics. This is a post hoc descriptive analysis, not an era-interaction test. See that folder's README for execution and interpretation.

Publication preparation scripts that modify author Word documents or Desktop delivery folders are intentionally excluded. The conceptual Figure S5 renderer is retained because it is independent of those author files.

## Verification and limitations

Release preparation checks Python syntax, permitted file types, hard-coded credential patterns, source hashes, and the existing synthetic AIPW formula fixtures. No patient-data reanalysis was run solely to publish this repository. Passing release checks does not independently validate the source databases, clinical phenotypes, or causal assumptions. Figure scripts require outputs from the documented preceding stages; no automatic end-to-end runner or bundled data is claimed.

## Citation and reuse

Until the manuscript has a verified publication record, cite the repository URL together with the exact Git commit used. No publication DOI, invented author list or registration identifier is supplied. Public visibility does not confer access to the governed databases. No open-source license has been selected by the authors for this initial release; public inspection should not be mistaken for an unrestricted reuse license.
