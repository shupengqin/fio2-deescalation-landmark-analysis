# Reproducing the analysis

## Environment

Use a private local copy of this repository. Python 3.14.5 and the statistical versions in `analysis/software_versions.json` were recorded for the study. Install `requirements.txt` into an isolated environment. Auxiliary plotting/database packages are not fully version-pinned; exact binary equivalence across platforms is not promised.

Configure these process environment variables before extraction:

| Variable | Purpose / default |
|---|---|
| `PGHOST` | PostgreSQL host; default `127.0.0.1` |
| `PGPORT` | Port; default `5432`; use your actual installation |
| `PGUSER` | Database user; default `postgres` |
| `PGPASSWORD` | Required local database password; never put a real value in Git |
| `MIMIC_DATABASE` | Local database name; default `mimiciv31` |
| `EICU_DATABASE` | Local database name; default `eicu` |
| `SICDB_DIR` | Directory containing authorized SICdb compressed CSV archives |
| `AGGREGATE_SOURCE_DIR` | Folder containing the three source balance CSVs and corresponding diagnostics JSONs for v7 Figure 2 |

MIMIC queries expect `mimiciv_hosp`, `mimiciv_icu` and the derived tables referenced in SQL, including ventilation and blood-gas definitions. eICU queries expect the schemas/table names in the extraction code. Inspect SQL against your installed release before running. The SICdb signal and disposition codes are source-specific.

Scripts write outputs beside themselves. Do not run them in a directory synchronized publicly. Keep all analysis stage folders as siblings under `analysis/`. Most scripts have top-level execution and are not safe to import as a library. Do not execute all `.py` files indiscriminately.

## Execution order

Run the following from each indicated stage folder. `<db>` means repeat for `mimic`, `eicu` and `sicdb`, except where explicitly limited. The v2 reference and later post hoc analyses have separate outputs; never choose a result by significance.

### Reference analysis: audited_reanalysis_v2

1. `python extract_verified.py <db>` for each source, then `python build_cohorts.py`.
2. `python analyze_verified.py <db> 300` for each source.
3. `python extra_analyses.py <db>` for each source.
4. `python source_exclusion_sensitivity.py extract`, then `python source_exclusion_sensitivity.py <db>`.
5. Run `quantify_residual_confounding.py`, `qa_and_tables.py`, and `test_estimand.py`.
6. Run `reviewer_priority_analyses.py`, `eicu_hospital_bootstrap.py`, and `primary_bootstrap_1000.py`. The last script supplies the reported reference 1,000-refit intervals.

### Recording and history audits: audited_revision_v3_20260916

1. Run `audit_candidates.py`.
2. Run `extract_timestamp_audit.py mimic` and `extract_timestamp_audit.py eicu`.
3. Run `extract_available_labs.py`, `build_revision_audits.py`, and `validate_revision_data.py`.
4. Run `revision_sensitivity.py <db> 1000` for each source.
5. `parallel_mimic_extension.py` is an alternative execution of the MIMIC laboratory/history configuration using its documented seed and draws. Do not count it as an additional independent analysis.
6. Run `aggregate_revision_evidence.py`, then `make_revision_figures.py` if figures are needed.

### Hospital and minute checks: audited_revision_v4_20260916

1. Run `hospital_and_exposure.py` for treatment magnitudes, all hospital omissions and the both-arm hospital analysis.
2. Run `extract_sic_minutes.py`, then `validate_sic_minutes.py`.
3. `make_figures.py` uses generated aggregate outputs plus earlier figure files. Historical numbering is remapped by later stages; do not treat every intermediate figure as a distinct manuscript figure.

### Paired definitions and model checks: audited_revision_v5_20260917

1. Run `extract_timing_support.py`, then `ventilator_corroboration.py`.
2. Run `support_models.py` and `additional_support_models.py`.
3. Run `paired_minute_models.py`, `rare_event_models.py`, and `crossfit_models.py`.
4. Run `minute_uncertainty_signal.py`, then `audit_new_results.py`.
5. Run `make_figures.py` to produce the paired-comparison figure and preserve the earlier measured-hypoxaemia figure as S6.

### Independent eligibility and repeated folds: audited_revision_v6_20260917

1. Run `eicu_regularization.py`. `eicu_split_sensitivity.py` supplies the earlier split diagnostic; its numeric-only balance measure does not replace all-term balance.
2. Run `audit_interval_fields.py`.
3. Run `sic_candidates.py`, then `sic_measurement_analysis.py`.
4. Run `make_figures.py` for Figures S7–S8 from aggregate outputs. The original optional clinical-review rendering is separately retained as `render_clinical_review.py`; it requires locally generated private arrays and must only be run in a governed environment. A blank review form is not completed clinical adjudication.

### Final display refinements

- `audited_revision_v7_20260917/make_figure2.py`: all numeric, categorical and missingness balance terms. Set `AGGREGATE_SOURCE_DIR` to the folder with the generated source diagnostic CSVs expected by this script.
- `editorial_revision_v9_20260918/redraw_s5.py`: standalone conceptual schematic; no patient data required.

## Verification boundary

The release preparation did not rerun extraction or patient-level bootstrap analyses. Only data-free source checks and mathematical fixtures were run. Governed datasets and intermediate outputs are not provided; those are needed to execute the statistical scripts. Bootstrap counts, seeds, estimand and source-specific definitions are retained in the scripts. Check the resulting sample sizes, missingness, support diagnostics, failures and risk bounds before comparing manuscript estimates.
