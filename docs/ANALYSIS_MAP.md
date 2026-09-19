# Analysis and manuscript mapping

| Analysis | Stage and key files |
|---|---|
| Source-specific eligibility, first classifiable opportunity, landmark and measured outcomes | v2 `extract_verified.py`, SQL files, `build_cohorts.py` |
| Reference spline AIPW, IPTW, propensity diagnostics and patient bootstrap | v2 `analyze_verified.py`, `primary_bootstrap_1000.py` |
| Linear/overlap/support and repeated-opportunity analyses | v2 `extra_analyses.py`, `reviewer_priority_analyses.py` |
| ECMO, support initiation and treatment limitation restrictions | v2 `source_exclusion_sensitivity.py` |
| E-values | v2 `quantify_residual_confounding.py` |
| First-physiological eligibility, recording availability, timing, prior history and labs | v3 extraction/audit scripts, `revision_sensitivity.py` |
| Hospital influence and treatment magnitude | v4 `hospital_and_exposure.py` |
| SICdb minute reconstruction | v4 `extract_sic_minutes.py`, `validate_sic_minutes.py` |
| Paired hourly/minute mortality and hypoxaemia comparisons | v5 `paired_minute_models.py`, `minute_uncertainty_signal.py` |
| Electronic-entry timing and eICU ventilator corroboration | v5 extraction/support scripts |
| Cross-fitting and sparse outcomes | v5 `crossfit_models.py`, `rare_event_models.py` |
| eICU repeated-fold/site-penalty diagnostics | v6 `eicu_regularization.py`, `eicu_split_sensitivity.py` |
| Independent SICdb eligibility and trajectory morphology | v6 `sic_candidates.py`, `sic_measurement_analysis.py` |
| All-term reference balance main figure | v7 `make_figure2.py` |
| Conceptual confounding/selection schematic | editorial v9 `redraw_s5.py` |

Full-cohort mortality is the reference endpoint in the revised series. Earlier exploratory work used recorded 28-day death; the change is disclosed in manuscript Methods and the analysis-development table. Folder timestamps are not evidence of preregistration. The separate private review renderer creates a local clinical signal packet, not completed clinical adjudication.

`source_provenance.json` maps public scripts to hashes of the retained author analysis package. Portability edits change input configuration, not estimator functions, thresholds, seeds or bootstrap counts. The optional private renderer was separated from the v6 aggregate renderer. Document-production and literature-audit programs are not statistical dependencies and were omitted.
