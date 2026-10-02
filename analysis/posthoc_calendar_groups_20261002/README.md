# Post hoc MIMIC-IV calendar-group analysis

This extension uses the already constructed first-classifiable landmark cohort. It stratifies the same records by the five existing de-identified calendar groups: 2008–2010, 2011–2013, 2014–2016, 2017–2019 and 2020–2022. These are mapped source groups, not exact admission dates. It introduces no new extraction, exposure threshold or outcome definition.

Run in a private copy with controlled-access cohort data:

```text
python calendar_group_analysis.py --cohort ../audited_reanalysis_v2/mimic_cohort.pkl --estimator ../audited_reanalysis_v2/analyze_verified.py --output PRIVATE_OUTPUT_DIRECTORY --bootstrap 300 --seed 20261002
```

The estimator refits the existing spline AIPW treatment and outcome models within each group. Calendar group is constant within strata and is removed from the feature matrix. Sex and unit indicators retain the reference cohort's category vocabulary. Each patient-cluster resample refits imputation, spline knots, scaling and both nuisance models. Probabilities are bounded at 0.01 and 0.99 as in the reference analysis.

Outputs are aggregate group counts, observed mortality, reference-model risk differences and percentile intervals, effective samples, propensity-bound counts, all-term residual balance, bootstrap convergence and arm-risk diagnostics. Replicate risk estimates are not clipped or selectively discarded for statistical significance.

The script writes an extension plan before estimating group results. That record documents this execution; it is not retrospective evidence of prospective registration. The analysis is exploratory, does not test an era interaction and does not establish absence of temporal confounding. Patients contributing stays to more than one calendar group mean that group results cannot be assumed independent for cross-group tests. No pooled re-estimate or revised primary estimand is created.

The public repository contains code only. Governed cohort files and generated private outputs must not be committed.
