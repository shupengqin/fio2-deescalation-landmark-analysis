"""Post hoc MIMIC-IV calendar-group mortality diagnostics.

This script requires controlled-access cohort data. It produces aggregate results,
not patient records. Calendar groups are de-identified admission-year groups,
not exact admission dates. No interaction or confirmatory era effect is tested.
"""
import argparse
import importlib.util
import json
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cohort', type=Path, required=True)
    parser.add_argument('--estimator', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--bootstrap', type=int, default=300)
    parser.add_argument('--seed', type=int, default=20261002)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location('reference_estimator', args.estimator)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    import numpy as np
    import pandas as pd
    import sklearn
    import scipy
    cohort = pd.read_pickle(args.cohort).reset_index(drop=True)
    assert len(cohort) == 18043, 'Unexpected reference-cohort size'
    assert cohort['calendar_group'].notna().all(), 'Calendar-group missingness requires explicit handling'
    num = [c for c in ['fio2', 'sat', 'age', 'icu_hour', 'peep', 'rr_vent', 'rr_vital', 'mbp', 'heart_rate', 'fio2_trend', 'sat_trend']
           if c in cohort and cohort[c].notna().any()]
    cat = ['sex', 'unit']
    levels = {c: sorted(cohort[c].fillna('Missing').astype(str).unique()) for c in cat}
    groups = sorted(cohort['calendar_group'].astype(str).unique())
    # Write the extension plan before inspecting stratum effect estimates.
    # This timestamp does not turn the extension into prospective registration.
    plan = dict(created_utc=datetime.now(timezone.utc).isoformat(), status='post hoc',
                calendar_groups=groups, cohort_rows=len(cohort), numeric=num, categorical=cat,
                categorical_levels=levels, outcome='hospital_mortality', model='reference cubic B-spline AIPW; logistic C=1',
                ps_bounds=[.01,.99], bootstrap_requested=args.bootstrap, seed=args.seed,
                nuisance_refit='imputation, spline knots, scaling, propensity, outcome in each resample',
                clustering='patient within calendar group; groups not assumed independent for comparisons',
                purpose='descriptive mortality and support diagnostics; no interaction testing',
                estimator_sha256=hashlib.sha256(args.estimator.read_bytes()).hexdigest(),
                cohort_sha256=hashlib.sha256(args.cohort.read_bytes()).hexdigest(),
                versions=dict(python=sys.version.split()[0],pandas=pd.__version__,numpy=np.__version__,sklearn=sklearn.__version__,scipy=scipy.__version__))
    (args.output/'analysis_plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')
    # Verify compatibility with the manuscript reference point estimate.
    full_cat = ['sex','calendar_group','unit']
    full_levels = {c:sorted(cohort[c].fillna('Missing').astype(str).unique()) for c in full_cat}
    full = module.estimate(cohort,'hospital_mortality','spline',num,full_cat,full_levels)
    assert abs(full['aipw_rd'] - (-.0149)) < .0001, full
    results, replicas, failures, balance_rows = [], [], [], []
    for index, group in enumerate(groups):
        q = cohort[cohort['calendar_group'].astype(str).eq(group)].reset_index(drop=True)
        point, raw_ps, bounded_ps, weights, _ = module.estimate(q,'hospital_mortality','spline',num,cat,levels,True)
        a = q.A.to_numpy()
        balance = []
        for c in num:
            x = pd.to_numeric(q[c],errors='coerce').to_numpy(float)
            for label,xx in [(c,x),(c+' missing',(~np.isfinite(x)).astype(float))]:
                val = module.smd(xx,a,weights)
                balance.append(val)
                balance_rows.append(dict(calendar_group=group,variable=label,smd_after=val))
        for c in cat:
            for level in levels[c]:
                val = module.smd(q[c].fillna('Missing').astype(str).eq(level).to_numpy(float),a,weights)
                balance.append(val)
                balance_rows.append(dict(calendar_group=group,variable=c+': '+level,smd_after=val))
        clusters = q.groupby('pid',dropna=False).indices
        ids = list(clusters)
        rng = np.random.default_rng(args.seed+index)
        draws = []
        for b in range(args.bootstrap):
            ix = np.concatenate([clusters[ids[j]] for j in rng.integers(0,len(ids),len(ids))])
            try:
                out = module.estimate(q.iloc[ix],'hospital_mortality','spline',num,cat,levels)
                if not all(np.isfinite(v) for v in out.values()): raise ValueError('nonfinite result')
                draws.append(out)
                replicas.append(dict(calendar_group=group,replicate=b,**out))
            except Exception as exc:
                failures.append(dict(calendar_group=group,replicate=b,type=type(exc).__name__,message=str(exc)))
            if (b+1)%50==0:
                print(group,b+1,'/',args.bootstrap,flush=True)
        row = dict(calendar_group=group,n=len(q),patients=q.pid.nunique(),
                   reduction_n=int(q.A.sum()),stable_n=int((q.A==0).sum()),
                   reduction_deaths=int(q.loc[q.A==1,'hospital_mortality'].sum()),
                   stable_deaths=int(q.loc[q.A==0,'hospital_mortality'].sum()),
                   baseline_fio2_100_reduction_pct=100*q.loc[q.A==1,'fio2'].eq(100).mean(),
                   baseline_fio2_100_stable_pct=100*q.loc[q.A==0,'fio2'].eq(100).mean(),
                   icu_hour_reduction_median=q.loc[q.A==1,'icu_hour'].median(),
                   icu_hour_stable_median=q.loc[q.A==0,'icu_hour'].median(),
                   rd=point['aipw_rd'],risk_reduction=point['risk1'],risk_stable=point['risk0'],
                   ps_below_01_n=int((raw_ps<.01).sum()),ps_above_99_n=int((raw_ps>.99).sum()),
                   ess_reduction=float(weights[a==1].sum()**2/(weights[a==1]@weights[a==1])),
                   ess_stable=float(weights[a==0].sum()**2/(weights[a==0]@weights[a==0])),
                   max_abs_smd=float(np.nanmax(np.abs(balance))),
                   bootstrap_requested=args.bootstrap,bootstrap_success=len(draws),bootstrap_failures=args.bootstrap-len(draws),
                   arm_risk_outside_01=sum(any(z[k]<0 or z[k]>1 for k in ['risk1','risk0']) for z in draws),
                   uncertainty_status='complete' if len(draws)>=.95*args.bootstrap else 'failed uncertainty gate')
        if draws:
            row['ci_low'],row['ci_high']=np.quantile([z['aipw_rd'] for z in draws],[.025,.975])
        else:
            row['ci_low']=row['ci_high']=np.nan
        results.append(row)
        pd.DataFrame(results).to_csv(args.output/'mimic_calendar_groups.csv',index=False)
        pd.DataFrame(replicas).to_csv(args.output/'mimic_calendar_group_bootstrap.csv',index=False)
        pd.DataFrame(failures,columns=['calendar_group','replicate','type','message']).to_csv(args.output/'mimic_calendar_group_failures.csv',index=False)
        pd.DataFrame(balance_rows).to_csv(args.output/'mimic_calendar_group_balance.csv',index=False)
        print(group,'complete',json.dumps(row),flush=True)
    assert sum(r['n'] for r in results)==len(cohort)
    assert all(r['uncertainty_status']=='complete' for r in results), 'Some intervals failed completion gate'
    print('All five calendar-group diagnostics complete.',flush=True)


if __name__=='__main__':
    main()
