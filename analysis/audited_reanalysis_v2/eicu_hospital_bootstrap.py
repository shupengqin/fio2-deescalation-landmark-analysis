"""Hospital-cluster bootstrap sensitivity for eICU primary hospital mortality."""
from pathlib import Path
import numpy as np,pandas as pd
from analyze_verified import estimate
R=Path(__file__).resolve().parent
d=pd.read_pickle(R/'eicu_cohort.pkl').reset_index(drop=True)
num=[c for c in ['fio2','sat','age','icu_hour','peep','rr_vent','rr_vital','mbp','heart_rate','fio2_trend','sat_trend'] if c in d and d[c].notna().any()]
cat=[c for c in ['sex','calendar_group','unit','site'] if c in d]
levels={c:sorted(d[c].fillna('Missing').astype(str).unique()) for c in cat}
point=estimate(d,'hospital_mortality','spline',num,cat,levels)
groups=d.groupby('site',dropna=False).indices
sites=list(groups); rng=np.random.default_rng(26109); vals=[]; fails=0
for b in range(1000):
    ix=np.concatenate([groups[sites[j]] for j in rng.integers(0,len(sites),len(sites))])
    q=d.iloc[ix].reset_index(drop=True)
    try: vals.append(estimate(q,'hospital_mortality','spline',num,cat,levels)['aipw_rd'])
    except Exception: fails+=1
    if (b+1)%100==0: print(b+1,flush=True)
out=pd.DataFrame([dict(database='eICU-CRD',analysis='hospital-cluster bootstrap primary',n=len(d),patients=d.pid.nunique(),hospitals=d.site.nunique(),rd=point['aipw_rd'],ci_low=np.quantile(vals,.025),ci_high=np.quantile(vals,.975),bootstrap_requested=1000,bootstrap_success=len(vals),bootstrap_failures=fails)])
out.to_csv(R/'eicu_hospital_bootstrap.csv',index=False)
print(out.to_string(index=False))
