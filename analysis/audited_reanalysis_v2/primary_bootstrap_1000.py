"""1,000-refit patient-cluster bootstrap for each database primary endpoint."""
from pathlib import Path
import sys,numpy as np,pandas as pd
from analyze_verified import estimate
R=Path(__file__).resolve().parent
def run(db):
 d=pd.read_pickle(R/f'{db}_cohort.pkl').reset_index(drop=True)
 num=[c for c in ['fio2','sat','age','icu_hour','peep','rr_vent','rr_vital','mbp','heart_rate','fio2_trend','sat_trend'] if c in d and d[c].notna().any()]
 cat=[c for c in ['sex','calendar_group']+(['unit'] if db!='sicdb' else [])+(['site'] if db=='eicu' else []) if c in d]
 levels={c:sorted(d[c].fillna('Missing').astype(str).unique()) for c in cat}
 point=estimate(d,'hospital_mortality','spline',num,cat,levels)
 groups=d.groupby('pid',dropna=False).indices; ids=list(groups); rng=np.random.default_rng(80317+len(db)); vals=[];fails=0
 for b in range(1000):
  ix=np.concatenate([groups[ids[j]] for j in rng.integers(0,len(ids),len(ids))]);q=d.iloc[ix]
  try: vals.append(estimate(q,'hospital_mortality','spline',num,cat,levels)['aipw_rd'])
  except Exception:fails+=1
  if (b+1)%200==0:print(db,b+1,flush=True)
 out=dict(database=db,analysis='primary hospital mortality, 1000 patient-cluster refits',n=len(d),patients=d.pid.nunique(),treated_n=int(d.A.sum()),deaths=int(d.hospital_mortality.sum()),rd=point['aipw_rd'],ci_low=np.quantile(vals,.025),ci_high=np.quantile(vals,.975),bootstrap_requested=1000,bootstrap_success=len(vals),bootstrap_failures=fails)
 pd.DataFrame([out]).to_csv(R/f'{db}_primary_bootstrap1000.csv',index=False);print(out,flush=True)
for db in sys.argv[1:] or ['mimic','eicu','sicdb']:run(db)
