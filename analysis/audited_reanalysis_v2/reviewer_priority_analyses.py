"""Reviewer-priority analyses: baseline-FiO2 support strata and delayed-opportunity sensitivity."""
from pathlib import Path
import sys, numpy as np, pandas as pd
from analyze_verified import estimate
R=Path(__file__).resolve().parent
DBS=['mimic','eicu','sicdb']
def run(db,B=300):
 d=pd.read_pickle(R/f'{db}_cohort.pkl').reset_index(drop=True)
 num=[c for c in ['fio2','sat','age','icu_hour','peep','rr_vent','rr_vital','mbp','heart_rate','fio2_trend','sat_trend'] if c in d and d[c].notna().any()]
 cat=[c for c in ['sex','calendar_group']+(['unit'] if db!='sicdb' else [])+(['site'] if db=='eicu' else []) if c in d]
 levels={c:sorted(d[c].fillna('Missing').astype(str).unique()) for c in cat}
 masks={
  'baseline_FiO2_40_50':(d.fio2>=40)&(d.fio2<=50),
  'baseline_FiO2_50_70':(d.fio2>50)&(d.fio2<=70),
  'baseline_FiO2_70_90':(d.fio2>70)&(d.fio2<=90),
  'baseline_FiO2_90_100':(d.fio2>90)&(d.fio2<=100),
  'ICU_hour_ge_6':d.icu_hour>=6,
 }
 out=[]
 for name,mask in masks.items():
  q=d[mask].reset_index(drop=True)
  if len(q)<100 or q.A.nunique()<2 or q.hospital_mortality.nunique()<2: continue
  point=estimate(q,'hospital_mortality','spline',num,cat,levels)
  groups=q.groupby('pid',dropna=False).indices;ids=list(groups);rng=np.random.default_rng(99173+len(out));bs=[]
  for b in range(B):
   ix=np.concatenate([groups[ids[j]] for j in rng.integers(0,len(ids),len(ids))])
   try:bs.append(estimate(q.iloc[ix],'hospital_mortality','spline',num,cat,levels)['aipw_rd'])
   except Exception: pass
  out.append(dict(database=db,analysis=name,n=len(q),patients=q.pid.nunique(),treated_n=int(q.A.sum()),deaths=int(q.hospital_mortality.sum()),rd=point['aipw_rd'],ci_low=np.quantile(bs,.025),ci_high=np.quantile(bs,.975),bootstrap_success=len(bs)))
 pd.DataFrame(out).to_csv(R/f'{db}_reviewer_priority.csv',index=False)
 print(db,pd.DataFrame(out).to_string(index=False))
if __name__=='__main__':
 for db in sys.argv[1:] or DBS:run(db)
