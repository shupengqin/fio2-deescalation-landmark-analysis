from pathlib import Path
import sys,json,pandas as pd,numpy as np
O=Path(__file__).resolve().parent;R=O.parent/'audited_reanalysis_v2';sys.path.insert(0,str(R))
from analyze_verified import smd
rows=[];bals=[];crude=[];strata=[]
for db in ['mimic','eicu','sicdb']:
 d=pd.read_pickle(R/f'{db}_cohort.pkl');z=pd.read_pickle(R/f'{db}_model_details.pkl');assert np.array_equal(d.id,z.id)
 sp=json.loads((R/f'{db}_model_spec.json').read_text());a=d.A.to_numpy();e=z.ps_raw.to_numpy();w=np.where(a==1,1-e,e)
 b=[]
 for c in sp['numeric']:
  x=pd.to_numeric(d[c]).to_numpy(float);b.append((c,smd(x,a,w)))
  if np.isnan(x).any():b.append((c+' missing',smd(np.isnan(x).astype(float),a,w)))
 for c in sp['categorical']:
  for lev in d[c].fillna('Missing').astype(str).unique():b.append((c+': '+lev,smd(d[c].fillna('Missing').astype(str).eq(lev).to_numpy(float),a,w)))
 bals.extend(dict(database=db,variable=k,smd_overlap=v) for k,v in b)
 for aa in [0,1]:
  ww=w[a==aa];rows.append(dict(database=db,A=aa,n=int((a==aa).sum()),ess=ww.sum()**2/(ww@ww),max_abs_smd=np.nanmax([abs(v) for k,v in b]),top1pct_weight_share=np.sort(ww)[-max(1,int(np.ceil(len(ww)*.01))):].sum()/ww.sum()))
 for out in ['hospital_mortality','death28_recorded']:
  for aa,g in d.groupby('A'):crude.append(dict(database=db,outcome=out,A=aa,n=len(g),deaths=int(g[out].sum()),pct=100*g[out].mean()))
 old=pd.read_csv(R/f'{db}_reviewer_priority.csv')
 for label,lo,hi in [('40_50',40,50),('50_70',50,70),('70_90',70,90),('90_100',90,100)]:
  g=d[d.fio2.ge(lo) if lo==40 else d.fio2.gt(lo)];g=g[g.fio2.le(hi)];match=old[old.analysis.eq('baseline_FiO2_'+label)]
  rr=dict(database=db,analysis=label,n=len(g),treated_n=int(g.A.sum()),stable_n=int((g.A==0).sum()))
  if len(match):rr.update(match.iloc[0][['rd','ci_low','ci_high','bootstrap_success']].to_dict());rr['status']='exploratory, unadjusted CI'
  else:rr['status']='not estimated: original subgroup minimum size rule (n<100 or absent class)'
  strata.append(rr)
pd.DataFrame(rows).to_csv(O/'overlap_full_diagnostics.csv',index=False);pd.DataFrame(bals).to_csv(O/'overlap_all_variable_balance.csv',index=False);pd.DataFrame(crude).to_csv(O/'crude_outcomes.csv',index=False);pd.DataFrame(strata).to_csv(O/'complete_baseline_fio2_strata.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
