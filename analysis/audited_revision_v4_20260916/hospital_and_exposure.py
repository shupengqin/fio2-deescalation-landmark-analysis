import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import sys,json,warnings
import numpy as np,pandas as pd
from sklearn.exceptions import ConvergenceWarning
O=Path(__file__).resolve().parent;R=O.parent/'audited_reanalysis_v2';sys.path.insert(0,str(R))
from analyze_verified import estimate,diagnostics,smd
rows=[];joint=[]
for db in ['mimic','eicu','sicdb']:
 d=pd.read_pickle(R/f'{db}_cohort.pkl')
 for a,g in d.groupby('A'):
  rows.append(dict(database=db,A=int(a),n=len(g),baseline100=int(g.fio2.eq(100).sum()),next_median=g.fio2_next.median(),delta_q25=g.delta.quantile(.25),delta_median=g.delta.median(),delta_q75=g.delta.quantile(.75),drop10_20=int((g.delta.le(-10)&g.delta.gt(-20)).sum()),drop20_30=int((g.delta.le(-20)&g.delta.gt(-30)).sum()),drop30plus=int(g.delta.le(-30).sum())))
  counts,x,y=np.histogram2d(g.fio2,g.fio2_next,bins=[np.arange(20,111,10),np.arange(20,111,10)])
  for i,j in np.ndindex(counts.shape):joint.append(dict(database=db,A=int(a),baseline_left=x[i],next_left=y[j],n=int(counts[i,j])))
pd.DataFrame(rows).to_csv(O/'treatment_versions.csv',index=False);pd.DataFrame(joint).to_csv(O/'treatment_joint.csv',index=False)
d=pd.read_pickle(R/'eicu_cohort.pkl').reset_index(drop=True);spec=json.loads((R/'eicu_model_spec.json').read_text());num=spec['numeric'];cat=spec['categorical'];levels=spec['levels']
s=d.groupby('site').agg(n=('A','size'),reduction_n=('A','sum'),deaths=('hospital_mortality','sum'));s['stable_n']=s.n-s.reduction_n;s['reduction_pct']=100*s.reduction_n/s.n
for a,name in [(1,'reduction'),(0,'stable')]:s[name+'_deaths']=d[d.A.eq(a)].groupby('site').hospital_mortality.sum().reindex(s.index,fill_value=0)
s=s.sort_index();s['hospital_label']=['H'+str(i+1).zfill(2) for i in range(len(s))]
s.drop(columns=[],errors='ignore').reset_index(drop=True).to_csv(O/'hospital_support.csv',index=False)
base=estimate(d,'hospital_mortality','spline',num,cat,levels)['aipw_rd'];loo=[]
for site,r in s.iterrows():
 q=d[d.site.ne(site)];z=estimate(q,'hospital_mortality','spline',num,cat,levels)
 loo.append(dict(hospital_label=r.hospital_label,excluded_n=int(r.n),excluded_reduction_n=int(r.reduction_n),remaining_n=len(q),rd=z['aipw_rd'],change_from_full=z['aipw_rd']-base))
pd.DataFrame(loo).to_csv(O/'hospital_leave_one_out.csv',index=False)
keep=s[(s.reduction_n>0)&(s.stable_n>0)].index;q=d[d.site.isin(keep)].reset_index(drop=True)
point,raw,e,w,dr=estimate(q,'hospital_mortality','spline',num,cat,levels,True)
bb=[]
for c in num:
 x=pd.to_numeric(q[c],errors='coerce').to_numpy(float);bb.append(smd(x,q.A.to_numpy(),w));bb.append(smd(np.isnan(x).astype(float),q.A.to_numpy(),w))
for c in cat:
 for lev in levels[c]:bb.append(smd(q[c].fillna('Missing').astype(str).eq(lev).to_numpy(float),q.A.to_numpy(),w))
groups=q.groupby('site').indices;ids=list(groups);rng=np.random.default_rng(916401);bs=[];fail=[]
for b in range(1000):
 ix=np.concatenate([groups[ids[j]] for j in rng.integers(0,len(ids),len(ids))])
 try:
  with warnings.catch_warnings():
   warnings.simplefilter('error',ConvergenceWarning);z=estimate(q.iloc[ix],'hospital_mortality','spline',num,cat,levels)
  bs.append(dict(replicate=b,**z))
 except (ValueError,ConvergenceWarning) as exc:fail.append(dict(replicate=b,error=str(exc)))
 if (b+1)%100==0:print('both-arm hospital bootstrap',b+1,flush=True)
pd.DataFrame(bs).to_csv(O/'hospital_both_arms_replicates.csv',index=False);pd.DataFrame(fail,columns=['replicate','error']).to_csv(O/'hospital_both_arms_failures.csv',index=False)
ww=w[q.A.to_numpy()==1]
result=dict(n=len(q),patients=q.pid.nunique(),hospitals=q.site.nunique(),treated_n=int(q.A.sum()),rd=point['aipw_rd'],ci_low=np.quantile([z['aipw_rd'] for z in bs],.025),ci_high=np.quantile([z['aipw_rd'] for z in bs],.975),max_abs_smd=np.nanmax(np.abs(bb)),treated_ess=ww.sum()**2/(ww@ww),bootstrap_requested=1000,bootstrap_success=len(bs),bootstrap_failures=len(fail))
pd.DataFrame([result]).to_csv(O/'hospital_both_arms_result.csv',index=False)
print(json.dumps(result),flush=True)
