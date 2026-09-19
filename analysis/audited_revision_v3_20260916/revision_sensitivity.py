import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
from pathlib import Path
import sys,json,warnings
import numpy as np,pandas as pd
from sklearn.exceptions import ConvergenceWarning
O=Path(__file__).resolve().parent;R=O.parent/'audited_reanalysis_v2';sys.path.insert(0,str(R))
from analyze_verified import estimate,smd

def diag(q,nums,cats,eraw,w):
 a=q.A.to_numpy();bb=[]
 for c in nums:
  x=pd.to_numeric(q[c],errors='coerce').to_numpy(float);bb.append((c,smd(x,a,w)))
  if np.isnan(x).any():bb.append((c+' missing',smd(np.isnan(x).astype(float),a,w)))
 for c in cats:
  for lev in q[c].fillna('Missing').astype(str).unique():bb.append((c+': '+lev,smd(q[c].fillna('Missing').astype(str).eq(lev).to_numpy(float),a,w)))
 out={'max_abs_smd':float(np.nanmax([abs(x[1]) for x in bb])),'ps_lt01_n':int((eraw<.01).sum()),'ps_gt99_n':int((eraw>.99).sum())}
 for aa in [0,1]:
  ww=w[a==aa];out[f'ess_{aa}']=float(ww.sum()**2/(ww@ww));out[f'top1pct_weight_share_{aa}']=float(np.sort(ww)[-max(1,int(np.ceil(len(ww)*.01))):].sum()/ww.sum())
 return out,bb

def run(db,B=1000):
 d=pd.read_pickle(O/f'{db}_augmented_cohort.pkl').reset_index(drop=True)
 core=[c for c in ['fio2','sat','age','icu_hour','peep','rr_vent','rr_vital','mbp','heart_rate'] if c in d and d[c].notna().any()]
 cats=[c for c in ['sex','calendar_group']+(['unit'] if db!='sicdb' else [])+(['site'] if db=='eicu' else []) if c in d]
 configs=[('core_no_exact_trends',np.ones(len(d),bool),core),('clinical40_70_established6h',d.fio2.between(40,70)&d.hours_since_support.ge(6),core),('first_physiological_opportunity',d.first_to_selected_delay_hr.eq(0),core)]
 if db!='sicdb':
  configs += [('timestamp30_90_alignment15',d.setting_gap_min.between(30,90)&d.baseline_alignment_min.le(15),core),('available_history_and_blood_gas',np.ones(len(d),bool),core+['history_change_per_hr','history_age_hr','po2','pco2','ph','lactate','po2_age_hr','pco2_age_hr','ph_age_hr','lactate_age_hr'])]
 out=[];balances=[];reps=[];failures=[]
 for j,(name,mask,nums) in enumerate(configs):
  q=d[mask].reset_index(drop=True);nums=[c for c in nums if q[c].notna().any()];levels={c:sorted(q[c].fillna('Missing').astype(str).unique()) for c in cats}
  meta=dict(database=db,analysis=name,n=len(q),patients=q.pid.nunique(),treated_n=int(q.A.sum()),deaths=int(q.hospital_mortality.sum()))
  if len(q)<100 or q.A.value_counts().min()<20 or q.hospital_mortality.nunique()<2:
   out.append(dict(**meta,status='insufficient support; no adjusted estimate',bootstrap_requested=0,bootstrap_success=0));continue
  point,raw,e,w,dr=estimate(q,'hospital_mortality','spline',nums,cats,levels,True);dg,bb=diag(q,nums,cats,raw,w)
  balances.extend(dict(database=db,analysis=name,variable=k,smd_after=v) for k,v in bb)
  groups=q.groupby('pid').indices;ids=list(groups);rng=np.random.default_rng(916260+j+100*['mimic','eicu','sicdb'].index(db));bs=[]
  (O/f'{db}_{name}_model.json').write_text(json.dumps(dict(numeric=nums,categorical=cats,bootstrap_requested=B,seed=916260+j+100*['mimic','eicu','sicdb'].index(db)),indent=2))
  print(db,name,'n',len(q),'treated',int(q.A.sum()),'starting',flush=True)
  for i in range(B):
   ix=np.concatenate([groups[ids[k]] for k in rng.integers(0,len(ids),len(ids))])
   try:
    z=estimate(q.iloc[ix],'hospital_mortality','spline',nums,cats,levels);bs.append(z['aipw_rd']);reps.append(dict(database=db,analysis=name,replicate=i,rd=z['aipw_rd'],risk1=z['risk1'],risk0=z['risk0']))
   except (ValueError,ConvergenceWarning) as exc:failures.append(dict(database=db,analysis=name,replicate=i,error=str(exc)))
   if (i+1)%250==0:print(db,name,i+1,flush=True)
  out.append(dict(**meta,**point,**dg,ci_low=np.quantile(bs,.025),ci_high=np.quantile(bs,.975),bootstrap_requested=B,bootstrap_success=len(bs),status='complete' if len(bs)>=.95*B else 'incomplete uncertainty'))
  pd.DataFrame(out).to_csv(O/f'{db}_revision_results.csv',index=False)
  pd.DataFrame(balances).to_csv(O/f'{db}_revision_balance.csv',index=False)
  pd.DataFrame(reps).to_csv(O/f'{db}_revision_replicates.csv',index=False)
  pd.DataFrame(failures,columns=['database','analysis','replicate','error']).to_csv(O/f'{db}_revision_failures.csv',index=False)
  print(db,name,'RD',point['aipw_rd'],'CI',out[-1]['ci_low'],out[-1]['ci_high'],'SMD',dg['max_abs_smd'],flush=True)
 pd.DataFrame(out).to_csv(O/f'{db}_revision_results.csv',index=False)
if __name__=='__main__':run(sys.argv[1],int(sys.argv[2]) if len(sys.argv)>2 else 1000)
