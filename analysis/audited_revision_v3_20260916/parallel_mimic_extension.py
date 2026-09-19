import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
from pathlib import Path
import sys,json
import numpy as np,pandas as pd
from concurrent.futures import ProcessPoolExecutor,as_completed
O=Path(__file__).resolve().parent;R=O.parent/'audited_reanalysis_v2';sys.path.insert(0,str(R))
from analyze_verified import estimate
from revision_sensitivity import diag

def initialize():
 global d,nums,cats,levels,groups,ids
 d=pd.read_pickle(O/'mimic_augmented_cohort.pkl').reset_index(drop=True)
 spec=json.loads((O/'mimic_available_history_and_blood_gas_model.json').read_text());nums=spec['numeric'];cats=spec['categorical'];levels={c:sorted(d[c].fillna('Missing').astype(str).unique()) for c in cats};groups=d.groupby('pid').indices;ids=list(groups)
def batch(items):
 out=[];fails=[]
 for b,chosen in items:
  ix=np.concatenate([groups[ids[k]] for k in chosen])
  try:
   z=estimate(d.iloc[ix],'hospital_mortality','spline',nums,cats,levels);out.append(dict(database='mimic',analysis='available_history_and_blood_gas',replicate=b,rd=z['aipw_rd'],risk1=z['risk1'],risk0=z['risk0']))
  except Exception as exc:fails.append(dict(database='mimic',analysis='available_history_and_blood_gas',replicate=b,error=str(exc)))
 return out,fails
if __name__=='__main__':
 initialize();name='available_history_and_blood_gas';rng=np.random.default_rng(916264)
 # Same 1,000 draws as the serial script, with execution order independent of draw generation.
 jobs=[(b,rng.integers(0,len(ids),len(ids))) for b in range(1000)]
 point,raw,e,w,dr=estimate(d,'hospital_mortality','spline',nums,cats,levels,True);dg,bb=diag(d,nums,cats,raw,w)
 outs=[];fails=[]
 with ProcessPoolExecutor(max_workers=6,initializer=initialize) as ex:
  fs=[ex.submit(batch,jobs[i:i+10]) for i in range(0,1000,10)]
  for f in as_completed(fs):
   z,err=f.result();outs+=z;fails+=err
   pd.DataFrame(outs).sort_values('replicate').to_csv(O/'mimic_extended_parallel_checkpoint.csv',index=False)
   if (len(outs)+len(fails))%50==0:print('extended complete',len(outs),'failed',len(fails),flush=True)
 bs=[z['rd'] for z in outs]
 rr=dict(database='mimic',analysis=name,n=len(d),patients=d.pid.nunique(),treated_n=int(d.A.sum()),deaths=int(d.hospital_mortality.sum()),**point,**dg,ci_low=np.quantile(bs,.025),ci_high=np.quantile(bs,.975),bootstrap_requested=1000,bootstrap_success=len(bs),status='complete' if len(bs)>=950 else 'incomplete uncertainty')
 for filename,extra in [('mimic_revision_results.csv',[rr]),('mimic_revision_balance.csv',[dict(database='mimic',analysis=name,variable=k,smd_after=v) for k,v in bb]),('mimic_revision_replicates.csv',outs),('mimic_revision_failures.csv',fails)]:
  prev=pd.read_csv(O/filename);prev=prev[prev.analysis.ne(name)];pd.concat([prev,pd.DataFrame(extra)],ignore_index=True).to_csv(O/filename,index=False)
 print(json.dumps(rr),flush=True)
