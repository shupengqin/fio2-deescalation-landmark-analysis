import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
from pathlib import Path
import sys,json,warnings
import numpy as np,pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
O=Path(__file__).resolve().parent;R=O.parent/'audited_reanalysis_v2';V=O.parent/'audited_revision_v3_20260916';W=O.parent/'audited_revision_v4_20260916'
sys.path.insert(0,str(R))
from analyze_verified import estimate,design,smd
def spec(db):
 s=json.loads((R/f'{db}_model_spec.json').read_text());return s['numeric'],s['categorical'],s['levels']
def diag(d,num,cat,levels):
 z,raw,e,w,dr=estimate(d,'hospital_mortality','spline',num,cat,levels,True);a=d.A.to_numpy();bal=[]
 for c in num:
  x=pd.to_numeric(d[c],errors='coerce').to_numpy(float);bal.append((c,smd(x,a,w)));bal.append((c+' missing',smd((~np.isfinite(x)).astype(float),a,w)))
 for c in cat:
  for lev in levels[c]:bal.append((c+': '+lev,smd(d[c].fillna('Missing').astype(str).eq(lev).to_numpy(float),a,w)))
 out={'n':len(d),'patients':d.pid.nunique(),'treated_n':int(d.A.sum()),'max_abs_smd':float(np.nanmax(np.abs([x[1] for x in bal]))),'ps_clip_n':int(((raw<.01)|(raw>.99)).sum())}
 for aa in [0,1]:
  ww=w[a==aa];out[f'ess_{aa}']=ww.sum()**2/(ww@ww)
 return out,bal
def sample_indices(d,rng):
 groups=d.groupby('pid').indices;ids=list(groups)
 return np.concatenate([groups[ids[k]] for k in rng.integers(0,len(ids),len(ids))])
def standardize(d,out,num,cat,levels,mode='spline'):
 a=d.A.to_numpy(float);y=d[out].to_numpy(float);X=design(d,num,cat,levels,mode)
 with warnings.catch_warnings():
  warnings.simplefilter('error',ConvergenceWarning)
  model=LogisticRegression(C=1,max_iter=2500,tol=1e-7).fit(np.column_stack([a,X]),y)
  r1=model.predict_proba(np.column_stack([np.ones(len(d)),X]))[:,1].mean();r0=model.predict_proba(np.column_stack([np.zeros(len(d)),X]))[:,1].mean()
 return {'rd':r1-r0,'risk1':r1,'risk0':r0}
