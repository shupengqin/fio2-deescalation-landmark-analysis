import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import ast,sys,warnings
import numpy as np,pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.exceptions import ConvergenceWarning
O=Path(__file__).resolve().parent;V=O.parent/'audited_revision_v5_20260917';R=O.parent/'audited_reanalysis_v2';sys.path.insert(0,str(V))
t=ast.parse((V/'crossfit_models.py').read_text());t.body=[x for x in t.body if isinstance(x,(ast.Import,ast.ImportFrom,ast.FunctionDef))];ns={};exec(compile(t,'defs','exec'),ns)
num,cat,levels=ns['spec']('eicu')
def predictions(d,seed,mode):
 ids=np.array(sorted(d.pid.unique()));perm=np.random.default_rng(seed).permutation(len(ids));fm={ids[i]:j%5 for j,i in enumerate(perm)};fold=d.pid.map(fm).to_numpy()
 a=d.A.to_numpy(float);y=d.hospital_mortality.to_numpy(float);e=np.empty(len(d));m1=e.copy();m0=e.copy()
 for k in range(5):
  test=fold==k;train=~test
  if not test.any():continue
  assert not set(d.loc[test,'pid']) & set(d.loc[train,'pid'])
  cats=cat if mode=='original' else [c for c in cat if c!='site']
  X,Z=ns['matrices'](d[train],d[test],num,cats,levels)
  if mode!='original':
   # Full, unscaled hospital indicators imply equal ridge coefficient penalties.
   site=levels['site'];X=np.column_stack([X,np.column_stack([d.loc[train,'site'].astype(str).eq(s) for s in site])]);Z=np.column_stack([Z,np.column_stack([d.loc[test,'site'].astype(str).eq(s) for s in site])])
  with warnings.catch_warnings():
   warnings.simplefilter('error',ConvergenceWarning)
   ps=LogisticRegression(C=1,max_iter=2500,tol=1e-7).fit(X,a[train]);om=LogisticRegression(C=1,max_iter=2500,tol=1e-7).fit(np.column_stack([a[train],X]),y[train])
  e[test]=ps.predict_proba(Z)[:,1];m1[test]=om.predict_proba(np.column_stack([np.ones(test.sum()),Z]))[:,1];m0[test]=om.predict_proba(np.column_stack([np.zeros(test.sum()),Z]))[:,1]
 return e,m1,m0
def estimate(d,seeds,mode,details=False):
 raw,m1,m0=np.mean([predictions(d,s,mode) for s in seeds],axis=0);e=np.clip(raw,.01,.99);a=d.A.to_numpy();y=d.hospital_mortality.to_numpy();h1=m1+a/e*(y-m1);h0=m0+(1-a)/(1-e)*(y-m0);w=a/e+(1-a)/(1-e)
 bal=[]
 for c in num:
  x=pd.to_numeric(d[c],errors='coerce').to_numpy(float);bal.extend([(c,ns['smd'](x,a,w)),(c+' missing',ns['smd']((~np.isfinite(x)).astype(float),a,w))])
 for c in cat:
  bal.extend([(c+': '+s,ns['smd'](d[c].fillna('Missing').astype(str).eq(s).to_numpy(float),a,w)) for s in levels[c]])
 wt=w[a==1];r=dict(rd=float((h1-h0).mean()),risk1=h1.mean(),risk0=h0.mean(),max_abs_smd=float(np.nanmax(np.abs([v for k,v in bal]))),treated_ess=wt.sum()**2/(wt@wt),ps_clip_n=int(((raw<.01)|(raw>.99)).sum()),weight_max=w.max(),top10_weight_share=np.sort(w)[-10:].sum()/w.sum())
 if details:return r,pd.DataFrame(bal,columns=['variable','smd']),pd.DataFrame(dict(pid=d.pid,site=d.site,A=a,weight=w,contribution=h1-h0))
 return r
if __name__=='__main__':
 d=pd.read_pickle(R/'eicu_cohort.pkl').reset_index(drop=True);rows=[]
 for mode in ['original','equal_site_penalty']:
  for seed in range(917501,917521):rows.append(dict(mode=mode,seed=seed,**estimate(d,[seed],mode)))
  print(mode,'20 splits complete',flush=True)
 pd.DataFrame(rows).to_csv(O/'eicu_site_penalty_splits.csv',index=False)
 seeds=list(range(917501,917506));point,bal,detail=estimate(d,seeds,'equal_site_penalty',True);bal.to_csv(O/'eicu_repeated_balance.csv',index=False)
 detail.groupby('site').agg(n=('A','size'),reductions=('A','sum'),weight_sum=('weight','sum'),weight_max=('weight','max'),contribution_sum=('contribution','sum')).reset_index().to_csv(O/'eicu_repeated_site_influence.csv',index=False)
 bs=[];failed=[];rng=np.random.default_rng(917601)
 for b in range(300):
  q=d.iloc[ns['sample_indices'](d,rng)].reset_index(drop=True)
  try:bs.append(dict(replicate=b,**estimate(q,seeds,'equal_site_penalty')))
  except (ValueError,ConvergenceWarning) as ex:failed.append(dict(replicate=b,error=str(ex)))
  if (b+1)%25==0:
   pd.DataFrame(bs).to_csv(O/'eicu_repeated_bootstrap.csv',index=False);print('repeat bootstrap',b+1,flush=True)
 z=pd.DataFrame(bs);lo,hi=z.rd.quantile([.025,.975]);point.update(ci_low=lo,ci_high=hi,bootstrap_success=len(z),bootstrap_failures=len(failed),risk_outside_n=int(((z.risk1<0)|(z.risk1>1)|(z.risk0<0)|(z.risk0>1)).sum()))
 pd.DataFrame([point]).to_csv(O/'eicu_repeated_result.csv',index=False);pd.DataFrame(failed,columns=['replicate','error']).to_csv(O/'eicu_repeated_failures.csv',index=False);print(point,flush=True)
