from analysis_utils import *
from sklearn.preprocessing import SplineTransformer,StandardScaler
def matrices(train,test,num,cat,levels):
 tr=[];te=[]
 for c in num:
  x=pd.to_numeric(train[c],errors='coerce').to_numpy(float);u=pd.to_numeric(test[c],errors='coerce').to_numpy(float)
  bad=~np.isfinite(x);ubad=~np.isfinite(u);med=np.median(x[~bad]) if (~bad).any() else 0
  x=np.where(bad,med,x);u=np.where(ubad,med,u)
  if c in ['fio2','sat','age','icu_hour'] and len(np.unique(x))>=5:
   s=SplineTransformer(n_knots=4,degree=3,knots='quantile',extrapolation='linear',include_bias=False);xt=s.fit_transform(x[:,None]);ut=s.transform(u[:,None])
  else:xt=x[:,None];ut=u[:,None]
  tr.extend([xt,bad.astype(float)[:,None]]);te.extend([ut,ubad.astype(float)[:,None]])
 for c in cat:
  lev=levels[c][1:]
  for target,df in [(tr,train),(te,test)]:
   z=df[c].fillna('Missing').astype(str).to_numpy();target.append(np.column_stack([z==k for k in lev]).astype(float) if lev else np.zeros((len(df),1)))
 scale=StandardScaler();return scale.fit_transform(np.column_stack(tr)),scale.transform(np.column_stack(te))
def fit(d,num,cat,levels):
 n=len(d);a=d.A.to_numpy(float);y=d.hospital_mortality.to_numpy(float);e=np.empty(n);m1=np.empty(n);m0=np.empty(n)
 for k in range(5):
  test=d.fold.eq(k).to_numpy();train=~test
  if not test.any():continue
  # Resampled copies of a patient cannot cross between training and held-out fold.
  assert not set(d.loc[test,'pid'])&set(d.loc[train,'pid'])
  X,Z=matrices(d[train],d[test],num,cat,levels)
  with warnings.catch_warnings():
   warnings.simplefilter('error',ConvergenceWarning)
   ps=LogisticRegression(C=1,max_iter=2500,tol=1e-7).fit(X,a[train]);e[test]=ps.predict_proba(Z)[:,1]
   om=LogisticRegression(C=1,max_iter=2500,tol=1e-7).fit(np.column_stack([a[train],X]),y[train])
   m1[test]=om.predict_proba(np.column_stack([np.ones(test.sum()),Z]))[:,1];m0[test]=om.predict_proba(np.column_stack([np.zeros(test.sum()),Z]))[:,1]
 raw=e.copy();e=np.clip(e,.01,.99);r1=(m1+a/e*(y-m1)).mean();r0=(m0+(1-a)/(1-e)*(y-m0)).mean()
 return dict(rd=r1-r0,risk1=r1,risk0=r0,ps_clip_n=int(((raw<.01)|(raw>.99)).sum())),a/e+(1-a)/(1-e)
for db in sys.argv[1:] or ['mimic','eicu','sicdb']:
 d=pd.read_pickle(R/f'{db}_cohort.pkl').reset_index(drop=True);num,cat,levels=spec(db)
 ids=np.array(sorted(d.pid.unique()));rng=np.random.default_rng(917501);perm=rng.permutation(len(ids));foldmap={ids[i]:int(j%5) for j,i in enumerate(perm)};d['fold']=d.pid.map(foldmap)
 point,w=fit(d,num,cat,levels);bal=[];a=d.A.to_numpy()
 for c in num:
  x=pd.to_numeric(d[c],errors='coerce').to_numpy(float);bal.append((c,smd(x,a,w)));bal.append((c+' missing',smd((~np.isfinite(x)).astype(float),a,w)))
 for c in cat:
  for k in levels[c]:bal.append((c+': '+k,smd(d[c].fillna('Missing').astype(str).eq(k).to_numpy(float),a,w)))
 bs=[];fail=[];rng=np.random.default_rng(917510+['mimic','eicu','sicdb'].index(db))
 for b in range(300):
  try:z,_=fit(d.iloc[sample_indices(d,rng)].reset_index(drop=True),num,cat,levels);bs.append(dict(replicate=b,**z))
  except (ValueError,ConvergenceWarning) as exc:fail.append(dict(replicate=b,error=str(exc)))
  if (b+1)%25==0:print(db,'crossfit',b+1,flush=True)
  pd.DataFrame(bs).to_csv(O/f'{db}_crossfit_replicates.csv',index=False)
 ci=np.quantile([x['rd'] for x in bs],[.025,.975]);ww=w[a==1]
 result=dict(database=db,n=len(d),treated_n=int(d.A.sum()),**point,ci_low=ci[0],ci_high=ci[1],max_abs_smd=float(np.nanmax(np.abs([x[1] for x in bal]))),treated_ess=ww.sum()**2/(ww@ww),bootstrap_requested=300,bootstrap_success=len(bs),bootstrap_failures=len(fail))
 pd.DataFrame([result]).to_csv(O/f'{db}_crossfit_result.csv',index=False);pd.DataFrame(bal,columns=['variable','smd_after']).to_csv(O/f'{db}_crossfit_balance.csv',index=False);pd.DataFrame(fail,columns=['replicate','error']).to_csv(O/f'{db}_crossfit_failures.csv',index=False)
 print(result,flush=True)
