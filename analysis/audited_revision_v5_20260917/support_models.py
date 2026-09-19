from analysis_utils import *
d=pd.read_pickle(O/'eicu_timing_support_private.pkl');num,cat,levels=spec('eicu');res=[];reps=[];fails=[];balance=[]
for j,(name,mask) in enumerate([('eicu_both_entries_available',d.baseline_entry_by_boundary&d.classification_entry_by_boundary),('eicu_explicit_vent_interval',d.explicit_vent_interval)]):
 q=d[mask].copy().reset_index(drop=True)
 if len(q)<100 or q.A.nunique()<2 or q.A.value_counts().min()<20:
  res.append(dict(analysis=name,n=len(q),treated_n=int(q.A.sum()),status='insufficient support; no model fitted'));continue
 point=estimate(q,'hospital_mortality','spline',num,cat,levels);dg,bb=diag(q,num,cat,levels);balance.extend(dict(analysis=name,variable=k,smd_after=v) for k,v in bb);bs=[];rng=np.random.default_rng(917300+j)
 for b in range(1000):
  try:
   z=estimate(q.iloc[sample_indices(q,rng)],'hospital_mortality','spline',num,cat,levels);bs.append(z);reps.append(dict(analysis=name,replicate=b,**z))
  except (ValueError,ConvergenceWarning) as exc:fails.append(dict(analysis=name,replicate=b,error=str(exc)))
  if (b+1)%200==0:print(name,b+1,flush=True)
 ci=np.quantile([z['aipw_rd'] for z in bs],[.025,.975]);res.append(dict(analysis=name,**dg,rd=point['aipw_rd'],ci_low=ci[0],ci_high=ci[1],bootstrap_success=len(bs),bootstrap_failures=1000-len(bs),status='complete'))
 pd.DataFrame(res).to_csv(O/'support_model_results.csv',index=False);pd.DataFrame(reps).to_csv(O/'support_model_replicates.csv',index=False);pd.DataFrame(fails,columns=['analysis','replicate','error']).to_csv(O/'support_model_failures.csv',index=False);pd.DataFrame(balance).to_csv(O/'support_model_balance.csv',index=False)
pd.DataFrame(res).to_csv(O/'support_model_results.csv',index=False)
print(pd.DataFrame(res).to_string(index=False),flush=True)
