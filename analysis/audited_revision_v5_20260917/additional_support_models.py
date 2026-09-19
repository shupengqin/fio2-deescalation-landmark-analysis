from analysis_utils import *
res=[];reps=[];fails=[];balances=[]
for j,db in enumerate(['mimic','eicu']):
 d=pd.read_pickle(O/('mimic_timing_private.pkl' if db=='mimic' else 'eicu_corroboration_private.pkl'));num,cat,levels=spec(db)
 mask=d.baseline_entry_by_boundary&d.classification_entry_by_boundary if db=='mimic' else d.corroborating_settings
 name=db+('_both_entries_available' if db=='mimic' else '_corroborating_settings');q=d[mask].copy().reset_index(drop=True)
 if len(q)<100 or q.A.nunique()<2 or q.A.value_counts().min()<20:
  res.append(dict(analysis=name,n=len(q),treated_n=int(q.A.sum()),status='insufficient support'));continue
 point=estimate(q,'hospital_mortality','spline',num,cat,levels);dg,bb=diag(q,num,cat,levels);balances.extend(dict(analysis=name,variable=k,smd_after=v) for k,v in bb);bs=[];rng=np.random.default_rng(917350+j)
 for b in range(1000):
  try:z=estimate(q.iloc[sample_indices(q,rng)],'hospital_mortality','spline',num,cat,levels);bs.append(z);reps.append(dict(analysis=name,replicate=b,**z))
  except (ValueError,ConvergenceWarning) as exc:fails.append(dict(analysis=name,replicate=b,error=str(exc)))
  if (b+1)%100==0:print(name,b+1,flush=True)
 ci=np.quantile([z['aipw_rd'] for z in bs],[.025,.975]);res.append(dict(analysis=name,**dg,rd=point['aipw_rd'],ci_low=ci[0],ci_high=ci[1],bootstrap_success=len(bs),bootstrap_failures=1000-len(bs),status='complete'))
 pd.DataFrame(res).to_csv(O/'additional_support_results.csv',index=False);pd.DataFrame(reps).to_csv(O/'additional_support_replicates.csv',index=False);pd.DataFrame(fails,columns=['analysis','replicate','error']).to_csv(O/'additional_support_failures.csv',index=False);pd.DataFrame(balances).to_csv(O/'additional_support_balance.csv',index=False)
print(pd.DataFrame(res).to_string(index=False),flush=True)
