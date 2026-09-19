from analysis_utils import *
d=pd.read_pickle(R/'sicdb_cohort.pkl');q=d[d.n_observed_bins6.eq(6)].copy();num,cat,levels=spec('sicdb');out='low88_any6';res=[];reps=[];fail=[]
configs=[('original_spline_standardization',num,cat,levels,'spline'),('compact_linear_standardization',['fio2','sat','age','icu_hour'],[],{},'linear')]
for j,(name,ns,cs,ls,mode) in enumerate(configs):
 point=standardize(q,out,ns,cs,ls,mode);rng=np.random.default_rng(917200+j);bs=[]
 for b in range(1000):
  try:
   z=standardize(q.iloc[sample_indices(q,rng)],out,ns,cs,ls,mode);bs.append(z);reps.append(dict(analysis=name,replicate=b,**z))
  except (ValueError,ConvergenceWarning) as exc:fail.append(dict(analysis=name,replicate=b,error=str(exc)))
  if (b+1)%200==0:print(name,b+1,flush=True)
 ci=np.quantile([z['rd'] for z in bs],[.025,.975]);res.append(dict(analysis=name,n=len(q),treated_n=int(q.A.sum()),events=int(q[out].sum()),**point,ci_low=ci[0],ci_high=ci[1],bootstrap_success=len(bs),bootstrap_failures=1000-len(bs)))
 pd.DataFrame(res).to_csv(O/'rare_event_results.csv',index=False);pd.DataFrame(reps).to_csv(O/'rare_event_replicates.csv',index=False);pd.DataFrame(fail,columns=['analysis','replicate','error']).to_csv(O/'rare_event_failures.csv',index=False)
print(pd.DataFrame(res).to_string(index=False),flush=True)
