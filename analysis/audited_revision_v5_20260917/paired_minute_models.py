from analysis_utils import *
d=pd.read_pickle(R/'sicdb_cohort.pkl');p=pd.read_pickle(W/'sic_minute_validation_private.pkl')
d=d.merge(p.drop(columns=['pid','A','hospital_mortality']),on='id',validate='one_to_one')
d['A_minute']=d.minute_class.map({'reduction':1,'stable':0})
num,cat,levels=spec('sicdb');results=[];reps=[];fail=[];balance=[]
common=d[d.A_minute.notna()].copy()
configs=[('common_hospital',common,'hospital_mortality'),('qualified_dense_hospital',common[common.baseline30_pass&common.setting0_n.ge(48)&common.setting1_n.ge(48)],'hospital_mortality'),('common_dense_follow_low88',common[common.follow_n.ge(288)],'consecutive5_low88')]
for j,(name,q,out) in enumerate(configs):
 q=q.copy().reset_index(drop=True);q[out]=q[out].astype(int)
 if len(q)<100 or min(q.A.value_counts().min(),q.A_minute.value_counts().min())<20:
  results.append(dict(analysis=name,status='insufficient support',n=len(q)));continue
 point=[];diagnostics=[]
 for field in ['A','A_minute']:
  z=q.copy();z.A=z[field].astype(int)
  point.append(estimate(z,out,'spline',num,cat,levels));dg,bb=diag(z,num,cat,levels);diagnostics.append(dg)
  balance.extend(dict(analysis=name,definition=field,variable=k,smd_after=v) for k,v in bb)
 rng=np.random.default_rng(917100+j);bs=[]
 for b in range(1000):
  ix=sample_indices(q,rng);sample=q.iloc[ix].copy();row={'analysis':name,'replicate':b}
  try:
   for field in ['A','A_minute']:
    z=sample.copy();z.A=z[field].astype(int);v=estimate(z,out,'spline',num,cat,levels)
    for key in ['aipw_rd','risk1','risk0']:row[field+'_'+key]=v[key]
   row['difference']=row['A_minute_aipw_rd']-row['A_aipw_rd'];bs.append(row);reps.append(row)
  except (ValueError,ConvergenceWarning) as exc:fail.append(dict(analysis=name,replicate=b,error=str(exc)))
  if (b+1)%100==0:print(name,b+1,flush=True)
 bb=pd.DataFrame(bs)
 for k,field in enumerate(['A','A_minute']):
  ci=np.quantile(bb[field+'_aipw_rd'],[.025,.975]);results.append(dict(analysis=name,definition=field,outcome=out,**diagnostics[k],events=int(q[out].sum()),rd=point[k]['aipw_rd'],risk1=point[k]['risk1'],risk0=point[k]['risk0'],ci_low=ci[0],ci_high=ci[1],bootstrap_success=len(bs),bootstrap_failures=1000-len(bs),risk_outside_n=int(((bb[field+'_risk1']<0)|(bb[field+'_risk1']>1)|(bb[field+'_risk0']<0)|(bb[field+'_risk0']>1)).sum()),status='complete'))
 ci=np.quantile(bb.difference,[.025,.975]);results.append(dict(analysis=name,definition='paired_difference',outcome=out,n=len(q),rd=point[1]['aipw_rd']-point[0]['aipw_rd'],ci_low=ci[0],ci_high=ci[1],bootstrap_success=len(bs),bootstrap_failures=1000-len(bs),status='complete'))
 pd.DataFrame(results).to_csv(O/'minute_model_results.csv',index=False);pd.DataFrame(reps).to_csv(O/'minute_model_replicates.csv',index=False);pd.DataFrame(balance).to_csv(O/'minute_model_balance.csv',index=False);pd.DataFrame(fail,columns=['analysis','replicate','error']).to_csv(O/'minute_model_failures.csv',index=False)
print(pd.DataFrame(results).to_string(index=False),flush=True)
