from analysis_utils import *
d=pd.read_pickle(W/'sic_minute_validation_private.pkl').reset_index(drop=True)
out=[];rep=[]
def summarize(q,value,name,a,coverage):
 q=q.copy();q['value']=value
 g=q.groupby('pid').agg(total=('value','sum'),n=('value','size'))
 rng=np.random.default_rng(917620+len(out));bs=[]
 for b in range(1000):
  x=g.iloc[rng.integers(0,len(g),len(g))];v=x.total.sum()/x.n.sum();bs.append(v)
  rep.append(dict(measure=name,A=a,coverage=coverage,replicate=b,value=v))
 ci=np.quantile(bs,[.025,.975]);out.append(dict(measure=name,A=a,coverage=coverage,n=len(q),patients=len(g),estimate=q.value.mean(),ci_low=ci[0],ci_high=ci[1],bootstrap_success=len(bs)))
for a,g in d.groupby('A'):
 expected='reduction' if a else 'stable'
 for coverage,q in [('paired',g[g.minute_class.ne('missing')]),('dense_pair',g[g.minute_class.ne('missing')&g.setting0_n.ge(48)&g.setting1_n.ge(48)])]:
  summarize(q,q.minute_class.eq(expected).astype(float),'same_treatment_class',a,coverage)
 for coverage,q in [('any_follow',g[g.follow_n.gt(0)]),('dense_follow',g[g.follow_n.ge(288)])]:
  for k in [88,90]:
   qk=q[q[f'hourly_low{k}'].notna()]
   for col in [f'minute_low{k}',f'consecutive5_low{k}']:
    summarize(qk,qk[col].astype(float)-qk[f'hourly_low{k}'].astype(float),col+'_minus_hourly',a,coverage)
pd.DataFrame(out).to_csv(O/'minute_uncertainty_results.csv',index=False)
pd.DataFrame(rep).to_csv(O/'minute_uncertainty_replicates.csv',index=False)
# Numeric flags describe potentially abrupt traces; they do not adjudicate artifacts.
raw=pd.read_pickle(W/'sic_embedded_raw_selected.pkl');signals={}
for r in raw[raw.DataID.eq(710)].itertuples():
 if pd.isna(r.rawdata):continue
 s=str(r.rawdata);b=bytes.fromhex(s[2:] if s.startswith('0x') else s)
 x=np.frombuffer(b,dtype='<f4').astype(float);present=np.frombuffer(b,dtype='<u4')!=0
 x[~present|~np.isfinite(x)|(x<1)|(x>100)]=np.nan
 signals[(r.CaseID,int(r.Offset/3600))]=x
c=pd.read_pickle(R/'sicdb_cohort.pkl').set_index('id');empty=np.full(60,np.nan);rows=[];traces=[]
for r in d.itertuples():
 h=int(c.loc[r.id,'hr']);x=np.concatenate([signals.get((r.id,h+j),empty) for j in range(2,8)])
 adj=np.isfinite(x[1:])&np.isfinite(x[:-1]);jump=adj&(np.abs(np.diff(x))>=10)
 isolated=(x[1:-1]<88)&(x[:-2]>=94)&(x[2:]>=94)
 rows.append(dict(id=r.id,pid=r.pid,A=r.A,follow_n=r.follow_n,jump10=bool(jump.any()),isolated88=bool(isolated.any()),jump_count=int(jump.sum()),isolated_count=int(isolated.sum())))
 if len(traces)<6 and isolated.any():traces.append(x)
pd.DataFrame(rows).to_pickle(O/'minute_signal_flags_private.pkl')
flags=pd.DataFrame(rows);agg=[]
for a,q in flags.groupby('A'):
 for cov,g in [('any_follow',q[q.follow_n.gt(0)]),('dense_follow',q[q.follow_n.ge(288)])]:
  agg.append(dict(A=a,coverage=cov,n=len(g),jump10_cases=int(g.jump10.sum()),isolated88_cases=int(g.isolated88.sum()),jump10_pairs=int(g.jump_count.sum()),isolated88_minutes=int(g.isolated_count.sum())))
pd.DataFrame(agg).to_csv(O/'minute_signal_flags.csv',index=False)
np.save(O/'flagged_traces_private.npy',np.array(traces))
print(pd.DataFrame(out).to_string(index=False));print(pd.DataFrame(agg).to_string(index=False))
