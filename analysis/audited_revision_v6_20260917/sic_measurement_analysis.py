from pathlib import Path
import json
import numpy as np,pandas as pd
O=Path(__file__).resolve().parent;R=O.parent/'audited_reanalysis_v2'
def run_exists(x,n):
 return bool(len(x)>=n and (np.convolve(x.astype(int),np.ones(n,dtype=int),'valid')>=n).any())
def morphology(f0,f1,n):
 b=f0[np.isfinite(f0)]
 if not len(b) or np.isfinite(f1).sum()<48:return 'insufficient'
 low=np.isfinite(f1)&(f1<=b[-1]-10);runs=np.where(np.convolve(low.astype(int),np.ones(n,dtype=int),'valid')>=n)[0]
 if not len(runs):return 'no_qualifying_run'
 tail=f1[runs[0]+n:];rebound=run_exists(np.isfinite(tail)&(tail>=b[-1]-5),n)
 if np.isfinite(f1[-15:]).all() and np.all(f1[-15:]<=b[-1]-10):return 'sustained_to_end'
 if rebound:return 'reduction_with_rebound'
 return 'run_without_confirmed_persistence'
def cls(delta):return 'missing' if not np.isfinite(delta) else 'reduction' if delta<=-10 else 'stable' if -5<=delta<=5 else 'other'
def main():
 keys=pd.read_pickle(O/'sic_signal_keys_private.pkl');ix=pd.Index(keys.key);arr=np.load(O/'sic_candidate_signals_private.npy',mmap_mode='r');last=np.load(O/'sic_source_row_ids_private.npy',mmap_mode='r');flags=pd.read_pickle(O/'sic_quality_flags_private.pkl');bad=set()
 for r in flags.itertuples():
  k=ix.get_indexer([r.case_id*100000+r.hr])[0];j=0 if r.signal==2283 else 1
  if last[k,j]==r.source_row_id:bad.add((int(k),j))
 pd.DataFrame([dict(signal=s,flagged_selected_rows=sum(j==(0 if s==2283 else 1) for k,j in bad)) for s in [2283,710]]).to_csv(O/'sic_selected_quality_flags.csv',index=False)
 h=pd.read_pickle(R/'sicdb_hours.pkl').set_index(['id','hr']);d=pd.read_pickle(O/'sic_candidates_private.pkl');empty=np.full(60,np.nan);rows=[]
 def get(case,hr,j):
  k=ix.get_indexer([int(case)*100000+int(hr)])[0]
  return np.array(arr[k,j],float) if k>=0 and (k,j) not in bad else empty.copy()
 for r in d.itertuples():
  f0=get(r.id,r.hr,0);f1=get(r.id,r.hr+1,0);s0=get(r.id,r.hr,1);b=f0[np.isfinite(f0)];f=f1[np.isfinite(f1)];tail=s0[30:];n=np.isfinite(tail).sum();high=float((tail>=96).sum()/n) if n else np.nan
  nxt=h.loc[(r.id,r.hr+1),'fio2'] if (r.id,r.hr+1) in h.index else np.nan
  hourly=bool(r.fio2>=40 and 96<=r.sat<=100);minute=bool(len(b) and b[-1]>=40 and n>=24 and high>=.8)
  hc=cls(nxt-r.fio2 if 21<=nxt<=100 else np.nan);mc=cls(f[-1]-b[-1] if len(b) and len(f) else np.nan)
  rows.append(dict(id=r.id,pid=r.pid,hr=r.hr,hourly_eligible=hourly,minute_eligible=minute,hour_class=hc,minute_class=mc,survives_landmark=r.survives_landmark,hospital_mortality=r.hospital_mortality,baseline_fio2_last=b[-1] if len(b) else np.nan,next_fio2_last=f[-1] if len(f) else np.nan,baseline_setting_n=len(b),next_setting_n=len(f),baseline_sat_tail_n=n,baseline_sat_high_fraction=high,morph5=morphology(f0,f1,5),morph10=morphology(f0,f1,10)))
 q=pd.DataFrame(rows);q.to_pickle(O/'sic_trajectory_private.pkl')
 flow=[];cohorts={}
 for name in ['hour','minute']:
  eligible=q.hourly_eligible if name=='hour' else q.minute_eligible
  z=q[eligible].copy();flow.append(dict(definition=name,stage='baseline eligible',hours=len(z),cases=z.id.nunique()))
  z=z[z[name+'_class'].ne('missing')&z.survives_landmark];flow.append(dict(definition=name,stage='classifiable and survived',hours=len(z),cases=z.id.nunique()))
  z=z.sort_values(['id','hr']).drop_duplicates('id');flow.append(dict(definition=name,stage='first classifiable',hours=len(z),cases=z.id.nunique()))
  z=z[z[name+'_class'].isin(['reduction','stable'])&z.hospital_mortality.notna()];cohorts[name]=z;z.to_pickle(O/f'sic_{name}_independent_cohort_private.pkl');flow.append(dict(definition=name,stage='binary and hospital status',hours=len(z),cases=z.id.nunique()))
 pd.DataFrame(flow).to_csv(O/'sic_independent_eligibility.csv',index=False)
 old=pd.read_pickle(R/'sicdb_cohort.pkl');a=cohorts['hour'];assert set(a.id)==set(old.id);assert np.array_equal(a.set_index('id').hr.sort_index(),old.set_index('id').hr.sort_index())
 cross=q.groupby(['hourly_eligible','minute_eligible'],dropna=False).agg(hours=('id','size'),cases=('id','nunique')).reset_index();cross.to_csv(O/'sic_eligibility_cross.csv',index=False)
 shared=q[q.hourly_eligible&q.minute_eligible&q.survives_landmark];shared.groupby(['hour_class','minute_class']).size().rename('n').reset_index().to_csv(O/'sic_shared_classification.csv',index=False)
 p=cohorts['hour'][['id','hr','hour_class']].merge(cohorts['minute'][['id','hr','minute_class']],on='id',how='outer',suffixes=('_hour','_minute'),indicator=True)
 out=dict(hour_cases=len(cohorts['hour']),minute_cases=len(cohorts['minute']),common_cases=int(p._merge.eq('both').sum()),hour_only=int(p._merge.eq('left_only').sum()),minute_only=int(p._merge.eq('right_only').sum()),common_same_hour=int((p._merge.eq('both')&p.hr_hour.eq(p.hr_minute)).sum()))
 (O/'sic_independent_overlap.json').write_text(json.dumps(out,indent=2))
 # patient-cluster bootstrap proportions; neither all hours nor repeated stays treated as independent people.
 stats=[];rng=np.random.default_rng(917602)
 configs=[('hour_first',cohorts['hour'],'hour_class'),('minute_first',cohorts['minute'],'minute_class')]
 for pop,z,classfield in configs:
  for treatment in ['reduction','stable']:
   g=z[z[classfield].eq(treatment)].reset_index(drop=True);ids=g.pid.unique();indices=g.groupby('pid').indices
   metrics={'last_value_class_agrees':g.minute_class.eq(g.hour_class).astype(float).to_numpy()}
   for n in [5,10]:
    for cat in ['sustained_to_end','reduction_with_rebound','run_without_confirmed_persistence','no_qualifying_run','insufficient']:metrics[f'morph{n}_{cat}']=g[f'morph{n}'].eq(cat).astype(float).to_numpy()
   bs={k:[] for k in metrics}
   for b in range(1000):
    take=np.concatenate([indices[i] for i in rng.choice(ids,len(ids),replace=True)])
    for k,v in metrics.items():bs[k].append(v[take].mean())
   for k,v in metrics.items():lo,hi=np.quantile(bs[k],[.025,.975]);stats.append(dict(population=pop,treatment=treatment,metric=k,n=len(g),events=int(v.sum()),estimate=v.mean(),ci_low=lo,ci_high=hi))
 pd.DataFrame(stats).to_csv(O/'sic_trajectory_summary.csv',index=False)
 # Quality review sampling is stratified by signal morphology, independent of mortality.
 sample=[]
 for cat,g in q[q.hourly_eligible|q.minute_eligible].groupby('morph5'):
  z=g.sample(min(12,len(g)),random_state=917603);sample.append(z)
 review=pd.concat(sample).sample(frac=1,random_state=917604).reset_index(drop=True);review.insert(0,'review_id',['T%03d'%(i+1) for i in range(len(review))]);review.to_pickle(O/'sic_blinded_review_key_private.pkl')
 pd.DataFrame(dict(review_id=review.review_id,signal_quality='',setting_change='',persistence='',confidence='',assessor='',notes='')).to_csv(O/'clinical_adjudication_blank.csv',index=False)
 print(pd.DataFrame(flow).to_string(index=False));print(out);print('bad selected signal rows',len(bad));print('clinical review sample',len(review))
if __name__=='__main__':main()
