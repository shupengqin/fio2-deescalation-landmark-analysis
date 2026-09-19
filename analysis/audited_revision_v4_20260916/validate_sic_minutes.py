from pathlib import Path
import pandas as pd,numpy as np,json
O=Path(__file__).resolve().parent;R=O.parent/'audited_reanalysis_v2'
d=pd.read_pickle(R/'sicdb_cohort.pkl').set_index('id');z=pd.read_pickle(O/'sic_embedded_raw_selected.pkl')
signals={};qa=[];errors=[];unavailable=[]
for r in z.itertuples():
 if pd.isna(r.rawdata):
  unavailable.append(dict(signal=int(r.DataID),declared_count=int(r.cnt)));continue
 try:
  raw=str(r.rawdata);b=bytes.fromhex(raw[2:] if raw.startswith('0x') else raw)
  if len(b)!=240:raise ValueError('raw length '+str(len(b)))
  x=np.frombuffer(b,dtype='<f4').astype(float);present=np.frombuffer(b,dtype='<u4')!=0;x[~present]=np.nan
  qa.append(dict(signal=int(r.DataID),declared_count=r.cnt,decoded_count=int(present.sum()),count_match=int(present.sum())==int(r.cnt),mean_error=abs(np.nanmean(x)-r.Val) if present.any() else np.nan))
  valid=np.isfinite(x)&((x>=1)&(x<=100) if r.DataID==710 else (x>=21)&(x<=100));x[~valid]=np.nan
  signals[(r.CaseID,r.DataID,int(r.Offset/3600))]=x
 except (ValueError,TypeError) as exc:errors.append(dict(signal=int(r.DataID),error=str(exc)))
pd.DataFrame(errors,columns=['signal','error']).to_csv(O/'minute_decode_failures.csv',index=False)
pd.DataFrame(unavailable).groupby(['signal','declared_count']).size().rename('hours_without_rawdata').reset_index().to_csv(O/'minute_raw_unavailable.csv',index=False)
qa=pd.DataFrame(qa);out=[]
for signal,g in qa.groupby('signal'):out.append(dict(signal=signal,hours=len(g),count_matches=int(g.count_match.sum()),max_abs_mean_error=g.mean_error.max(),mean_error_p99=g.mean_error.quantile(.99),mean_within_001_n=int(g.mean_error.le(.01).sum())))
pd.DataFrame(out).to_csv(O/'minute_decode_checks.csv',index=False)
empty=np.full(60,np.nan);rows=[]
for id,r in d.iterrows():
 h=int(r.hr);sat=signals.get((id,710,h),empty);f0=signals.get((id,2283,h),empty);f1=signals.get((id,2283,h+1),empty)
 tail=sat[30:];n=np.isfinite(tail).sum();meet=(tail>=96).sum();last0=f0[np.isfinite(f0)];last1=f1[np.isfinite(f1)]
 delta=last1[-1]-last0[-1] if len(last0) and len(last1) else np.nan
 cls='missing' if not np.isfinite(delta) else 'reduction' if delta<=-10 else 'stable' if -5<=delta<=5 else 'other'
 follow=np.concatenate([signals.get((id,710,h+j),empty) for j in range(2,8)])
 row=dict(id=id,pid=r.pid,A=int(r.A),n_baseline30=int(n),baseline30_meets=int(meet),baseline30_pass=bool(n>=24 and meet/n>=.8),minute_class=cls,minute_delta=delta,setting0_n=len(last0),setting1_n=len(last1),follow_n=int(np.isfinite(follow).sum()),hospital_mortality=int(r.hospital_mortality))
 for k in [88,90]:
  low=follow<k;row[f'minute_low{k}']=bool(low.any()) if np.isfinite(follow).any() else np.nan;row[f'hourly_low{k}']=r[f'low{k}_any6'];row[f'low{k}_minutes']=int(low.sum());row[f'consecutive5_low{k}']=bool((np.convolve(low.astype(int),np.ones(5,dtype=int),'valid')>=5).any())
 rows.append(row)
p=pd.DataFrame(rows);p.to_pickle(O/'sic_minute_validation_private.pkl')
summary=[];agreement=[];safety=[]
for a,g in p.groupby('A'):
 paired=g[g.minute_class.ne('missing')];expected='reduction' if a else 'stable'
 dense=paired[paired.setting0_n.ge(48)&paired.setting1_n.ge(48)]
 summary.append(dict(A=a,n=len(g),baseline30_atleast24_n=int(g.n_baseline30.ge(24).sum()),baseline30_pass_n=int(g.baseline30_pass.sum()),paired_settings_n=len(paired),same_class_n=int(paired.minute_class.eq(expected).sum()),dense_pair_n=len(dense),dense_same_n=int(dense.minute_class.eq(expected).sum()),follow_any_n=int(g.follow_n.gt(0).sum()),follow_ge288_n=int(g.follow_n.ge(288).sum()),follow_median=g.follow_n.median(),follow_q25=g.follow_n.quantile(.25),follow_q75=g.follow_n.quantile(.75)))
 for cls in ['reduction','stable','other','missing']:agreement.append(dict(hourly_class=expected,minute_class=cls,n=int(g.minute_class.eq(cls).sum())))
 for restrict,mask in [('at_least_one',g.follow_n.gt(0)),('at_least_288',g.follow_n.ge(288))]:
  q=g[mask]
  for k in [88,90]:
   both=q[q[f'hourly_low{k}'].notna()];mn=both[f'minute_low{k}'].astype(bool);hr=both[f'hourly_low{k}'].astype(bool)
   safety.append(dict(A=a,coverage=restrict,threshold=k,n=len(both),minute_events=int(mn.sum()),consecutive5_events=int(both[f'consecutive5_low{k}'].sum()),hourly_events=int(hr.sum()),minute_only=int((mn&~hr).sum()),hourly_only=int((~mn&hr).sum()),low_minutes_median=both[f'low{k}_minutes'].median(),low_minutes_q75=both[f'low{k}_minutes'].quantile(.75)))
pd.DataFrame(summary).to_csv(O/'minute_validation_summary.csv',index=False);pd.DataFrame(agreement).to_csv(O/'minute_exposure_agreement.csv',index=False);pd.DataFrame(safety).to_csv(O/'minute_hypoxaemia_comparison.csv',index=False)
print(pd.DataFrame(out).to_string(index=False));print(pd.DataFrame(summary).to_string(index=False));print(pd.DataFrame(safety).to_string(index=False))
