import os
from pathlib import Path
import sys,time,json
import numpy as np,pandas as pd
O=Path(__file__).resolve().parent;R=O.parent/'audited_reanalysis_v2'
def prepare():
 h=pd.read_pickle(R/'sicdb_hours.pkl');m=pd.read_csv((Path(os.environ['SICDB_DIR']) / 'cases.csv.gz')).rename(columns={'CaseID':'id','PatientID':'pid','AgeOnAdmission':'age'})
 h=h.merge(m[['id','pid','age','ICUOffset','TimeOfStay','OffsetOfDeath','HospitalDischargeDay','HospitalDischargeType']],on='id',validate='many_to_one')
 h['icu_hour']=h.hr-h.ICUOffset/3600
 h=h[h.age.ge(18)&h.icu_hour.ge(2)&h.icu_hour.lt(168)&h.fio2.between(21,100)].copy().reset_index(drop=True)
 v=pd.read_csv((Path(os.environ['SICDB_DIR']) / 'data_range.csv.gz'));v=v[v.DataID.isin([720,3041])].rename(columns={'CaseID':'case_id'})
 groups=v.groupby('case_id');valid=np.zeros(len(h),bool)
 for case,idx in h.groupby('id').indices.items():
  if case not in groups.indices:continue
  z=groups.get_group(case);t=(h.iloc[idx].hr.to_numpy()+1)*3600
  valid[idx]=((t[:,None]>=z.Offset.to_numpy())&(t[:,None]<z.OffsetEnd.to_numpy())).any(axis=1)
 h=h[valid].sort_values(['id','hr']).reset_index(drop=True)
 h['landmark_hr']=h.hr+2;h['survives_landmark']=(h.TimeOfStay/3600>=h.landmark_hr)&(h.OffsetOfDeath.isna()|(h.OffsetOfDeath/3600>h.landmark_hr))&(h.HospitalDischargeDay.isna()|(h.HospitalDischargeDay*24+24>=h.landmark_hr))
 h['hospital_mortality']=h.HospitalDischargeType.map({2026:0,2028:1,3129:0,3130:1,3131:0,3132:0})
 h.to_pickle(O/'sic_candidates_private.pkl')
 keys=pd.concat([h[['id','hr']],h[['id','hr']].assign(hr=h.hr+1)]).drop_duplicates().sort_values(['id','hr']).reset_index(drop=True)
 assert keys.hr.between(0,99999).all();keys['key']=keys.id.astype('int64')*100000+keys.hr.astype('int64');keys.to_pickle(O/'sic_signal_keys_private.pkl')
 print('candidate hours',len(h),'cases',h.id.nunique(),'signal hour keys',len(keys),flush=True)
def extract():
 keys=pd.read_pickle(O/'sic_signal_keys_private.pkl');ix=pd.Index(keys.key)
 arr=np.lib.format.open_memmap(O/'sic_candidate_signals_private.npy',mode='w+',dtype='float32',shape=(len(keys),2,60));arr[:]=np.nan
 last=np.full((len(keys),2),-1,dtype='int64');stats={'scanned':0,'decoded':0,'count_mismatch':0,'mean_error_max':0.,'missing_array':0,'decode_failure':0,'invalid_minute_values':0}
 errors=[];quality=[];start=time.time()
 for c in pd.read_csv((Path(os.environ['SICDB_DIR']) / 'data_float_h.csv.gz'),chunksize=150000,usecols=['id','CaseID','DataID','Offset','Val','cnt','rawdata']):
  stats['scanned']+=len(c);c=c[c.DataID.isin([710,2283])].copy()
  if len(c):
   c['row_index']=ix.get_indexer(c.CaseID.astype('int64')*100000+(c.Offset/3600).astype('int64'));c=c[c.row_index.ge(0)]
   for row in c.itertuples():
    j=0 if row.DataID==2283 else 1;k=row.row_index
    if row.id<=last[k,j]:continue
    last[k,j]=row.id;arr[k,j,:]=np.nan
    if pd.isna(row.rawdata):stats['missing_array']+=1;continue
    try:
     s=str(row.rawdata);b=bytes.fromhex(s[2:] if s.startswith('0x') else s)
     if len(b)!=240:raise ValueError('wrong byte length')
     x=np.frombuffer(b,dtype='<f4').copy();present=np.frombuffer(b,dtype='<u4')!=0;x[~present]=np.nan
     stats['decoded']+=1;stats['count_mismatch']+=int(present.sum()!=row.cnt)
     err=abs(float(np.nanmean(x.astype(float)))-row.Val) if present.any() else 0
     stats['mean_error_max']=max(stats['mean_error_max'],err)
     if present.sum()!=row.cnt or err>.01:quality.append(dict(case_id=row.CaseID,hr=int(row.Offset/3600),signal=row.DataID,source_row_id=row.id,count_expected=row.cnt,count_decoded=int(present.sum()),mean_error=err))
     good=np.isfinite(x)&(x>=(21 if j==0 else 1))&(x<=100);stats['invalid_minute_values']+=int((present&~good).sum());x[~good]=np.nan;arr[k,j,:]=x
    except (ValueError,TypeError) as exc:stats['decode_failure']+=1;errors.append({'signal':row.DataID,'error':str(exc)})
  if stats['scanned']%3000000==0:print('scanned',stats['scanned'],'seconds',round(time.time()-start),flush=True)
 arr.flush();pd.DataFrame(quality,columns=['case_id','hr','signal','source_row_id','count_expected','count_decoded','mean_error']).to_pickle(O/'sic_quality_flags_private.pkl');np.save(O/'sic_source_row_ids_private.npy',last);stats['signal_keys']=len(keys);stats['present_fio2_hours']=int((last[:,0]>=0).sum());stats['present_sat_hours']=int((last[:,1]>=0).sum())
 (O/'sic_decode_audit.json').write_text(json.dumps(stats,indent=2));pd.DataFrame(errors,columns=['signal','error']).to_csv(O/'sic_decode_failures.csv',index=False);print(stats,flush=True)
if __name__=='__main__':
 prepare();extract()
