import os
from pathlib import Path
import pandas as pd,numpy as np,json
O=Path(__file__).resolve().parent;R=O.parent/'audited_reanalysis_v2'
selection=[];delays=[];flows=[];events=[];phase=[];times=[];treatments=[];availability=[]
for db in ['mimic','eicu','sicdb']:
 d=pd.read_pickle(R/f'{db}_cohort.pkl').copy();q=pd.read_pickle(O/f'{db}_preobservation_candidates.pkl').sort_values(['id','hr'])
 first=q.drop_duplicates('id',keep='first').copy();first['next_available']=first.fio2_next.notna()
 first['survived_landmark']=first.icu_end_hr.ge(first.landmark_hr)&(first.death_hr.isna()|first.death_hr.gt(first.landmark_hr))&((first.hosp_end_hr+(24 if db=='sicdb' else 0)).ge(first.landmark_hr)|first.hosp_end_hr.isna())
 for label,g in first.groupby('next_available'):
  for c in ['fio2','sat','age','icu_hour','peep','mbp','heart_rate','hosp_death']:
   if c not in g:continue
   z=pd.to_numeric(g[c],errors='coerce');selection.append(dict(database=db,next_available=int(label),variable=c,n=len(g),observed_n=z.notna().sum(),mean=z.mean(),median=z.median(),q25=z.quantile(.25),q75=z.quantile(.75)))
 event=dict(database=db,physiological_first_stays=len(first),first_next_available=int(first.next_available.sum()),first_next_missing=int((~first.next_available).sum()),deaths_during_classification=int(((first.death_hr>first.decision_hr)&(first.death_hr<=first.landmark_hr)).sum()),left_icu_during_classification=int(((first.icu_end_hr>=first.decision_hr)&(first.icu_end_hr<first.landmark_hr)).sum()))
 events.append(event)
 base=first.set_index('id');d['first_physiological_hr']=d.id.map(base.hr);d['first_to_selected_delay_hr']=d.hr-d.first_physiological_hr
 delays.append(dict(database=db,selected_n=len(d),delayed_n=int(d.first_to_selected_delay_hr.gt(0).sum()),delay_median=d.first_to_selected_delay_hr.median(),delay_q75=d.first_to_selected_delay_hr.quantile(.75),delay_max=d.first_to_selected_delay_hr.max()))
 census=pd.read_csv(O/f'{db}_source_census.csv').iloc[0];oldflow=pd.read_csv(R/f'{db}_flow.csv')
 for stage,n in [('All source ICU stays/cases',census.all_stays),('Adult stays/cases',census.adult_stays),('Adult with documented invasive-support/airway evidence at any time',census.adult_with_airway_evidence)]:flows.append(dict(database=db,stage=stage,stays=int(n)))
 for ix in [5,6,7,8,9,10]:flows.append(dict(database=db,stage=oldflow.iloc[ix].stage,stays=int(oldflow.iloc[ix].stays)))
 if db=='mimic':
  v=pd.read_csv(R/'mimic_vent.csv');v=v[v.ventilation_status.eq('InvasiveVent')]
 elif db=='eicu':v=pd.read_csv(R/'eicu_airway.csv');v=v[v.airwaytype.isin(['Oral ETT','Nasal ETT','Tracheostomy','Double-Lumen Tube','Cricothyrotomy'])]
 else:
  v=pd.read_csv((Path(os.environ['SICDB_DIR']) / 'data_range.csv.gz'));v=v[v.DataID.isin([720,3041])].drop(columns='id').rename(columns={'CaseID':'id'});v['start_hr']=v.Offset/3600;v['end_hr']=v.OffsetEnd/3600
 if db=='eicu':d['support_start_hr']=d.start_hr.where(d.start_hr.le(d.decision_hr))
 else:
  gg=v.groupby('id');starts=[]
  for z in d.itertuples():
   a=gg.get_group(z.id) if z.id in gg.indices else pd.DataFrame(columns=v.columns)
   a=a[a.start_hr.le(z.decision_hr)&a.end_hr.gt(z.decision_hr)]
   starts.append(a.start_hr.max() if len(a) else np.nan)
  d['support_start_hr']=starts
 d['hours_since_support']=d.decision_hr-d.support_start_hr
 for a,g in d.groupby('A'):
  phase.append(dict(database=db,A=a,n=len(g),fio2_100_n=int(g.fio2.eq(100).sum()),support_start_known=int(g.support_start_hr.notna().sum()),support_within6h_n=int(g.hours_since_support.lt(6).sum()),delta_median=g.delta.median(),delta_q25=g.delta.quantile(.25),delta_q75=g.delta.quantile(.75)))
  treatments.append(dict(database=db,A=a,n=len(g),deaths=int(g.hospital_mortality.sum()),mortality_pct=100*g.hospital_mortality.mean()))
 if db in ['mimic','eicu']:
  f=pd.read_pickle(O/f'{db}_fio2_timestamps.pkl');s=pd.read_pickle(O/f'{db}_sat_timestamps.pkl')
  for t,cols in [(f,['minute','fio2']),(s,['minute','sat'])]:
   for c in cols:t[c]=pd.to_numeric(t[c])
  f0=f[(f.minute>=f.hr*60)&(f.minute<(f.hr+1)*60)].groupby('id',sort=False).tail(1).set_index('id')
  f1=f[(f.minute>=(f.hr+1)*60)&(f.minute<(f.hr+2)*60)].groupby('id',sort=False).tail(1).set_index('id')
  s0=s[(s.minute>=s.hr*60)&(s.minute<(s.hr+1)*60)].groupby('id',sort=False).tail(1).set_index('id')
  assert np.allclose(d.fio2,d.id.map(f0.fio2)) and np.allclose(d.fio2_next,d.id.map(f1.fio2)),db
  d['setting_gap_min']=d.id.map(f1.minute)-d.id.map(f0.minute)
  d['baseline_alignment_min']=(d.id.map(f0.minute)-d.id.map(s0.minute)).abs()
  d['baseline_oxygen_age_min']=d.decision_hr*60-d.id.map(f0.minute)
  exp=s[(s.minute>=(s.hr+1)*60)&(s.minute<(s.hr+2)*60)]
  count=exp.groupby('id').size();low=exp.groupby('id').sat.min()
  d['classification_sat_n']=d.id.map(count).fillna(0);d['classification_low88']=d.id.map(low).lt(88).where(d.id.map(count).notna())
  # History is strictly before the baseline bin, at most four hours old; no carrying across stays.
  prev=f[(f.minute<f.hr*60)&(f.minute>=(f.hr-4)*60)].groupby('id',sort=False).tail(1).set_index('id')
  d['history_age_hr']=(d.id.map(f0.minute)-d.id.map(prev.minute))/60
  d['history_change_per_hr']=(d.fio2-d.id.map(prev.fio2))/d.history_age_hr
  d['history_age_hr']=d.history_age_hr.where(d.history_age_hr.gt(0))
  if db=='mimic':
   b=pd.read_pickle(O/'mimic_available_arterial_labs.pkl');b['variable']=b.itemid.map({50821:'po2',50818:'pco2',50820:'ph',50813:'lactate'});b['value']=b.valuenum
  else:
   b=pd.read_pickle(O/'eicu_predecision_bg.pkl');b=b[b.revised_minute.le((b.hr+1)*60)&b.revised_minute.notna()].copy();b['variable']=b.labname.map({'paO2':'po2','paCO2':'pco2','pH':'ph','lactate':'lactate'});b['value']=b.labresult
  b['minute']=pd.to_numeric(b.minute);b['value']=pd.to_numeric(b.value,errors='coerce')
  limits={'po2':(20,700),'pco2':(5,200),'ph':(6.5,8),'lactate':(.1,30)}
  for col,(lo,hi) in limits.items():
   z=b[b.variable.eq(col)&pd.to_numeric(b.value,errors='coerce').between(lo,hi)].groupby('id',sort=False).tail(1).set_index('id')
   d[col]=pd.to_numeric(d.id.map(z.value),errors='coerce');d[col+'_age_hr']=((d.decision_hr*60)-d.id.map(z.minute))/60
  for a,g in d.groupby('A'):
   for col in ['setting_gap_min','baseline_alignment_min','baseline_oxygen_age_min']:
    z=g[col];times.append(dict(database=db,A=a,variable=col,n=len(z),min=z.min(),q25=z.quantile(.25),median=z.median(),q75=z.quantile(.75),max=z.max(),over60_n=int(z.gt(60).sum())))
   events.append(dict(database=db,group='primary treatment '+str(a),n=len(g),classification_measured_n=int(g.classification_sat_n.gt(0).sum()),classification_low88_n=int(g.classification_low88.fillna(0).sum())))
  for c in ['history_age_hr','history_change_per_hr','po2','pco2','ph','lactate']:
   availability.append(dict(database=db,variable=c,n=len(d),observed_n=int(d[c].notna().sum()),missing_pct=100*d[c].isna().mean()))
 else:
  hours=pd.read_pickle(R/'sicdb_hours.pkl').set_index(['id','hr'])
  key=pd.MultiIndex.from_arrays([d.id,d.hr+1]);zz=hours.reindex(key);obs=zz.sat.notna().to_numpy()
  for a in [0,1]:
   mask=d.A.eq(a).to_numpy();events.append(dict(database=db,group='primary treatment '+str(a),n=int(mask.sum()),classification_measured_n=int((obs&mask).sum()),classification_low88_n=int((zz.sat.lt(88).to_numpy()&mask).sum())))
 d.to_pickle(O/f'{db}_augmented_cohort.pkl')
for name,rows in [('selection_comparison',selection),('selection_delay',delays),('stay_flow',flows),('classification_events',events),('care_phase',phase),('timestamp_summary',times),('raw_mortality',treatments),('new_covariate_availability',availability)]:
 pd.DataFrame(rows).to_csv(O/(name+'.csv'),index=False)
print(pd.DataFrame(delays).to_string(index=False));print(pd.DataFrame(times).to_string(index=False));print(pd.DataFrame(availability).to_string(index=False))
