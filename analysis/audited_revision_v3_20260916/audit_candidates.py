import os
from pathlib import Path
import pandas as pd, numpy as np, json, sys
R=Path(__file__).resolve().parent.parent/'audited_reanalysis_v2'
OUT=Path(__file__).resolve().parent
flow=[]

def record(db,stage,d,previous=None):
    if stage.startswith('Current SpO2'): d.to_pickle(OUT/f'{db}_preobservation_candidates.pkl')
    if stage.startswith('Alive and'): d.to_pickle(OUT/f'{db}_precontrast_candidates.pkl')
    flow.append(dict(database=db,stage=stage,rows=len(d),stays=d.id.nunique(),patients=d.pid.nunique() if 'pid' in d else np.nan,
                     excluded_rows=np.nan if previous is None else previous-len(d)))

def interval_mask(d,v):
    valid=np.zeros(len(d),dtype=bool)
    groups=v.groupby('id')
    for id,idx in d.groupby('id').indices.items():
        if id not in groups.indices: continue
        z=groups.get_group(id); t=d.iloc[idx].decision_hr.to_numpy()
        valid[idx]=((t[:,None]>=z.start_hr.to_numpy())&(t[:,None]<z.end_hr.to_numpy())).any(axis=1)
    return valid

def prepare(db):
    if db=='sicdb':
        h=pd.read_pickle(R/'sicdb_hours.pkl')
        m=pd.read_csv((Path(os.environ['SICDB_DIR']) / 'cases.csv.gz')).rename(columns={'CaseID':'id','PatientID':'pid','AgeOnAdmission':'age','Sex':'sex','AdmissionYear':'calendar_group'})
        m['icu_offset_hr']=m.ICUOffset/3600
        m['hosp_end_hr']=m.HospitalDischargeDay*24 # day resolution, not exact time
        m['icu_end_hr']=m.TimeOfStay/3600
        m['death_hr']=m.OffsetOfDeath/3600
        m['hosp_death']=m.HospitalDischargeType.map({2026:0,2028:1,3129:0,3130:1,3131:0,3132:0})
        m['transfer']=m.HospitalDischargeType.eq(3131)
        m['unit']=m.HospitalUnit.astype(str)
        v=pd.read_csv((Path(os.environ['SICDB_DIR']) / 'data_range.csv.gz'))
        v=v[v.DataID.isin([720,3041])].drop(columns='id').rename(columns={'CaseID':'id'})
        v['start_hr']=v.Offset/3600; v['end_hr']=v.OffsetEnd/3600
    else:
        h=pd.read_csv(R/f'{db}_hours.csv')
        m=pd.read_csv(R/f'{db}_meta.csv'); m['icu_offset_hr']=0
        if db=='mimic':
            v=pd.read_csv(R/'mimic_vent.csv')
            v=v[v.ventilation_status.eq('InvasiveVent')]
        else:
            v=pd.read_csv(R/'eicu_airway.csv')
    adult=m[m.age.ge(18)&m.age.notna()]
    anyair=v[v.id.isin(adult.id)]
    if db=='eicu': anyair=anyair[anyair.airwaytype.isin(['Oral ETT','Nasal ETT','Tracheostomy','Double-Lumen Tube','Cricothyrotomy'])]
    pd.DataFrame([dict(database=db,all_stays=len(m),adult_stays=len(adult),adult_with_airway_evidence=anyair.id.nunique())]).to_csv(OUT/f'{db}_source_census.csv',index=False)
    assert not h.duplicated(['id','hr']).any(),db
    h=h.sort_values(['id','hr'])
    for c,lo,hi in [('fio2',21,100),('sat',1,100),('mbp',10,250),('heart_rate',10,300),('rr_vital',1,100),('rr_vent',1,100),('peep',0,40)]:
        if c in h: h.loc[~h[c].between(lo,hi),c]=np.nan
    h=h.merge(m,on='id',validate='many_to_one')
    h['icu_hour']=h.hr-h.icu_offset_hr
    h['decision_hr']=(h.hr+1).astype(float)
    h['landmark_hr']=h.hr+2
    # Observed history must be from the exact preceding bin, not previous available record.
    prev=h[['id','hr','fio2','sat']].copy(); prev.hr+=1
    prev=prev.rename(columns={'fio2':'fio2_prev','sat':'sat_prev'})
    h=h.merge(prev,on=['id','hr'],how='left',validate='one_to_one')
    h['fio2_trend']=h.fio2-h.fio2_prev; h['sat_trend']=h.sat-h.sat_prev
    nxt=h[['id','hr','fio2','sat']].copy(); nxt.hr-=1
    nxt=nxt.rename(columns={'fio2':'fio2_next','sat':'sat_next'})
    h=h.merge(nxt,on=['id','hr'],how='left',validate='one_to_one')
    # Start flow at candidate bins with recorded valid FiO2; never imply all source admissions screened.
    d=h[h.fio2.notna()].copy(); record(db,'Valid FiO2 bins in extract window',d)
    for label,mask in [
        ('Adult age >=18 (SICdb: lowest accepted age bin 20)',d.age.ge(18)&d.age.notna())]:
        n=len(d); d=d[mask].copy(); record(db,label,d,n)
    n=len(d); d=d[d.icu_hour.ge(2)&d.icu_hour.lt(168)].copy(); record(db,'Current bin starts 2 to <168 h after ICU admission',d,n)
    if db=='eicu':
        # Latest airway assessment observed before decision, within 24 h.
        a=v[v.status_hr.notna()].sort_values('status_hr').drop_duplicates(['id','status_hr'],keep='last')
        d=pd.merge_asof(d.sort_values('decision_hr'),a.sort_values('status_hr'),left_on='decision_hr',right_on='status_hr',by='id',direction='backward',tolerance=24)
        vent=d.airwaytype.isin(['Oral ETT','Nasal ETT','Tracheostomy','Double-Lumen Tube','Cricothyrotomy'])
        vent &= (d.start_hr.isna()|d.start_hr.le(d.decision_hr))
        vent &= (d.end_hr.isna()|d.end_hr.le(d.start_hr)|d.end_hr.gt(d.decision_hr))
    else: vent=interval_mask(d,v)
    n=len(d); d=d[vent].copy(); record(db,'Documented invasive ventilation/airway at decision',d,n)
    for label,condition in [
        ('Current FiO2 >=40%',lambda q:q.fio2.ge(40)),
        ('Current SpO2 >=96% (hourly mean in SICdb)',lambda q:q.sat.ge(96)&q.sat.le(100)),
        ('Next-bin FiO2 available; next-bin saturation not required',lambda q:q.fio2_next.notna()),
        ('Alive and in observed ICU/hospital stay through landmark',lambda q:q.icu_end_hr.ge(q.landmark_hr)&(q.death_hr.isna()|q.death_hr.gt(q.landmark_hr))&((q.hosp_end_hr+(24 if db=='sicdb' else 0)).ge(q.landmark_hr)|q.hosp_end_hr.isna()))]:
        n=len(d); d=d[condition(d)].copy(); record(db,label,d,n)
    d=d.sort_values(['id','hr'])
    repeated=d.copy()
    n=len(d); d=d.drop_duplicates('id',keep='first').copy(); record(db,'First eligible decision per stay/case',d,n)
    d['delta']=d.fio2_next-d.fio2
    d['treatment']=np.select([d.delta.le(-10),d.delta.between(-5,5)],['deescalation','stable'],default='other')
    d.groupby('treatment').size().rename('n').to_csv(OUT/f'{db}_first_decision_categories.csv')
    n=len(d); d=d[d.treatment.ne('other')].copy(); record(db,'Retain >=10 pp reduction or stable +/-5 pp',d,n)
    d['A']=d.treatment.eq('deescalation').astype(int)
    # Observed in-hospital death, not all known subsequent deaths. For SICdb hospital discharge is day-resolution.
    n=len(d); d=d[d.hosp_death.notna()].copy(); record(db,'Known hospital disposition',d,n)
    d['hospital_mortality']=d.hosp_death.astype(int)
    d['death28_recorded']=((d.death_hr>d.landmark_hr)&(d.death_hr<=d.landmark_hr+672)).astype(int)
    # Descriptive post-discharge ascertainment categories, do not call discharges censoring for in-hospital endpoint.
    d['alive_discharge_before28']=d.hosp_death.eq(0)&d.hosp_end_hr.le(d.landmark_hr+672)
    if db=='sicdb':
        d['hospital_before_landmark']=d.hosp_end_hr+24<d.landmark_hr
    else: d['hospital_before_landmark']=d.hosp_end_hr<d.landmark_hr
    assert not d.hospital_before_landmark.any(),db
    # Measurement-available outcomes. Start strictly after end of exposure bin (h+2 through h+7).
    for c in ['n_observed_bins6','n_measurements6','n_low88_bins6','n_low90_bins6']: d[c]=0
    idx=h.set_index(['id','hr'])
    for shift in range(2,8):
        key=pd.MultiIndex.from_arrays([d.id,d.hr+shift])
        z=idx.reindex(key)
        obs=z.sat.notna().to_numpy()
        d['n_observed_bins6']+=obs.astype(int)
        d['n_measurements6']+=np.where(obs,z.sat_n.fillna(0).to_numpy(),0)
        for k in [88,90]:
            # SICdb sat_min is intentionally the hourly mean, with explicit source-specific labeling.
            d[f'n_low{k}_bins6']+=(obs&z.sat_min.lt(k).to_numpy()).astype(int)
    for k in [88,90]:
        d[f'low{k}_any6']=(d[f'n_low{k}_bins6']>0).astype(float)
        d.loc[d.n_observed_bins6.eq(0),f'low{k}_any6']=np.nan
    d.to_pickle(OUT/f'{db}_cohort.pkl')
    # Local patient-level data only; never copied into shareable submission package.
    repeated.to_pickle(OUT/f'{db}_repeated_candidates.pkl')
    pd.DataFrame([r for r in flow if r['database']==db]).to_csv(OUT/f'{db}_flow.csv',index=False)
    print(db,'final',len(d),'treated',int(d.A.sum()),'patients',d.pid.nunique(),'hospital deaths',int(d.hospital_mortality.sum()),flush=True)

if __name__=='__main__':
    for db in sys.argv[1:] or ['mimic','eicu','sicdb']:
        prepare(db)
    pd.DataFrame(flow).to_csv(OUT/('flow_'+'_'.join(sys.argv[1:] or ['all'])+'.csv'),index=False)
