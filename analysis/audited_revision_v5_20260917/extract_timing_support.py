import os
from analysis_utils import *
import psycopg2
def connection(db):
 c=psycopg2.connect(host=os.environ.get('PGHOST','127.0.0.1'),port=int(os.environ.get('PGPORT','5432')),user=os.environ.get('PGUSER','postgres'),password=os.environ['PGPASSWORD'],dbname={'mimiciv31':os.environ.get('MIMIC_DATABASE','mimiciv31'),'eicu':os.environ.get('EICU_DATABASE','eicu')}.get((db), (db)));c.set_session(readonly=True);return c
def fetch(c,sql,params=()):
 with c.cursor() as q:
  q.execute("SET statement_timeout='8min'");q.execute(sql,params);return pd.DataFrame(q.fetchall(),columns=[x[0] for x in q.description])
e=pd.read_pickle(R/'eicu_cohort.pkl');f=pd.read_pickle(V/'eicu_fio2_timestamps.pkl');rows=[]
for shift,label,field in [(0,'baseline','fio2'),(1,'classification','fio2_next')]:
 q=f[f.minute.ge((f.hr+shift)*60)&f.minute.lt((f.hr+shift+1)*60)].sort_values(['id','minute','entry_minute']).drop_duplicates('id',keep='last').set_index('id')
 assert np.allclose(e[field],e.id.map(q.fio2))
 e[label+'_entry_minute']=e.id.map(q.entry_minute);e[label+'_entry_by_boundary']=e[label+'_entry_minute'].le((e.hr+shift+1)*60)
 for a,g in e.groupby('A'):rows.append(dict(database='eicu',phase=label,A=a,n=len(g),matched_n=len(g),late_n=int((~g[label+'_entry_by_boundary']).sum()),entry_missing_n=int(g[label+'_entry_minute'].isna().sum())))
e['explicit_vent_interval']=e.start_hr.notna()&e.end_hr.notna()&e.start_hr.le(e.decision_hr)&e.end_hr.gt(e.decision_hr)
c=connection('eicu');m=fetch(c,"""WITH c AS (SELECT unnest(%s::bigint[]) AS id,unnest(%s::integer[]) AS hr)
 SELECT c.id,c.hr,r.respchartoffset,r.respchartentryoffset,r.respchartvaluelabel,r.respchartvalue
 FROM c JOIN respiratorycharting r ON r.patientunitstayid=c.id
 AND r.respchartoffset>=(c.hr+1)*60-120 AND r.respchartoffset<(c.hr+1)*60
 WHERE r.respchartvaluelabel IN ('Mechanical Ventilator Mode','Ventilator Support Mode','Non-invasive Ventilation Mode','Bipap Delivery Mode')
 ORDER BY c.id,r.respchartoffset,r.respchartentryoffset""",(e.id.astype(int).tolist(),e.hr.astype(int).tolist()));c.close()
m.to_pickle(O/'eicu_mode_records_private.pkl')
m.groupby(['respchartvaluelabel','respchartvalue'],dropna=False).agg(records=('id','size'),stays=('id','nunique')).reset_index().to_csv(O/'eicu_mode_label_counts.csv',index=False)
e['mode_record_within2h']=e.id.isin(m.id);e['mechanical_mode_label_within2h']=e.id.isin(m[m.respchartvaluelabel.eq('Mechanical Ventilator Mode')].id)
e.to_pickle(O/'eicu_timing_support_private.pkl')
rr=[]
for a,g in e.groupby('A'):
 rr.append(dict(A=a,n=len(g),explicit_interval_n=int(g.explicit_vent_interval.sum()),mode_record_n=int(g.mode_record_within2h.sum()),mechanical_mode_label_n=int(g.mechanical_mode_label_within2h.sum()),both_entries_by_boundary_n=int((g.baseline_entry_by_boundary&g.classification_entry_by_boundary).sum())))
pd.DataFrame(rr).to_csv(O/'eicu_support_audit.csv',index=False)
print('eicu support',rr,flush=True)
# Exact original FiO2 charttime/value match against raw chartevents; earliest storage among identical records.
d=pd.read_pickle(R/'mimic_cohort.pkl');f=pd.read_pickle(V/'mimic_fio2_timestamps.pkl');selected=[]
for shift,label in [(0,'baseline'),(1,'classification')]:
 q=f[f.minute.ge((f.hr+shift)*60)&f.minute.lt((f.hr+shift+1)*60)].groupby('id',sort=False).tail(1).copy();q['phase']=label;selected.append(q)
s=pd.concat(selected);c=connection('mimiciv31')
items=fetch(c,"SELECT itemid,label FROM mimiciv_icu.d_items WHERE itemid IN (223835,3420,3422,190)")
items.to_csv(O/'mimic_fio2_item_dictionary.csv',index=False)
raw=fetch(c,"""WITH c AS (SELECT unnest(%s::bigint[]) AS id, unnest(%s::integer[]) AS hr)
 SELECT c.id,c.hr,extract(epoch from(ch.charttime-i.intime))/60 AS minute,
 extract(epoch from(ch.storetime-i.intime))/60 AS entry_minute,ch.valuenum,ch.itemid
 FROM c JOIN mimiciv_icu.icustays i ON i.stay_id=c.id JOIN mimiciv_icu.chartevents ch ON ch.stay_id=c.id
 AND ch.charttime>=i.intime+c.hr*interval '1 hour' AND ch.charttime<i.intime+(c.hr+2)*interval '1 hour'
 WHERE ch.itemid IN (223835,3420,3422,190) AND ch.valuenum IS NOT NULL""",(d.id.astype(int).tolist(),d.hr.astype(int).tolist()));c.close()
raw.to_pickle(O/'mimic_fio2_raw_times_private.pkl');raw['minute']=pd.to_numeric(raw.minute);raw['entry_minute']=pd.to_numeric(raw.entry_minute);raw['value']=pd.to_numeric(raw.valuenum)
raw.loc[raw.value.between(.21,1),'value']*=100
s['minute']=pd.to_numeric(s.minute);s['fio2']=pd.to_numeric(s.fio2)
joined=s.merge(raw,on=['id','hr','minute'],how='left');joined=joined[np.isclose(joined.fio2,joined.value,equal_nan=False)]
matched=joined.groupby(['id','phase']).entry_minute.min()
for shift,label in [(0,'baseline'),(1,'classification')]:
 d[label+'_raw_entry']=d.id.map(matched.xs(label,level='phase'));d[label+'_entry_by_boundary']=d[label+'_raw_entry'].le((d.hr+shift+1)*60)
 for a,g in d.groupby('A'):rows.append(dict(database='mimic',phase=label,A=a,n=len(g),matched_n=int(g[label+'_raw_entry'].notna().sum()),late_n=int((g[label+'_raw_entry']>(g.hr+shift+1)*60).sum()),entry_missing_n=int(g[label+'_raw_entry'].isna().sum())))
d.to_pickle(O/'mimic_timing_private.pkl');pd.DataFrame(rows).to_csv(O/'record_availability_audit.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False),flush=True)
