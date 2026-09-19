import os
from analysis_utils import *
import psycopg2,re
d=pd.read_pickle(O/'eicu_timing_support_private.pkl');c=psycopg2.connect(host=os.environ.get('PGHOST','127.0.0.1'),port=int(os.environ.get('PGPORT','5432')),user=os.environ.get('PGUSER','postgres'),password=os.environ['PGPASSWORD'],dbname={'mimiciv31':os.environ.get('MIMIC_DATABASE','mimiciv31'),'eicu':os.environ.get('EICU_DATABASE','eicu')}.get(('eicu'), ('eicu')));c.set_session(readonly=True)
rate=['Vent Rate','Adult Con Setting Set RR'];vol=['Tidal Volume (set)','Adult Con Setting Set Vt','Set Vt (Servo,LTV)'];pressure=['Pressure Support','Pressure Control','Inspiratory Pressure, Set','PS above PEEP'];peep=['PEEP','PEEP/CPAP']
with c.cursor() as q:
 q.execute("""WITH c AS (SELECT unnest(%s::bigint[]) AS id,unnest(%s::integer[]) AS hr)
 SELECT c.id,r.respchartoffset,r.respchartentryoffset,r.respchartvaluelabel AS label,r.respchartvalue AS value
 FROM c JOIN respiratorycharting r ON r.patientunitstayid=c.id AND r.respchartoffset>=(c.hr+1)*60-120 AND r.respchartoffset<(c.hr+1)*60
 WHERE r.respchartvaluelabel=ANY(%s) OR r.respchartvaluelabel LIKE 'NIV%%' OR r.respchartvaluelabel IN('Non-invasive Ventilation Mode','Bipap Delivery Mode')""",(d.id.astype(int).tolist(),d.hr.astype(int).tolist(),rate+vol+pressure+peep))
 z=pd.DataFrame(q.fetchall(),columns=[x[0] for x in q.description])
c.close();z.to_pickle(O/'eicu_ventilator_settings_private.pkl');z['numeric']=pd.to_numeric(z.value,errors='coerce');summ=[]
for name,labels,lo,hi in [('set_rate',rate,1,60),('set_volume',vol,50,2000),('pressure_assist',pressure,1,60),('peep',peep,0,40)]:
 d[name+'_evidence']=d.id.isin(z[z.label.isin(labels)&z.numeric.between(lo,hi)].id)
d['niv_label']=d.id.isin(z[z.label.str.startswith('NIV')|z.label.isin(['Non-invasive Ventilation Mode','Bipap Delivery Mode'])].id)
d['corroborating_settings']=((d.set_rate_evidence&d.set_volume_evidence)|(d.pressure_assist_evidence&d.peep_evidence))&~d.niv_label
for a,g in d.groupby('A'):
 row={'A':a,'n':len(g)}
 for col in ['set_rate_evidence','set_volume_evidence','pressure_assist_evidence','peep_evidence','niv_label','corroborating_settings']:row[col+'_n']=int(g[col].sum())
 summ.append(row)
pd.DataFrame(summ).to_csv(O/'eicu_corroboration_counts.csv',index=False);d.to_pickle(O/'eicu_corroboration_private.pkl');print(pd.DataFrame(summ).to_string(index=False),flush=True)
