from pathlib import Path
import os,sys,json
import pandas as pd,numpy as np,psycopg2
O=Path(__file__).resolve().parent;R=O.parent/'audited_reanalysis_v2'
db=sys.argv[1];d=pd.read_pickle(R/f'{db}_cohort.pkl')
con=psycopg2.connect(host=os.environ.get('PGHOST','127.0.0.1'),port=int(os.environ.get('PGPORT','5432')),user=os.environ.get('PGUSER','postgres'),password=os.environ['PGPASSWORD'],dbname={'mimiciv31':os.environ.get('MIMIC_DATABASE','mimiciv31'),'eicu':os.environ.get('EICU_DATABASE','eicu')}.get(('mimiciv31' if db=='mimic' else 'eicu'), ('mimiciv31' if db=='mimic' else 'eicu')));con.set_session(readonly=True)
def fetch(sql,params):
 with con.cursor() as c:
  c.execute("SET statement_timeout='12min'");c.execute(sql,params);return pd.DataFrame(c.fetchall(),columns=[x[0] for x in c.description])
ids=d.id.astype(int).tolist();hrs=d.hr.astype(int).tolist()
cte='WITH c AS (SELECT unnest(%s::bigint[]) AS id, unnest(%s::integer[]) AS hr) '
if db=='mimic':
 sql=cte+'''SELECT c.id,c.hr,extract(epoch from(v.charttime-i.intime))/60 AS minute,v.fio2
 FROM c JOIN mimiciv_icu.icustays i ON i.stay_id=c.id JOIN mimiciv_derived.ventilator_setting v ON v.stay_id=c.id
 AND v.charttime>=i.intime+greatest(0,c.hr-4)*interval '1 hour' AND v.charttime<i.intime+(c.hr+2)*interval '1 hour'
 WHERE v.fio2 BETWEEN 21 AND 100 ORDER BY c.id,v.charttime'''
 f=fetch(sql,(ids,hrs));f.to_pickle(O/f'{db}_fio2_timestamps.pkl');print(db,'setting timestamps',len(f),flush=True)
 sql=cte+'''SELECT c.id,c.hr,extract(epoch from(v.charttime-i.intime))/60 AS minute,v.spo2 AS sat
 FROM c JOIN mimiciv_icu.icustays i ON i.stay_id=c.id JOIN mimiciv_derived.vitalsign v ON v.stay_id=c.id
 AND v.charttime>=i.intime+c.hr*interval '1 hour' AND v.charttime<i.intime+(c.hr+2)*interval '1 hour'
 WHERE v.spo2 BETWEEN 1 AND 100 ORDER BY c.id,v.charttime'''
 s=fetch(sql,(ids,hrs));s.to_pickle(O/f'{db}_sat_timestamps.pkl');print(db,'saturation timestamps',len(s),flush=True)
 sql=cte+'''SELECT c.id,c.hr,extract(epoch from(b.charttime-i.intime))/60 AS minute,b.po2,b.pco2,b.ph,b.lactate
 FROM c JOIN mimiciv_icu.icustays i ON i.stay_id=c.id JOIN mimiciv_derived.bg b ON b.hadm_id=i.hadm_id
 AND b.charttime>=i.intime+(c.hr+1-6)*interval '1 hour' AND b.charttime<i.intime+(c.hr+1)*interval '1 hour'
 WHERE b.specimen='ART.' ORDER BY c.id,b.charttime'''
 b=fetch(sql,(ids,hrs));b.to_pickle(O/f'{db}_predecision_bg.pkl');print(db,'arterial blood gas rows',len(b),flush=True)
else:
 sql=cte+'''SELECT c.id,c.hr,r.respchartoffset AS minute,r.respchartentryoffset AS entry_minute,
 CASE WHEN r.respchartvalue ~ '^ *[0-9]+([.][0-9]+)? *%%? *$' THEN replace(trim(r.respchartvalue),'%%','')::numeric END AS fio2
 FROM c JOIN respiratorycharting r ON r.patientunitstayid=c.id AND r.respchartoffset>=greatest(0,c.hr-4)*60 AND r.respchartoffset<(c.hr+2)*60
 WHERE r.respchartvaluelabel IN('FiO2','FIO2 (%%)','Set Fraction of Inspired Oxygen (FIO2)') ORDER BY c.id,r.respchartoffset,r.respchartentryoffset'''
 f=fetch(sql,(ids,hrs));f['fio2']=pd.to_numeric(f.fio2);f.loc[f.fio2.between(.21,1),'fio2']*=100;f=f[f.fio2.between(21,100)];f.to_pickle(O/f'{db}_fio2_timestamps.pkl');print(db,'setting timestamps',len(f),flush=True)
 sql=cte+'''SELECT c.id,c.hr,v.observationoffset AS minute,v.sao2 AS sat FROM c JOIN vitalperiodic v ON v.patientunitstayid=c.id
 AND v.observationoffset>=c.hr*60 AND v.observationoffset<(c.hr+2)*60 WHERE v.sao2 BETWEEN 1 AND 100 ORDER BY c.id,v.observationoffset'''
 s=fetch(sql,(ids,hrs));s.to_pickle(O/f'{db}_sat_timestamps.pkl');print(db,'saturation timestamps',len(s),flush=True)
 sql=cte+'''SELECT c.id,c.hr,l.labresultoffset AS minute,l.labname,l.labresult,l.labresultrevisedoffset AS revised_minute
 FROM c JOIN lab l ON l.patientunitstayid=c.id AND l.labresultoffset>=(c.hr+1-6)*60 AND l.labresultoffset<(c.hr+1)*60
 WHERE lower(l.labname) IN('pao2','paco2','ph','lactate') ORDER BY c.id,l.labresultoffset,l.labresultrevisedoffset'''
 b=fetch(sql,(ids,hrs));b.to_pickle(O/f'{db}_predecision_bg.pkl');print(db,'blood gas and lactate rows',len(b),'labels',b.labname.unique().tolist(),flush=True)
con.close()
(O/f'{db}_timestamp_extraction_done.json').write_text(json.dumps(dict(setting_rows=len(f),sat_rows=len(s),lab_rows=len(b))))
