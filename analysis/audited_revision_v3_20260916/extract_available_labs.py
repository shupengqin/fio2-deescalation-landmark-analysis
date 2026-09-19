from pathlib import Path
import os,psycopg2,pandas as pd
O=Path(__file__).resolve().parent;R=O.parent/'audited_reanalysis_v2';d=pd.read_pickle(R/'mimic_cohort.pkl')
c=psycopg2.connect(host=os.environ.get('PGHOST','127.0.0.1'),port=int(os.environ.get('PGPORT','5432')),user=os.environ.get('PGUSER','postgres'),password=os.environ['PGPASSWORD'],dbname={'mimiciv31':os.environ.get('MIMIC_DATABASE','mimiciv31'),'eicu':os.environ.get('EICU_DATABASE','eicu')}.get(('mimiciv31'), ('mimiciv31')));c.set_session(readonly=True)
with c.cursor() as x:
 x.execute('select itemid,label,fluid from mimiciv_hosp.d_labitems where itemid in (50821,50818,50820,50813)');print(x.fetchall())
 sql='''WITH c AS (SELECT unnest(%s::bigint[]) id,unnest(%s::int[]) hr)
 SELECT c.id,c.hr,extract(epoch from(l.charttime-i.intime))/60 AS minute,extract(epoch from(l.storetime-i.intime))/60 AS available_minute,
 l.itemid,l.valuenum,l.valueuom
 FROM c JOIN mimiciv_icu.icustays i ON i.stay_id=c.id JOIN mimiciv_hosp.labevents l ON l.hadm_id=i.hadm_id
 AND l.charttime>=i.intime+(c.hr-5)*interval '1 hour' AND l.charttime<i.intime+(c.hr+1)*interval '1 hour'
 AND l.storetime<=i.intime+(c.hr+1)*interval '1 hour'
 WHERE l.itemid IN(50821,50818,50820,50813) AND l.valuenum IS NOT NULL
 AND EXISTS(SELECT 1 FROM mimiciv_derived.bg b WHERE b.hadm_id=l.hadm_id AND b.charttime=l.charttime AND b.specimen='ART.')
 ORDER BY c.id,l.charttime,l.storetime,l.labevent_id'''
 x.execute("SET statement_timeout='10min'");x.execute(sql,(d.id.astype(int).tolist(),d.hr.astype(int).tolist()));z=pd.DataFrame(x.fetchall(),columns=[v[0] for v in x.description]);z.to_pickle(O/'mimic_available_arterial_labs.pkl');print('available arterial labs',len(z))
c.close()
