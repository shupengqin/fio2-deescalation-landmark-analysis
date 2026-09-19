"""Read-only extraction. All output goes to this new audit directory."""
from pathlib import Path
import os, sys, json, time, gzip, shutil
import pandas as pd
import psycopg2

R=Path(__file__).resolve().parent

def copy(db,name,query):
    dest=R/(name+'.csv')
    if dest.exists() and dest.with_suffix('.done').exists(): return
    (R/(name+'.sql')).write_text(query,encoding='utf8')
    con=psycopg2.connect(host=os.environ.get('PGHOST','127.0.0.1'),port=int(os.environ.get('PGPORT','5432')),user=os.environ.get('PGUSER','postgres'),password=os.environ['PGPASSWORD'],dbname={'mimiciv31':os.environ.get('MIMIC_DATABASE','mimiciv31'),'eicu':os.environ.get('EICU_DATABASE','eicu')}.get((db), (db)))
    con.set_session(readonly=True)
    with con.cursor() as c, dest.open('w',encoding='utf8',newline='') as f:
        c.execute("SET statement_timeout='15min'")
        c.copy_expert('COPY ('+query+') TO STDOUT WITH CSV HEADER', f)
    con.close(); dest.with_suffix('.done').write_text(str(time.time()))
    print(name,dest.stat().st_size,flush=True)

if sys.argv[1]=='mimic':
    copy('mimiciv31','mimic_meta',"""SELECT i.stay_id AS id,i.subject_id AS pid,i.hadm_id,
    p.anchor_age+extract(year from a.admittime)-p.anchor_year AS age,p.gender AS sex,
    p.anchor_year_group AS calendar_group,i.first_careunit AS unit,
    extract(epoch from (a.deathtime-i.intime))/3600 AS death_hr,
    extract(epoch from (a.dischtime-i.intime))/3600 AS hosp_end_hr,
    extract(epoch from (i.outtime-i.intime))/3600 AS icu_end_hr,
    a.hospital_expire_flag AS hosp_death
    FROM mimiciv_icu.icustays i JOIN mimiciv_hosp.admissions a USING(hadm_id)
    JOIN mimiciv_hosp.patients p ON p.subject_id=i.subject_id""")
    copy('mimiciv31','mimic_vent',"""SELECT v.stay_id AS id,
    extract(epoch from (v.starttime-i.intime))/3600 AS start_hr,
    extract(epoch from (v.endtime-i.intime))/3600 AS end_hr,ventilation_status
    FROM mimiciv_derived.ventilation v JOIN mimiciv_icu.icustays i USING(stay_id)""")
    copy('mimiciv31','mimic_hours',"""WITH f AS (
      SELECT v.stay_id AS id,floor(extract(epoch from(v.charttime-i.intime))/3600)::int AS hr,
      (array_agg(fio2 ORDER BY charttime DESC) FILTER(WHERE fio2 BETWEEN 21 AND 100))[1] AS fio2,
      (array_agg(peep ORDER BY charttime DESC) FILTER(WHERE peep BETWEEN 0 AND 40))[1] AS peep,
      (array_agg(respiratory_rate_total ORDER BY charttime DESC) FILTER(WHERE respiratory_rate_total BETWEEN 1 AND 100))[1] AS rr_vent
      FROM mimiciv_derived.ventilator_setting v JOIN mimiciv_icu.icustays i USING(stay_id)
      WHERE v.charttime>=i.intime AND v.charttime<i.outtime AND v.charttime<i.intime+interval '8 days'
      GROUP BY v.stay_id,hr
    ), s AS (
      SELECT v.stay_id AS id,floor(extract(epoch from(v.charttime-i.intime))/3600)::int AS hr,
      (array_agg(spo2 ORDER BY charttime DESC) FILTER(WHERE spo2 BETWEEN 1 AND 100))[1] AS sat,
      count(spo2) FILTER(WHERE spo2 BETWEEN 1 AND 100) AS sat_n,
      min(spo2) FILTER(WHERE spo2 BETWEEN 1 AND 100) AS sat_min,
      (array_agg(coalesce(mbp,mbp_ni) ORDER BY charttime DESC) FILTER(WHERE coalesce(mbp,mbp_ni) BETWEEN 10 AND 250))[1] AS mbp,
      (array_agg(heart_rate ORDER BY charttime DESC) FILTER(WHERE heart_rate BETWEEN 10 AND 300))[1] AS heart_rate,
      (array_agg(resp_rate ORDER BY charttime DESC) FILTER(WHERE resp_rate BETWEEN 1 AND 100))[1] AS rr_vital
      FROM mimiciv_derived.vitalsign v JOIN mimiciv_icu.icustays i USING(stay_id)
      WHERE v.charttime>=i.intime AND v.charttime<i.outtime AND v.charttime<i.intime+interval '8 days'
      GROUP BY v.stay_id,hr
    ) SELECT coalesce(f.id,s.id) AS id,coalesce(f.hr,s.hr) AS hr,fio2,peep,rr_vent,sat,sat_n,sat_min,mbp,heart_rate,rr_vital
    FROM f FULL JOIN s USING(id,hr)""")

elif sys.argv[1]=='eicu':
    copy('eicu','eicu_meta',"""SELECT patientunitstayid AS id,uniquepid AS pid,patienthealthsystemstayid AS hadm_id,
    CASE WHEN age ~ '^[0-9]+$' THEN age::numeric WHEN age='> 89' THEN 90 WHEN age='< 18' THEN 17 ELSE NULL END AS age,
    gender AS sex,hospitalid AS site,hospitaldischargeyear AS calendar_group,unittype AS unit,
    hospitaldischargeoffset/60.0 AS hosp_end_hr,unitdischargeoffset/60.0 AS icu_end_hr,
    CASE WHEN hospitaldischargestatus='Expired' THEN 1 WHEN hospitaldischargestatus='Alive' THEN 0 ELSE NULL END AS hosp_death,
    CASE WHEN hospitaldischargestatus='Expired' THEN hospitaldischargeoffset/60.0 ELSE NULL END AS death_hr
    FROM patient""")
    copy('eicu','eicu_airway',"""SELECT patientunitstayid AS id,respcarestatusoffset/60.0 AS status_hr,airwaytype,
    ventstartoffset/60.0 AS start_hr,ventendoffset/60.0 AS end_hr FROM respiratorycare
    WHERE airwaytype IS NOT NULL AND airwaytype<>''""")
    copy('eicu','eicu_hours',"""WITH fr AS (
      SELECT patientunitstayid AS id,respchartoffset,respchartentryoffset,
      CASE WHEN respchartvalue ~ '^ *[0-9]+([.][0-9]+)? *%? *$'
      THEN replace(trim(respchartvalue),'%','')::numeric ELSE NULL END AS val
      FROM respiratorycharting WHERE respchartvaluelabel IN('FiO2','FIO2 (%)','Set Fraction of Inspired Oxygen (FIO2)')
      AND respchartoffset>=0 AND respchartoffset<11520
    ), f AS (
      SELECT id,floor(respchartoffset/60.0)::int AS hr,
      (array_agg(CASE WHEN val BETWEEN .21 AND 1 THEN val*100 ELSE val END ORDER BY respchartoffset DESC,respchartentryoffset DESC)
      FILTER(WHERE val BETWEEN 21 AND 100 OR val BETWEEN .21 AND 1))[1] AS fio2
      FROM fr GROUP BY id,hr
    ), s AS (
      SELECT patientunitstayid AS id,floor(observationoffset/60.0)::int AS hr,
      (array_agg(sao2 ORDER BY observationoffset DESC) FILTER(WHERE sao2 BETWEEN 1 AND 100))[1] AS sat,
      count(sao2) FILTER(WHERE sao2 BETWEEN 1 AND 100) AS sat_n,
      min(sao2) FILTER(WHERE sao2 BETWEEN 1 AND 100) AS sat_min,
      (array_agg(systemicmean ORDER BY observationoffset DESC) FILTER(WHERE systemicmean BETWEEN 10 AND 250))[1] AS mbp,
      (array_agg(heartrate ORDER BY observationoffset DESC) FILTER(WHERE heartrate BETWEEN 10 AND 300))[1] AS heart_rate,
      (array_agg(respiration ORDER BY observationoffset DESC) FILTER(WHERE respiration BETWEEN 1 AND 100))[1] AS rr_vital
      FROM vitalperiodic WHERE observationoffset>=0 AND observationoffset<11520
      AND patientunitstayid IN(SELECT DISTINCT id FROM fr)
      GROUP BY patientunitstayid,hr
    ) SELECT coalesce(f.id,s.id) AS id,coalesce(f.hr,s.hr) AS hr,fio2,sat,sat_n,sat_min,mbp,heart_rate,rr_vital
    FROM f FULL JOIN s USING(id,hr)""")

elif sys.argv[1]=='sicdb':
    # SICdb offsets include preceding surgery. Bin boundaries remain source-hour boundaries;
    # eligibility uses actual ICUOffset, not the case-open time.
    dest=R/'sicdb_hours.pkl'
    if not dest.exists():
        cases=pd.read_csv((Path(os.environ['SICDB_DIR']) / 'cases.csv.gz')).set_index('CaseID')
        ids={710:'sat',2283:'fio2',2278:'peep',2282:'rr_vent',719:'rr_vital',703:'mbp',707:'heart_rate'}
        chunks=[]; n=0
        for c in pd.read_csv((Path(os.environ['SICDB_DIR']) / 'data_float_h.csv.gz'),usecols=['id','CaseID','DataID','Offset','Val','cnt'],chunksize=750000):
            n+=len(c)
            c=c[c.DataID.isin(ids)].copy()
            offset=c.CaseID.map(cases.ICUOffset)
            c=c[(c.Offset>=offset-3600)&(c.Offset<offset+8*86400)]
            chunks.append(c)
            if n%7500000==0: print('sicdb scanned',n,flush=True)
        c=pd.concat(chunks,ignore_index=True).sort_values('id').drop_duplicates(['CaseID','DataID','Offset'],keep='last')
        s=c[c.DataID==710].set_index(['CaseID','Offset'])[['cnt']].rename(columns={'cnt':'sat_n'})
        c['variable']=c.DataID.map(ids)
        h=c.pivot(index=['CaseID','Offset'],columns='variable',values='Val').join(s).reset_index().rename(columns={'CaseID':'id'})
        h['hr']=h.Offset/3600; h['sat_min']=h.sat # hourly mean endpoint, not minimum or duration
        h.to_pickle(dest)
        print('sicdb_hours',h.shape,flush=True)
