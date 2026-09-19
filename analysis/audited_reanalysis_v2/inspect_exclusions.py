from pathlib import Path
import os,psycopg2,pandas as pd
R=Path(__file__).resolve().parent
queries={
 'mimiciv31':["SELECT itemid,label,linksto FROM mimiciv_icu.d_items WHERE lower(label) ~ '(ecmo|extracorp|code status|intubat)' ORDER BY itemid"],
 'eicu':["SELECT column_name FROM information_schema.columns WHERE table_name='careplangeneral' ORDER BY ordinal_position", "SELECT DISTINCT cplgroup,cplitemvalue FROM public.careplangeneral WHERE lower(cplgroup) ~ '(resusc|code status|care limit)' LIMIT 100", "SELECT DISTINCT treatmentstring FROM public.treatment WHERE lower(treatmentstring) ~ '(ecmo|extracorp|comfort|withdrawal|do not resusc)' LIMIT 70"]}
for db,qs in queries.items():
 c=psycopg2.connect(host=os.environ.get('PGHOST','127.0.0.1'),port=int(os.environ.get('PGPORT','5432')),user=os.environ.get('PGUSER','postgres'),password=os.environ['PGPASSWORD'],dbname={'mimiciv31':os.environ.get('MIMIC_DATABASE','mimiciv31'),'eicu':os.environ.get('EICU_DATABASE','eicu')}.get((db), (db)));c.set_session(readonly=True)
 for q in qs:
  try:
   with c.cursor() as cur:
    cur.execute(q);print(db,cur.fetchall())
  except Exception as exc:c.rollback();print(type(exc).__name__,str(exc)[:200])
 c.close()
d=pd.read_csv((Path(os.environ['SICDB_DIR']) / 'd_references.csv.gz'));mask=d.astype(str).apply(lambda x:x.str.contains('ECMO|extracorp|limitation|resusc|comfort|intubat',case=False,regex=True)).any(axis=1)
print(d[mask].to_string(index=False))
