"""Source-specific exclusions use only evidence at/before the decision; refit inference."""
from analyze_verified import *
import os,psycopg2

def extract():
    queries={
      'mimic':('mimiciv31',"""SELECT c.stay_id AS id,extract(epoch FROM(c.charttime-i.intime))/3600 AS hr,c.itemid,c.value,c.valuenum
      FROM mimiciv_icu.chartevents c JOIN mimiciv_icu.icustays i USING(stay_id)
      WHERE c.itemid IN(223758,228687,229784,229270,229277,229280)
      AND c.charttime<=i.intime+interval '8 days'"""),
      'eicu':('eicu',"""SELECT patientunitstayid AS id,cplitemoffset/60.0 AS hr,cplitemvalue AS value
      FROM public.careplangeneral WHERE cplgroup='Care Limitation' AND cplitemoffset<=8*1440""")}
    for db,(database,query) in queries.items():
      dest=R/f'{db}_exclusion_events.csv'
      if dest.exists():continue
      (R/f'{db}_exclusion_events.sql').write_text(query)
      c=psycopg2.connect(host=os.environ.get('PGHOST','127.0.0.1'),port=int(os.environ.get('PGPORT','5432')),user=os.environ.get('PGUSER','postgres'),password=os.environ['PGPASSWORD'],dbname={'mimiciv31':os.environ.get('MIMIC_DATABASE','mimiciv31'),'eicu':os.environ.get('EICU_DATABASE','eicu')}.get((database), (database)));c.set_session(readonly=True)
      with c.cursor() as cur,dest.open('w',newline='',encoding='utf8') as f:cur.copy_expert('COPY ('+query+') TO STDOUT WITH CSV HEADER',f)
      c.close();print(db,'exclusion event extraction complete',flush=True)
    dest=R/'sicdb_ecmo_events.csv'
    if not dest.exists():
      chunks=[]
      for i,z in enumerate(pd.read_csv((Path(os.environ['SICDB_DIR']) / 'data_float_h.csv.gz'),chunksize=1000000,usecols=['CaseID','DataID','Offset','Val'])):
        q=z[z.DataID.isin([2023,2024,2025])&z.Val.gt(0)].copy()
        if len(q):chunks.append(q)
        if (i+1)%10==0:print('SICdb ECMO source scan',i+1,'million rows',flush=True)
      pd.concat(chunks).rename(columns={'CaseID':'id'}).to_csv(dest,index=False)
      print('SICdb exclusion event extraction complete',flush=True)

def flags(db,d):
    f=pd.DataFrame(index=d.index)
    f['recent_support']=False;f['ecmo_before_decision']=False;f['limitation_before_decision']=False
    if db=='mimic':
      v=pd.read_csv(R/'mimic_vent.csv');v=v[v.ventilation_status.eq('InvasiveVent')]
    elif db=='sicdb':
      v=pd.read_csv((Path(os.environ['SICDB_DIR']) / 'data_range.csv.gz'));v=v[v.DataID.isin([720,3041])].drop(columns='id').rename(columns={'CaseID':'id'})
      v['start_hr']=v.Offset/3600;v['end_hr']=v.OffsetEnd/3600
    else:v=d[['id','start_hr','end_hr']].copy()
    vg=v.groupby('id')
    for ix,z in d.iterrows():
      if z.id in vg.indices:
        u=vg.get_group(z.id)
        f.loc[ix,'recent_support']=bool(((u.start_hr<=z.decision_hr)&(u.start_hr>z.decision_hr-2)).any())
    if db in ['mimic','eicu']:
      e=pd.read_csv(R/f'{db}_exclusion_events.csv').sort_values('hr');g=e.groupby('id')
      for ix,z in d.iterrows():
        if z.id not in g.indices:continue
        u=g.get_group(z.id);u=u[u.hr<=z.decision_hr]
        if db=='mimic':
          f.loc[ix,'ecmo_before_decision']=bool((u.itemid.isin([229270,229277,229280])&u.valuenum.gt(0)).any())
          q=u[u.itemid.isin([223758,228687,229784])]
          if len(q):f.loc[ix,'limitation_before_decision']=bool(__import__('re').search(r'dnr|dni|comfort|no cpr|do not|cmo',str(q.iloc[-1]['value']),__import__('re').I))
        else:
          # Conservative exclusion if any specific limitation was documented before the decision.
          f.loc[ix,'limitation_before_decision']=bool(u.value.isin(['Comfort measures only','Do not resuscitate','No augmentation of care','No blood draws','No blood products','No cardioversion','No CPR','No intubation','No vasopressors/inotropes']).any())
    else:
      # Hourly means become baseline information only at the end of the source bin.
      e=pd.read_csv(R/'sicdb_ecmo_events.csv');e['hr']=e.Offset/3600+1;g=e.groupby('id')
      for ix,z in d.iterrows():
        if z.id in g.indices:f.loc[ix,'ecmo_before_decision']=bool(g.get_group(z.id).hr.le(z.decision_hr).any())
    return f

def run_exclusion(db,B=300):
    d=pd.read_pickle(R/f'{db}_cohort.pkl').reset_index(drop=True)
    f=flags(db,d);d['source_exclusion']=f.any(axis=1)
    spec=json.loads((R/f'{db}_model_spec.json').read_text());num=spec['numeric'];cat=spec['categorical'];levels=spec['levels']
    audit=[dict(database=db,criterion=c,n=int(f[c].sum()),source_available=not ((db=='eicu' and c=='ecmo_before_decision') or (db=='sicdb' and c=='limitation_before_decision'))) for c in f]
    audit.append(dict(database=db,criterion='any_source_specific_exclusion',n=int(d.source_exclusion.sum()),source_available=True))
    pd.DataFrame(audit).to_csv(R/f'{db}_exclusion_counts.csv',index=False)
    rows=[]
    configs=[('hospital_source_exclusions','hospital_mortality',~d.source_exclusion),
             ('low88_source_exclusions','low88_any6',~d.source_exclusion&d.low88_any6.notna()),
             ('hospital_complete_numeric','hospital_mortality',d[num].notna().all(axis=1))]
    for no,(name,out,mask) in enumerate(configs):
      q=d[mask].reset_index(drop=True)
      if len(q)<100 or q.A.nunique()<2 or q[out].nunique()<2 or min(q.A.value_counts().min(),q[out].value_counts().min())<5:
        rows.append(dict(database=db,analysis=name,n=len(q),status='insufficient sample/events'));continue
      point=estimate(q,out,'spline',num,cat,levels)
      groups=q.groupby('pid').indices;ids=list(groups);rng=np.random.default_rng(79411+no);bs=[]
      for b in range(B):
        ix=np.concatenate([groups[ids[j]] for j in rng.integers(0,len(ids),len(ids))])
        bs.append(estimate(q.iloc[ix],out,'spline',num,cat,levels)['aipw_rd'])
        if (b+1)%100==0:print(db,name,b+1,flush=True)
      rows.append(dict(database=db,analysis=name,n=len(q),patients=q.pid.nunique(),treated_n=int(q.A.sum()),events=int(q[out].sum()),aipw_rd=point['aipw_rd'],aipw_rd_low=np.quantile(bs,.025),aipw_rd_high=np.quantile(bs,.975),bootstrap_success=len(bs),status='complete'))
      pd.DataFrame(rows).to_csv(R/f'{db}_source_sensitivity.csv',index=False)
    pd.DataFrame(rows).to_csv(R/f'{db}_source_sensitivity.csv',index=False)

if __name__=='__main__':
  if sys.argv[1]=='extract':extract()
  else:run_exclusion(sys.argv[1])
