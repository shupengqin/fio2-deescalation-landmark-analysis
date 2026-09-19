from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,sys,platform
R=Path(__file__).resolve().parent
dbs=['mimic','eicu','sicdb'];display=dict(mimic='MIMIC-IV',eicu='eICU-CRD',sicdb='SICdb')
assertions=[];summaries=[];follow=[];safety=[]
for db in dbs:
    d=pd.read_pickle(R/f'{db}_cohort.pkl')
    checks={
      'one row per stay':not d.id.duplicated().any(),
      'stable comparator never >5 pp increase':d.loc[d.A.eq(0),'delta'].between(-5,5).all(),
      'de-escalation threshold':d.loc[d.A.eq(1),'delta'].le(-10).all(),
      'alive through landmark':(d.death_hr.isna()|d.death_hr.gt(d.landmark_hr)).all(),
      'records with missing next-bin saturation retained':d.sat_next.isna().any(),
      'death coding binary':d.hospital_mortality.isin([0,1]).all(),
      'hour grid':((d.landmark_hr-d.decision_hr)==1).all(),
    }
    for key,val in checks.items():assertions.append(dict(database=db,check=key,passed=bool(val)))
    assert all(checks.values()),db
    summ=dict(database=db,n=len(d),patients=d.pid.nunique(),treated=int(d.A.sum()),treatment_pct=d.A.mean()*100,
              mortality_n=int(d.hospital_mortality.sum()),mortality_pct=100*d.hospital_mortality.mean(),
              recorded28_n=int(d.death28_recorded.sum()),recorded28_pct=100*d.death28_recorded.mean(),
              next_sat_missing_n=int(d.sat_next.isna().sum()),sites=d.site.nunique() if 'site' in d else 1,
              calendar_groups='; '.join(map(str,sorted(d.calendar_group.dropna().unique()))))
    summaries.append(summ)
    for a,z in d.groupby('A'):
        for outcome in ['low88_any6','low90_any6']:
            q=z[z[outcome].notna()]
            safety.append(dict(database=db,A=int(a),outcome=outcome,n=len(q),events=int(q[outcome].sum()),event_pct=q[outcome].mean()*100))
        follow.append(dict(database=db,A=int(a),n=len(z),hospital_death_n=int(z.hospital_mortality.sum()),
                           recorded28_death_n=int(z.death28_recorded.sum()),discharged_alive_before28_n=int(z.alive_discharge_before28.sum()),
                           ascertainment='1-year/6-month recorded mortality, distinct from hospital discharge' if db=='sicdb' else 'Hospital-disposition-based mortality; no post-discharge ascertainment used'))
pd.DataFrame(summaries).to_csv(R/'cohort_summary.csv',index=False)
pd.DataFrame(follow).to_csv(R/'outcome_ascertainment.csv',index=False)
pd.DataFrame(safety).to_csv(R/'safety_counts_by_treatment.csv',index=False)
pd.DataFrame(assertions).to_csv(R/'cohort_integrity_checks.csv',index=False)
pd.concat([pd.read_csv(R/f'{db}_flow.csv') for db in dbs]).to_csv(R/'flow_all.csv',index=False)
for typ in ['balance','support_comparison','results','measurement_density','extra_results','source_sensitivity','exclusion_counts']:
    fs=[R/f'{db}_{typ}.csv' for db in dbs]
    if all(f.exists() for f in fs):pd.concat([pd.read_csv(f) for f in fs]).to_csv(R/f'all_{typ}.csv',index=False)
# A reproducibility manifest contains hashes and file names only, never patient data.
fs=list(R.glob('*.py'))+list(R.glob('*.sql'))+list(R.glob('*_cohort.pkl'))
man=[]
for f in fs:man.append(dict(filename=f.name,sha256=hashlib.file_digest(f.open('rb'),'sha256').hexdigest(),bytes=f.stat().st_size))
pd.DataFrame(man).to_csv(R/'analysis_manifest.csv',index=False)
import sklearn, scipy
(R/'software_versions.json').write_text(json.dumps(dict(python=sys.version,pandas=pd.__version__,numpy=np.__version__,sklearn=sklearn.__version__,scipy=scipy.__version__),indent=2))
print(pd.DataFrame(summaries).to_string(index=False))
