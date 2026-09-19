from pathlib import Path
import pandas as pd,numpy as np,json
O=Path(__file__).resolve().parent;R=O.parent/'audited_reanalysis_v2'
checks=[]
for db in ['mimic','eicu','sicdb']:
 a=pd.read_pickle(R/f'{db}_cohort.pkl');b=pd.read_pickle(O/f'{db}_cohort.pkl');checks.append(dict(database=db,check='reconstructed cohort identities and treatment match',passed=a[['id','pid','hr','A','hospital_mortality']].reset_index(drop=True).equals(b[['id','pid','hr','A','hospital_mortality']].reset_index(drop=True))))
 d=pd.read_pickle(O/f'{db}_augmented_cohort.pkl');checks.append(dict(database=db,check='all selected times at or after first physiology',passed=bool(d.first_to_selected_delay_hr.ge(0).all())))
 if db!='sicdb':
  checks.append(dict(database=db,check='record intervals positive and below 120 min',passed=bool(d.setting_gap_min.gt(0).all()&d.setting_gap_min.lt(120).all())))
  checks.append(dict(database=db,check='new lab sample ages between zero and six hours',passed=bool(all(d[c].dropna().between(0,6).all() for c in ['po2_age_hr','pco2_age_hr','ph_age_hr','lactate_age_hr']))))
pd.DataFrame(checks).to_csv(O/'revision_data_checks.csv',index=False)
print(pd.DataFrame(checks).to_string(index=False));assert all(c['passed'] for c in checks)
