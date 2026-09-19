from pathlib import Path
import pandas as pd
O=Path(__file__).resolve().parent;R=O.parent/'audited_reanalysis_v2';rows=[]
for name,z in [('all_airway_records',pd.read_csv(R/'eicu_airway.csv')),('selected_records',pd.read_pickle(R/'eicu_cohort.pkl'))]:
 s=pd.to_numeric(z.start_hr,errors='coerce');e=pd.to_numeric(z.end_hr,errors='coerce')
 rows.append(dict(population=name,n=len(z),start_missing=int(s.isna().sum()),end_missing=int(e.isna().sum()),start_zero=int(s.eq(0).sum()),end_zero=int(e.eq(0).sum()),end_le_start=int((e<=s).sum()),valid_positive_interval=int((e>s).sum())))
pd.DataFrame(rows).to_csv(O/'eicu_interval_field_audit.csv',index=False)
