import os
from pathlib import Path
import pandas as pd,numpy as np,time,json
O=Path(__file__).resolve().parent;R=O.parent/'audited_reanalysis_v2'
d=pd.read_pickle(R/'sicdb_cohort.pkl');hrs=d.set_index('id').hr
pieces=[];n=0;start=time.time()
for c in pd.read_csv((Path(os.environ['SICDB_DIR']) / 'data_float_h.csv.gz'),chunksize=150000,usecols=['id','CaseID','DataID','Offset','Val','cnt','rawdata']):
 n+=len(c);c=c[c.DataID.isin([710,2283])&c.CaseID.isin(hrs.index)].copy()
 if len(c):
  b=c.CaseID.map(hrs)*3600;c=c[c.Offset.ge(b)&c.Offset.lt(b+8*3600)];pieces.append(c)
 if n%3000000==0:print('rawdata scanned',n,'elapsed',round(time.time()-start),flush=True)
z=pd.concat(pieces,ignore_index=True).sort_values('id').drop_duplicates(['CaseID','DataID','Offset'],keep='last')
z.to_pickle(O/'sic_embedded_raw_selected.pkl')
print('selected signal hours',len(z),'elapsed',time.time()-start,flush=True)
