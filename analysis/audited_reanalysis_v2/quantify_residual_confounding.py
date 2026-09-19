from pathlib import Path
import pandas as pd,numpy as np
R=Path(__file__).resolve().parent
def evalue(rr):
    rr=max(rr,1/rr)
    return rr+np.sqrt(rr*(rr-1))
rows=[]
for db in ['mimic','eicu','sicdb']:
    r=pd.read_csv(R/f'{db}_results.csv').query("analysis=='hospital_primary'").iloc[0]
    b=pd.read_csv(R/f'{db}_bootstrap.csv').query("analysis=='hospital_primary'")
    assert len(b)==300 and b.risk0.gt(0).all() and b.risk1.gt(0).all()
    rr=r.risk1/r.risk0;lo,hi=np.quantile(b.risk1/b.risk0,[.025,.975])
    ci_e=1 if lo<=1<=hi else evalue(hi if hi<1 else lo)
    rows.append(dict(database=db,rr=rr,rr_low=lo,rr_high=hi,evalue_point=evalue(rr),evalue_ci=ci_e))
pd.DataFrame(rows).to_csv(R/'evalue_sensitivity.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
