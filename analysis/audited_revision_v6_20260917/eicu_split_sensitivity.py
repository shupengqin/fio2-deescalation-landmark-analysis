from pathlib import Path
import ast,sys,json
import numpy as np,pandas as pd
V=Path(__file__).resolve().parent.parent/'audited_revision_v5_20260917';O=Path(__file__).resolve().parent;R=V.parent/'audited_reanalysis_v2';sys.path.insert(0,str(V))
tree=ast.parse((V/'crossfit_models.py').read_text(encoding='utf8'));tree.body=[x for x in tree.body if isinstance(x,(ast.Import,ast.ImportFrom,ast.FunctionDef))];ns={};exec(compile(tree,'defs','exec'),ns)
d=pd.read_pickle(R/'eicu_cohort.pkl').reset_index(drop=True);num,cat,levels=ns['spec']('eicu');ids=np.array(sorted(d.pid.unique()))
rows=[]
for seed in range(917501,917521):
 rng=np.random.default_rng(seed);perm=rng.permutation(len(ids));fmap={ids[i]:int(j%5) for j,i in enumerate(perm)};d['fold']=d.pid.map(fmap)
 point,w=ns['fit'](d,num,cat,levels);a=d.A.to_numpy();rows.append(dict(seed=seed,**point,max_abs_smd_numeric=float(np.nanmax(np.abs([ns['smd'](pd.to_numeric(d[c],errors='coerce').to_numpy(float),a,w) for c in num]))),treated_ess=float((w[a==1].sum()**2)/(w[a==1]@w[a==1]))))
pd.DataFrame(rows).to_csv(O/'eicu_20_split_sensitivity.csv',index=False)
q=pd.DataFrame(rows);summary=pd.DataFrame([dict(n_splits=len(q),rd_mean=q.rd.mean(),rd_sd=q.rd.std(ddof=1),rd_min=q.rd.min(),rd_max=q.rd.max(),risk1_mean=q.risk1.mean(),risk0_mean=q.risk0.mean(),ess_median=q.treated_ess.median())]);summary.to_csv(O/'eicu_20_split_summary.csv',index=False);print(summary.to_string(index=False))
