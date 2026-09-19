"""Correct AIPW with patient-cluster bootstrap that refits BOTH nuisance models.
All reanalyses are post-review; none is described as preregistered.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'; os.environ['OMP_NUM_THREADS']='1'; os.environ['MKL_NUM_THREADS']='1'
from pathlib import Path
import json, sys, warnings, time
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import SplineTransformer, StandardScaler
from sklearn.exceptions import ConvergenceWarning
from threadpoolctl import threadpool_limits
threadpool_limits(1)
R=Path(__file__).resolve().parent

def aipw(a,y,e,m1,m0):
    h1=m1+a/e*(y-m1)
    h0=m0+(1-a)/(1-e)*(y-m0)
    return h1,h0

def design(d,num,cat,catlevels,mode):
    arrays=[]
    for c in num:
        x=pd.to_numeric(d[c],errors='coerce').to_numpy(float)
        bad=~np.isfinite(x)
        if bad.all(): x=np.zeros_like(x)
        else: x=np.where(bad,np.nanmedian(x[~bad]),x)
        if mode=='spline' and c in ['fio2','sat','age','icu_hour'] and len(np.unique(x))>=5:
            z=SplineTransformer(n_knots=4,degree=3,knots='quantile',extrapolation='linear',include_bias=False).fit_transform(x[:,None])
        else: z=x[:,None]
        arrays.extend([z,bad.astype(float)[:,None]])
    for c in cat:
        levels=catlevels[c]
        x=d[c].fillna('Missing').astype(str).to_numpy()
        arrays.append(np.column_stack([x==k for k in levels[1:]]).astype(float) if len(levels)>1 else np.zeros((len(d),1)))
    return StandardScaler().fit_transform(np.column_stack(arrays))

def estimate(d,outcome,mode,num,cat,levels,return_detail=False):
    a=d.A.to_numpy(float); y=d[outcome].to_numpy(float)
    if min(np.unique(a).size,np.unique(y).size)<2: raise ValueError('single class')
    X=design(d,num,cat,levels,mode)
    with warnings.catch_warnings():
        warnings.simplefilter('error',ConvergenceWarning)
        ps=LogisticRegression(C=1,max_iter=2500,tol=1e-7).fit(X,a)
        eraw=ps.predict_proba(X)[:,1]
        e=np.clip(eraw,.01,.99)
        om=LogisticRegression(C=1,max_iter=2500,tol=1e-7).fit(np.column_stack([a,X]),y)
        m1=om.predict_proba(np.column_stack([np.ones(len(d)),X]))[:,1]
        m0=om.predict_proba(np.column_stack([np.zeros(len(d)),X]))[:,1]
    h1,h0=aipw(a,y,e,m1,m0)
    w=a/e+(1-a)/(1-e)
    # IPTW is Hájek-normalized, no unreported percentile truncation. e bound is explicit.
    ip1=np.average(y[a==1],weights=w[a==1]); ip0=np.average(y[a==0],weights=w[a==0])
    est=dict(aipw_rd=h1.mean()-h0.mean(),risk1=h1.mean(),risk0=h0.mean(),iptw_rd=ip1-ip0,iptw_risk1=ip1,iptw_risk0=ip0)
    if return_detail: return est,eraw,e,w,h1-h0
    return est

def smd(x,a,w):
    ok=np.isfinite(x); x=x[ok]; aa=a[ok]; ww=w[ok]
    if not ((aa==0).any() and (aa==1).any()): return np.nan
    # Fixed unweighted pooled variance denominator for before/after comparability.
    den=np.sqrt((np.var(x[aa==1])+np.var(x[aa==0]))/2)
    diff=np.average(x[aa==1],weights=ww[aa==1])-np.average(x[aa==0],weights=ww[aa==0])
    return diff/den if den>0 else (0 if abs(diff)<1e-10 else np.nan)

def diagnostics(db,d,num,cat,e,w,eraw):
    rows=[];a=d.A.to_numpy()
    for c in num:
        x=pd.to_numeric(d[c],errors='coerce').to_numpy(float)
        row=dict(database=db,variable=c,missing_n=int((~np.isfinite(x)).sum()),missing_pct=100*np.mean(~np.isfinite(x)),
                 smd_before=smd(x,a,np.ones(len(d))),smd_after=smd(x,a,w))
        for val,lab in [(1,'deescalation'),(0,'stable')]:
            z=x[a==val]; z=z[np.isfinite(z)]
            row[lab+'_n']=int((a==val).sum()); row[lab+'_missing_n']=int(((a==val)&~np.isfinite(x)).sum())
            for q,label in [(25,'q25'),(50,'median'),(75,'q75')]: row[lab+'_'+label]=np.percentile(z,q) if len(z) else np.nan
        rows.append(row)
        if (~np.isfinite(x)).any():
            b=(~np.isfinite(x)).astype(float)
            rows.append(dict(database=db,variable=c+' missing',smd_before=smd(b,a,np.ones(len(d))),smd_after=smd(b,a,w)))
    for c in cat:
        for level in sorted(d[c].fillna('Missing').astype(str).unique()):
            b=d[c].fillna('Missing').astype(str).eq(level).to_numpy(float)
            rows.append(dict(database=db,variable=c+': '+level,smd_before=smd(b,a,np.ones(len(d))),smd_after=smd(b,a,w),
                deescalation_pct=100*b[a==1].mean(),stable_pct=100*b[a==0].mean()))
    pd.DataFrame(rows).to_csv(R/f'{db}_balance.csv',index=False)
    q=d[['id','pid','A']+num+cat].copy();q['ps_raw']=eraw;q['ps_bounded']=e;q['iptw']=w;q.to_pickle(R/f'{db}_model_details.pkl')
    # Source data for plots: fixed bins, no patient identifiers.
    hist=[]
    for aa in [0,1]:
        counts,edges=np.histogram(eraw[a==aa],bins=np.linspace(0,1,41))
        hist.extend(dict(database=db,A=aa,left=edges[i],right=edges[i+1],n=int(v)) for i,v in enumerate(counts))
    pd.DataFrame(hist).to_csv(R/f'{db}_ps_histogram.csv',index=False)
    diag=dict(database=db,n=len(d),patients=d.pid.nunique(),treated_n=int(a.sum()),ps_clip_low=int((eraw<.01).sum()),ps_clip_high=int((eraw>.99).sum()),
              overlap_retained=int(((eraw>=.10)&(eraw<=.90)).sum()),max_abs_smd=float(pd.DataFrame(rows).smd_after.abs().max()))
    for aa in [0,1]:
        ww=w[a==aa];diag[f'ess_{aa}']=float(ww.sum()**2/(ww@ww));diag[f'weight_max_{aa}']=float(ww.max());diag[f'weight_p99_{aa}']=float(np.quantile(ww,.99))
    (R/f'{db}_diagnostics.json').write_text(json.dumps(diag,indent=2))
    # Retained/excluded characterization, uses actual model support and includes treatment prevalence.
    support=[]
    for retained,mask in [(True,(eraw>=.10)&(eraw<=.90)),(False,(eraw<.10)|(eraw>.90))]:
        z=d[mask]
        for c in num+['A','hospital_mortality']:
            xx=pd.to_numeric(z[c],errors='coerce')
            support.append(dict(database=db,retained=retained,n=len(z),variable=c,mean=xx.mean(),median=xx.median(),q25=xx.quantile(.25),q75=xx.quantile(.75)))
    pd.DataFrame(support).to_csv(R/f'{db}_support_comparison.csv',index=False)

def run(db,B):
    d=pd.read_pickle(R/f'{db}_cohort.pkl').reset_index(drop=True)
    num=[c for c in ['fio2','sat','age','icu_hour','peep','rr_vent','rr_vital','mbp','heart_rate','fio2_trend','sat_trend'] if c in d and d[c].notna().any()]
    # SICdb HospitalUnit is the last unit using the case, so exclude it as potentially post-treatment.
    cat=[c for c in ['sex','calendar_group']+(['unit'] if db!='sicdb' else [])+(['site'] if db=='eicu' else []) if c in d]
    levels={c:sorted(d[c].fillna('Missing').astype(str).unique()) for c in cat}
    (R/f'{db}_model_spec.json').write_text(json.dumps(dict(numeric=num,categorical=cat,levels=levels,primary='cubic B-splines, 4 quantile knots; logistic C=1; no treatment interactions',ps_guard=[.01,.99],bootstrap_refits='median, spline, scale, propensity and outcome'),indent=2))
    base,eraw,e,w,dr=estimate(d,'hospital_mortality','spline',num,cat,levels,True)
    diagnostics(db,d,num,cat,e,w,eraw)
    configs=[('hospital_primary','hospital_mortality','spline',np.ones(len(d),bool)),
             ('hospital_linear','hospital_mortality','linear',np.ones(len(d),bool)),
             ('recorded28_secondary','death28_recorded','spline',np.ones(len(d),bool)),
             ('low88_observed','low88_any6','spline',d.low88_any6.notna().to_numpy()),
             ('low90_observed','low90_any6','spline',d.low90_any6.notna().to_numpy()),
             ('low88_six_bins','low88_any6','spline',d.n_observed_bins6.eq(6).to_numpy()),
             ('hospital_overlap','hospital_mortality','spline',(eraw>=.1)&(eraw<=.9))]
    if db=='sicdb': configs.append(('hospital_no_transfer','hospital_mortality','spline',~d.transfer.to_numpy()))
    allres=[]; allbs=[]
    for n,(name,out,mode,mask) in enumerate(configs):
        q=d[mask].reset_index(drop=True)
        if len(q)<100 or min(q.A.value_counts().min(),q[out].value_counts().min())<5:
            allres.append(dict(database=db,analysis=name,n=len(q),status='insufficient events'));continue
        point=estimate(q,out,mode,num,cat,levels)
        groups=q.groupby('pid',dropna=False).indices; ids=list(groups);rng=np.random.default_rng(71401+n)
        bs=[];fail=[]
        for b in range(B):
            chosen=rng.integers(0,len(ids),len(ids));ix=np.concatenate([groups[ids[j]] for j in chosen])
            try: z=estimate(q.iloc[ix],out,mode,num,cat,levels);bs.append(z);allbs.append(dict(database=db,analysis=name,replicate=b,**z))
            except (ValueError,ConvergenceWarning) as exc: fail.append(str(exc))
            if (b+1)%100==0: print(db,name,b+1,flush=True)
        rr=dict(database=db,analysis=name,n=len(q),patients=q.pid.nunique(),treated_n=int(q.A.sum()),events=int(q[out].sum()),**point,
                bootstrap_requested=B,bootstrap_success=len(bs),bootstrap_failures=len(fail),status='complete' if len(bs)>=.95*B else 'failed uncertainty gate')
        for key in ['aipw_rd','iptw_rd','risk1','risk0']:
            rr[key+'_low'],rr[key+'_high']=np.quantile([x[key] for x in bs],[.025,.975])
        allres.append(rr)
        pd.DataFrame(allres).to_csv(R/f'{db}_results.csv',index=False)
        pd.DataFrame(allbs).to_csv(R/f'{db}_bootstrap.csv',index=False)
        print(db,name,'RD',rr['aipw_rd'],'CI',rr['aipw_rd_low'],rr['aipw_rd_high'],flush=True)
    density=[]
    for a,z in d.groupby('A'):
        density.append(dict(database=db,A=int(a),n=len(z),no_measurement_n=int(z.n_observed_bins6.eq(0).sum()),at_least_one_pct=100*z.n_observed_bins6.gt(0).mean(),six_bins_pct=100*z.n_observed_bins6.eq(6).mean(),
                measurement_median=z.n_measurements6.median(),measurement_q25=z.n_measurements6.quantile(.25),measurement_q75=z.n_measurements6.quantile(.75),
                bins_median=z.n_observed_bins6.median(),bins_q25=z.n_observed_bins6.quantile(.25),bins_q75=z.n_observed_bins6.quantile(.75)))
    pd.DataFrame(density).to_csv(R/f'{db}_measurement_density.csv',index=False)

if __name__=='__main__': run(sys.argv[1],int(sys.argv[2]) if len(sys.argv)>2 else 300)
