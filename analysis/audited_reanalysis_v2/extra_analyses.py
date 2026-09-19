from analyze_verified import *

def overlap(q,out,num,cat,levels,mode='spline'):
    a=q.A.to_numpy(float);y=q[out].to_numpy(float)
    X=design(q,num,cat,levels,mode)
    ps=LogisticRegression(C=1,max_iter=2500,tol=1e-7).fit(X,a).predict_proba(X)[:,1]
    w=np.where(a==1,1-ps,ps)
    r1=np.average(y[a==1],weights=w[a==1]); r0=np.average(y[a==0],weights=w[a==0])
    return dict(rd=r1-r0,risk1=r1,risk0=r0),w

def run_extra(db,B=300):
    d=pd.read_pickle(R/f'{db}_cohort.pkl')
    spec=json.loads((R/f'{db}_model_spec.json').read_text());num=spec['numeric'];cat=spec['categorical'];levels=spec['levels']
    rr=[]
    # New overlap-population estimand chosen because positivity audit failed; do not relabel as ATE.
    point,w=overlap(d,'hospital_mortality',num,cat,levels)
    a=d.A.to_numpy()
    bal=[]
    for c in num:
        x=pd.to_numeric(d[c],errors='coerce').to_numpy(float)
        bal.append(dict(database=db,variable=c,smd_overlap=smd(x,a,w)))
    pd.DataFrame(bal).to_csv(R/f'{db}_overlap_weight_balance.csv',index=False)
    for name,q,out in [('hospital_overlap_weighted',d,'hospital_mortality')]:
        groups=q.groupby('pid').indices;ids=list(groups);rng=np.random.default_rng(77141);bs=[]
        for b in range(B):
            ix=np.concatenate([groups[ids[k]] for k in rng.integers(0,len(ids),len(ids))])
            z,_=overlap(q.iloc[ix],out,num,cat,levels);bs.append(z['rd'])
        rr.append(dict(database=db,analysis=name,n=len(q),rd=point['rd'],ci_low=np.quantile(bs,.025),ci_high=np.quantile(bs,.975),bootstrap=B))
    z=pd.read_pickle(R/f'{db}_repeated_candidates.pkl')
    z['delta']=z.fio2_next-z.fio2
    z=z[(z.delta.le(-10)|z.delta.between(-5,5))&z.hosp_death.notna()].copy()
    z['A']=z.delta.le(-10).astype(int);z['hospital_mortality']=z.hosp_death.astype(int)
    # For repeated rows, include the same declared covariates and fit models on all opportunities;
    # resample unique patients with all of their rows, preserving within-patient dependence.
    levels2={c:sorted(z[c].fillna('Missing').astype(str).unique()) for c in cat}
    point=estimate(z,'hospital_mortality','spline',num,cat,levels2)
    groups=z.groupby('pid').indices;ids=list(groups);rng=np.random.default_rng(78141);bs=[]
    for b in range(B):
        ix=np.concatenate([groups[ids[k]] for k in rng.integers(0,len(ids),len(ids))])
        bs.append(estimate(z.iloc[ix],'hospital_mortality','spline',num,cat,levels2)['aipw_rd'])
        if (b+1)%100==0: print(db,'repeated bootstrap',b+1,flush=True)
    rr.append(dict(database=db,analysis='hospital_repeated_secondary',n=len(z),rd=point['aipw_rd'],ci_low=np.quantile(bs,.025),ci_high=np.quantile(bs,.975),bootstrap=B))
    pd.DataFrame(rr).to_csv(R/f'{db}_extra_results.csv',index=False)
    print(db,rr,flush=True)

if __name__=='__main__':run_extra(sys.argv[1])
