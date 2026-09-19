from analysis_utils import *
checks=[]
for stem in ['minute_model','rare_event','support_model','additional_support']:
 r=pd.read_csv(O/(stem+'_results.csv'));b=pd.read_csv(O/(stem+'_replicates.csv'))
 keys=['analysis','definition'] if stem=='minute_model' else ['analysis']
 for _,row in r.iterrows():
  if row.get('status','complete')!='complete':continue
  q=b.copy()
  q=q[q.analysis.eq(row.analysis)]
  value=('difference' if row.definition=='paired_difference' else row.definition+'_aipw_rd') if stem=='minute_model' else ('rd' if 'rd' in q else 'aipw_rd')
  assert len(q)==1000 and q.replicate.nunique()==1000,(stem,row[keys].to_dict())
  assert np.allclose(np.quantile(q[value],[.025,.975]),[row.ci_low,row.ci_high])
  risk_cols=([row.definition+'_risk1',row.definition+'_risk0'] if stem=='minute_model' and row.definition!='paired_difference' else ['risk1','risk0'] if 'risk1' in q else [])
  if stem=='minute_model':assert np.allclose(q.difference,q.A_minute_aipw_rd-q.A_aipw_rd)
  checks.append(dict(analysis=str(row['analysis']),definition=row.get('definition',''),n_replicates=len(q),risk_outside_n=int(((q[risk_cols]<0)|(q[risk_cols]>1)).any(axis=1).sum()) if risk_cols else None))
cf=[]
for db in ['mimic','eicu','sicdb']:
 r=pd.read_csv(O/f'{db}_crossfit_result.csv').iloc[0];b=pd.read_csv(O/f'{db}_crossfit_replicates.csv')
 assert len(b)==300 and b.replicate.nunique()==300 and np.allclose(np.quantile(b.rd,[.025,.975]),[r.ci_low,r.ci_high])
 outside=((b[['risk1','risk0']]<0)|(b[['risk1','risk0']]>1)).any(axis=1)
 cf.append(dict(database=db,replicates=300,risk_outside_n=int(outside.sum()),minimum_risk=b[['risk1','risk0']].min().min(),maximum_risk=b[['risk1','risk0']].max().max()))
pd.DataFrame(cf).to_csv(O/'crossfit_risk_range_audit.csv',index=False)
# Original rare-event risk behavior, reported without clipping or dropping replicates.
b=pd.read_csv(R/'sicdb_bootstrap.csv');q=b[b.analysis.eq('low88_six_bins')]
rr=dict(analysis='original_sicdb_six_bin_low88',replicates=len(q),risk_outside_n=int(((q[['risk1','risk0']]<0)|(q[['risk1','risk0']]>1)).any(axis=1).sum()),minimum_risk=q[['risk1','risk0']].min().min(),maximum_risk=q[['risk1','risk0']].max().max())
pd.DataFrame([rr]).to_csv(O/'rare_original_risk_range_audit.csv',index=False)
pd.DataFrame(checks).to_csv(O/'bootstrap_consistency_checks.csv',index=False)
print(pd.DataFrame(cf).to_string(index=False));print(rr);print('All new aggregate confidence intervals reproduced from saved replicates.')
