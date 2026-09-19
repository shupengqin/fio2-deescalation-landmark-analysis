"""Quantitative figures use audited CSV outputs only, with no hand-entered study arrays.
Contract: one flow figure and three 4-panel figures show uncertainty, measurement coverage,
and support limitations. They do not claim a mortality benefit or an individualized rule.
Python/matplotlib, 183 mm width, editable PDF/SVG and 600 dpi TIFF plus preview PNG.
"""
from pathlib import Path
import pandas as pd,numpy as np,json,textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
R=Path(__file__).resolve().parent;OUT=R/'figures';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':['Arial','DejaVu Sans'],'font.size':8,'axes.titlesize':9,'axes.labelsize':8,
    'xtick.labelsize':7,'ytick.labelsize':8,'pdf.fonttype':42,'svg.fonttype':'none',
    'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'legend.frameon':False})
dbs=['mimic','eicu','sicdb'];names={'mimic':'MIMIC-IV','eicu':'eICU-CRD','sicdb':'SICdb'}
colors={'mimic':'#355B78','eicu':'#BB763D','sicdb':'#568B83'}
res=pd.read_csv(R/'all_results.csv');extra=pd.read_csv(R/'all_extra_results.csv')
summary=pd.read_csv(R/'cohort_summary.csv').set_index('database')

def save(fig,name):
    fig.savefig(OUT/(name+'.pdf'),facecolor='white')
    fig.savefig(OUT/(name+'.svg'),facecolor='white')
    fig.savefig(OUT/(name+'.tiff'),dpi=600,facecolor='white')
    fig.savefig(OUT/(name+'.png'),dpi=300,facecolor='white')
    plt.close(fig)

def panel(ax,letter,title):
    ax.set_title(title,loc='left',pad=12)
    ax.text(-.10,1.08,letter,transform=ax.transAxes,fontweight='bold',fontsize=11)

def forest(ax,analysis,letter,title,extras=False):
    t=(extra if extras else res);t=t[t.analysis.eq(analysis)].set_index('database')
    assert all(db in t.index for db in dbs),analysis
    for i,db in enumerate(dbs):
        z=t.loc[db];m=z['rd' if extras else 'aipw_rd']*100
        lo=z['ci_low' if extras else 'aipw_rd_low']*100;hi=z['ci_high' if extras else 'aipw_rd_high']*100
        ax.plot([lo,hi],[i,i],color=colors[db],lw=1.5);ax.scatter(m,i,color=colors[db],s=25,zorder=3)
    ax.axvline(0,color='#7D858C',lw=.8,ls='--');ax.set_yticks(range(3),[names[x] for x in dbs]);ax.set_ylim(2.7,-.7)
    ax.set_xlabel('Risk difference (percentage points)');ax.grid(axis='x',color='#E9ECEF',lw=.5)
    panel(ax,letter,title)

# Figure 1: a single connected flow, with denominators explicitly labeled.
flow=pd.read_csv(R/'flow_all.csv');fig,ax=plt.subplots(figsize=(7.2,8.15));ax.axis('off');ax.set_xlim(0,1);ax.set_ylim(0,1)
ax.text(.5,.98,'Construction of the first-opportunity cohort',ha='center',fontsize=11,fontweight='bold')
labels=['Valid FiO₂ bins in source window','Adult age criterion','2 to <168 hours after ICU admission','Invasive ventilation / airway documented','Current FiO₂ ≥40%','Current saturation ≥96%','Next-bin FiO₂ available','Observed alive through landmark','First eligible opportunity per stay','Reduction ≥10 pp or stable ±5 pp','Known hospital disposition']
xs=[.405,.64,.875]
for x,db in zip(xs,dbs):ax.text(x,.925,names[db],ha='center',color=colors[db],fontweight='bold')
for i,lab in enumerate(labels):
    y=.865-i*.073
    ax.text(.015,y,textwrap.fill(lab,28),ha='left',va='center',fontsize=7.3)
    for x,db in zip(xs,dbs):
        z=flow[flow.database.eq(db)].iloc[i]
        box=FancyBboxPatch((x-.103,y-.025),.206,.05,boxstyle='round,pad=0.003,rounding_size=0.004',facecolor='#F4F7F9',edgecolor=colors[db],lw=.7)
        ax.add_patch(box)
        detail='recorded bins' if i==0 else f"removed {int(z.excluded_rows):,}"
        ax.text(x,y+.007,f'{int(z.rows):,}',ha='center',va='center',fontsize=8.4)
        ax.text(x,y-.011,detail,ha='center',va='center',fontsize=6.9,color='#53616B')
        if i<10:ax.annotate('',xy=(x,y-.048),xytext=(x,y-.028),arrowprops=dict(arrowstyle='->',lw=.6,color='#7D858C'))
ax.text(.01,.048,'Rows are recorded hourly opportunities until the first-opportunity step. Final counts are ICU stays/cases.',fontsize=7)
ax.text(.01,.030,'Flow begins with valid FiO₂ records, not with all database admissions. Full counts are in Supplementary Table S1.',fontsize=7)
fig.subplots_adjust(left=.03,right=.99,top=.99,bottom=.01);save(fig,'Figure_1')

fig,axes=plt.subplots(2,2,figsize=(7.2,6.15));fig.subplots_adjust(left=.115,right=.985,bottom=.10,top=.91,wspace=.55,hspace=.80)
forest(axes[0,0],'hospital_primary','A','In-hospital mortality: primary')
forest(axes[0,1],'recorded28_secondary','B','Recorded 28-day death: secondary')
forest(axes[1,0],'hospital_linear','C','Linear nuisance models')
forest(axes[1,1],'hospital_overlap_weighted','D','Hospital death: overlap population',True)
save(fig,'Figure_2')

density=pd.read_csv(R/'all_measurement_density.csv')
fig,axes=plt.subplots(2,2,figsize=(7.2,6.15));fig.subplots_adjust(left=.115,right=.985,bottom=.17,top=.91,wspace=.55,hspace=.85)
forest(axes[0,0],'low88_observed','A','Recorded hypoxaemia <88%')
forest(axes[0,1],'low90_observed','B','Recorded hypoxaemia <90%')
ax=axes[1,0];panel(ax,'C','Six observed hourly bins (%)')
for a,col,label in [(1,'#355B78','Reduction ≥10 pp'),(0,'#ACB6BE','Stable ±5 pp')]:
    yy=np.arange(3)+(a-.5)*.28
    v=[density[density.database.eq(db)&density.A.eq(a)].iloc[0].six_bins_pct for db in dbs]
    ax.barh(yy,v,height=.25,color=col,label=label)
ax.set_yticks(range(3),[names[x] for x in dbs]);ax.invert_yaxis();ax.set_xlim(0,100);ax.set_xlabel('Percentage of primary cohort')
fig.legend(*ax.get_legend_handles_labels(),fontsize=7,loc='lower center',ncol=2,bbox_to_anchor=(.55,.025))
ax=axes[1,1];panel(ax,'D','Six-hour measurement counts')
for a,col,label in [(1,'#355B78','Reduction ≥10 pp'),(0,'#ACB6BE','Stable ±5 pp')]:
    for i,db in enumerate(dbs):
        z=density[density.database.eq(db)&density.A.eq(a)].iloc[0];y=i+(a-.5)*.25
        ax.plot([z.measurement_q25,z.measurement_q75],[y,y],color=col,lw=1.5)
        ax.scatter(z.measurement_median,y,s=16,color=col)
ax.set_yticks(range(3),[names[x] for x in dbs]);ax.invert_yaxis();ax.set_xlabel('Count: median and IQR (log scale)');ax.set_xscale('log')
assert density.measurement_q25.min()>0
ax.set_xlim(density.measurement_q25.min()*.65,density.measurement_q75.max()*1.4)
save(fig,'Figure_3')

fig,axes=plt.subplots(2,2,figsize=(7.2,6.5));fig.subplots_adjust(left=.22,right=.985,bottom=.16,top=.91,wspace=.83,hspace=.78)
friendly={'fio2':'FiO₂','sat':'SpO₂','age':'Age','icu_hour':'ICU hour','mbp':'MAP','heart_rate':'Heart rate','peep':'PEEP','rr_vent':'Ventilator RR','rr_vital':'Monitor RR','fio2_trend':'Prior FiO₂ change','sat_trend':'Prior SpO₂ change'}
for ax,db,letter in zip(axes.flat,dbs,'ABC'):
    b=pd.read_csv(R/f'{db}_balance.csv');b=b[b.variable.isin(friendly)].copy();b['rank']=b.smd_after.abs();b=b.sort_values('rank')
    y=np.arange(len(b));ax.scatter(b.smd_before.abs(),y,color='#ACB6BE',s=15,label='Before');ax.scatter(b.smd_after.abs(),y,color=colors[db],s=18,label='After IPTW')
    ax.set_yticks(y,[friendly[v] for v in b.variable],fontsize=7);ax.axvline(.1,color='#7D858C',ls='--',lw=.7);ax.set_xlabel('Absolute SMD');panel(ax,letter,names[db])
    if db=='mimic':ax.legend(loc='lower right',fontsize=6.5)
ax=axes[1,1];panel(ax,'D','Support and treated ESS')
for i,db in enumerate(dbs):
    z=json.loads((R/f'{db}_diagnostics.json').read_text())
    ax.barh(i-.14,100*z['overlap_retained']/z['n'],height=.25,color=colors[db])
    ax.barh(i+.14,100*z['ess_1']/z['treated_n'],height=.25,color='#ACB6BE')
ax.set_yticks(range(3),[names[x] for x in dbs]);ax.invert_yaxis();ax.set_xlim(0,100);ax.set_xlabel('Percentage')
fig.text(.715,.028,'Colour: retained with PS 0.10–0.90\nGrey: treated ESS / treated count',ha='center',fontsize=7)
save(fig,'Figure_4')
# Supplementary positivity distributions and balance under an alternative weighting target.
fig,axes=plt.subplots(2,2,figsize=(7.2,6.4));fig.subplots_adjust(left=.13,right=.985,bottom=.16,top=.91,wspace=.50,hspace=.80)
for ax,db,letter in zip(axes.flat,dbs,'ABC'):
    hist=pd.read_csv(R/f'{db}_ps_histogram.csv')
    for a,col,label in [(1,'#355B78','Reduction ≥10 pp'),(0,'#ACB6BE','Stable ±5 pp')]:
        z=hist[hist.A.eq(a)];heights=z.n/z.n.sum()*100
        ax.stairs(heights.to_numpy(),np.r_[z.left.to_numpy(),z.right.iloc[-1]],color=col,lw=1.1,label=label)
    ax.axvspan(.10,.90,color='#EDF1F3',zorder=-1)
    ax.set_xlim(0,1);ax.set_xlabel('Estimated propensity score');ax.set_ylabel('Percentage within treatment arm');panel(ax,letter,names[db])
fig.legend(*axes[0,0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.55,.025),ncol=2,fontsize=7)
ax=axes[1,1];panel(ax,'D','Numeric balance with overlap weights')
for i,db in enumerate(dbs):
    b=pd.read_csv(R/f'{db}_overlap_weight_balance.csv')
    ax.scatter(b.smd_overlap.abs(),np.full(len(b),i),color=colors[db],s=17)
ax.set_yticks(range(3),[names[x] for x in dbs]);ax.set_ylim(2.6,-.6);ax.axvline(.10,color='#7D858C',ls='--',lw=.7);ax.set_xlabel('Absolute SMD across numeric variables')
save(fig,'Figure_S1')
print('Four figures exported from audited results.')
