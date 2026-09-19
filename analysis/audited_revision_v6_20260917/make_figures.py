from pathlib import Path
import json
import numpy as np,pandas as pd,matplotlib as mpl,matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
O=Path(__file__).resolve().parent;F=O/'figures';F.mkdir(exist_ok=True)
# Contract: quantitative A-D grids, complete source diagnostics, no causal or accuracy claim.
# Python backend; 7.2 inch double-column figures; editable PDF/SVG and 600dpi TIFF.
mpl.rcParams.update({'font.family':'Arial','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
def export(fig,name):
 fig.savefig(F/(name+'.pdf'),bbox_inches='tight')
 fig.savefig(F/(name+'.svg'),bbox_inches='tight')
 fig.savefig(F/(name+'.png'),dpi=600,bbox_inches='tight')
 fig.savefig(F/(name+'.tiff'),dpi=600,bbox_inches='tight')
 plt.close(fig)
q=pd.read_csv(O/'eicu_site_penalty_splits.csv');fig,axs=plt.subplots(2,2,figsize=(7.2,6.2),layout='constrained');colors=['#b67940','#477a95']
for mode,color in zip(['original','equal_site_penalty'],colors):
 z=q[q['mode'].eq(mode)];label='Original' if mode=='original' else 'Equal hospital penalty'
 axs[0,0].plot(range(1,21),z.rd*100,'o-',ms=3,lw=.8,color=color,label=label)
 axs[0,1].plot(range(1,21),z.max_abs_smd,'o-',ms=3,lw=.8,color=color)
 axs[1,0].plot(range(1,21),z.treated_ess,'o-',ms=3,lw=.8,color=color)
 axs[1,1].plot(range(1,21),z.top10_weight_share*100,'o-',ms=3,lw=.8,color=color)
titles=['Split-dependent risk differences','Residual imbalance','Reduction effective sample size','Top 10 records: weight share'];yl=['Risk difference (pp)','Maximum absolute SMD','Effective sample size','Percentage of total weight']
for i,ax in enumerate(axs.flat):
 ax.set_title(titles[i],loc='left',fontsize=9,pad=16);ax.text(-.16,1.12,'ABCD'[i],transform=ax.transAxes,fontweight='bold',fontsize=12);ax.set_xlabel('Patient split');ax.set_ylabel(yl[i]);ax.set_xticks([1,5,10,15,20])
axs[0,0].axhline(0,color='.5',ls='--',lw=.7);axs[0,0].legend(fontsize=7,loc='best');axs[0,1].axhline(.1,color='.5',ls='--',lw=.7);export(fig,'Figure_S7')
flow=pd.read_csv(O/'sic_independent_eligibility.csv');over=json.loads((O/'sic_independent_overlap.json').read_text());m=pd.read_csv(O/'sic_trajectory_summary.csv');fig,axs=plt.subplots(2,2,figsize=(7.2,6.4),layout='constrained')
stages=['baseline eligible','classifiable and survived','first classifiable','binary and hospital status'];names=['Baseline eligible','Classifiable + survived','First classifiable','Final binary cohort']
for j,(mode,c) in enumerate(zip(['hour','minute'],colors)):
 z=flow[flow.definition.eq(mode)].set_index('stage').loc[stages];axs[0,0].barh(np.arange(4)+(j-.5)*.32,z.cases,height=.3,color=c,label='Hourly' if mode=='hour' else 'Minute')
axs[0,0].set_yticks(range(4),names,fontsize=7);axs[0,0].invert_yaxis();axs[0,0].set_xlabel('Cases');axs[0,0].legend(fontsize=7)
vv=[over['hour_only'],over['common_cases'],over['minute_only']];axs[0,1].bar(['Hourly\nonly','Shared','Minute\nonly'],vv,color=[colors[0],'#8b969d',colors[1]])
for i,v in enumerate(vv):axs[0,1].text(i,v+35,str(v),ha='center',fontsize=8)
axs[0,1].set_ylim(0,2600);axs[0,1].set_ylabel('Cases')
cats=['sustained_to_end','reduction_with_rebound','run_without_confirmed_persistence','no_qualifying_run','insufficient'];labs=['Persists to bin end','Run then rebound','Run, end not confirmed','No qualifying run','Insufficient coverage']
for ax,pop in zip(axs[1],['hour_first','minute_first']):
 z=m[m.population.eq(pop)&m.treatment.eq('reduction')].set_index('metric').loc[['morph5_'+x for x in cats]];v=z.estimate.to_numpy()*100;ax.errorbar(v,range(5),xerr=np.array([v-z.ci_low.to_numpy()*100,z.ci_high.to_numpy()*100-v]),fmt='o',color=colors[0] if pop=='hour_first' else colors[1],capsize=2,ms=4);ax.set_yticks(range(5),labs,fontsize=7);ax.invert_yaxis();ax.set_xlim(0,100);ax.set_xlabel('Percent of reductions (95% CI)')
for i,ax in enumerate(axs.flat):ax.set_title(['Independent eligibility','Final cohort overlap','Hourly-defined reductions','Minute-defined reductions'][i],loc='left',fontsize=9,pad=16);ax.text(-.16,1.12,'ABCD'[i],transform=ax.transAxes,fontweight='bold',fontsize=12)
export(fig,'Figure_S8')
