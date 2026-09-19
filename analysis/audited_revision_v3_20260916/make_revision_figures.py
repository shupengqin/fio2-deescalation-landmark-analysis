from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parent;R=O.parent/'audited_reanalysis_v2';OUT=O/'figures';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':['Arial','DejaVu Sans'],'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'legend.frameon':False})
dbs=['mimic','eicu','sicdb'];names=dict(mimic='MIMIC-IV',eicu='eICU-CRD',sicdb='SICdb');colors=dict(mimic='#355B78',eicu='#BB763D',sicdb='#568B83')
def save(fig,name):
 fig.canvas.draw()
 fig.savefig(OUT/(name+'.pdf'),facecolor='white');fig.savefig(OUT/(name+'.svg'),facecolor='white');fig.savefig(OUT/(name+'.png'),dpi=300,facecolor='white');fig.savefig(OUT/(name+'.tiff'),dpi=600,facecolor='white',pil_kwargs={'compression':'tiff_lzw'});plt.close(fig)
def panel(ax,letter,title):
 ax.set_title(title,loc='left',pad=12);ax.text(-.13,1.08,letter,transform=ax.transAxes,fontweight='bold',fontsize=11)
res=pd.read_csv(R/'all_results.csv');extra=pd.read_csv(R/'all_extra_results.csv');primary=pd.concat([pd.read_csv(R/f'{db}_primary_bootstrap1000.csv') for db in dbs])
fig,axes=plt.subplots(2,2,figsize=(7.2,6.15));fig.subplots_adjust(left=.13,right=.985,bottom=.10,top=.91,wspace=.60,hspace=.80)
for ax,analysis,letter,title in zip(axes.flat,['hospital_primary','recorded28_secondary','hospital_linear','hospital_overlap_weighted'],'ABCD',['In-hospital mortality: full cohort','Recorded 28-day death: secondary','Linear nuisance models','Hospital death: overlap population']):
 for i,db in enumerate(dbs):
  if analysis=='hospital_primary':z=primary[primary.database.eq(db)].iloc[0];m,lo,hi=[z[c]*100 for c in ['rd','ci_low','ci_high']]
  elif analysis=='hospital_overlap_weighted':z=extra[extra.database.eq(db)&extra.analysis.eq(analysis)].iloc[0];m,lo,hi=[z[c]*100 for c in ['rd','ci_low','ci_high']]
  else:z=res[res.database.eq(db)&res.analysis.eq(analysis)].iloc[0];m,lo,hi=[z[c]*100 for c in ['aipw_rd','aipw_rd_low','aipw_rd_high']]
  ax.plot([lo,hi],[i,i],color=colors[db],lw=1.5);ax.scatter(m,i,color=colors[db],s=25,zorder=3)
 ax.axvline(0,color='#7D858C',lw=.8,ls='--');ax.set_yticks(range(3),[names[x] for x in dbs]);ax.set_ylim(2.7,-.7);ax.set_xlabel('Risk difference (percentage points)');ax.grid(axis='x',color='#E9ECEF',lw=.5);panel(ax,letter,title)
save(fig,'Figure_2')
# All panels are descriptive and use complete aggregate counts or treatment-specific quantiles.
event=pd.read_csv(O/'classification_events.csv');delay=pd.read_csv(O/'selection_delay.csv');phase=pd.read_csv(O/'care_phase.csv');time=pd.read_csv(O/'timestamp_summary.csv')
fig,axs=plt.subplots(2,2,figsize=(7.2,6.55));fig.subplots_adjust(left=.13,right=.985,bottom=.14,top=.90,wspace=.60,hspace=.90)
ax=axs[0,0];panel(ax,'A','Next-bin FiO₂ at first eligibility')
for i,db in enumerate(dbs):
 z=event[event.database.eq(db)&event.physiological_first_stays.notna()].iloc[0];v=100*z.first_next_available/z.physiological_first_stays
 ax.barh(i,v,color=colors[db],height=.45);ax.text(v+1,i,f'{v:.1f}%',va='center',fontsize=7)
ax.set_yticks(range(3),[names[d] for d in dbs]);ax.invert_yaxis();ax.set_xlim(0,115);ax.set_xticks([0,50,100]);ax.set_xlabel('Stays with next-bin recording (%)')
ax=axs[0,1];panel(ax,'B','Selected later than first eligibility')
for i,db in enumerate(dbs):
 z=delay[delay.database.eq(db)].iloc[0];v=100*z.delayed_n/z.selected_n
 ax.barh(i,v,color=colors[db],height=.45);ax.text(v+1,i,f'{v:.1f}%',va='center',fontsize=7)
ax.set_yticks(range(3),[names[d] for d in dbs]);ax.invert_yaxis();ax.set_xlim(0,100);ax.set_xlabel('Final cohort with delayed selection (%)')
ax=axs[1,0];panel(ax,'C','Baseline FiO₂ equal to 100%')
for a,col,label in [(1,'#355B78','Reduction ≥10 pp'),(0,'#ACB6BE','Stable ±5 pp')]:
 z=phase[phase.A.eq(a)].set_index('database');v=[100*z.loc[d,'fio2_100_n']/z.loc[d,'n'] for d in dbs]
 ax.barh(np.arange(3)+(a-.5)*.28,v,color=col,height=.25,label=label)
ax.set_yticks(range(3),[names[d] for d in dbs]);ax.invert_yaxis();ax.set_xlim(0,100);ax.set_xlabel('Percentage within treatment arm')
fig.legend(*ax.get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.54,.022),ncol=2,fontsize=7)
ax=axs[1,1];panel(ax,'D','FiO₂–SpO₂ timestamp separation')
for a,col in [(1,'#355B78'),(0,'#ACB6BE')]:
 for i,db in enumerate(['mimic','eicu']):
  z=time[time.database.eq(db)&time.A.eq(a)&time.variable.eq('baseline_alignment_min')].iloc[0];y=i+(a-.5)*.28
  ax.plot([z.q25,z.q75],[y,y],color=col,lw=2);ax.scatter(z['median'],y,s=24,color=col,zorder=3)
ax.set_yticks(range(2),[names[d] for d in ['mimic','eicu']]);ax.set_ylim(1.5,-.6);ax.set_xlim(-2,60);ax.set_xlabel('Minutes, median and IQR')
save(fig,'Figure_S2')
print('Updated Figure 2 and new four-panel Figure S2 from aggregate source tables.')
