from pathlib import Path
import shutil
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
from matplotlib.path import Path as MplPath
O=Path(__file__).resolve().parent;V=O.parent/'audited_revision_v3_20260916';R=O.parent/'audited_reanalysis_v2';P=V/'delivery'/'figures';F=O/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':['Arial','DejaVu Sans'],'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})
colors=['#355B78','#BB763D','#568B83'];dbs=['mimic','eicu','sicdb'];names=['MIMIC-IV','eICU-CRD','SICdb']
def save(fig,name):
 fig.canvas.draw()
 fig.savefig(F/(name+'.pdf'),facecolor='white')
 fig.savefig(F/(name+'.svg'),facecolor='white')
 fig.savefig(F/(name+'.png'),dpi=300,facecolor='white')
 fig.savefig(F/(name+'.tiff'),dpi=600,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
 plt.close(fig)
def panel(ax,l,t):ax.set_title(t,loc='left',pad=11);ax.text(-.13,1.08,l,transform=ax.transAxes,fontweight='bold',fontsize=11)
for old,new in [(4,2),(2,3),(3,4)]:
 for p in P.glob(f'Figure_{old}.*'):shutil.copy2(p,F/p.name.replace(f'Figure_{old}',f'Figure_{new}'))
for p in P.glob('Figure_S*'):shutil.copy2(p,F/p.name)
# Figure1: all stages are stays/cases. Source-evidence denominator is not true ventilation prevalence.
flow=pd.read_csv(V/'stay_flow.csv');fig,ax=plt.subplots(figsize=(7.2,8.1));ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
labels=['All source ICU stays/cases','Adult stays/cases','Documented invasive-support /\nairway evidence at any time','Physiologically eligible at\nleast once in analysis window','Eligible occasion with\nnext-bin FiO₂ recorded','Classifiable occasion with\nsurvival through landmark','First classifiable occasion\nselected within stay/case','Reduction ≥10 pp or\nstability within ±5 pp','Known hospital disposition']
xs=[.45,.665,.88];ys=np.linspace(.85,.12,9);w=.18;hh=.059
ax.text(.02,.98,'Source selection and construction of the analysis cohorts',fontsize=11,fontweight='bold',va='top')
for x,n,c in zip(xs,names,colors):ax.text(x,.925,n,ha='center',color=c,fontweight='bold',fontsize=9)
for i,(y,lab) in enumerate(zip(ys,labels)):
 ax.text(.015,y,lab,va='center',fontsize=8)
 for x,db,c in zip(xs,dbs,colors):
  counts=flow[flow.database.eq(db)].stays.to_numpy();n=int(counts[i]);ax.add_patch(FancyBboxPatch((x-w/2,y-hh/2),w,hh,boxstyle='round,pad=.005',lw=.8,edgecolor=c,facecolor='#F5F8FA'))
  ax.text(x,y,f'{n:,}',ha='center',va='center',fontsize=9)
  if i:ax.add_patch(FancyArrowPatch((x,ys[i-1]-hh/2-.006),(x,y+hh/2+.006),arrowstyle='->',mutation_scale=8,lw=.7,color='#919AA3'))
ax.text(.015,.046,'All counts are stays/cases. First-occasion selection does not remove additional stays.\nDocumentation rules differ by source; incomplete documentation is not absence of ventilation.\nHourly exclusion counts and distinct patients are provided in Supplementary Table S1.',fontsize=7,va='center')
fig.subplots_adjust(left=.02,right=.995,top=.99,bottom=.02);save(fig,'Figure_1')
v=pd.read_csv(O/'treatment_versions.csv');h=pd.read_csv(O/'hospital_support.csv');loo=pd.read_csv(O/'hospital_leave_one_out.csv');b=pd.read_csv(O/'hospital_both_arms_result.csv').iloc[0]
assert (h.n > 0).all(), 'Logarithmic hospital-size axis requires positive counts'
fig,axs=plt.subplots(2,2,figsize=(7.2,6.5));fig.subplots_adjust(left=.12,right=.98,bottom=.12,top=.91,wspace=.55,hspace=.80)
ax=axs[0,0];panel(ax,'A','Observed reduction magnitude')
for j,db in enumerate(dbs):
 z=v[v.database.eq(db)&v.A.eq(1)].iloc[0];ax.bar(np.arange(3)+(j-1)*.24,[100*z[c]/z.n for c in ['drop10_20','drop20_30','drop30plus']],width=.22,color=colors[j],label=names[j])
ax.set_xticks(range(3),['10–<20','20–<30','≥30']);ax.set_xlabel('Decrease (percentage points)');ax.set_ylabel('Within reduction group (%)');ax.legend(fontsize=7)
ax=axs[0,1];panel(ax,'B','eICU hospital treatment support');ax.scatter(h.n,h.reduction_pct,s=18,c=colors[1],alpha=.75);ax.set_xlabel('Records per hospital');ax.set_ylabel('Reduction records (%)');ax.set_xscale('log');ax.set_ylim(-4,104)
ax=axs[1,0];panel(ax,'C','Omit one eICU hospital at a time');z=loo.sort_values('excluded_n',ascending=False);ax.scatter(np.arange(1,len(z)+1),z.rd*100,s=12,c=colors[1]);ax.axhline(-4.147077266248236,color='#68737D',ls='--',lw=.8);ax.set_xlabel('Omitted hospital, ordered by size');ax.set_ylabel('Risk difference (pp)')
ax=axs[1,1];panel(ax,'D','eICU hospital-level uncertainty');p=pd.read_csv(R/'eicu_primary_bootstrap1000.csv').iloc[0];hosp=pd.read_csv(R/'eicu_hospital_bootstrap.csv').iloc[0]
for i,z in enumerate([p,hosp,b]):ax.plot([z.ci_low*100,z.ci_high*100],[i,i],lw=1.5,color=colors[1]);ax.scatter(z.rd*100,i,s=22,color=colors[1])
ax.axvline(0,color='#68737D',ls='--',lw=.8);ax.set_yticks(range(3),['Patient clusters','Hospital clusters','Both-arm hospitals']);ax.set_ylim(2.6,-.6);ax.set_xlabel('Risk difference (pp)');save(fig,'Figure_S3')
m=pd.read_csv(O/'minute_validation_summary.csv').set_index('A');agr=pd.read_csv(O/'minute_exposure_agreement.csv');low=pd.read_csv(O/'minute_hypoxaemia_comparison.csv')
fig,axs=plt.subplots(2,2,figsize=(7.2,6.5));fig.subplots_adjust(left=.12,right=.98,bottom=.12,top=.91,wspace=.60,hspace=.86)
ax=axs[0,0];panel(ax,'A','Baseline final 30 minutes');bottom=np.zeros(2)
for name,col,c in [('Meets rule','pass','#568B83'),('Does not meet','fail','#B3BFC8'),('Insufficient coverage','missing','#E2E6EA')]:
 vals=np.array([m.loc[a,'baseline30_pass_n'] if col=='pass' else m.loc[a,'baseline30_atleast24_n']-m.loc[a,'baseline30_pass_n'] if col=='fail' else m.loc[a,'n']-m.loc[a,'baseline30_atleast24_n'] for a in [1,0]])/np.array([m.loc[a,'n'] for a in [1,0]])*100
 ax.bar([0,1],vals,bottom=bottom,width=.55,color=c,label=name);bottom+=vals
ax.set_xticks([0,1],['Reduction','Stable']);ax.set_ylabel('Within hourly-defined group (%)');ax.legend(fontsize=6.5,loc='lower left',bbox_to_anchor=(-.08,-.48))
ax=axs[0,1];panel(ax,'B','FiO₂ summary classification')
mat=np.array([[agr[(agr.hourly_class.eq(r))&(agr.minute_class.eq(c))].n.iloc[0] for c in ['reduction','stable','other']] for r in ['reduction','stable']]);pct=100*mat/mat.sum(axis=1)[:,None]
ax.imshow(pct,cmap='Blues',vmin=0,vmax=100,aspect='auto')
for (i,j),v0 in np.ndenumerate(mat):ax.text(j,i,f'{v0}\n({pct[i,j]:.1f}%)',ha='center',va='center',color='white' if pct[i,j]>65 else '#243745',fontsize=8)
ax.set_xticks(range(3),['Reduction','Stable','Other']);ax.set_yticks([0,1],['Reduction','Stable']);ax.set_xlabel('Last-minute-value classification');ax.set_ylabel('Hourly-mean classification')
for ax,k,letter in [(axs[1,0],88,'C'),(axs[1,1],90,'D')]:
 panel(ax,letter,f'Saturation below {k}%')
 for j,a in enumerate([1,0]):
  z=low[low.A.eq(a)&low.threshold.eq(k)&low.coverage.eq('at_least_one')].iloc[0]
  ax.bar(np.arange(3)+(j-.5)*.32,[100*z[c]/z.n for c in ['hourly_events','minute_events','consecutive5_events']],width=.30,color=['#355B78','#B0BCC6'][j],label=['Reduction','Stable'][j])
 ax.set_xticks(range(3),['Low hourly\nmean','Any low\nminute','≥5 consecutive\nlow minutes']);ax.set_ylabel('Observed records (%)');ax.legend(fontsize=7)
save(fig,'Figure_S4')
# Conceptual DAG: assumed paths, not estimated relations; inclusion is conditioned on.
fig,ax=plt.subplots(figsize=(7.2,6.0));ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
pos={'L':(.20,.86),'U':(.80,.86),'A':(.25,.56),'Y':(.75,.56),'R':(.20,.26),'V':(.80,.26),'S':(.50,.08)}
texts={'L':'Measured pre-decision\nstate and site practices','U':'Unmeasured trajectory\nand events before adjustment','A':'Actual oxygen\nadjustment','Y':'Subsequent\nhospital death','R':'Classifiable FiO₂\nrecord available (R)','V':'Alive and observed\nthrough landmark (V)','S':'Included\nR = 1 and V = 1'}
for key,(x,y) in pos.items():
 w=.32 if key!='S' else .27;hh=.115 if key!='S' else .10
 ax.add_patch(FancyBboxPatch((x-w/2,y-hh/2),w,hh,boxstyle='round,pad=.007',facecolor='#F4F7F9' if key!='S' else '#E3ECF0',edgecolor='#536E80',lw=1));ax.text(x,y,texts[key],ha='center',va='center',fontsize=8)
edges=[('L','A',0),('L','Y',.1),('U','A',-.1),('U','Y',0),('A','Y',0),('A','R',0),('A','V',-.1),('R','S',0),('V','S',0)]
for a,b,rad in edges:
 x1,y1=pos[a];x2,y2=pos[b];dx=x2-x1;dy=y2-y1
 if abs(dy)<.03:st=(x1+.175*np.sign(dx),y1);en=(x2-.175*np.sign(dx),y2)
 else:st=(x1,y1-.065 if dy<0 else y1+.065);en=(x2,y2+.065 if dy<0 else y2-.065)
 ax.add_patch(FancyArrowPatch(st,en,connectionstyle=f'arc3,rad={rad}',arrowstyle='-|>',mutation_scale=9,lw=.8,color='#78858F',zorder=0))
for pts in [[(.033,.86),(.012,.86),(.012,.26),(.033,.26)],[(.967,.86),(.988,.86),(.988,.26),(.967,.26)],[(.70,.795),(.50,.71),(.50,.37),(.28,.325)]]:
 ax.add_patch(FancyArrowPatch(path=MplPath(pts,[MplPath.MOVETO]+[MplPath.LINETO]*(len(pts)-1)),arrowstyle='-|>',mutation_scale=9,lw=.8,color='#78858F',zorder=0))
ax.text(.5,.98,'Proposed confounding and selection structure',ha='center',va='top',fontsize=11,fontweight='bold');fig.subplots_adjust(left=.03,right=.97,bottom=.03,top=.98);save(fig,'Figure_S5')
print('Figure1 rebuilt, Figure2–4 reordered, S3–S5 added.')
