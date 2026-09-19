from pathlib import Path
import shutil
import numpy as np,pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parent;W=O.parent/'audited_revision_v4_20260916';F=O/'figures';F.mkdir(exist_ok=True)
for p in (W/'figures').iterdir():
 if p.is_file():shutil.copy2(p,F/p.name)
# Retain measured-hypoxaemia figure, now supplementary, without changing its data.
for ext in ['png','pdf','svg','tiff']:
 p=F/f'Figure_4.{ext}'
 if p.exists():shutil.copy2(p,F/f'Figure_S6.{ext}')
mpl.rcParams.update({'font.family':'Arial','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'legend.frameon':False})
blue='#396A8B';orange='#B67C40';grey='#565B61'
u=pd.read_csv(O/'minute_uncertainty_results.csv');m=pd.read_csv(O/'minute_model_results.csv')
fig,axes=plt.subplots(2,2,figsize=(7.2,6.35));fig.subplots_adjust(left=.22,right=.98,bottom=.11,top=.93,wspace=.87,hspace=.72)
def panel(ax,letter,title):
 ax.set_title(title,loc='left',pad=12);ax.text(-.38,1.10,letter,transform=ax.transAxes,fontweight='bold',fontsize=12)
 ax.grid(axis='x',color='#E5E7E9',lw=.6);ax.set_axisbelow(True)
def point(ax,y,r,color):
 v,lo,hi=np.array([r.estimate,r.ci_low,r.ci_high])*100 if 'estimate' in r else np.array([r.rd,r.ci_low,r.ci_high])*100
 ax.errorbar(v,y,xerr=[[v-lo],[hi-v]],fmt='o',color=color,ms=5,capsize=3,lw=1.3)
ax=axes[0,0];panel(ax,'A','Treatment-class retention')
for y,a,c in [(1,1,orange),(0,0,blue)]:point(ax,y,u[u.measure.eq('same_treatment_class')&u.coverage.eq('paired')&u.A.eq(a)].iloc[0],c)
ax.set_yticks([1,0],['Hourly reduction\nn = 494','Hourly stable\nn = 1,949']);ax.set_xlim(0,100);ax.set_ylim(-.55,1.55);ax.set_xlabel('Same last-value class (%)')
ax=axes[0,1];panel(ax,'B','In-hospital mortality')
q=m[m.analysis.eq('common_hospital')]
for y,definition,c in [(2,'A',blue),(1,'A_minute',orange),(0,'paired_difference',grey)]:point(ax,y,q[q.definition.eq(definition)].iloc[0],c)
ax.set_yticks([2,1,0],['Hourly mean','Last minute','Paired difference']);ax.set_xlim(-12,10);ax.set_ylim(-.7,2.7);ax.axvline(0,color='#999999',lw=.8,ls='--');ax.set_xlabel('Risk difference (pp)');ax.text(.02,.97,'Same 2,265 cases',transform=ax.transAxes,va='top',fontsize=8)
ax=axes[1,0];panel(ax,'C','Detection beyond hourly means')
for y,a,measure,c in [(3,1,'minute_low88_minus_hourly',orange),(2,0,'minute_low88_minus_hourly',blue),(1,1,'consecutive5_low88_minus_hourly',orange),(0,0,'consecutive5_low88_minus_hourly',blue)]:
 point(ax,y,u[u.measure.eq(measure)&u.coverage.eq('any_follow')&u.A.eq(a)].iloc[0],c)
ax.set_yticks([3,2,1,0],['Any minute: reduce','Any minute: stable','5-min run: reduce','5-min run: stable']);ax.set_ylim(-.6,3.6);ax.set_xlim(-2,40);ax.axvline(0,color='#999999',lw=.8,ls='--');ax.set_xlabel('Paired prevalence difference (pp)')
ax=axes[1,1];panel(ax,'D','Five-minute runs below 88%')
q=m[m.analysis.eq('common_dense_follow_low88')]
for y,definition,c in [(2,'A',blue),(1,'A_minute',orange),(0,'paired_difference',grey)]:point(ax,y,q[q.definition.eq(definition)].iloc[0],c)
ax.set_yticks([2,1,0],['Hourly mean','Last minute','Paired difference']);ax.set_ylim(-.7,2.7);ax.set_xlim(-5,5);ax.axvline(0,color='#999999',lw=.8,ls='--');ax.set_xlabel('Risk difference (pp)');ax.text(.02,.97,'Same 2,061 cases',transform=ax.transAxes,va='top',fontsize=8)
fig.savefig(F/'Figure_4.pdf');fig.savefig(F/'Figure_4.svg')
fig.savefig(F/'Figure_4.png',dpi=300);fig.savefig(F/'Figure_4.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'});plt.close(fig)
print('Figure 4 exported; original Figure 4 retained as S6.')
