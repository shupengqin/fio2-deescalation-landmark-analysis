import os
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

O = Path(__file__).resolve().parent
SOURCE = Path(os.environ['AGGREGATE_SOURCE_DIR'])
OUT = O / 'figures'
OUT.mkdir(parents=True, exist_ok=True)
# Quantitative grid: every saved model term contributes to panels A-C.
# Panel D shows the information retained by the reduction-arm weights.
plt.rcParams.update({'font.family':'Times New Roman', 'font.size':10,
    'axes.titlesize':12, 'axes.labelsize':10, 'pdf.fonttype':42,
    'svg.fonttype':'none', 'axes.spines.top':False, 'axes.spines.right':False})
names = {'mimic':'MIMIC-IV', 'eicu':'eICU-CRD', 'sicdb':'SICdb'}
pretty = {'sat_trend missing':'Saturation trend missing',
          'unit: Cardiac Vascular Intensive Care Unit (CVICU)':'CVICU indicator',
          'fio2':'FiO₂', 'mbp missing':'MAP missing', 'sat':'Saturation',
          'site: 443':'Hospital indicator', 'peep missing':'PEEP missing',
          'fio2_trend':'FiO₂ trend', 'calendar_group: 2018':'2018 indicator'}
groups = ['Numeric', 'Categorical', 'Missingness']
colors = ['#376B86', '#8A7092', '#C07849']
fig, axes = plt.subplots(2, 2, figsize=(7.5, 6.2))
fig.subplots_adjust(left=.13, right=.98, bottom=.10, top=.92, hspace=.57, wspace=.43)
source_rows=[]; diagnostics=[]
for j, db in enumerate(names):
    ax = axes.flat[j]
    frame = pd.read_csv(SOURCE / f'{db}_balance.csv')
    frame['group'] = frame.variable.map(lambda s:'Missingness' if s.endswith(' missing') else 'Categorical' if ': ' in s else 'Numeric')
    frame['absolute_smd'] = frame.smd_after.abs()
    assert frame.absolute_smd.notna().all()
    d = json.loads((SOURCE / f'{db}_diagnostics.json').read_text())
    assert abs(frame.absolute_smd.max()-d['max_abs_smd']) < 1e-10
    diagnostics.append(d)
    ax.axvspan(.10, .33, color='#F5EEE7', zorder=0)
    ax.axvline(.10, ls='--', color='#6B6B6B', lw=.9)
    labels=[]
    for k, (group, color) in enumerate(zip(groups, colors)):
        sub=frame[frame.group.eq(group)].sort_values('variable')
        y=2-k
        jitter=np.linspace(-.10, .10, len(sub))
        ax.scatter(sub.absolute_smd, y+jitter, s=17, color=color, alpha=.72, edgecolors='none', zorder=3)
        peak=sub.loc[sub.absolute_smd.idxmax()]
        ax.scatter([peak.absolute_smd], [y+jitter[list(sub.index).index(peak.name)]], marker='D', s=37,
                   facecolor=color, edgecolor='#222222', linewidth=.7, zorder=4)
        peak_label=pretty.get(peak.variable, peak.variable)
        ax.text(.325, y+.32, f'{peak_label}: {peak.absolute_smd:.3f}', ha='right', va='center', fontsize=9)
        labels.append(f'{group}\n({len(sub)} terms)')
    ax.set_yticks([2,1,0], labels)
    ax.set_ylim(-.4,2.65); ax.set_xlim(-.005,.335)
    ax.set_xticks([0,.1,.2,.3]);ax.set_xlabel('Absolute SMD after IPTW')
    ax.tick_params(axis='y',length=0,labelsize=9)
    ax.set_title(names[db], loc='left', pad=10)
    ax.text(-.23,1.085,'ABC'[j], transform=ax.transAxes, fontsize=13, fontweight='bold')
    ax.spines['left'].set_visible(False)
    for row in frame.itertuples():
        source_rows.append(dict(database=db,variable=row.variable,group=row.group,smd_before=row.smd_before,smd_after=row.smd_after))
ax=axes.flat[3]
pct=[100*d['ess_1']/d['treated_n'] for d in diagnostics]
ax.barh([2,1,0],pct,height=.40,color='#376B86')
for y, p, d in zip([2,1,0],pct,diagnostics):
    ax.text(p+.8,y,f'{p:.1f}%',va='center',fontsize=10)
    ax.text(.4,y+.34,f"ESS {d['ess_1']:.1f} / {d['treated_n']:,} reduction records",fontsize=9,va='center')
ax.set_yticks([2,1,0],list(names.values()));ax.tick_params(axis='y',length=0,labelsize=9)
ax.set_xlim(0,42);ax.set_ylim(-.4,2.65);ax.set_xticks([0,10,20,30,40])
ax.set_xlabel('Reduction ESS / reduction records (%)')
ax.set_title('Weighted treatment information',loc='left',pad=10)
ax.text(-.23,1.085,'D',transform=ax.transAxes,fontsize=13,fontweight='bold')
ax.spines['left'].set_visible(False)
fig.canvas.draw()
# Boundary check for every visible text object; no axis labels outside canvas.
bounds=fig.bbox
outside=[]
for t in fig.findobj(matplotlib.text.Text):
    if t.get_visible() and t.get_text():
        box=t.get_window_extent(fig.canvas.get_renderer())
        if box.x0<-.5 or box.y0<-.5 or box.x1>bounds.x1+.5 or box.y1>bounds.y1+.5:
            outside.append(t.get_text())
assert not outside, outside
fig.savefig(OUT/'Figure_2.pdf')
fig.savefig(OUT/'Figure_2.svg')
fig.savefig(OUT/'Figure_2.png',dpi=600)
fig.savefig(OUT/'Figure_2.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
fig.savefig(OUT/'Figure_2_preview.png',dpi=300)
pd.DataFrame(source_rows).to_csv(O/'figure2_all_terms.csv',index=False,encoding='utf-8-sig')
(O/'figure2_qa.json').write_text(json.dumps({'terms':len(source_rows),'all_saved_terms_included':True,
    'maxima_match_primary_diagnostics':True,'out_of_canvas_text':outside},indent=2))
plt.close(fig)
