from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from PIL import Image
O=Path(__file__).resolve().parent/'new_figures';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Times New Roman','font.size':12,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(10,8.2),facecolor='white')
ax=fig.add_axes([.03,.02,.94,.96]);ax.set(xlim=(0,10),ylim=(0,10));ax.axis('off')
blue='#315D76'; teal='#467C76';grey='#647580'
def box(x,y,w,h,t,color=blue,fill='#F0F5F8',fs=12):
 p=FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0.035,rounding_size=0.06',lw=1.05,edgecolor=color,facecolor=fill,zorder=3);ax.add_patch(p)
 ax.text(x,y,t,ha='center',va='center',fontsize=fs,zorder=4)
 return p
def arrow(a,b):
 ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=13,lw=1.2,color=grey,zorder=2))
ax.text(.05,9.75,'A',fontweight='bold',fontsize=17);ax.text(.48,9.75,'Confounding of the adjustment–outcome association',fontsize=14)
box(5,8.83,4.5,.72,'Measured pre-decision state\nand site practices')
box(1.9,7.50,2.75,.8,'Actual oxygen\nadjustment')
box(8.1,7.50,2.75,.8,'Subsequent\nhospital death')
box(5,6.18,4.5,.72,'Unmeasured clinical trajectory\nbefore adjustment',teal,'#F0F6F3')
for a,b in [((3.3,8.43),(1.9,7.94)),((6.7,8.43),(8.1,7.94)),((3.3,6.58),(1.9,7.06)),((6.7,6.58),(8.1,7.06)),((3.32,7.5),(6.68,7.5))]:arrow(a,b)
ax.plot([0,10],[5.48,5.48],color='#D9E1E6',lw=.8)
ax.text(.05,5.07,'B',fontweight='bold',fontsize=17);ax.text(.48,5.07,'Selection into the recorded comparison',fontsize=14)
for y,target in [(4.1,'Classifiable next-bin\nrecord available (R)'),(2.65,'Alive and observed\nthrough landmark (V)')]:
 box(1.53,y,2.72,.86,'Actual oxygen\nadjustment')
 box(5,y,2.95,.86,target,teal,'#EAF3F1')
 box(8.48,y,2.72,.86,'Clinical state and\nrecording practices')
 arrow((2.94,y),(3.47,y));arrow((7.06,y),(6.53,y))
box(5,.9,5.5,.72,'Analysis requires R = 1 and V = 1',teal,'#DFEEEA')
# Inclusion is a rule, not an additional biological causal link.
ax.text(5,1.63,'Both conditions required',ha='center',fontsize=11,color=grey)
ax.plot([3.48,6.52],[1.94,1.94],color=grey,lw=1)
for suffix in ['pdf','svg','png']:
 fig.savefig(O/f'Figure_S5.{suffix}',dpi=600,facecolor='white',bbox_inches=None)
Image.open(O/'Figure_S5.png').save(O/'Figure_S5.tiff',compression='tiff_lzw',dpi=(600,600))
fig.savefig(O/'Figure_S5_preview.png',dpi=150)
print('S5 exported: PDF SVG PNG TIFF; conceptual diagram, no data altered.')
