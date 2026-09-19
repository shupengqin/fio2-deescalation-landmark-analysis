"""Optional local-only clinical review packet. Outputs contain controlled data."""
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
# Local blinded signal packet: no outcome, original group, algorithmic classification or source ID.
private=O/'private_clinical_review';private.mkdir(exist_ok=True);review=pd.read_pickle(O/'sic_blinded_review_key_private.pkl');keys=pd.read_pickle(O/'sic_signal_keys_private.pkl');idx=pd.Index(keys.key);arr=np.load(O/'sic_candidate_signals_private.npy',mmap_mode='r')
with PdfPages(private/'blinded_signal_review.pdf') as pdf:
 for r in review.itertuples():
  i0=idx.get_indexer([r.id*100000+r.hr])[0];i1=idx.get_indexer([r.id*100000+r.hr+1])[0];fig,ax=plt.subplots(2,1,figsize=(8,5),sharex=True,layout='constrained')
  for j,label in [(0,'FiO2 (%)'),(1,'SpO2 (%)')]:
   vals=np.r_[arr[i0,j],arr[i1,j] if i1>=0 else np.full(60,np.nan)];ax[j].plot(np.arange(-60,60),vals,lw=1,color='#477a95');ax[j].axvline(0,color='.5',ls='--');ax[j].set_ylabel(label);ax[j].set_ylim(15 if j==0 else 50,105)
  ax[0].set_title('Signal review '+r.review_id,loc='left');ax[1].set_xlabel('Minutes relative to bin boundary');pdf.savefig(fig);plt.close(fig)
print('Two supplemental figures and 60-record local blinded packet exported')
