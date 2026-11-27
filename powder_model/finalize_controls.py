import json,csv
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import nnls
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from powder_fit import *

controls=json.loads((ROOT/'central_control_results.json').read_text())
dirs=directions(160,320)
fig,axs=plt.subplots(2,1,figsize=(9,6.5),sharex=True)
validation=[]
for rec,color in zip(controls,['#087b94','#cc703a']):
    p=np.array(rec['parameters']);total=0;vr={'source_seed':rec['source_seed'],'grid':'160x320 octant'}
    for ax,series,rho in zip(axs,['500','20'],p[9:]):
        fit=Fit(np.empty(0),series,.12,False,True);aa=np.zeros((len(fit.x),4))
        for site in range(4):
            for iso in [0,1]:
                nu_l=p[4+site]*(GAMMA_RATIO if iso else 1)
                raw=stick_profile(p[site],NQR_SETS['thesis2009'][site,iso],nu_l,dirs,central_only=True)
                smooth=gaussian_filter1d(raw,p[8]/STEP,mode='constant')
                aa[:,site]+=(rho if iso else 1)*np.interp(fit.x-nu_l,REL,smooth)
        t=(fit.x-fit.x.min())/(fit.x.max()-fit.x.min());base=np.column_stack([1-t,t])
        design=np.column_stack([aa,base]);coefs=nnls(design[fit.mask],fit.yy[fit.mask])[0]
        pred=design@coefs*fit.yscale;res=((pred-fit.y)/fit.yscale)[fit.mask];sse=float(res@res);total+=sse
        vr[series]={'SSE':sse,'baseline_endpoints':(coefs[4:]*fit.yscale).tolist(),'amplitudes':(coefs[:4]*fit.yscale).tolist()}
        np.savetxt(ROOT/f'central_exact_seed{rec["source_seed"]}_{series}ms.csv',np.column_stack([fit.x,fit.y,pred,fit.mask]),delimiter=',',header='frequency_MHz,observed,predicted,used_in_fit',comments='')
        ax.plot(fit.x,pred,color=color,lw=1.5,label=f'Conditional solution: eta1 = {p[0]:.3f}')
        if rec is controls[0]:ax.plot(fit.x,np.where(fit.mask,fit.y,np.nan),'k.',ms=3,label='Data outside coil mask')
        ax.set_title(f'{series} ms repetition time',loc='left',fontsize=11)
    vr['SSE_joint']=total;validation.append(vr);print(vr,flush=True)
for ax in axs:
    for f in [75.941,81.346]:ax.axvspan(f-.12,f+.12,color='0.92',zorder=-1)
    ax.grid(alpha=.2);ax.set_ylim(bottom=0);ax.set_ylabel('Echo integral (a.u.)')
axs[0].legend(fontsize=8,loc='upper right');axs[1].set_xlabel('Frequency (MHz)');axs[1].set_xlim(72,84.8)
fig.suptitle('Central-transition powder model with nonnegative empirical background',fontsize=11)
fig.tight_layout();fig.savefig(ROOT/'powder_model_comparison.png',dpi=240,bbox_inches='tight');fig.savefig(ROOT/'powder_model_comparison.pdf',bbox_inches='tight')
(ROOT/'central_exact_validation.json').write_text(json.dumps(validation,indent=2))

rows=[]
for r in json.loads((ROOT/'exact_joint_refinements.json').read_text()):
    rows.append({'model':'All transitions; free linear correction','branch':r['source_seed'],'calculation':'Direct exact-H refinement, 80x160 octant','SSE_joint':r['SSE_joint'],'eta':r['eta'],'nuQ':r['nuQ63_MHz'],'nuZ':r['nuZ63_MHz'],'sigma':r['sigma_MHz'],'baseline':'Negative compensation, not physical echo background'})
for r,val in zip(controls,validation):
    rows.append({'model':'Central transitions; nonnegative linear background','branch':r['source_seed'],'calculation':'Library fit; direct exact-H evaluation on 160x320 octant','SSE_joint':val['SSE_joint'],'eta':r['eta'],'nuQ':r['nuQ63_MHz'],'nuZ':r['nuZ63_MHz'],'sigma':r['sigma_MHz'],'baseline':'Nonnegative, approximately 600 a.u. / 130 a.u.'})
for r in json.loads((ROOT/'joint_literature_sensitivity.json').read_text()):
    s=r['long'];rows.append({'model':'2013 NQR set; all transitions; free linear correction','branch':r['warm_start_seed'],'calculation':'Warm-start library fit','SSE_joint':r['SSE_joint'],'eta':s['eta'],'nuQ':s['nuQ63_MHz'],'nuZ':s['effective_63_zeeman_MHz'],'sigma':s['sigma_MHz'],'baseline':'Negative compensation, not physical echo background'})
fields=['model','branch','calculation','SSE_joint']+[f'eta_s{i}' for i in range(1,5)]+[f'nuQ63_s{i}_MHz' for i in range(1,5)]+[f'nuZ63_s{i}_MHz' for i in range(1,5)]+['sigma_MHz','FWHM_MHz','baseline']
with (ROOT/'FINAL_MODEL_COMPARISON.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
    for row in rows:
        out={k:row[k] for k in ['model','branch','calculation','SSE_joint','baseline']}
        for k in range(4):out.update({f'eta_s{k+1}':row['eta'][k],f'nuQ63_s{k+1}_MHz':row['nuQ'][k],f'nuZ63_s{k+1}_MHz':row['nuZ'][k]})
        out.update(sigma_MHz=row['sigma'],FWHM_MHz=row['sigma']*np.sqrt(8*np.log(2)))
        w.writerow(out)
