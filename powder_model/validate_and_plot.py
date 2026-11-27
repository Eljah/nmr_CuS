from pathlib import Path
import json
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import nnls
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from powder_fit import *

def exact_result(record,nu=80,nphi=160):
    dirs=directions(nu,nphi);fit=Fit(np.empty(0),record['series_ms'],record['mask_half_width_MHz'],record['common_eta'],record['free_centres'])
    pars=record['params'];eta,cen,sigma,rho=fit.unpack(pars)
    a=np.empty((len(fit.x),4))
    for site in range(4):
        a[:,site]=0
        for iso,nqr in enumerate(NQR_SETS[record['dataset']][site]):
            nu_l=cen[site]*(GAMMA_RATIO if iso else 1)
            raw=stick_profile(eta[site],nqr,nu_l,dirs)
            prof=gaussian_filter1d(raw,sigma/STEP,mode='constant')
            a[:,site]+=(rho if iso else 1)*np.interp(fit.x-nu_l,REL,prof,left=0,right=0)
    ma=a[fit.mask];pa=ma-fit.Q@(fit.Q.T@ma);amps=nnls(pa,fit.py)[0]
    bg=np.linalg.lstsq(fit.base[fit.mask],fit.yy[fit.mask]-ma@amps,rcond=None)[0]
    pred=(a@amps+fit.base@bg)*fit.yscale
    sse=float(np.sum(((pred-fit.y)[fit.mask]/fit.yscale)**2))
    return fit,pred,{'grid':f'{nu}x{nphi}','normalized_SSE':sse,'amplitudes':(amps*fit.yscale).tolist(),'baseline':(bg*fit.yscale).tolist()}


def main():
    records=[]
    fig,axs=plt.subplots(3,1,figsize=(8,9),sharex=True,gridspec_kw={'height_ratios':[2,2,1]})
    for tag,color,label in [('500ms_main','#c56e2b','Common eta, common shift'),('500ms_four_eta','#007a92','Four eta, common shift'),('500ms_four_eta_four_centres','#8a598c','Four eta, four shifts')]:
        rec=json.loads((ROOT/f'fit_thesis2009_{tag}.json').read_text())
        for n in [40,80,160]:
            fit,pred,validation=exact_result(rec,n,2*n)
            records.append({'tag':tag,**validation})
            if n==160:
                np.savetxt(ROOT/f'exact_converged_{tag}.csv',np.column_stack([fit.x,fit.y,pred,fit.mask]),delimiter=',',header='frequency_MHz,observed,predicted,used_in_fit',comments='')
                axs[0].plot(fit.x,pred,color=color,label=label,lw=1.2)
                axs[1].plot(fit.x,pred,color=color,lw=1.2)
        print(tag,records[-3:],flush=True)
    axs[0].plot(fit.x,np.where(fit.mask,fit.y,np.nan),'k.',ms=3,label='500 ms: coil regions omitted')
    axs[1].plot(fit.x,np.where(fit.mask,fit.y,np.nan),'k.',ms=3)
    axs[0].set_ylim(0,11000);axs[1].set_ylim(0,3200)
    axs[0].legend(fontsize=8)
    for ax in axs[:2]:
        ax.set_ylabel('Stored echo integral (a.u.)');ax.grid(alpha=.2)
        for f in [75.941,81.346]:ax.axvspan(f-.12,f+.12,color='0.92',zorder=-1)
    profile=np.genfromtxt(ROOT/'eta_profile_thesis2009.csv',delimiter=',',names=True)
    axs[2].plot(profile['eta'],profile['SSE']/profile['SSE'].min(),color='#c56e2b')
    # Third panel needs its own eta axis.
    axs[2].remove();ax=fig.add_subplot(3,1,3)
    ax.plot(profile['eta'],profile['SSE']/profile['SSE'].min(),color='#c56e2b')
    ax.set_xlabel('Common eta (restricted model)');ax.set_ylabel('SSE / minimum');ax.grid(alpha=.2)
    axs[1].set_xlabel('Frequency (MHz)')
    axs[0].set_xlim(72,84.8)
    fig.tight_layout()
    fig.savefig(ROOT/'powder_fit_conditional.png',dpi=220,bbox_inches='tight')
    fig.savefig(ROOT/'powder_fit_conditional.pdf',bbox_inches='tight')
    (ROOT/'numerical_validation.json').write_text(json.dumps(records,indent=2))

if __name__ == '__main__':
    main()
