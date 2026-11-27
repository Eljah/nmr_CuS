"""Local exact-Hamiltonian refinement; no translated frequency library."""
import json,functools,time
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import nnls,least_squares
from powder_fit import *

dirs=directions(80,160)
source=json.loads((ROOT/'joint_results.json').read_text())
long=Fit(np.empty(0),'500',.12,False,True);short=Fit(np.empty(0),'20',.12,False,True)
LOW=np.array([0]*4+[74.8]*4+[.025,.05,.05],float)
HIGH=np.array([1]*4+[76.8]*4+[.25,.8,.8],float)
STEPS=np.array([.002]*4+[.0015]*4+[.001,.002,.002])
name='thesis2009'

@functools.lru_cache(maxsize=300)
def profiles(site,eta,centre):
    return np.array([stick_profile(eta,NQR_SETS[name][site,iso],centre*(GAMMA_RATIO if iso else 1),dirs) for iso in [0,1]])

def design(p,fit,rho):
    aa=np.empty((len(fit.x),4));sigma=p[8]
    for site in range(4):
        eta,centre=p[site],p[4+site]
        raw=profiles(site,float(eta),float(centre))
        smooth=gaussian_filter1d(raw,sigma/STEP,axis=1,mode='constant')
        aa[:,site]=np.interp(fit.x-centre,REL,smooth[0])+rho*np.interp(fit.x-centre*GAMMA_RATIO,REL,smooth[1])
    return aa

def series_result(p,fit,rho):
    a=design(p,fit,rho);ma=a[fit.mask];pa=ma-fit.Q@(fit.Q.T@ma)
    amps=nnls(pa,fit.py)[0];bg=np.linalg.lstsq(fit.base[fit.mask],fit.yy[fit.mask]-ma@amps,rcond=None)[0]
    pred=a@amps+fit.base@bg
    return (pred-fit.yy)[fit.mask],pred*fit.yscale,amps*fit.yscale,bg*fit.yscale

def residual(p):return np.r_[series_result(p,long,p[9])[0],series_result(p,short,p[10])[0]]

def jac(p):
    cols=[]
    for i,step in enumerate(STEPS):
        lo=p.copy();hi=p.copy();lo[i]=max(LOW[i],p[i]-step);hi[i]=min(HIGH[i],p[i]+step)
        cols.append((residual(hi)-residual(lo))/(hi[i]-lo[i]))
    return np.asarray(cols).T

out=[]
for src in [min(source,key=lambda r:r['SSE_joint']),[r for r in source if r['seed']==67][0]]:
    t=time.time();p=np.array(src['params']);before=float(residual(p)@residual(p))
    opt=least_squares(residual,p,jac=jac,bounds=(LOW,HIGH),max_nfev=85,ftol=2e-6,xtol=2e-6,gtol=2e-6,x_scale='jac')
    p=opt.x
    r={'source_seed':src['seed'],'exact_grid':'80x160 octant','initial_exact_SSE':before,'SSE_joint':float(opt.fun@opt.fun),'parameters':p.tolist(),'eta':p[:4].tolist(),'nuZ63_MHz':p[4:8].tolist(),'sigma_MHz':float(p[8]),'rho500':float(p[9]),'rho20':float(p[10]),'nuQ63_MHz':(NQR_SETS[name][:,0]/np.sqrt(1+p[:4]**2/3)).tolist(),'nfev':opt.nfev,'message':opt.message,'success':bool(opt.success),'elapsed_seconds':time.time()-t}
    for fit,rho,key in [(long,p[9],'long'),(short,p[10],'short')]:
        res,pred,amps,bg=series_result(p,fit,rho)
        r[key]={'normalized_SSE':float(res@res),'amplitudes':amps.tolist(),'baseline':bg.tolist()}
        np.savetxt(ROOT/f'exact_refined_seed{src["seed"]}_{fit.series}ms.csv',np.column_stack([fit.x,fit.y,pred,fit.mask]),delimiter=',',header='frequency_MHz,observed,predicted,used_in_fit',comments='')
    out.append(r);print(json.dumps(r),flush=True)
(ROOT/'exact_joint_refinements.json').write_text(json.dumps(out,indent=2))
