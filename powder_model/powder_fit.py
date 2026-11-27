"""Exploratory exact-Hamiltonian powder fit for the uploaded Cu NMR spectra.

Four *spectroscopic groups*, not a claim of four crystallographic Cu sites.
No physical K tensor is inferred from the fitted effective Zeeman frequency.
"""
from pathlib import Path
import json, csv, argparse, time
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import nnls, differential_evolution, minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
DATA=ROOT.parent/'data_analysis/export'
GAMMA_RATIO=1.071210656566704
NQR_SETS={
    'thesis2009':np.array([[11.97,11.07],[15.27,14.13],[16.18,14.97],[16.96,15.69]]),
    'publication2013':np.array([[11.90,11.02],[15.21,14.10],[16.26,15.05],[16.81,15.55]])}
ETAS=np.linspace(0,1,41)
SIGMAS=np.linspace(.025,.25,24)
REL=np.linspace(-18,18,3601)
STEP=REL[1]-REL[0]
REFERENCE=75.8
I=1.5
M=np.arange(I,-I-1,-1)
JP=np.zeros((4,4),complex)
for j in range(1,4):JP[j-1,j]=np.sqrt(I*(I+1)-M[j]*(M[j]+1))
JX=(JP+JP.T)/2
JY=(JP-JP.T)/(2j)
JZ=np.diag(M)
SPIN=np.array([JX,JY,JZ])
QUAD0=3*JZ@JZ-I*(I+1)*np.eye(4)
QUADE=JX@JX-JY@JY

def directions(nu=40,nphi=80):
    # The D2 and time-reversal symmetries allow one octant. Uniform cos(theta), phi.
    u=(np.arange(nu)+.5)/nu
    phi=(np.arange(nphi)+.5)/nphi*np.pi/2
    U,P=np.meshgrid(u,phi,indexing='ij')
    return np.array([np.sqrt(1-U.ravel()**2)*np.cos(P.ravel()),
                     np.sqrt(1-U.ravel()**2)*np.sin(P.ravel()),U.ravel()]).T

def stick_profile(eta,nqr,nu_l,dirs,central_only=False):
    nq=nqr/np.sqrt(1+eta**2/3)
    H=-nu_l*np.einsum('na,aij->nij',dirs,SPIN)+nq/6*(QUAD0+eta*QUADE)
    vals,vecs=np.linalg.eigh(H)
    # A fixed transverse RF field averaged over the azimuth about B0:
    # 1/2*(sum_a |<u|Ia|v>|^2 - |<u|n.I|v>|^2).
    mat=np.einsum('nki,akl,nlj->naij',vecs.conj(),SPIN,vecs,optimize=True)
    par=np.einsum('na,naij->nij',dirs,mat)
    weights=(np.sum(np.abs(mat)**2,axis=1)-np.abs(par)**2)/2
    result=np.zeros(len(REL))
    pairs=[(1,2)] if central_only else [(a,b) for a in range(4) for b in range(a+1,4)]
    for a,b in pairs:
        offsets=vals[:,b]-vals[:,a]-nu_l
        w=weights[:,a,b]/len(dirs)
        # Linear deposit reduces discretization error before convolution.
        z=(offsets-REL[0])/STEP;lo=np.floor(z).astype(int);f=z-lo
        good=(lo>=0)&(lo<len(REL)-1)
        np.add.at(result,lo[good],w[good]*(1-f[good]))
        np.add.at(result,lo[good]+1,w[good]*f[good])
    return result/STEP

def build_library(name,nu=40,nphi=80):
    path=ROOT/f'library_{name}_{nu}x{nphi}.npz'
    if path.exists():return np.load(path)['templates']
    print('building',path.name,flush=True)
    dirs=directions(nu,nphi)
    templates=np.empty((len(ETAS),len(SIGMAS),4,2,len(REL)),np.float32)
    for ie,eta in enumerate(ETAS):
        for site,pair in enumerate(NQR_SETS[name]):
            for iso,nqr in enumerate(pair):
                profile=stick_profile(eta,nqr,REFERENCE*(GAMMA_RATIO if iso else 1),dirs)
                for si,sigma in enumerate(SIGMAS):
                    templates[ie,si,site,iso]=gaussian_filter1d(profile,sigma/STEP,mode='constant')
        if ie%10==0:print('eta',eta,flush=True)
    np.savez_compressed(path,templates=templates)
    return templates

def brackets(grid,v):
    k=int(np.clip(np.searchsorted(grid,v)-1,0,len(grid)-2));a=(v-grid[k])/(grid[k+1]-grid[k])
    return k,a

class Fit:
    def __init__(self,lib,series='500',mask_width=.12,common_eta=True,free_centres=False):
        raw=np.genfromtxt(DATA/'opj_Spec05122016.csv',delimiter=',',names=True)
        x=raw['Freq'];y=raw['IntQ1' if series=='500' else 'IntQ']
        valid=np.isfinite(x)&np.isfinite(y)&(~np.isclose(x,81.4,atol=1e-7))
        self.x=x[valid];self.y=y[valid]
        self.mask=(abs(self.x-75.941)>mask_width)&(abs(self.x-81.346)>mask_width)
        self.yscale=max(self.y[self.mask]);self.yy=self.y/self.yscale
        self.lib=lib;self.series=series;self.mask_width=mask_width
        self.common_eta=common_eta;self.free_centres=free_centres
        self.base=np.column_stack([np.ones(len(self.x)),(self.x-78)/6])
        self.Q=np.linalg.qr(self.base[self.mask])[0]
        self.py=self.yy[self.mask]-self.Q@(self.Q.T@self.yy[self.mask])
    def component(self,eta,sigma,site):
        ie,a=brackets(ETAS,eta);si,b=brackets(SIGMAS,sigma)
        return ((1-a)*(1-b)*self.lib[ie,si,site]+a*(1-b)*self.lib[ie+1,si,site]
               +(1-a)*b*self.lib[ie,si+1,site]+a*b*self.lib[ie+1,si+1,site])
    def unpack(self,p):
        n_eta=1 if self.common_eta else 4
        eta=np.repeat(p[0],4) if self.common_eta else np.array(p[:4])
        nc=4 if self.free_centres else 1
        cen=np.array(p[n_eta:n_eta+nc]) if nc==4 else np.repeat(p[n_eta],4)
        sigma,rho=p[-2:]
        return eta,cen,sigma,rho
    def design(self,p,xx=None):
        eta,cen,sigma,rho=self.unpack(p);xx=self.x if xx is None else xx
        A=np.empty((len(xx),4))
        for j in range(4):
            tmp=self.component(eta[j],sigma,j)
            A[:,j]=np.interp(xx-cen[j],REL,tmp[0],left=0,right=0)+rho*np.interp(xx-cen[j]*GAMMA_RATIO,REL,tmp[1],left=0,right=0)
        return A
    def solve(self,p,ret=False):
        A=self.design(p);a=A[self.mask]
        pa=a-self.Q@(self.Q.T@a)
        amps=nnls(pa,self.py)[0]
        bg=np.linalg.lstsq(self.base[self.mask],self.yy[self.mask]-a@amps,rcond=None)[0]
        pred=A@amps+self.base@bg
        res=(pred-self.yy)[self.mask]
        if ret:return float(res@res),pred*self.yscale,amps*self.yscale,bg*self.yscale
        return float(res@res)
    def bounds(self):
        return [(0,1)]*(1 if self.common_eta else 4)+[(75.2,76.3)]*(4 if self.free_centres else 1)+[(SIGMAS[0],SIGMAS[-1]),(.15,.8)]
    def optimize(self,seed=42,popsize=8,maxiter=110):
        de=differential_evolution(self.solve,self.bounds(),seed=seed,popsize=popsize,maxiter=maxiter,tol=2e-5,polish=False,workers=1)
        fit=minimize(self.solve,de.x,method='Nelder-Mead',bounds=self.bounds(),options={'maxiter':3000,'xatol':1e-7,'fatol':1e-10})
        return fit
    def profile_eta(self,p,grid=np.linspace(0,1,41)):
        assert self.common_eta
        rows=[]
        for eta in grid:
            def fun(v):return self.solve(np.r_[eta,v])
            v=np.array(p[1:])
            opt=minimize(fun,v,method='Powell',bounds=self.bounds()[1:],options={'maxiter':80,'xtol':1e-5,'ftol':1e-7})
            rows.append((eta,opt.fun,*opt.x))
        return np.asarray(rows)

def summary(fit,opt,name,tag):
    sse,pred,amps,bg=fit.solve(opt.x,True)
    eta,cen,sigma,rho=fit.unpack(opt.x)
    dic={'dataset':name,'tag':tag,'series_ms':fit.series,'mask_half_width_MHz':fit.mask_width,
         'common_eta':fit.common_eta,'free_centres':fit.free_centres,
         'eta':eta.tolist(),'effective_63_zeeman_MHz':cen.tolist(),'sigma_MHz':sigma,'isotope65_scale':rho,
         'amplitudes':amps.tolist(),'baseline':bg.tolist(),'normalized_SSE':sse,
         'masked_n':int(fit.mask.sum()),'normalized_RMS':float(np.sqrt(sse/fit.mask.sum())),
         'nuQ63_MHz':(NQR_SETS[name][:,0]/np.sqrt(1+eta**2/3)).tolist(),
         'params':opt.x.tolist(),'optimizer_success':bool(opt.success),
         'normalization':fit.yscale,
         'note':'Conditional exploratory fit. All transitions, isotropic magnetic shift, common Gaussian broadening. No confidence interval or unique site assignment is asserted.'}
    (ROOT/f'fit_{name}_{tag}.json').write_text(json.dumps(dic,indent=2))
    np.savetxt(ROOT/f'fit_{name}_{tag}.csv',np.column_stack([fit.x,fit.y,pred,fit.mask]),delimiter=',',header='frequency_MHz,observed,predicted,used_in_fit',comments='')
    return dic

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--set',default='thesis2009');ap.add_argument('--extended',action='store_true');args=ap.parse_args()
    name=args.set;lib=build_library(name)
    all_results=[]
    for series,width,tag in [('500',.12,'500ms_main'),('20',.12,'20ms_main'),('500',.08,'500ms_narrow_mask'),('500',.18,'500ms_wide_mask')]:
        fit=Fit(lib,series,width);opt=fit.optimize();dic=summary(fit,opt,name,tag);all_results.append(dic)
        print(tag,json.dumps(dic),flush=True)
        if tag=='500ms_main':
            prof=fit.profile_eta(opt.x);np.savetxt(ROOT/f'eta_profile_{name}.csv',prof,delimiter=',',header='eta,SSE,effective63_MHz,sigma_MHz,isotope65_scale',comments='')
    if args.extended:
        for common,centres,tag in [(False,False,'500ms_four_eta'),(True,True,'500ms_four_centres'),(False,True,'500ms_four_eta_four_centres')]:
            fit=Fit(lib,'500',.12,common,centres);opt=fit.optimize(maxiter=180,popsize=10);dic=summary(fit,opt,name,tag);all_results.append(dic);print(tag,json.dumps(dic),flush=True)
    (ROOT/f'summary_{name}.json').write_text(json.dumps(all_results,indent=2))

if __name__=='__main__':main()
