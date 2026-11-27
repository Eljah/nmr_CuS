"""Control model: central transitions and nonnegative empirical baseline."""
import json
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import nnls,minimize
from powder_fit import *

def central_library():
    file=ROOT/'central_library_thesis2009_40x80.npz'
    if file.exists():return np.load(file)['templates']
    dirs=directions(40,80)
    a=np.empty((len(ETAS),len(SIGMAS),4,2,len(REL)),np.float32)
    for ie,eta in enumerate(ETAS):
        for site,pair in enumerate(NQR_SETS['thesis2009']):
            for iso,nqr in enumerate(pair):
                raw=stick_profile(eta,nqr,REFERENCE*(GAMMA_RATIO if iso else 1),dirs,central_only=True)
                for si,sigma in enumerate(SIGMAS):a[ie,si,site,iso]=gaussian_filter1d(raw,sigma/STEP,mode='constant')
    np.savez_compressed(file,templates=a);return a

class CentralFit(Fit):
    def solve(self,p,ret=False):
        a=self.design(p)
        t=(self.x-self.x.min())/(self.x.max()-self.x.min())
        b=np.column_stack([1-t,t])
        design=np.column_stack([a,b])
        coefs=nnls(design[self.mask],self.yy[self.mask])[0]
        pred=design@coefs;res=(pred-self.yy)[self.mask]
        if ret:return float(res@res),pred*self.yscale,coefs[:4]*self.yscale,coefs[4:]*self.yscale
        return float(res@res)

def main():
    lib=central_library();long=CentralFit(lib,'500',.12,False,True);short=CentralFit(lib,'20',.12,False,True)
    bounds=[(0,1)]*4+[(74.8,76.8)]*4+[(.025,.25),(.05,.8),(.05,.8)]
    def objective(p):return long.solve(p[:10])+short.solve(np.r_[p[:9],p[10]])
    source=json.loads((ROOT/'exact_joint_refinements.json').read_text());out=[]
    for orig in source:
        p=orig['parameters'];initial=objective(p)
        opt=minimize(objective,p,method='Nelder-Mead',bounds=bounds,options={'maxiter':6000,'xatol':1e-7,'fatol':1e-10})
        p=opt.x;rec={'source_seed':orig['source_seed'],'model':'Central transitions with nonnegative linear baseline at both frequency-range endpoints','initial_SSE':initial,'SSE_joint':float(opt.fun),'parameters':p.tolist(),'eta':p[:4].tolist(),'nuZ63_MHz':p[4:8].tolist(),'sigma_MHz':float(p[8]),'rho500':float(p[9]),'rho20':float(p[10]),'nuQ63_MHz':(NQR_SETS['thesis2009'][:,0]/np.sqrt(1+p[:4]**2/3)).tolist(),'success':bool(opt.success),'message':opt.message}
        for fit,params,key in [(long,p[:10],'long'),(short,np.r_[p[:9],p[10]],'short')]:
            sse,pred,amps,base=fit.solve(params,True)
            rec[key]={'normalized_SSE':sse,'amplitudes':amps.tolist(),'nonnegative_baseline_endpoints':base.tolist()}
            np.savetxt(ROOT/f'central_control_seed{orig["source_seed"]}_{fit.series}ms.csv',np.column_stack([fit.x,fit.y,pred,fit.mask]),delimiter=',',header='frequency_MHz,observed,predicted,used_in_fit',comments='')
        out.append(rec);print(json.dumps(rec),flush=True)
    (ROOT/'central_control_results.json').write_text(json.dumps(out,indent=2))

if __name__=='__main__':main()
