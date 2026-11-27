"""Reproduce the central-transition and satellite extent counterexample."""
import json
import numpy as np
from scipy.optimize import differential_evolution
from powder_fit import ROOT,SPIN,QUAD0,QUADE

rows=[]
for nqr in [1.5,16.96]:
    for eta in [0,1]:
        nuq=nqr/np.sqrt(1+eta*eta/3)
        def freq(p):
            u,phi=p
            n=np.array([np.sqrt(1-u*u)*np.cos(phi),np.sqrt(1-u*u)*np.sin(phi),u])
            h=-76*np.einsum('a,aij->ij',n,SPIN)+nuq/6*(QUAD0+eta*QUADE)
            return np.diff(np.linalg.eigvalsh(h))
        extrema=[]
        for t in range(3):
            lo=differential_evolution(lambda p:freq(p)[t],[(0,1),(0,np.pi/2)],tol=1e-12,seed=42)
            hi=differential_evolution(lambda p:-freq(p)[t],[(0,1),(0,np.pi/2)],tol=1e-12,seed=43)
            extrema.append({'transition':t,'min_MHz':float(lo.fun),'max_MHz':float(-hi.fun),'min_offset_kHz':float((lo.fun-76)*1000),'max_offset_kHz':float((-hi.fun-76)*1000)})
        rows.append({'nuNQR_MHz':nqr,'eta':eta,'nuQ_MHz':nuq,'nuL_MHz':76,'transitions':extrema})
(ROOT/'central_extent_examples.json').write_text(json.dumps(rows,indent=2))
