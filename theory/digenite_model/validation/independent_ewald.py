"""Independent NumPy/SciPy Ewald and neutral-cell reference calculation.

Run from any directory. --full-direct reproduces the saved n=40,80 and
40/80/120-A spherical crosschecks and is more memory/time intensive.
G is in A^-3 for dimensionless charges q/e. The variable eta inside ewald()
is its splitting parameter kappa, not the EFG asymmetry called asym below.
The independent SI conversion differs from the Java input by <1e-9 relative.
"""
from pathlib import Path
import sys
import numpy as np
from scipy.special import erfc
import json, time

AR = 16.16
ALPHA = np.deg2rad(13 + 56/60)
c = np.cos(ALPHA)
az = AR*np.sqrt((1+2*c)/3)
ap = AR*np.sqrt(2*(1-c)/3)
avec = np.array([[ap,0,az],[-ap/2,ap*np.sqrt(3)/2,az],[-ap/2,-ap*np.sqrt(3)/2,az]])
bvec = 2*np.pi*np.linalg.inv(avec).T
volume = abs(np.linalg.det(avec))
xs = np.array([.5,.060,-.060,.133,-.133,.250,-.250,.350,-.350,0,.2,-.2,.4,-.4])
frac = np.repeat(xs[:,None],3,axis=1)
pos = frac@avec
q = np.array([10/9]*9+[-2]*5)
probe_indices = [0,1,3,5,7]
k_e = 8.9875517923e9
e = 1.602176634e-19
h = 6.62607015e-34
unit = k_e*e*1e30

def grid(n):
    r = np.arange(-n,n+1)
    return np.stack(np.meshgrid(r,r,r,indexing='ij'),axis=-1).reshape(-1,3)

def ewald(eta=0.3, exponent=32):
    rc = np.sqrt(exponent)/eta
    gc = 2*eta*np.sqrt(exponent)
    nr = int(np.ceil(rc*np.linalg.norm(bvec,axis=1).max()/(2*np.pi)+1.5))
    ng = int(np.ceil(gc*AR/(2*np.pi)))
    trans = grid(nr)@avec
    g = grid(ng)@bvec
    g2 = np.sum(g*g,axis=1)
    g = g[(g2>1e-14)&(g2<=gc*gc)]
    g2 = np.sum(g*g,axis=1)
    rho = np.exp(-1j*(g@pos.T))@q
    result=[]
    for pidx in probe_indices:
        origin = pos[pidx]
        vreal = np.zeros((3,3))
        for rj,qj in zip(pos,q):
            d = trans+rj-origin
            r2 = np.sum(d*d,axis=1)
            d = d[(r2>1e-18)&(r2<=rc*rc)]
            r2 = np.sum(d*d,axis=1)
            r = np.sqrt(r2)
            exp = np.exp(-eta*eta*r2)
            er = erfc(eta*r)
            a = 3*er/r**5 + (6*eta/r**4+4*eta**3/r**2)/np.sqrt(np.pi)*exp
            b = er/r**3 + 2*eta/np.sqrt(np.pi)*exp/r**2
            vreal += qj*(np.einsum('n,ni,nj->ij',a,d,d)-np.eye(3)*b.sum())
        rg = np.real(rho*np.exp(1j*(g@origin)))
        weight = -(4*np.pi/volume)*np.exp(-g2/(4*eta*eta))*rg/g2
        vrecip = np.einsum('n,ni,nj->ij',weight,g,g)
        vself = np.eye(3)*4*q[pidx]*eta**3/(3*np.sqrt(np.pi))
        total=vreal+vrecip+vself
        eig=np.linalg.eigvalsh(total)
        eig=eig[np.argsort(np.abs(eig))]
        vzz=eig[-1]*unit
        asym=abs((eig[0]-eig[1])/eig[-1])
        nu63=abs(e*(-.211e-28)*vzz)/(2*h)*np.sqrt(1+asym**2/3)/1e6
        nu65=abs(e*(-.195e-28)*vzz)/(2*h)*np.sqrt(1+asym**2/3)/1e6
        result.append(dict(x=xs[pidx],tensor_raw=total.tolist(),eig_SI=(eig*unit).tolist(),Vzz=vzz,eta=asym,nu63_MHz=nu63,nu65_MHz=nu65,real=vreal.tolist(),reciprocal=vrecip.tolist(),self=vself.tolist()))
    return dict(method='ewald',eta=eta,exponent=exponent,real_cutoff=rc,reciprocal_cutoff=gc,nr=nr,ng=ng,result=result)

def neutral_cell_sum(n):
    trans = grid(n)@avec
    ans=[]
    for pidx in probe_indices:
        total=np.zeros((3,3))
        for rj,qj in zip(pos,q):
            d = trans+rj-pos[pidx]
            r2 = np.sum(d*d,axis=1)
            d=d[r2>1e-18]
            r2=np.sum(d*d,axis=1)
            total += qj*(3*np.einsum('n,ni,nj->ij',r2**-2.5,d,d)-np.eye(3)*(r2**-1.5).sum())
        ans.append((total*unit).tolist())
    return dict(method='cell_cube',n=n,tensors=ans)

def spherical_sum(rc):
    nr=int(np.ceil(rc*np.linalg.norm(bvec,axis=1).max()/(2*np.pi)+1.5))
    trans=grid(nr)@avec
    ans=[]
    for pidx in probe_indices:
        total=np.zeros((3,3))
        for rj,qj in zip(pos,q):
            d=trans+rj-pos[pidx]
            r2=np.sum(d*d,axis=1)
            d=d[(r2>1e-18)&(r2<=rc*rc)]
            r2=np.sum(d*d,axis=1)
            total += qj*(3*np.einsum('n,ni,nj->ij',r2**-2.5,d,d)-np.eye(3)*(r2**-1.5).sum())
        ans.append((total*unit).tolist())
    return dict(method='sphere',R=rc,nr=nr,tensors=ans)

if __name__=='__main__':
    print(json.dumps(dict(AR=AR,alpha_deg=np.rad2deg(ALPHA),volume=volume,primitive=avec.tolist(),reciprocal=bvec.tolist(),netq=q.sum())))
    allres=[]
    for eta in [.15,.25,.4]:
        r=ewald(eta,36)
        allres.append(r)
        print(json.dumps(r))
    for n in ([10,20,40,80] if '--full-direct' in sys.argv else [10,20]):
        r=neutral_cell_sum(n)
        allres.append(r)
        print(json.dumps(r))
    for radius in ([40,80,120] if '--full-direct' in sys.argv else []):
        r=spherical_sum(radius)
        allres.append(r)
        print(json.dumps(r))
    with (Path(__file__).resolve().parent / 'independent_ewald_recomputed.json').open('w') as f:
        json.dump(allres,f,indent=2)
