"""Geometrical EFG proxies for anilite. NOT an electronic-structure prediction.

Uses the complete, separately preserved transcription of AMCSD 0009376.
Units of the tensor are elementary-charge / angstrom**3, omitting 1/(4*pi*eps0).
Only eigenvalue ratios, relative tensors and their model sensitivity are used.
No antishielding, covalence, on-site electronic EFG, or NQR frequencies are fitted.
"""
from pathlib import Path
from itertools import product
import csv
import json
import re
import numpy as np
from scipy.special import erfc

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parents[1]
CIF = BASE / 'research_literature/anilite_amcsd0009376_transcribed.cif'
text = CIF.read_text()
cell = np.array([float(re.search(rf'_cell_length_{a}\s+(\S+)', text).group(1)) for a in 'abc'])
atoms = {m.group(1): np.array([float(m.group(i)) for i in (2,3,4)])
         for m in re.finditer(r'^(Cu\d+|S\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)', text, re.M)}

def orbit(p):
    x,y,z=p
    return np.array(sorted({tuple(round(v%1,8) for v in op) for op in [
        (x,y,z),(x,.5-y,z),(-x,.5+y,-z),(.5-x,.5+y,.5+z),
        (.5+x,.5-y,.5-z),(.5+x,y,.5-z),(.5-x,-y,.5+z),(-x,-y,-z)]}))

orbits = {label:orbit(p) for label,p in atoms.items()}
full_labels=[]; full_positions=[]
for label,positions in orbits.items():
    full_labels += [label]*len(positions)
    full_positions.extend(positions)
full_labels=np.array(full_labels); full_positions=np.array(full_positions)
assert len(full_labels)==44
assert sum(x.startswith('Cu') for x in full_labels)==28
assert sum(x.startswith('S') for x in full_labels)==16

def tensor_kernel(r):
    r2=np.einsum('ni,ni->n',r,r)
    return (3*r[:,:,None]*r[:,None,:]-r2[:,None,None]*np.eye(3))/r2[:,None,None]**2.5

def tensor_properties(t):
    t=(t+t.T)/2
    w,v=np.linalg.eigh(t)
    order=np.argsort(np.abs(w))
    xx,yy,zz=w[order]
    eta=float((xx-yy)/zz) if abs(zz)>1e-15 else None
    return {'eta':eta,'eigenvalues_xx_yy_zz':w[order].tolist(),
            'abs_Vzz_proxy':float(abs(zz)),
            'principal_axes_columns_xx_yy_zz':v[:,order].tolist(),
            'tensor_abc':t.tolist(),'trace':float(np.trace(t))}

def ewald_tensor(position, charges_by_label, alpha):
    """Periodic Coulomb Hessian, self charge removed; alpha in inverse angstrom."""
    sitecart=position*cell
    rc=6/alpha
    nc=np.ceil(rc/cell).astype(int)+1
    trans=np.array(list(product(*(range(-n,n+1) for n in nc))))
    pos=(trans[:,None,:]+full_positions[None,:,:]).reshape(-1,3)*cell
    q0=np.array([charges_by_label[l] for l in full_labels])
    q=np.tile(q0,len(trans))
    r=pos-sitecart
    norm=np.linalg.norm(r,axis=1)
    select=(norm>1e-8)&(norm<rc)
    r=r[select];norm=norm[select];q=q[select]
    ex=np.exp(-(alpha*norm)**2)
    a=3*erfc(alpha*norm)/norm**5+6*alpha/np.sqrt(np.pi)*ex/norm**4+4*alpha**3/np.sqrt(np.pi)*ex/norm**2
    b=-erfc(alpha*norm)/norm**3-2*alpha/np.sqrt(np.pi)*ex/norm**2
    real=np.einsum('n,ni,nj->ij',q*a,r,r)+np.eye(3)*np.sum(q*b)
    kc=12*alpha
    nk=np.ceil(kc*cell/(2*np.pi)).astype(int)
    k=np.array(list(product(*(range(-n,n+1) for n in nk))))*2*np.pi/cell
    k2=np.einsum('ni,ni->n',k,k)
    select=(k2>0)&(k2<kc**2)
    k=k[select];k2=k2[select]
    charge_cos=np.cos(k@(sitecart[None,:]-full_positions*cell).T)@q0
    reciprocal=-4*np.pi/np.prod(cell)*np.einsum('n,ni,nj->ij',np.exp(-k2/(4*alpha**2))*charge_cos/k2,k,k)
    # The self Gaussian has an isotropic Hessian. Remove it explicitly.
    same=np.linalg.norm((full_positions-position)*cell,axis=1)<1e-8
    self_q=float(q0[same].sum())
    total=real+reciprocal+4*alpha**3/(3*np.sqrt(np.pi))*self_q*np.eye(3)
    trace_before=float(np.trace(total))
    assert abs(trace_before)<1e-9,trace_before
    total-=np.trace(total)/3*np.eye(3)
    return total,trace_before

max_cutoff=100.0
# A sphere centred on each Cu site is sampled, then radial cutoffs are applied.
nmax=np.ceil(max_cutoff/cell).astype(int)+1
translations=np.array(list(product(*(range(-n,n+1) for n in nmax))))
all_positions=(translations[:,None,:]+full_positions[None,:,:]).reshape(-1,3)*cell
all_labels=np.tile(full_labels,len(translations))
is_copper=np.array([x.startswith('Cu') for x in all_labels])

scenarios={'uniform_Cu_8over7': {l:(8/7 if l.startswith('Cu') else -2.) for l in atoms}}
for target in ['Cu1','Cu2','Cu3','Cu4','Cu5']:
    q={l:(1. if l.startswith('Cu') else -2.) for l in atoms}
    q[target]+=4/len(orbits[target])
    scenarios[f'excess_charge_on_{target}']=q
for label,q in scenarios.items():
    total=sum(len(orbits[l])*v for l,v in q.items())
    assert abs(total)<1e-12,(label,total)

sites=[]; convergence=[]; proxy_rows=[]
cutoffs=[15.,20.,30.,40.,60.,80.,100.]
for site,p in atoms.items():
    if not site.startswith('Cu'): continue
    delta=all_positions-p*cell
    radii=np.linalg.norm(delta,axis=1)
    use=(radii>1e-8)&(radii<=max_cutoff)
    r=delta[use]; radius=radii[use]; labels=all_labels[use]
    kernel=tensor_kernel(r)
    sulfur=np.array([l.startswith('S') for l in labels])
    inds=np.where(sulfur)[0][np.argsort(radius[sulfur])[:6]]
    geometric={'site':site,'multiplicity':len(orbits[site]),
               'wyckoff':'4c' if len(orbits[site])==4 else '8d',
               'site_symmetry':'.m.' if len(orbits[site])==4 else '1',
               'nearest_S_A':radius[inds].tolist(),
               'nearest_S_vectors_A':r[inds].tolist()}
    first=np.where(sulfur&(radius<2.6))[0]
    vec=r[first]
    geometric['short_S_count']=len(first)
    geometric['short_S_angles_deg']=[float(np.degrees(np.arccos(np.clip(np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b),-1,1))))
                                      for ia,a in enumerate(vec) for b in vec[ia+1:]]
    if len(first)==3:
        normal=np.cross(vec[1]-vec[0],vec[2]-vec[0]);normal/=np.linalg.norm(normal)
        geometric['Cu_distance_from_three_S_plane_A']=float(abs(np.dot(normal,vec[0])))
        geometric['three_S_plane_normal_abc']=normal.tolist()
    for cutoff in [2.6,3.3]:
        sel=sulfur&(radius<cutoff)
        props=tensor_properties(-2*kernel[sel].sum(axis=0))
        geometric[f'S_shell_below_{cutoff}A']={'n_ligands':int(sel.sum()),**props}
        proxy_rows.append({'site':site,'model':f'S_only_shell_{cutoff}A','cutoff_A':cutoff,
                           'eta':props['eta'],'abs_Vzz_proxy':props['abs_Vzz_proxy']})
    for name,q in scenarios.items():
        charges=np.array([q[l] for l in labels])
        weighted=kernel*charges[:,None,None]
        percut={}
        for cutoff in cutoffs:
            props=tensor_properties(weighted[radius<=cutoff].sum(axis=0))
            percut[str(cutoff)]=props
            convergence.append({'site':site,'charge_model':name,'cutoff_A':cutoff,
                                'eta':props['eta'],'abs_Vzz_proxy':props['abs_Vzz_proxy']})
        geometric[name]=percut
        proxy_rows.append({'site':site,'model':name,'cutoff_A':100.,
                          'eta':percut['100.0']['eta'],'abs_Vzz_proxy':percut['100.0']['abs_Vzz_proxy']})
    geometric['periodic_ewald']={}
    for name,q in scenarios.items():
        vals={}
        for alpha in [.2,.3,.4]:
            tensor,tr=ewald_tensor(p,q,alpha)
            vals[str(alpha)]={**tensor_properties(tensor),'trace_before_projection':tr}
        reference=np.array(vals['0.3']['tensor_abc'])
        error=max(float(np.max(np.abs(np.array(v['tensor_abc'])-reference))) for v in vals.values())
        assert error<1e-9,error
        geometric['periodic_ewald'][name]={'alpha_checks':vals,'maximum_tensor_change':error,**vals['0.3']}
        proxy_rows.append({'site':site,'model':'periodic_Ewald_'+name,'cutoff_A':'periodic',
                           'eta':vals['0.3']['eta'],'abs_Vzz_proxy':vals['0.3']['abs_Vzz_proxy']})
    sites.append(geometric)

output={'source':str(CIF.relative_to(BASE)),
        'model_status':'Unscreened geometry-only point-charge proxy; not a physical prediction of Cu EFG or eta.',
        'cell_A':cell.tolist(),'formal_charge_models':scenarios,
        'max_cutoff_A':max_cutoff,'sites':sites,
        'limitations':['Room-temperature average atomic positions, not a 5 K refinement.',
                       'All point charges are assumptions; no orbital charge density, covalence, screening or antishielding.',
                       'Formal charge redistribution examples test sensitivity; they are not proposed Cu valence assignments.',
                       'No conversion to nuclear quadrupole frequency or absolute physical EFG is made.']}
(ROOT/'geometry_efg_proxy.json').write_text(json.dumps(output,indent=2)+'\n')
for filename,rows in [('geometry_efg_proxy_summary.csv',proxy_rows),('geometry_efg_proxy_convergence.csv',convergence)]:
    with (ROOT/filename).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print(json.dumps({'source':output['source'],'models':list(scenarios),
                  'sites':[{'site':s['site'],'symmetry':s['site_symmetry'],
                            'short_S_count':s['short_S_count'],
                            'uniform_eta_Ewald':s['periodic_ewald']['uniform_Cu_8over7']['eta'],
                            'abs_Vzz_Ewald':s['periodic_ewald']['uniform_Cu_8over7']['abs_Vzz_proxy'],
                            'uniform_eta_spherical_100A':s['uniform_Cu_8over7']['100.0']['eta'],
                            'eta_80A':s['uniform_Cu_8over7']['80.0']['eta'],
                            'charge_model_eta_range':[min(s[m]['100.0']['eta'] for m in scenarios),max(s[m]['100.0']['eta'] for m in scenarios)],
                            'S_shell_eta_2p6A':s['S_shell_below_2.6A']['eta'],
                            'S_shell_eta_3p3A':s['S_shell_below_3.3A']['eta']} for s in sites]},indent=2))
