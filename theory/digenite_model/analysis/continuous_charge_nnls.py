#!/usr/bin/env python3
"""Global fit over the five-vertex convex family of neutral Cu charges.

Independent Python audit; the underlying EFG calculation remains Java.
Reads the Java-generated full-precision CSV relative to this script.
Run from any directory: python /path/to/continuous_charge_nnls.py
Dependencies: Python 3, numpy, scipy.

Domain: fixed Donnay geometry, five Cu-orbit extra-charge vertices,
w >= 0 and sum(w)=1, all sulfur charges -2, one additional positive
charge beyond all Cu(+1), and one common positive frequency scale a.
No valence EFG tensor, independent site scales, or symmetry breaking
is introduced.  The Cu multiplicities are (1,2,2,2,2), so the charges
are q_i=1+w_i/m_i.  All charge combinations remain neutral.

Exact equivalence used in the exhaustive numerical audit:
 min_{t>=0} ||abs(Mt)-y||^2
 = min_{s in {-1,+1}^4} min_{t>=0} ||diag(s)Mt-y||^2.
The equality follows by choosing the sign of each component of Mt.
The t=0 case includes the a=0 boundary; the reported optimum has a>0.
The five choices of omitted site and all 24 assignments are enumerated,
for 5*24*16 = 1920 convex NNLS problems. KKT residuals and the actual
absolute-value predictions are verified. Outputs are written to the
model results directory, and the input CSV SHA-256 is recorded.
For a>0 and w>=0 with sum(w)=1, use t=a*w>=0.  Every assignment
and every sign branch is linear nonnegative least squares.  Their minimum
equals min_t ||abs(M_selected @ t)-observations||_2, since the sign choice
can be minimized independently for every component.  Actual absolute
frequencies and KKT residuals are checked after each solve.
"""
import csv
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.optimize import nnls

MODEL_ROOT=Path(__file__).resolve().parents[1]
SOURCE=MODEL_ROOT/'results/tensors_and_frequencies_30.csv'
OUTPUT=MODEL_ROOT/'results/continuous_charge_optimum.json'
CSV_OUTPUT=MODEL_ROOT/'results/continuous_charge_best_by_omitted_site.csv'
SITES=['DgA','DgB','DgC','DgD','DgE']
VERTICES=['extra_charge_'+site for site in SITES]
OBS=np.array([11.97,15.27,16.18,16.96])
rows=list(csv.DictReader(SOURCE.open()))
lookup={(r['scenario'],r['site']):r for r in rows}
M=np.array([[float(lookup[v,s]['principal_Gzz_signed']) for v in VERTICES] for s in SITES])
uniform=np.array([float(lookup['uniform_Cu10over9',s]['principal_Gzz_signed']) for s in SITES])
uniform_w=np.array([1,2,2,2,2])/9
assert np.max(np.abs(M@uniform_w-uniform))<1e-12
assert np.max([abs(float(r['eta'])) for r in rows])<1e-10

records=[]
max_kkt=0.0
nnls_count=0
for selected in itertools.permutations(range(5),4):
    selected=np.array(selected)
    matrix=M[selected,:]
    omitted=next(i for i in range(5) if i not in selected)
    for signs_tuple in itertools.product([-1.,1.],repeat=4):
        signs=np.array(signs_tuple)
        A=signs[:,None]*matrix
        t,rnorm=nnls(A,OBS,maxiter=1000)
        grad=A.T@(A@t-OBS)
        active=t>1e-10
        kkt=max(np.max(np.abs(grad[active])) if np.any(active) else 0.,
                np.max(np.maximum(-grad[~active],0.)) if np.any(~active) else 0.)
        max_kkt=max(max_kkt,kkt)
        assert kkt<1e-9
        all_signed=M@t
        all_freq=np.abs(all_signed)
        actual_residual=all_freq[selected]-OBS
        actual_rmse=float(np.sqrt(np.mean(actual_residual**2)))
        scale=float(t.sum())
        record={
          'selected_sites_by_observed_order':[SITES[i] for i in selected],
          'omitted_site':SITES[omitted],
          'sign_branch':signs.tolist(),
          't':t.tolist(),
          'scale_MHz_per_raw_G':scale,
          'weights_extra_A_B_C_D_E':(t/scale).tolist() if scale else None,
          'charges_CuA_B_C_D_E':(1+(t/scale)/np.array([1,2,2,2,2])).tolist() if scale else None,
          'branch_rmse_MHz':float(rnorm/2),
          'actual_rmse_MHz':actual_rmse,
          'actual_selected_predictions_MHz':all_freq[selected].tolist(),
          'actual_selected_residuals_MHz':actual_residual.tolist(),
          'all_signed_scale_times_G_A_B_C_D_E':all_signed.tolist(),
          'all_frequencies_MHz_A_B_C_D_E':all_freq.tolist(),
          'omitted_frequency_MHz':float(all_freq[omitted]),
          'sign_consistent':bool(np.all(signs*all_signed[selected]>=-1e-10)),
          'kkt_violation':float(kkt),
        }
        records.append(record)
        nnls_count+=1

best=min(records,key=lambda r:r['branch_rmse_MHz'])
best_actual=min(records,key=lambda r:r['actual_rmse_MHz'])
assert abs(best['actual_rmse_MHz']-best['branch_rmse_MHz'])<1e-10
assert abs(best_actual['actual_rmse_MHz']-best['branch_rmse_MHz'])<1e-10
assert best['sign_consistent']
best_omitted={s:min((r for r in records if r['omitted_site']==s),key=lambda r:r['branch_rmse_MHz']) for s in SITES}
result={
  'source_csv':str(SOURCE.relative_to(MODEL_ROOT)),
  'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
  'matrix_M_rows_sites_columns_extra_charge_orbits_raw_G':M.tolist(),
  'site_order':SITES,'vertex_order':VERTICES,'observations_MHz':OBS.tolist(),
  'uniform_weights':uniform_w.tolist(),
  'uniform_reconstruction_max_abs_error':float(np.max(np.abs(M@uniform_w-uniform))),
  'NNLS_problems':nnls_count,'max_NNLS_KKT_violation':max_kkt,
  'global_best':best,'best_by_omitted_site':best_omitted,
  'component_min_raw_G':M.min(axis=1).tolist(),
  'component_max_raw_G':M.max(axis=1).tolist(),
  'analytic_frequency_ratio_bound':{
    'observed_max_over_min':float(OBS.max()/OBS.min()),
    'if_C_selected_lower_bound':float(M[2].min()/max(np.max(np.abs(M[3])),np.max(np.abs(M[4])))),
    'if_C_omitted_lower_bound':float(M[1].min()/np.max(np.abs(M[4]))),
    'explanation':'If C is selected among four sites, at least one of D or E is also selected. C>=min(M_C)>0, while |D|<=max|M_D| and |E|<=max|M_E|. If C is omitted, B and E are both selected; B>=min(M_B)>0 and |E|<=max|M_E|. A common positive scale cancels from both bounds.'
  },
}
OUTPUT.write_text(json.dumps(result,indent=2)+'\n')
with CSV_OUTPUT.open('w',newline='') as fh:
    fields=['omitted_site','rmse_MHz','scale_MHz_per_raw_G','assignment_by_observed_order','omitted_frequency_MHz','sign_consistent']
    fields += ['weight_extra_'+s for s in SITES]
    fields += ['charge_'+s for s in SITES]
    fields += ['frequency_MHz_'+s for s in SITES]
    writer=csv.DictWriter(fh,fieldnames=fields,lineterminator='\n')
    writer.writeheader()
    for site in SITES:
        r=best_omitted[site]
        row=dict(omitted_site=site,rmse_MHz=r['actual_rmse_MHz'],scale_MHz_per_raw_G=r['scale_MHz_per_raw_G'],assignment_by_observed_order=';'.join(r['selected_sites_by_observed_order']),omitted_frequency_MHz=r['omitted_frequency_MHz'],sign_consistent=r['sign_consistent'])
        row.update({'weight_extra_'+s:v for s,v in zip(SITES,r['weights_extra_A_B_C_D_E'])})
        row.update({'charge_'+s:v for s,v in zip(SITES,r['charges_CuA_B_C_D_E'])})
        row.update({'frequency_MHz_'+s:v for s,v in zip(SITES,r['all_frequencies_MHz_A_B_C_D_E'])})
        writer.writerow(row)
print(json.dumps(dict(NNLS_problems=nnls_count,max_NNLS_KKT_violation=max_kkt,global_best=best,analytic_frequency_ratio_bound=result['analytic_frequency_ratio_bound']),indent=2))
