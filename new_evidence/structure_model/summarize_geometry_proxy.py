"""Make a concise structural-prior table from the explicit proxy calculation."""
from pathlib import Path
import json
import csv
import numpy as np
ROOT=Path(__file__).resolve().parent
r=json.loads((ROOT/'geometry_efg_proxy.json').read_text())
models=list(r['formal_charge_models'])
ref=next(s for s in r['sites'] if s['site']=='Cu2')['periodic_ewald']['uniform_Cu_8over7']['abs_Vzz_proxy']
rows=[]
for s in r['sites']:
    d=s['periodic_ewald']; u=d['uniform_Cu_8over7']
    rows.append({'site':s['site'],'Wyckoff':s['wyckoff'],'site_symmetry':s['site_symmetry'],
                 'multiplicity':s['multiplicity'],'coordination':'distorted_tetrahedral' if s['short_S_count']==4 else 'triangular',
                 'eta_uniform_charge_Ewald':u['eta'],'relative_abs_Vzz_uniform_charge_Cu2_1':u['abs_Vzz_proxy']/ref,
                 'eta_min_six_charge_models':min(d[m]['eta'] for m in models),
                 'eta_max_six_charge_models':max(d[m]['eta'] for m in models),
                 'abs_Vzz_min_six_charge_models':min(d[m]['abs_Vzz_proxy'] for m in models),
                 'abs_Vzz_max_six_charge_models':max(d[m]['abs_Vzz_proxy'] for m in models),
                 'maximum_Ewald_alpha_change':max(d[m]['maximum_tensor_change'] for m in models)})
with (ROOT/'structural_priors_table.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
vertices=[m for m in models if m!='uniform_Cu_8over7']
directional=[]
for s in r['sites']:
    d=s['periodic_ewald']
    axis=np.array(d['uniform_Cu_8over7']['principal_axes_columns_xx_yy_zz'])[:,2]
    values=[float(axis@np.array(d[m]['tensor_abc'])@axis) for m in vertices]
    directional.append({'site':s['site'],'signed_projection_min':min(values),'signed_projection_max':max(values)})
result={'table':rows,'directional_bounds':directional,
        'note':'Values are exploratory point-charge proxies. Six-charge-model ranges are model sensitivity, not uncertainty intervals.',
        'Cu3_rank_argument':'The spectral norm of Cu3 is at most 0.08088 for any convex combination of the five charge-redistribution vertices. At each other site a projection on its uniform-model principal axis has one sign at all five vertices and a magnitude at least 0.28488 (Cu1), 0.40150 (Cu2), 0.34181 (Cu4), or 0.38664 (Cu5). Hence Cu3 remains the smallest tensor norm throughout that restricted charge simplex, not merely at its vertices.'}
(ROOT/'structural_priors_summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
