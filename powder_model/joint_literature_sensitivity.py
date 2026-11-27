import json
import numpy as np
from scipy.optimize import minimize
from powder_fit import *

original=json.loads((ROOT/'joint_results.json').read_text())
name='publication2013';lib=build_library(name)
long=Fit(lib,'500',.12,False,True);short=Fit(lib,'20',.12,False,True)
bounds=[(0,1)]*4+[(74.8,76.8)]*4+[(.025,.25),(.05,.8),(.05,.8)]
def objective(p):return long.solve(p[:10])+short.solve(np.r_[p[:9],p[10]])
results=[]
for orig in [min(original,key=lambda r:r['SSE_joint']),[r for r in original if r['seed']==67][0]]:
    p=np.array(orig['params'])
    fit=minimize(objective,p,method='Nelder-Mead',bounds=bounds,options={'maxiter':6000,'xatol':1e-7,'fatol':1e-11})
    r1=type('Result',(),{'x':fit.x[:10],'success':fit.success})();r2=type('Result',(),{'x':np.r_[fit.x[:9],fit.x[10]],'success':fit.success})()
    a=summary(long,r1,name,f'joint_warm{orig["seed"]}_500ms');b=summary(short,r2,name,f'joint_warm{orig["seed"]}_20ms')
    record={'warm_start_seed':orig['seed'],'SSE_joint':float(fit.fun),'params':fit.x.tolist(),'long':a,'short':b}
    results.append(record);print(json.dumps(record),flush=True)
(ROOT/'joint_literature_sensitivity.json').write_text(json.dumps(results,indent=2))
