import json
import numpy as np
from scipy.optimize import differential_evolution,minimize
from powder_fit import *

name='thesis2009'
lib=build_library(name)
long=Fit(lib,'500',.12,False,True)
short=Fit(lib,'20',.12,False,True)
bounds=[(0,1)]*4+[(74.8,76.8)]*4+[(.025,.25),(.05,.8),(.05,.8)]
def objective(p):
    return long.solve(p[:10])+short.solve(np.r_[p[:9],p[10]])
records=[]
for seed in [11,29,67]:
    de=differential_evolution(objective,bounds,seed=seed,maxiter=200,popsize=9,tol=5e-6,polish=False)
    opt=minimize(objective,de.x,method='Nelder-Mead',bounds=bounds,options={'maxiter':4500,'xatol':1e-7,'fatol':1e-10})
    a=type('Result',(),{'x':opt.x[:10],'success':opt.success})()
    b=type('Result',(),{'x':np.r_[opt.x[:9],opt.x[10]],'success':opt.success})()
    rl=summary(long,a,name,f'joint_seed{seed}_500ms')
    rs=summary(short,b,name,f'joint_seed{seed}_20ms')
    record={'seed':seed,'SSE_joint':float(opt.fun),'params':opt.x.tolist(),'long':rl,'short':rs}
    records.append(record)
    print(json.dumps(record),flush=True)
(ROOT/'joint_results.json').write_text(json.dumps(records,indent=2))
