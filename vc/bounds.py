"""Whole-graph MILP bounds (incumbent, dual bound) for the planted rho=0.8 cases,
where HiGHS cannot prove optimality: results/bounds_<rho>_<c>.json"""
import json, sys, time
import numpy as np, scipy.sparse as sp
from scipy.optimize import milp, LinearConstraint, Bounds
sys.path.insert(0, '.')
from synthetic import planted
rho, c, tl = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
n, e, b = planted(150, rho, c, seed=1)
m = len(e); A = sp.csr_matrix((np.ones(2 * m), (np.repeat(np.arange(m), 2), e.ravel())), shape=(m, n))
t = time.time()
res = milp(c=np.ones(n), constraints=LinearConstraint(A, lb=1, ub=np.inf), integrality=np.ones(n),
           bounds=Bounds(0, 1), options=dict(time_limit=tl, disp=False))
out = dict(rho=rho, c=c, n=n, m=m, incumbent=res.fun, dual_bound=res.mip_dual_bound, gap=res.mip_gap,
           status=int(res.status), secs=time.time() - t)
print(out, flush=True)
json.dump(out, open(f'results/bounds_{rho}_{c}.json', 'w'))
