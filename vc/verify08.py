import sys, time, json
import numpy as np, networkx as nx, scipy.sparse as sp
from scipy.optimize import milp, linprog, LinearConstraint, Bounds
sys.path.insert(0, '.')
from synthetic import planted
from chybp import ChyBP
rho, c = 0.8, 0.5
n, e, b = planted(150, rho, c, seed=1)
print('instance hash', hash(e.tobytes()) % 10**8, 'n', n, 'm', len(e), flush=True)
bp = ChyBP(n, e, b, smax=64); bp.iterate(iters=400); cnt, cov = bp.decimate()
G = nx.Graph(e.tolist())
uncovered = [v for v in range(n) if not cov[v]]
bad = sum(1 for a, b_ in G.edges() if not cov[a] and not cov[b_])
indep = nx.is_independent_set(G, set(uncovered)) if hasattr(nx, 'is_independent_set') else (bad == 0)
print('complex cover', cnt, 'violated edges', bad, 'uncovered independent', indep, flush=True)
np.save('results/cover_0.8_0.5.npy', cov)
m = len(e); A = sp.csr_matrix((np.ones(2 * m), (np.repeat(np.arange(m), 2), e.ravel())), shape=(m, n))
lp = linprog(np.ones(n), A_ub=-A, b_ub=-np.ones(m), bounds=(0, 1), method='highs')
print('LP relaxation bound', lp.fun, lp.status, flush=True)
for presolve in (True, False):
    t = time.time()
    res = milp(c=np.ones(n), constraints=LinearConstraint(A, lb=1, ub=np.inf), integrality=np.ones(n),
               bounds=Bounds(0, 1), options=dict(time_limit=600, presolve=presolve, disp=False))
    x = res.x
    ok = None if x is None else bool(np.all((x[e[:, 0]] > 0.5) | (x[e[:, 1]] > 0.5)))
    print(f'milp presolve={presolve}: fun={res.fun} status={res.status} gap={res.mip_gap} dual={res.mip_dual_bound} '
          f'solution valid cover={ok} sum(x)={None if x is None else x.sum():} [{time.time()-t:.0f}s]', flush=True)
