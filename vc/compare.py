"""Exact minimum cover (MILP) vs leaf removal vs node-level min-sum vs
community-complex min-sum, on the degree-correlated ensemble."""
import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.optimize import Bounds, LinearConstraint, milp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from chybp import ChyBP          # noqa: E402
from ensemble import graph, leaf_removal  # noqa: E402

RES = Path(__file__).resolve().parent / 'results'
OUT = RES / 'compare.json'


def exact_vc(n, e, time_limit):
    """Minimum cover by MILP, one connected component at a time (the r = 1
    graphs are unions of degree classes, and HiGHS on the whole thing at once
    runs out of memory).  Returns an object with .fun, .status, .mip_gap, .x."""
    import networkx as nx
    from types import SimpleNamespace
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(e.tolist())
    total, status, gap = 0.0, 0, 0.0
    for comp in nx.connected_components(G):
        if len(comp) < 2:
            continue
        comp = sorted(comp)
        idx = {v: k for k, v in enumerate(comp)}
        sub = np.array([(idx[a], idx[b]) for a, b in G.subgraph(comp).edges()])
        if len(comp) == 2:
            total += 1
            continue
        m, nn = len(sub), len(comp)
        A = sp.csr_matrix((np.ones(2 * m), (np.repeat(np.arange(m), 2), sub.ravel())), shape=(m, nn))
        res = milp(c=np.ones(nn), constraints=LinearConstraint(A, lb=1, ub=np.inf),
                   integrality=np.ones(nn), bounds=Bounds(0, 1),
                   options=dict(time_limit=time_limit, disp=False))
        if res.x is None:
            return SimpleNamespace(fun=None, status=res.status, mip_gap=None, x=None)
        total += round(res.fun)
        status = max(status, res.status)
        gap = max(gap, res.mip_gap or 0.0)
    return SimpleNamespace(fun=total, status=status, mip_gap=gap, x=True)


def run_bp(n, e, blocks, smax, iters, damping):
    t = time.time()
    bp = ChyBP(n, e, blocks, smax=smax)
    it, diff = bp.iterate(iters=iters, damping=damping)
    E = bp.bethe_energy()
    dec, cov = bp.decimate()
    ee = np.asarray(e)
    assert bool(np.all(cov[ee[:, 0]] | cov[ee[:, 1]])), 'decimated cover is not a cover'
    return dict(smax=smax, regions=bp.n_regions, singletons=bp.n_singletons,
                ext_edges=len(bp.ext), iters=it, resid=float(diff), bethe=float(E),
                cover=dec, secs=time.time() - t)


def one(tau, r, n, seed, smaxes, time_limit, iters=400, damping=0.5):
    e = graph(tau, r, n, seed)
    blocks = np.load(RES / f'blocks_tau{tau}_r{r}_n{n}_s{seed}.npy')
    t = time.time()
    lr, core = leaf_removal(n, e)
    row = dict(tau=tau, r=r, n=n, seed=seed, m=int(len(e)), leaf=lr, core=len(core),
               isolated=int(n - len(np.unique(e))))
    t = time.time()
    res = exact_vc(n, e, time_limit)
    row.update(exact=float(res.fun) if res.x is not None else None,
               exact_status=int(res.status), exact_gap=float(res.mip_gap) if res.x is not None else None,
               exact_secs=time.time() - t)
    row['bp'] = [run_bp(n, e, blocks, s, iters, damping) for s in smaxes]
    return row


def main():
    ns = [int(x) for x in sys.argv[1:]] or [2000, 5000, 10000, 20000]
    rows = json.loads(OUT.read_text()) if OUT.exists() else []
    done = {(x['tau'], x['r'], x['n']) for x in rows}
    for n in ns:
        for r in [0.0, 0.5, 0.8, 1.0]:
            if (2.5, r, n) in done:
                continue
            row = one(2.5, r, n, 0, smaxes=[0, 128], time_limit=1800)
            rows.append(row)
            OUT.write_text(json.dumps(rows, indent=1))
            ex = row['exact']
            print(f"r={r} n={n:6d} m={row['m']} leaf={row['leaf']} core={row['core']} "
                  f"exact={ex} (status {row['exact_status']}, gap {row['exact_gap']}, {row['exact_secs']:.0f}s)", flush=True)
            for b in row['bp']:
                print(f"    smax={b['smax']:4d} regions={b['regions']:6d} singletons={b['singletons']:6d} "
                      f"ext={b['ext_edges']:6d} iters={b['iters']:4d} resid={b['resid']:.1e} "
                      f"bethe={b['bethe']:.2f} cover={b['cover']} [{b['secs']:.0f}s]", flush=True)


if __name__ == '__main__':
    main()
