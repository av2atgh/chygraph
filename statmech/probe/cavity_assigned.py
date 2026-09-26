"""Chapter 14's measurement repeated with each bond assigned to one clique.

Sec. 24.3 found that most of the chygraph recursion's error on overlapping
cliques is the bond counted twice, and that giving each bond to one of the
cliques containing it -- the factor graph then represents the pairwise model
exactly, and the recursion is ordinary belief propagation on it -- recovers
most of the error for nothing.  This runs that repair on the instances of
Chapter 14: the same Karrer--Newman graphs and the same 120 real
neighbourhoods (same cache, same seeds), and fresh draws of the hyperbolic
ensemble at the same n, kbar and tau, since the generator behind the original
thirty is no longer on disk.  For every instance and coupling it records the
recursion's error with the double count and with the assignment, from the
same solver (`cavity_clique.ChygraphBP`, `assign='all'` and `'once'`), and
the exact ln Z by enumeration.  Results in ``results/cavity_assigned.json``;
``book/figs/overlap.py`` reads them.

    python probe/cavity_assigned.py
"""

import json
import sys
import time
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'percolation' / 'src'))
from statmech.gbp import exact_log_Z, ising_factors  # noqa: E402
# the figures directory holds a script called statmech.py, so it goes on the
# path only after the package is imported
sys.path.append(str(Path(__file__).resolve().parents[2] / 'book' / 'figs'))
sys.path.append(str(Path(__file__).resolve().parent))
from cavity_clique import DAMPING, TOL, ChygraphBP  # noqa: E402
from loopseries import hyperbolic_graph  # noqa: E402

RESULTS = Path(__file__).parent / 'results'
OUT = RESULTS / 'cavity_assigned.json'
COUPLINGS = (0.3, 0.8)
SIZES = ((14, 3.5), (18, 5.0), (20, 7.0))
SEEDS = range(1, 11)


def paramagnetic(cx, edges, bJ, assign):
    """The symmetric fixed point, which the symmetric iteration keeps, and
    the Perron root of the recursion linearised there (stable below one)."""
    from statmech.bethehessian import (linearised_operator, spectral_radius,
                                       trivial_blocks)
    from statmech.loopseries import BinaryFactorGraph
    bp = ChygraphBP(cx, bJ, damping=0.0, edges=edges, assign=assign).run()
    fg = BinaryFactorGraph.promoted(cx, edges, bJ, assign)
    rho = spectral_radius(linearised_operator(fg, trivial_blocks(fg)))
    return bp, rho


def polarised(cx, edges, bJ, assign, init=0.5):
    """The fixed point reached from a polarised start, over the damping
    ladder from one half; the first to settle is kept."""
    best = None
    for d in DAMPING[1:]:
        bp = ChygraphBP(cx, bJ, damping=d, edges=edges, assign=assign)
        for k in bp.m_va:
            bp.m_va[k] = np.array([init, -init])
        bp.run(20000)
        if best is None or bp.residual < best.residual:
            best = bp
        if bp.residual < 1e-11:
            break
    return best


def magnetisation(bp):
    from scipy.special import logsumexp
    out = []
    for v in bp.nodes:
        acc = np.zeros(2)
        for a, c in enumerate(bp.A):
            if v in c:
                acc = acc + bp.m_av[(a, v)]
        p = np.exp(acc - logsumexp(acc))
        out.append(float(p[0] - p[1]))
    return float(np.median(np.abs(out)))


def solve(G, bJ):
    n = G.number_of_nodes()
    cx = [sorted(c) for c in nx.find_cliques(G)]
    edges = sorted({tuple(sorted(e)) for e in G.edges()})
    exact = exact_log_Z(ising_factors(edges, bJ), range(n))
    cov = {}
    for c in cx:
        for e in combinations(sorted(c), 2):
            cov[e] = cov.get(e, 0) + 1
    out = dict(n=n, beta_J=bJ, n_cliques=len(cx), n_bonds=len(edges),
               doubled_bonds=sum(1 for k in cov.values() if k > 1), exact=exact)
    for assign in ('all', 'once'):
        key = 'cavity' if assign == 'all' else 'cavity_once'
        bp, rho = paramagnetic(cx, edges, bJ, assign)
        out[key + '_para'] = bp.log_Z() - exact
        out[key + '_rho'] = rho
        bp = polarised(cx, edges, bJ, assign)
        out[key + '_pol'] = bp.log_Z() - exact
        out[key + '_pol_m'] = magnetisation(bp)
        out[key + '_pol_residual'] = bp.residual
        # the fixed point to report: the symmetric one while it is stable,
        # the polarised one once it is not
        out[key] = out[key + '_para'] if rho < 1 else out[key + '_pol']
        out[key + '_stable_para'] = bool(rho < 1)
    return out


def main():
    import types
    sys.modules.setdefault('hrg', types.SimpleNamespace(hrg_calibrated=None))
    from gbp_real import load
    from merge import karrer_graph

    rows = []
    t0 = time.time()
    R = {}
    for n, kbar in SIZES:
        for seed in SEEDS:
            G, R[n] = hyperbolic_graph(n, kbar, 2.5, np.random.default_rng(seed), R.get(n))
            G.remove_nodes_from([v for v in list(G) if G.degree(v) == 0])
            G = nx.convert_node_labels_to_integers(G)
            for bJ in COUPLINGS:
                rows.append(dict(ensemble='hyperbolic', seed=seed, **solve(G, bJ)))
    print(f'hyperbolic done, {time.time() - t0:.0f} s', flush=True)

    seen = set()
    for r in json.load(open(RESULTS / 'gbp_karrer.json'))['runs']:
        if (r['n'], r['seed']) in seen:
            continue
        seen.add((r['n'], r['seed']))
        G, _ = karrer_graph(r['n'], r['s_mean'], {3: r['triangles']}, r['seed'])
        G.remove_edges_from(nx.selfloop_edges(G))
        for bJ in COUPLINGS:
            rows.append(dict(ensemble='karrer', seed=r['seed'], **solve(G, bJ)))
    print(f'karrer done, {time.time() - t0:.0f} s', flush=True)

    cache, seen = {}, set()
    for r in json.load(open(RESULTS / 'gbp_real.json')):
        if (r['network'], r['ego']) in seen:
            continue
        seen.add((r['network'], r['ego']))
        G = cache.setdefault(r['key'], load(r['key']))
        ego = sorted([r['ego']] + list(G[r['ego']]), key=str)
        H = nx.convert_node_labels_to_integers(G.subgraph(ego))
        for bJ in COUPLINGS:
            rows.append(dict(ensemble='real', network=r['network'],
                             ego=str(r['ego']), **solve(H, bJ)))
    print(f'real done, {time.time() - t0:.0f} s', flush=True)

    json.dump(rows, OUT.open('w'), indent=1)
    summarise(rows)
    print('\nwrote', OUT)


def summarise(rows):
    for ens in ('hyperbolic', 'karrer', 'real'):
        for bJ in COUPLINGS:
            sel = [r for r in rows if r['ensemble'] == ens and r['beta_J'] == bJ]
            unstable = sum(not r['cavity_stable_para'] for r in sel)
            unstable_o = sum(not r['cavity_once_stable_para'] for r in sel)
            conv = sum(r['cavity_once_pol_residual'] < 1e-9 and
                       r['cavity_pol_residual'] < 1e-9 for r in sel)
            para = np.array([abs(r['cavity_para']) for r in sel])
            e = np.array([abs(r['cavity']) for r in sel])
            o = np.array([abs(r['cavity_once']) for r in sel])
            print(f'{ens:<11} bJ={bJ}: {len(sel)} runs; paramagnetic point unstable on '
                  f'{unstable} (double count) / {unstable_o} (assigned); polarised runs settled on {conv};\n'
                  f'    double count: paramagnetic median {np.median(para):.3f}, stable-point median {np.median(e):.3f} (max {e.max():.1f});\n'
                  f'    assigned: stable-point median {np.median(o):.3f} (max {o.max():.2f}), within 0.05 of ln 2 on '
                  f'{int(np.sum(np.abs(o - np.log(2)) < 0.05))}, below ln 2 + 0.05 on {int(np.sum(o < np.log(2) + 0.05))}, '
                  f'smaller than the double count on {int(np.sum(o < e))}')


if __name__ == '__main__':
    main()
