"""Is the chygraph-BP fixed point Vorob'ev-consistent on shared pairs?

Sec. 23.2 says the messages of the chygraph recursion enforce agreement on
single atoms only, so two complexes sharing a pair of atoms need not agree
on the pair, whereas the region-graph fixed point of Ch. 15 agrees on every
overlap by construction.  This measures the first claim on the instances of
``cavity_assigned.py`` at the stable fixed point: for every pair of cliques
sharing two or more atoms, the largest difference between the two clique
beliefs marginalised to the shared atoms, with the double count and with
each bond in one clique; and, as the control, the agreement on single atoms,
which a fixed point enforces.  Results in ``results/pair_consistency.json``.

    python probe/pair_consistency.py
"""

import json
import sys
import types
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.special import logsumexp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'percolation' / 'src'))
from statmech.bethehessian import (linearised_operator, spectral_radius,  # noqa: E402
                                   trivial_blocks)
from statmech.loopseries import BinaryFactorGraph  # noqa: E402
sys.modules.setdefault('hrg', types.SimpleNamespace(hrg_calibrated=None))
sys.path.append(str(Path(__file__).resolve().parents[2] / 'book' / 'figs'))
sys.path.append(str(Path(__file__).resolve().parent))
from cavity_clique import ChygraphBP  # noqa: E402
from cavity_assigned import COUPLINGS, SEEDS, SIZES, polarised  # noqa: E402
from loopseries import hyperbolic_graph  # noqa: E402

RESULTS = Path(__file__).parent / 'results'
OUT = RESULTS / 'pair_consistency.json'


def beliefs(bp):
    """Normalised belief of every complex, as an array over its members."""
    out = []
    for a, c in enumerate(bp.A):
        acc = bp.logf[a].copy()
        for i, u in enumerate(c):
            sh = [1] * len(c)
            sh[i] = 2
            acc = acc + bp.m_va[(u, a)].reshape(sh)
        out.append(np.exp(acc - logsumexp(acc)))
    return out


def marginal(b, c, keep):
    axes = tuple(i for i, v in enumerate(c) if v not in keep)
    m = b.sum(axis=axes) if axes else b
    # order the surviving axes as `keep`
    order = [v for v in c if v in keep]
    return np.moveaxis(m, range(len(order)), [keep.index(v) for v in order])


def disagreement(bp):
    b = beliefs(bp)
    pair, single = 0.0, 0.0
    for (i, ci), (j, cj) in combinations(enumerate(bp.A), 2):
        shared = sorted(set(ci) & set(cj))
        if not shared:
            continue
        for v in shared:
            single = max(single, float(np.abs(marginal(b[i], ci, [v]) - marginal(b[j], cj, [v])).max()))
        if len(shared) >= 2:
            pair = max(pair, float(np.abs(marginal(b[i], ci, shared) - marginal(b[j], cj, shared)).max()))
    return pair, single


def solve(G, bJ):
    cx = [sorted(c) for c in nx.find_cliques(G)]
    edges = sorted({tuple(sorted(e)) for e in G.edges()})
    out = dict(n=G.number_of_nodes(), beta_J=bJ)
    for assign in ('all', 'once'):
        fg = BinaryFactorGraph.promoted(cx, edges, bJ, assign)
        rho = spectral_radius(linearised_operator(fg, trivial_blocks(fg)))
        if rho < 1:
            bp = ChygraphBP(cx, bJ, damping=0.0, edges=edges, assign=assign).run()
        else:
            bp = polarised(cx, edges, bJ, assign)
        pair, single = disagreement(bp)
        out[assign] = dict(pair=pair, single=single, residual=bp.residual, rho=rho)
    return out


def main():
    from gbp_real import load
    from merge import karrer_graph
    rows = []
    R = {}
    for n, kbar in SIZES:
        for seed in SEEDS:
            G, R[n] = hyperbolic_graph(n, kbar, 2.5, np.random.default_rng(seed), R.get(n))
            G.remove_nodes_from([v for v in list(G) if G.degree(v) == 0])
            G = nx.convert_node_labels_to_integers(G)
            for bJ in COUPLINGS:
                rows.append(dict(ensemble='hyperbolic', **solve(G, bJ)))
    seen = set()
    for r in json.load(open(RESULTS / 'gbp_karrer.json'))['runs']:
        if (r['n'], r['seed']) in seen:
            continue
        seen.add((r['n'], r['seed']))
        G, _ = karrer_graph(r['n'], r['s_mean'], {3: r['triangles']}, r['seed'])
        G.remove_edges_from(nx.selfloop_edges(G))
        for bJ in COUPLINGS:
            rows.append(dict(ensemble='karrer', **solve(G, bJ)))
    cache, seen = {}, set()
    for r in json.load(open(RESULTS / 'gbp_real.json')):
        if (r['network'], r['ego']) in seen:
            continue
        seen.add((r['network'], r['ego']))
        G = cache.setdefault(r['key'], load(r['key']))
        ego = sorted([r['ego']] + list(G[r['ego']]), key=str)
        H = nx.convert_node_labels_to_integers(G.subgraph(ego))
        for bJ in COUPLINGS:
            rows.append(dict(ensemble='real', network=r['network'], **solve(H, bJ)))
    json.dump(rows, OUT.open('w'), indent=1)
    for assign in ('all', 'once'):
        p = np.array([r[assign]['pair'] for r in rows])
        s = np.array([r[assign]['single'] for r in rows])
        print(f"{assign}: single-atom disagreement max {s.max():.1e}; pair disagreement "
              f"min {p.min():.2e} median {np.median(p):.3f} max {p.max():.3f}; "
              f"above 1e-6 on {int(np.sum(p > 1e-6))}/{len(rows)}, above 0.01 on {int(np.sum(p > 0.01))}")
        for ens in ('hyperbolic', 'karrer', 'real'):
            for bJ in COUPLINGS:
                sel = [r[assign]['pair'] for r in rows if r['ensemble'] == ens and r['beta_J'] == bJ]
                print(f"   {ens:<11} bJ={bJ}: pair disagreement median {np.median(sel):.3f} max {max(sel):.3f}")
    print('wrote', OUT)


if __name__ == '__main__':
    main()
