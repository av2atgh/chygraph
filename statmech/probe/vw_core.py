"""Leaf-removal core on graphs drawn from VW03 Eq. (18), against r.

Section 11.8 locates the transition analytically: the fixed point of Eq. (11.11)
loses stability at r_c, and Sec. 11.4 says that point has three names, one of
which is the appearance of the leaf-removal core.  That identification was
established on the chy-degree axis.  Here it is asked along the *correlation*
axis instead, where nothing about it is automatic -- the map being destabilised
is resolved by degree class, not by layer, and no chygraph is involved.

So: draw graphs from the ensemble by VW03's own generator (probe/vw_clustering
wire()), run pure leaf removal on the simple graph, and report

    P_C(VC, G) = core fraction

the same quantity Sec. 11.9's first column reports for real networks.  The
generalised-leaf-removal error bound dx of [VW03] is recorded alongside, since
it is their own diagnostic and costs one extra pass.

Two structural coordinates of the same graph are recorded with it, transitivity
and the degree assortativity coefficient, so the transition can be located
against what the ensemble *is* and not only against the parameter r that was
turned.  Neither is a control variable: both are outcomes of r, and
probe/vw_clustering shows transitivity is the one of the two that does not
vanish with n.

Caching: probe/results/vw_core.json, read by book/figs/cover.py.
"""

import json
import sys
from pathlib import Path

import networkx as nx
import numpy as np

CC = Path.home() / 'av2atg' / 'computational_complexity' / 'code'
sys.path.insert(0, str(CC))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import leafremoval as lr                                        # noqa: E402
from hrg import to_csr                                          # noqa: E402
from vw_clustering import degrees, wire                         # noqa: E402

OUT = Path(__file__).resolve().parent / 'results' / 'vw_core.json'

TAUS = (2.5, 3.0)
RS = tuple(np.round(np.linspace(0.0, 1.0, 21), 3))
N = 200000
SEEDS = (0, 1, 2)


def simple(edges, n):
    """Drop self-loops and repeated pairs, then CSR."""
    e = edges[edges[:, 0] != edges[:, 1]]
    e = np.sort(e, axis=1)
    e = np.unique(e, axis=0)
    return e, to_csr(n, e[:, 0], e[:, 1])


def one(tau, r, n, seed):
    rng = np.random.default_rng(hash((tau, n, r, seed)) % 2**32)
    d = degrees(tau, n, rng)
    e, (indptr, indices) = simple(wire(d, r, rng), n)
    core, _ = lr.core(indptr, indices)
    g = lr.generalized(indptr, indices)
    G = nx.Graph()                       # for the two structural coordinates
    G.add_nodes_from(range(n))
    G.add_edges_from(map(tuple, e.tolist()))
    return dict(core=core / n, kbar=2 * len(e) / n, dx=g['dx'], xc=g['xc'],
                trans=nx.transitivity(G),
                assort=nx.degree_assortativity_coefficient(G))


def run():
    rows = []
    for tau in TAUS:
        print(f'\ntau = {tau}, n = {N}, {len(SEEDS)} seeds')
        print(f"{'r':>7}{'kbar':>8}{'P_C(VC,G)':>12}{'dx':>10}{'x_c':>9}"
              f"{'T':>10}{'a':>9}")
        for r in RS:
            out = [one(tau, float(r), N, s) for s in SEEDS]
            av = {k: float(np.mean([o[k] for o in out])) for k in out[0]}
            sd = float(np.std([o['core'] for o in out]))
            tsd = float(np.std([o['trans'] for o in out]))
            rows.append(dict(tau=tau, r=float(r), n=N, seeds=len(SEEDS),
                             core_sd=sd, trans_sd=tsd, **av))
            print(f"{r:>7.2f}{av['kbar']:>8.3f}{av['core']:>12.4f}"
                  f"{av['dx']:>10.5f}{av['xc']:>9.4f}{av['trans']:>10.5f}"
                  f"{av['assort']:>9.3f}")
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(rows, indent=1) + '\n')
    print(f'\nwrote {OUT}')


if __name__ == '__main__':
    run()
