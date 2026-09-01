"""Check the chygraph core map against pure leaf removal on real graphs.

Builds the graph the chygraph describes -- a bipartite configuration model of
nodes and cliques, expanded to a simple graph -- and runs the O(n+m) leaf
removal of ``~/av2atg/computational_complexity``.  Cardinality 2 is Erdos-Renyi
and checks the map against Bauer & Golinelli.
"""

import sys
from pathlib import Path

import numpy as np

CC = Path.home() / 'av2atg' / 'computational_complexity' / 'code'
sys.path.insert(0, str(CC))
import leafremoval as lr  # noqa: E402
from hrg import to_csr  # noqa: E402

import statmech.core as cp  # noqa: E402


def clique_graph(n, c, k, rng=1):
    """n nodes, Poisson(k) clique memberships each, cliques of cardinality c."""
    rng = np.random.default_rng(rng)
    stubs = np.repeat(np.arange(n), rng.poisson(k, n))
    rng.shuffle(stubs)
    stubs = stubs[:(stubs.size // c) * c].reshape(-1, c)
    src = np.concatenate([stubs[:, a] for a in range(c)
                          for b in range(a + 1, c)])
    dst = np.concatenate([stubs[:, b] for a in range(c)
                          for b in range(a + 1, c)])
    keep = src != dst
    lo = np.minimum(src[keep], dst[keep])
    hi = np.maximum(src[keep], dst[keep])
    _, u = np.unique(lo * n + hi, return_index=True)
    return lo[u], hi[u]


if __name__ == '__main__':
    n = 400_000
    print(f"n = {n}")
    print(f"{'c':>3}{'k':>7}{'kbar':>8}{'theory':>11}{'simulation':>12}{'err':>10}")
    for c, ks in ((2, (2.0, 3.0, 4.0, 6.0)),
                  (3, (0.1, 0.3, 1.0, 2.0)),
                  (4, (0.1, 0.5, 1.5))):
        for k in ks:
            s, d = clique_graph(n, c, k)
            sim = lr.core(*to_csr(n, s, d))[0] / n
            th = cp.clique_network(c, k).core_fraction()
            print(f"{c:>3}{k:>7.2f}{2*s.size/n:>8.3f}{th:>11.6f}"
                  f"{sim:>12.6f}{abs(th-sim):>10.2e}")
