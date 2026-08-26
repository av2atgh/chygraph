"""Erdos-Renyi control for the clique-moment probe.

ER has finite ``<k^2>`` and no geometry, so if a chygraph complex ensemble of
maximal cliques is ever well defined, it is here.  It is: ``sbar`` is flat to
three or four decimals across a 30x range in ``n``, and the clique number does
not move off 3.  This is the baseline against which the HRG and its
degree-matched configuration model are read in :mod:`clique_moments`.
"""

import sys
from collections import Counter
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path.home() / 'av2atg' / 'computational_complexity' / 'code'))
from hrg import erdos_renyi  # noqa: E402


def moments(n, src, dst):
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(zip(src.tolist(), dst.tolist()))
    cnt = Counter(len(c) for c in nx.find_cliques(G))
    c = np.array(sorted(cnt), float)
    w = np.array([cnt[int(i)] for i in c], float)
    w /= w.sum()
    m1, m2 = (w * c).sum(), (w * c * c).sum()
    return int(c.max()), float(m1), float(m2), float(m2 / m1 - 1)


if __name__ == '__main__':
    print("Erdos-Renyi: finite <k^2>, no geometry, no heavy tail")
    print(f"{'kbar':>6}{'n':>9}{'c_max':>7}{'<c>':>8}{'m2':>9}{'sbar':>8}")
    for kbar in (2.0, 4.0, 8.0):
        for n in (10_000, 30_000, 100_000, 300_000):
            cm, m1, m2, sbar = moments(n, *erdos_renyi(n, kbar, rng=1))
            print(f"{kbar:>6.0f}{n:>9}{cm:>7}{m1:>8.3f}{m2:>9.3f}{sbar:>8.3f}")
