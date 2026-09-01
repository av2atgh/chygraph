"""What ignoring clique overlap costs (WP5), measured on the HRG.

`core.py` predicts the leaf-removal core of a chygraph whose complexes are
treated as independent -- meeting in at most one node.  Real maximal cliques do
not: in a clustered graph they share edges, and `region.py` shows the Bethe
counting the chygraph uses is then wrong.

This measures the size of that error, on the graphs the whole programme is
about.  For each hyperbolic random graph:

  * enumerate the maximal cliques and record, per node, how many of each
    cardinality it belongs to -- the chygraph ensemble the graph actually has;
  * run `CorePercolation.from_samples` on that ensemble, which is the chygraph
    prediction with overlap ignored;
  * run pure leaf removal on the graph itself for the truth;
  * report `overlap_profile`, which says how far from treelike the clique
    family is.

The degree-matched configuration model is carried alongside as the control that
prediction 4 is about.
"""

import sys
from collections import Counter
from pathlib import Path

import networkx as nx
import numpy as np

CC = Path.home() / 'av2atg' / 'computational_complexity' / 'code'
sys.path.insert(0, str(CC))
import leafremoval as lr  # noqa: E402
from hrg import erased_configuration_model, hrg_calibrated, to_csr  # noqa: E402

import chygraph_statmech.core as cp  # noqa: E402
from chygraph_statmech.region import overlap_profile  # noqa: E402

CMAX = 8          # cardinalities 2..CMAX get their own layer; larger are pooled


def analyse(n, src, dst):
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(zip(src.tolist(), dst.tolist()))
    cliques = [frozenset(c) for c in nx.find_cliques(G) if len(c) >= 2]

    cards = list(range(2, CMAX + 1))
    idx = {c: i for i, c in enumerate(cards)}
    K = np.zeros((n, len(cards)))
    for a in cliques:
        j = idx.get(len(a), len(cards) - 1)
        for v in a:
            K[v, j] += 1
    keep = [i for i in range(len(cards)) if K[:, i].mean() > 0]
    K, cards = K[:, keep], [cards[i] for i in keep]

    predicted = cp.CorePercolation.from_samples(cards, K).core_fraction()
    measured = lr.core(*to_csr(n, src, dst))[0] / n
    prof = overlap_profile(cliques[:4000])       # profile a sample: O(pairs)
    return predicted, measured, prof, Counter(len(a) for a in cliques)


if __name__ == '__main__':
    n = 30_000
    print(f"n = {n}, cardinality layers 2..{CMAX} (larger pooled)\n")
    print(f"{'family':>7}{'tau':>5}{'kbar':>6}{'chygraph':>10}{'measured':>10}"
          f"{'ratio':>8}{'2+ shared':>11}{'edge cover':>12}")
    for tau in (2.5, 2.9):
        for kbar in (1.0, 2.0, 4.0):
            src, dst, r, th, R, C = hrg_calibrated(n, tau=tau, kbar=kbar,
                                                   rng=1, tol=0.01, max_iter=25)
            deg = np.bincount(np.concatenate((src, dst)), minlength=n)
            for fam, e in (('hrg', (src, dst)),
                           ('config', erased_configuration_model(deg, rng=1))):
                pred, meas, prof, _ = analyse(n, *e)
                ratio = pred / meas if meas > 1e-9 else float('inf')
                print(f"{fam:>7}{tau:>5.1f}{2*e[0].size/n:>6.2f}{pred:>10.4f}"
                      f"{meas:>10.4f}{ratio:>8.2f}{prof['shared_2plus']:>11.3f}"
                      f"{prof['edge_cover_mean']:>12.3f}", flush=True)
