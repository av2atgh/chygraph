"""Complexes as regions: the Kikuchi reading of a chygraph (WP5).

Everything up to here treats a node's complexes as independent.  That is the
Bethe approximation on the chygraph, and it is exact when the chygraph is
locally treelike -- which, at the node level, means **complexes meet in at most
one node**.  Two cliques sharing an edge are a loop the chygraph cannot see,
just as a triangle is a loop an ordinary graph cannot see.  Chygraphs move the
problem up one level; they do not remove it.

The region-graph construction is where it goes next.  Take the complexes as
regions, close the family under intersection, and give each region a counting
number by Mobius inversion

    c_R = 1 - sum over regions strictly containing R of c_{R'}

so that every node is counted exactly once, ``sum_{R containing v} c_R = 1``.
The Kikuchi free energy is then the region-weighted sum of exact free energies
inside each region.

**The reduction that makes this the right generalisation.** When complexes meet
in at most one node the family closes on the complexes plus the single nodes,
and a node in ``k`` complexes gets ``c_v = 1 - k``.  That is exactly the
counting of the Bethe free energy on a factor graph -- one term per complex,
``(1 - k_v)`` per node -- so the Kikuchi construction contains everything WP1
and WP4 compute, and differs from them only where complexes genuinely overlap.
:func:`RegionGraph.is_bethe` decides which case a given chygraph is in.

What this module does *not* do is run generalised belief propagation on the
region graph.  It builds the region graph, verifies the counting, and measures
how far a given chygraph is from the treelike case -- which is what decides
whether the independent-complex answer can be trusted at all.
"""

from collections import defaultdict
from itertools import combinations

import numpy as np


class RegionGraph:
    """Regions and Mobius counting numbers for a set of complexes.

    Args:
        complexes: iterable of iterables of node labels.
        max_rounds: how many times to close the family under intersection.
            Each round adds the intersections of the regions found so far;
            two rounds suffice unless complexes overlap in elaborate ways.
    """

    def __init__(self, complexes, max_rounds=2):
        self.complexes = [frozenset(a) for a in complexes]
        self.regions = self._close(self.complexes, max_rounds)
        self.counting = self._mobius(self.regions)

    # -- construction -------------------------------------------------------

    @staticmethod
    def _close(complexes, max_rounds):
        """Close the family under pairwise intersection, then add singletons."""
        regions = set(a for a in complexes if a)
        for _ in range(max_rounds):
            new = set()
            for a, b in combinations(regions, 2):
                s = a & b
                if s and s not in regions:
                    new.add(s)
            if not new:
                break
            regions |= new
        for a in complexes:
            for v in a:
                regions.add(frozenset({v}))
        return sorted(regions, key=lambda r: (-len(r), sorted(r)))

    @staticmethod
    def _mobius(regions):
        """``c_R = 1 - sum_{R' strictly containing R} c_{R'}``, largest first."""
        c = {}
        for r in regions:                      # already sorted by size, desc
            c[r] = 1 - sum(c[s] for s in c if r < s)
        return c

    # -- checks -------------------------------------------------------------

    def node_counts(self):
        """``sum_{R containing v} c_R`` per node; must be 1 everywhere."""
        out = defaultdict(int)
        for r, cr in self.counting.items():
            for v in r:
                out[v] += cr
        return dict(out)

    def counting_is_valid(self):
        return all(v == 1 for v in self.node_counts().values())

    def is_bethe(self):
        """True when complexes meet in at most one node.

        Then the region family is the complexes plus the single nodes, the
        counting numbers are ``1`` and ``1 - k_v``, and the Kikuchi free energy
        *is* the Bethe free energy the rest of this package computes.
        """
        return all(len(a & b) <= 1
                   for a, b in combinations(self.complexes, 2))

    def bethe_counting(self):
        """The counting numbers the Bethe treatment assumes: ``1 - k_v``."""
        k = defaultdict(int)
        for a in self.complexes:
            for v in a:
                k[v] += 1
        out = {a: 1 for a in self.complexes}
        for v, kv in k.items():
            out[frozenset({v})] = 1 - kv
        return out

    def bethe_error(self):
        """Total absolute discrepancy between Kikuchi and Bethe counting.

        Zero exactly when :meth:`is_bethe`; grows with the amount of
        higher-order overlap the Bethe treatment miscounts.
        """
        b = self.bethe_counting()
        keys = set(self.counting) | set(b)
        return sum(abs(self.counting.get(k, 0) - b.get(k, 0)) for k in keys)


# ---------------------------------------------------------------------------
# How far from treelike is a given graph's clique structure?
# ---------------------------------------------------------------------------

def overlap_profile(complexes):
    """Statistics of how the complexes intersect.

    ``shared_2plus`` is the fraction of intersecting complex pairs that share
    two or more nodes -- the pairs the chygraph mapping cannot represent.
    ``edge_cover_excess`` is how many times the average edge is covered by a
    complex; 1.0 means every edge sits in exactly one complex, which is the
    treelike case.
    """
    cs = [frozenset(a) for a in complexes]
    by_node = defaultdict(list)
    for i, a in enumerate(cs):
        for v in a:
            by_node[v].append(i)

    pairs, deep = set(), 0
    for members in by_node.values():
        for i, j in combinations(sorted(set(members)), 2):
            pairs.add((i, j))
    for i, j in pairs:
        if len(cs[i] & cs[j]) >= 2:
            deep += 1

    edge_cov = defaultdict(int)
    for a in cs:
        for u, v in combinations(sorted(a), 2):
            edge_cov[(u, v)] += 1
    cov = np.array(list(edge_cov.values()), dtype=float) if edge_cov else np.array([1.0])

    return {
        'n_complexes': len(cs),
        'n_intersecting_pairs': len(pairs),
        'shared_2plus': (deep / len(pairs)) if pairs else 0.0,
        'edge_cover_mean': float(cov.mean()),
        'edge_cover_max': int(cov.max()),
        'treelike': deep == 0,
    }
