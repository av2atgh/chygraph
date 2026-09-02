"""Does the assortativity knob of VW03 Eq. (18) also turn on clustering?

Section 11.8 draws its transition along `r` in

    e_dd' = q_d [ r delta_dd' + (1-r) q_d' ],

and the argument it feeds -- that the hyperbolic graph and its degree-matched
control sit at the same place on the {p_d, e_dd'} surface and are still in
different phases -- assumes the axis is a *pure* correlation axis.  If turning
`r` up also manufactured triangles, the section would be comparing a clustered
ensemble with an unclustered one and calling the difference correlation.

So measure it.  Graphs are drawn by VW03's own generator, a modified
Molloy-Reed:

  (i)   give node i a degree d_i ~ p_d and put d_i stubs in a list L;
  (ii)  pick a stub i in L at random;
  (iii) with probability r pick a second stub j in L with d_j = d_i, otherwise
        pick j in L at random;
  (iv)  join i and j, delete both.  Repeat until L is empty.

and then reported on the *simple* graph -- self-loops and repeated pairs
dropped -- because that is the object a clustering coefficient is defined on
and the object leaf removal would be run on.

Two things are asked of each ensemble.  Whether C rises with r at fixed n, and
whether whatever C it has survives n.  A configuration model has C ~ <k>/n, an
O(1/n) artefact of finite size and not clustering; the question is whether the
r = 1 corner is any different in kind.

Reported per (tau, r, n): mean degree of the simple graph, the degree
assortativity coefficient (the knob, measured), average local clustering
(Chapter 3's C), transitivity, and the fraction of stub pairings lost to
self-loops and multi-edges, since that fraction is what the r = 1 corner does
to the degree sequence it was asked to realise.

No caching: the whole scan is a couple of minutes.
"""

import sys
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'percolation' / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))

from statmech import vertexcover as vc                          # noqa: E402

TAUS = (2.5, 3.0)
RS = (0.0, 0.2, 0.4, 0.6, 0.7042, 0.8, 0.9, 1.0)
SIZES = (20000, 50000, 200000)
SEEDS = (0, 1)


# ---------------------------------------------------------------- the ensemble
def degrees(tau, n, rng):
    """d_i ~ d^-tau on 1..dmax with the structural cutoff dmax ~ n^(1/tau)."""
    dmax = max(2, int(n ** (1.0 / tau)))
    p = vc.scale_free(tau, dmax)
    d = rng.choice(len(p), size=n, p=p)
    if d.sum() % 2:                       # one spare stub, give it away
        d[rng.integers(n)] += 1
    return d


def wire(d, r, rng):
    """VW03's modified Molloy-Reed.  Returns the edge list, multi-edges kept.

    Two swap-and-pop pools are kept over the live stubs, one global and one per
    degree, so a step is O(1): drawing a partner of equal degree is a draw from
    that degree's own pool rather than a rejection loop, which matters at r = 1
    where rejection would almost never hit.
    """
    n = len(d)
    stub = np.repeat(np.arange(n), d)          # stub -> node
    m = len(stub)
    sdeg = d[stub]                             # stub -> its node's degree

    pool = np.arange(m)                        # global pool of live stubs
    at = np.arange(m)                          # stub -> position in pool
    npool = m

    order = np.argsort(sdeg, kind='stable')    # per-degree pools, packed
    start = np.searchsorted(sdeg[order], np.arange(sdeg.max() + 2))
    bpool = order.copy()
    bat = np.empty(m, dtype=np.int64)
    bat[bpool] = np.arange(m)
    nb = np.diff(start)                        # live count per degree

    def drop(s):
        nonlocal npool
        i, last = at[s], pool[npool - 1]       # from the global pool
        pool[i], at[last] = last, i
        npool -= 1
        g = sdeg[s]
        j = bat[s] - start[g]
        last = bpool[start[g] + nb[g] - 1]     # and from its degree's pool
        bpool[start[g] + j], bat[last] = last, start[g] + j
        nb[g] -= 1

    edges = np.empty((m // 2, 2), dtype=np.int64)
    assort = rng.random(m // 2) < r
    for e in range(m // 2):
        a = pool[rng.integers(npool)]
        drop(a)
        g = sdeg[a]
        if assort[e] and nb[g] > 0:
            b = bpool[start[g] + rng.integers(nb[g])]
        else:
            b = pool[rng.integers(npool)]
        drop(b)
        edges[e] = (stub[a], stub[b])
    return edges


def measure(edges, n):
    """Simple-graph statistics, plus the fraction of pairings that did not land."""
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from((int(a), int(b)) for a, b in edges if a != b)
    kept = G.number_of_edges()
    return dict(kbar=2 * kept / n,
                assort=nx.degree_assortativity_coefficient(G),
                C=nx.average_clustering(G),
                trans=nx.transitivity(G),
                lost=1 - kept / len(edges))


def run():
    print('VW03 Eq. (18): does r manufacture clustering?')
    print('  kbar  mean degree of the simple graph')
    print('  a     degree assortativity coefficient (the knob, measured)')
    print('  C     average local clustering, Chapter 3\'s C')
    print('  T     transitivity, 3 x triangles / triples')
    print('  lost  fraction of stub pairings dropped as self-loop or repeat')
    for tau in TAUS:
        for n in SIZES:
            print(f'\ntau = {tau}, n = {n}, {len(SEEDS)} seeds')
            print(f"{'r':>7}{'kbar':>8}{'a':>9}{'C':>10}{'T':>10}{'lost':>8}"
                  f"{'C n/kbar':>10}")
            for r in RS:
                out = []
                for s in SEEDS:
                    rng = np.random.default_rng((hash((tau, n, r, s)) % 2**32))
                    d = degrees(tau, n, rng)
                    out.append(measure(wire(d, r, rng), n))
                av = {k: float(np.mean([o[k] for o in out])) for k in out[0]}
                print(f"{r:>7.4f}{av['kbar']:>8.3f}{av['assort']:>9.3f}"
                      f"{av['C']:>10.5f}{av['trans']:>10.5f}{av['lost']:>8.3f}"
                      f"{av['C'] * n / av['kbar']:>10.2f}")
    print('\nC n/kbar is the finite-size scale: constant in n means C ~ kbar/n,')
    print('the configuration-model artefact, and growing means something else.')


if __name__ == '__main__':
    run()
