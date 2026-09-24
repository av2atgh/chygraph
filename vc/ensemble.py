"""Graphs from the degree-correlated ensemble of Sec. 11.8 (VW03 Eq. 18).

Thin wrapper over statmech/probe/vw_clustering.{degrees,wire}: same degrees
(p_d ~ d^-tau, cutoff n^(1/tau)), same modified Molloy-Reed wiring with
assortativity r, reduced to a simple graph.  Also a plain leaf removal, since
the ~/computational_complexity module the probes import is no longer present.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'statmech' / 'probe'))
sys.path.insert(0, str(HERE.parent / 'statmech' / 'src'))
from vw_clustering import degrees, wire  # noqa: E402


def graph(tau, r, n, seed):
    """Simple edge list (m, 2) of one draw; nodes 0..n-1."""
    rng = np.random.default_rng(hash((tau, n, r, seed)) % 2**32)
    d = degrees(tau, n, rng)
    e = wire(d, r, rng)
    e = e[e[:, 0] != e[:, 1]]
    e = np.sort(e, axis=1)
    return np.unique(e, axis=0)


def adjacency(n, edges):
    adj = [[] for _ in range(n)]
    for a, b in edges.tolist():
        adj[a].append(b)
        adj[b].append(a)
    return adj


def leaf_removal(n, edges):
    """Bauer-Golinelli leaf removal.  Returns (cover size, core node set).

    A leaf's neighbour goes into the cover, both are deleted.  What is left
    with no leaves is the core; the cover is a minimum one iff the core is
    empty (isolated vertices excluded).
    """
    adj = [set(a) for a in adjacency(n, edges)]
    deg = np.array([len(a) for a in adj])
    alive = np.ones(n, bool)
    cover = 0
    stack = [v for v in range(n) if deg[v] == 1]
    while stack:
        v = stack.pop()
        if not alive[v] or deg[v] != 1:
            continue
        (u,) = [w for w in adj[v] if alive[w]]
        cover += 1                                  # u into the cover
        for w in (u, v):
            alive[w] = False
        for w in adj[u]:
            if alive[w]:
                deg[w] -= 1
                if deg[w] == 1:
                    stack.append(w)
                elif deg[w] == 0:
                    pass
        deg[u] = deg[v] = 0
    core = [v for v in range(n) if alive[v] and deg[v] > 0]
    return cover, core
