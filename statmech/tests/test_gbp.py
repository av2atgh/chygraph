"""Generalised belief propagation on the region graph.

WP5 built the region graph and measured the counting; this is the algorithm the
manuscript's conclusions name as the standing item.  Three things have to hold
before any number it produces is worth quoting.

1. It must *contain* belief propagation.  On the two-layer region graph of a
   pairwise model --- one region per edge plus one per node --- the
   parent-to-child update has an empty denominator and reduces to ordinary BP,
   so GBP must reproduce a textbook BP implementation exactly, loopy error and
   all.
2. It must be exact where the region graph is a junction tree.  Two triangles
   sharing an edge is the manuscript's own Table IV example and is one.
3. Where it is not exact its fixed point must still *be* a fixed point: the
   parent-to-child update is derived from ``sum_{x_P \\ x_R} b_P = b_R``, so a
   converged run that violates that condition is a bug in the solver rather
   than a limit of the approximation.
"""

import itertools

import numpy as np
import pytest
from scipy.special import logsumexp

from statmech.gbp import (GBP, clique_edges, exact_log_Z,
                                   ising_factors, static_log_Z)
from statmech.region import RegionGraph


# ---------------------------------------------------------------------------
# an independent BP, written without reference to the region-graph code
# ---------------------------------------------------------------------------

def _bp_log_Z(edges, beta_J, nodes, iters=40_000, damping=0.5, tol=1e-14):
    """Textbook belief propagation on a pairwise Ising model."""
    nb = {v: [] for v in nodes}
    for u, v in edges:
        nb[u].append(v)
        nb[v].append(u)
    s = np.array([1.0, -1.0])
    K = beta_J * np.outer(s, s)
    m = {}
    for u, v in edges:
        m[(u, v)] = np.zeros(2)
        m[(v, u)] = np.zeros(2)

    def incoming(i, skip):
        acc = np.zeros(2)
        for k in nb[i]:
            if k != skip:
                acc = acc + m[(k, i)]
        return acc

    for _ in range(iters):
        worst = 0.0
        for (i, j) in list(m):
            new = logsumexp(K + incoming(i, j)[:, None], axis=0)
            new -= logsumexp(new)
            mix = damping * m[(i, j)] + (1 - damping) * new
            mix -= logsumexp(mix)
            worst = max(worst, float(np.abs(mix - m[(i, j)]).max()))
            m[(i, j)] = mix
        if worst < tol:
            break

    lnZ = 0.0
    for v in nodes:                                   # node terms, weight 1-k
        b = incoming(v, None)
        b = b - logsumexp(b)
        p = np.exp(b)
        lnZ += (1 - len(nb[v])) * float(-np.sum(p * b))
    for (u, v) in edges:                              # edge terms
        arr = K + incoming(u, v)[:, None] + incoming(v, u)[None, :]
        arr = arr - logsumexp(arr)
        p = np.exp(arr)
        lnZ += float(np.sum(p * (K - arr)))
    return lnZ


# ---------------------------------------------------------------------------
# 1.  GBP contains BP
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('name,edges', [
    ('chain', [(0, 1), (1, 2), (2, 3)]),
    ('star', [(0, 1), (0, 2), (0, 3)]),
    ('four-cycle', [(0, 1), (1, 2), (2, 3), (0, 3)]),
    ('theta', [(0, 1), (1, 2), (2, 3), (0, 3), (0, 2)]),
])
def test_gbp_on_pairwise_regions_is_belief_propagation(name, edges):
    """The reduction that makes GBP the right generalisation.

    Regions = one per edge plus one per node, which is the Bethe region graph:
    ``D(P,R)`` is empty, ``N(P,R)`` is the other incoming messages, and the
    update is BP.  Trees come out exact; the loops come out with BP's own error,
    which is the point --- the two agree on the mistake as well as the answer.
    """
    nodes = sorted({v for e in edges for v in e})
    rg = RegionGraph([list(e) for e in edges])
    for bJ in (0.3, 0.8):
        f = ising_factors(edges, bJ)
        g = GBP(rg, f, damping=0.5).run(40_000)
        assert g.converged()
        assert g.log_Z() == pytest.approx(_bp_log_Z(edges, bJ, nodes), abs=1e-10)
    # and on the trees BP is exact, so GBP is too
    if name in ('chain', 'star'):
        f = ising_factors(edges, 0.8)
        g = GBP(rg, f, damping=0.5).run(40_000)
        assert g.log_Z() == pytest.approx(exact_log_Z(f, nodes), abs=1e-12)


# ---------------------------------------------------------------------------
# 2.  exact where the region graph is a junction tree
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('complexes', [
    [[0, 1, 2], [1, 2, 3]],                     # Table IV
    [[0, 1, 2], [1, 2, 3], [2, 3, 4]],
    [[0, 1, 2, 3], [1, 2, 3, 4]],
    [[0, 1, 2], [2, 3, 4], [4, 5, 6]],          # treelike: meets in one node
])
def test_gbp_is_exact_on_junction_trees(complexes):
    """Where the closed complex family is a junction tree, GBP is exact."""
    nodes = sorted({v for a in complexes for v in a})
    rg = RegionGraph(complexes)
    for bJ in (0.2, 0.5, 1.0, 2.0):
        f = ising_factors(clique_edges(complexes), bJ)
        g = GBP(rg, f, damping=0.0).run(500)
        assert g.converged()
        assert g.log_Z() == pytest.approx(exact_log_Z(f, nodes), abs=1e-9)
        assert g.consistency() < 1e-9


def test_gbp_recovers_what_table_four_leaves_on_the_table():
    """Two triangles sharing an edge: the manuscript's own numbers, closed.

    Table IV reports the Mobius counting evaluated on *isolated* regions, which
    is closer than the Bethe counting but not exact.  Running the messages
    closes the gap entirely, and the two static errors come out at the values
    the table quotes.
    """
    cx = [[0, 1, 2], [1, 2, 3]]
    rg = RegionGraph(cx)
    quoted = {0.2: (-1.4e-3, +1.8e-2), 0.5: (-2.9e-2, +9.1e-2),
              1.0: (-6.6e-2, +3.7e-1), 2.0: (-1.7e-2, +1.3)}
    for bJ, (kik_q, bet_q) in quoted.items():
        f = ising_factors(clique_edges(cx), bJ)
        ex = exact_log_Z(f, range(4))
        kik = static_log_Z(rg.counting, f) - ex
        bet = static_log_Z(rg.bethe_counting(), f) - ex
        assert kik == pytest.approx(kik_q, rel=0.05)
        assert bet == pytest.approx(bet_q, rel=0.05)
        assert abs(kik) < abs(bet)                       # Kikuchi is closer
        assert bet > 0                                   # Bethe overestimates
        g = GBP(rg, f, damping=0.0).run(500)
        assert g.log_Z() - ex == pytest.approx(0.0, abs=1e-10)


# ---------------------------------------------------------------------------
# 3.  a converged fixed point really is one
# ---------------------------------------------------------------------------

def test_converged_fixed_point_is_marginal_consistent():
    """On a region graph that is *not* a junction tree, GBP is approximate --
    but its fixed point still satisfies the condition it was derived from, so
    the residual error is the approximation's and not the solver's."""
    cx = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]     # K4 by its triangles
    rg = RegionGraph(cx)
    f = ising_factors(clique_edges(cx), 0.2)
    g = GBP(rg, f, damping=0.99).run(20_000)
    assert g.converged(tol=1e-9)
    assert g.consistency() < 1e-8
    assert abs(g.log_Z() - exact_log_Z(f, range(4))) < 0.05


def test_beliefs_are_normalised_and_marginals_agree():
    cx = [[0, 1, 2], [1, 2, 3]]
    rg = RegionGraph(cx)
    g = GBP(rg, ising_factors(clique_edges(cx), 0.7), damping=0.0).run(500)
    for r in g.regions:
        b = g.belief(r)
        assert b.sum() == pytest.approx(1.0)
        assert (b >= 0).all()
    m = g.magnetisation()
    assert set(m) == {0, 1, 2, 3}
    assert all(abs(v) < 1e-9 for v in m.values())        # zero field, symmetric


# ---------------------------------------------------------------------------
# validity of the region graph itself
# ---------------------------------------------------------------------------

def test_region_family_must_cover_the_interaction():
    """A factor counted twice is refused rather than silently approximated.

    The Bethe region graph of two edge-sharing triangles counts the shared bond
    in both complexes, which is exactly the miscounting Sec. VII describes.  GBP
    will not run on it.
    """
    cx = [[0, 1, 2], [1, 2, 3]]

    class Bare:
        counting = RegionGraph(cx).bethe_counting()

    with pytest.raises(ValueError, match='counted 2 times'):
        GBP(Bare(), ising_factors(clique_edges(cx), 0.5))


def test_static_log_Z_matches_a_hand_enumeration():
    """The static counting used for Table IV, against a direct sum."""
    cx = [[0, 1, 2], [1, 2, 3]]
    rg = RegionGraph(cx)
    edges = clique_edges(cx)
    f = ising_factors(edges, 0.5)

    def lnZ(region):
        R = sorted(region)
        e = [(i, j) for i, j in edges if i in R and j in R]
        tot = 0.0
        for s in itertools.product((1, -1), repeat=len(R)):
            d = dict(zip(R, s))
            tot += np.exp(0.5 * sum(d[i] * d[j] for i, j in e))
        return np.log(tot)

    want = sum(c * lnZ(r) for r, c in rg.counting.items() if c)
    assert static_log_Z(rg.counting, f) == pytest.approx(want, rel=1e-12)
