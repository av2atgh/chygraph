"""Tests for the loop series (statmech.loopseries).

The Chertkov--Chernyak identity Z = Z_BP (1 + sum_C r_C) is exact on any
finite factor graph, so every test is the same test: enumerate the
generalised loops, evaluate the terms at the BP fixed point, and compare
ln(1 + sum r_C) with ln Z - ln Z_BP from exact enumeration.  The factor
graphs are the pairwise one and the promoted one of Chapter 19, which are
two representations of the same model with different loops.
"""

import math

import networkx as nx
import numpy as np

from statmech.loopseries import (
    BinaryFactorGraph, generalised_loops, loop_series, partial_sums,
)


def _identity(fg, tol=1e-12):
    fg.bp(damping=0.5)
    assert fg.converged(), fg.residual
    terms = loop_series(fg)
    total = sum(r for _, _, r in terms)
    assert abs(fg.exact_log_Z() - fg.log_Z_bethe() - math.log1p(total)) < tol
    return terms


def _ring(k):
    G = nx.Graph()
    hub = [3 * i for i in range(k)]
    for i in range(k):
        a, b, c = hub[i], 3 * i + 1, hub[(i + 1) % k]
        G.add_edges_from([(a, b), (b, c), (a, c)])
    return nx.convert_node_labels_to_integers(G)


def test_single_cycle_is_tanh_cubed():
    fg = BinaryFactorGraph.pairwise([(0, 1), (1, 2), (0, 2)], 0.5)
    terms = _identity(fg)
    assert len(terms) == 1
    assert abs(terms[0][2] - math.tanh(0.5) ** 3) < 1e-14


def test_two_triangles_promoted_has_one_loop_and_pairwise_four():
    cx = [(0, 1, 2), (1, 2, 3)]
    e = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]
    for bJ in (0.5, 1.0):
        assert len(_identity(BinaryFactorGraph.promoted(cx, e, bJ))) == 1
        assert len(_identity(BinaryFactorGraph.pairwise(e, bJ))) == 4
        # the double-count model is a different model, and the identity
        # holds for it too
        _identity(BinaryFactorGraph.promoted(cx, e, bJ, assign='all'))


def test_promoted_once_represents_the_pairwise_model():
    cx = [(0, 1, 2), (1, 2, 3)]
    e = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]
    a = BinaryFactorGraph.promoted(cx, e, 0.7).exact_log_Z()
    b = BinaryFactorGraph.pairwise(e, 0.7).exact_log_Z()
    c = BinaryFactorGraph.promoted(cx, e, 0.7, assign='all').exact_log_Z()
    assert abs(a - b) < 1e-12 and abs(c - b) > 1e-3


def test_ring_of_triangles_promoted_is_one_loop():
    for k in (4, 5):
        G = _ring(k)
        cx = [tuple(sorted(c)) for c in nx.find_cliques(G)]
        e = sorted(tuple(sorted(x)) for x in G.edges())
        for bJ in (0.3, 0.8):
            tp = _identity(BinaryFactorGraph.promoted(cx, e, bJ))
            assert len(tp) == 1 and sum(len(l) for _, l in tp[0][0]) == 2 * k
            tq = _identity(BinaryFactorGraph.pairwise(e, bJ))
            assert len(tq) > 100
            ps = partial_sums(tq)
            assert ps[0][0] == 6 and ps[-1][1] > ps[0][1] > 0


def test_random_graph_identity_and_loop_count():
    G = nx.gnm_random_graph(9, 14, seed=3)
    e = sorted(tuple(sorted(x)) for x in G.edges())
    fg = BinaryFactorGraph.pairwise(e, 0.4)
    terms = _identity(fg)
    # every generalised loop has minimum degree two
    for loop, size, _ in terms:
        deg = {}
        for a, legs in loop:
            for v in legs:
                deg[v] = deg.get(v, 0) + 1
        assert min(deg.values()) >= 2 and size == sum(len(l) for _, l in loop)


def test_max_edges_cutoff_is_a_prefix():
    G = nx.gnm_random_graph(8, 12, seed=1)
    e = sorted(tuple(sorted(x)) for x in G.edges())
    fg = BinaryFactorGraph.pairwise(e, 0.3).bp()
    full = {loop for loop in generalised_loops(fg)}
    cut = {loop for loop in generalised_loops(fg, max_edges=8)}
    assert cut == {l for l in full if sum(len(x) for _, x in l) <= 8}


def test_explicit_assignment_fills_the_unshared_bonds():
    # two triangles sharing the bond (1, 2): naming an owner for that bond
    # alone must leave every other bond in place, so the exact model does
    # not change and the bond moves from one table to the other
    cx, e = [(0, 1, 2), (1, 2, 3)], [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]
    fgs = [BinaryFactorGraph.promoted(cx, e, 0.5, {(1, 2): a}) for a in (0, 1)]
    exact = [fg.exact_log_Z() for fg in fgs]
    assert abs(exact[0] - exact[1]) < 1e-12
    assert abs(exact[0] - BinaryFactorGraph.pairwise(e, 0.5).exact_log_Z()) < 1e-12
    weight = [[float(np.abs(t).sum()) for _, t in fg.factors] for fg in fgs]
    assert weight[0][0] > weight[1][0] and weight[0][1] < weight[1][1]
    assert abs(sum(weight[0]) - sum(weight[1])) < 1e-12
    once = BinaryFactorGraph.promoted(cx, e, 0.5, 'once').bp()
    assert abs(fgs[0].bp().log_Z_bethe() - once.log_Z_bethe()) < 1e-12
