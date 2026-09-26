"""Tests for the Bethe Hessian of a chygraph (statmech.bethehessian).

The weighted Ihara--Bass identity is exact on any finite graph, so the
tests are the identity on small instances: at the trivial fixed point, at
polarised and random-field fixed points where every complex has at most
three members, and its failure where a complex of four members sits at a
random-field fixed point.  Around it: the analytic Jacobian against finite
differences, the block at the trivial point against ``u'(c)``, and the
Schur complement against Saade's Bethe Hessian on a graph.
"""

import networkx as nx
import numpy as np

from statmech.bethehessian import (
    bethe_hessian, belief_covariance, clique_factor_graph, edge_weights,
    factorise, instance_threshold, jacobian_block, linearised_operator,
    random_chygraph, spectral_radius, trivial_blocks, trivial_vertex_hessian,
    vertex_hessian,
)
from statmech.ising import clique_derivative, critical_coupling

SPIN = np.array([1.0, -1.0])
MIXED = [(0, 1, 2), (1, 2, 3), (3, 4), (2, 4, 5, 6), (5, 6, 7)]
SMALL = [(0, 1, 2), (1, 2, 3), (3, 4), (2, 4, 5), (5, 6, 7), (0, 7)]


def _identity(fg):
    T = linearised_operator(fg)
    alpha, beta, worst = edge_weights(fg)
    H, pref = bethe_hessian(fg, alpha, beta)
    S, det_cc = vertex_hessian(fg, alpha, beta)
    lhs = np.linalg.det(np.eye(T.shape[0]) - T.toarray())
    return lhs, pref * np.linalg.det(H.toarray()), \
        pref * det_cc * np.linalg.det(S.toarray()), worst


def test_block_is_uprime_at_trivial_point():
    fg = clique_factor_graph(MIXED, 0.4)
    for (sc, _), M in zip(fg.factors, trivial_blocks(fg)):
        c = len(sc)
        u = clique_derivative(c, 0.4)
        off = M[~np.eye(c, dtype=bool)]
        assert np.allclose(off, u, atol=1e-14)
        assert np.all(np.diag(M) == 0)


def test_jacobian_against_finite_differences():
    rng = np.random.default_rng(0)
    fg = clique_factor_graph(MIXED, 0.6, fields=rng.normal(0, 0.7, 8)).bp()
    assert fg.converged()
    eps = 1e-6
    for a, (sc, _) in enumerate(fg.factors):
        M = jacobian_block(fg, a)
        C, m = belief_covariance(fg, a)
        Mc = C / (1 - m ** 2)[:, None]
        np.fill_diagonal(Mc, 0.0)
        assert np.abs(Mc - M).max() < 1e-12
        for j, u in enumerate(sc):
            base = fg.m_va[(u, a)].copy()
            out = []
            for sgn in (+1, -1):
                fg.m_va[(u, a)] = base + sgn * eps * SPIN
                out.append(np.array([0.5 * float(fg._factor_to_var(a, v) @ SPIN)
                                     for v in sc]))
            fg.m_va[(u, a)] = base
            assert np.abs((out[0] - out[1]) / (2 * eps) - M[:, j]).max() < 1e-8


def test_identity_at_trivial_point_and_closed_form():
    fg = clique_factor_graph(MIXED, 0.4)
    blocks = trivial_blocks(fg)
    T = linearised_operator(fg, blocks)
    alpha, beta, worst = edge_weights(fg, blocks)
    assert worst < 1e-14
    H, pref = bethe_hessian(fg, alpha, beta)
    lhs = np.linalg.det(np.eye(T.shape[0]) - T.toarray())
    assert abs(lhs - pref * np.linalg.det(H.toarray())) < 1e-12
    S, pref2 = trivial_vertex_hessian(fg, [M[0, 1] for M in blocks])
    assert abs(lhs - pref2 * np.linalg.det(S.toarray())) < 1e-12
    Hs, _ = bethe_hessian(fg, alpha, beta, symmetric=True)
    assert np.allclose(Hs.toarray(), Hs.toarray().T)


def test_identity_at_nontrivial_points_up_to_three_members():
    rng = np.random.default_rng(1)
    fg = clique_factor_graph(SMALL, 0.6, fields=rng.normal(0, 0.7, 8)).bp()
    assert fg.converged()
    lhs, rhs, schur, worst = _identity(fg)
    assert worst < 1e-13
    assert abs(lhs - rhs) < 1e-12 and abs(lhs - schur) < 1e-12
    fg = clique_factor_graph(SMALL, 1.2)
    for k in fg.m_va:
        fg.m_va[k] = np.array([0.3, -0.3])
    fg.bp()
    assert fg.converged() and fg.magnetisation(0) > 0.9
    lhs, rhs, schur, worst = _identity(fg)
    assert worst < 1e-13
    assert abs(lhs - rhs) < 1e-12 and abs(lhs - schur) < 1e-12


def test_four_members_do_not_factorise_off_the_trivial_point():
    rng = np.random.default_rng(0)
    fg = clique_factor_graph(MIXED, 0.6, fields=rng.normal(0, 0.7, 8)).bp()
    res = [factorise(jacobian_block(fg, a))[2] for a in range(len(fg.factors))]
    cards = [len(sc) for sc, _ in fg.factors]
    for c, r in zip(cards, res):
        assert (r > 1e-3) if c == 4 else (r < 1e-13)
    lhs, rhs, _, _ = _identity(fg)
    assert abs(lhs - rhs) > 1e-4


def test_schur_complement_is_saade_on_a_graph():
    G = nx.random_regular_graph(3, 20, seed=1)
    E = [tuple(e) for e in G.edges()]
    fg = clique_factor_graph(E, 0.4)
    t = np.tanh(0.4)
    S, _ = trivial_vertex_hessian(fg, [t] * len(E))
    A = nx.to_numpy_array(G, nodelist=fg.nodes)
    r = 1.0 / t
    Hs = (r * r - 1) * np.eye(20) - r * A + np.diag(A.sum(1))
    assert np.abs(S.toarray() - (t * t / (1 - t * t)) * Hs).max() < 1e-12


def test_instance_threshold_by_either_route_and_near_the_ensemble():
    rng = np.random.default_rng(3)
    cx = random_chygraph(400, [2, 3], [1.5, 0.7], rng)
    bT = instance_threshold(cx, via='T')
    bH = instance_threshold(cx, via='H')
    assert abs(bT - bH) < 1e-7
    bc = critical_coupling([2, 3], [1.5, 0.7])
    assert abs(bT - bc) / bc < 0.1
