"""Tests for the finite-component size distribution (percolation.components).

Three routes to the same coefficients -- the Cauchy-integral iteration, the
exact marked-map series and Good's multivariate inversion -- are checked
against one another, and the first against the Borel distribution on an
Erdos-Renyi graph, which is the closed form the whole construction has to
reproduce.
"""

import numpy as np
from sympy import Rational

from percolation.giant import Chygraph, _tables, finite_pgf, hypergraph_giant
from percolation.applications import household_epidemic
from percolation.components import (
    ComponentDistribution, borel, symbolic_series, good_by_layer,
)


def _finite_two_layer():
    """Nodes with degree 1, 2, 3 in hyperedges of cardinality 2 or 3."""
    phi, phibar, g, gbar = _tables(2)
    deg = {1: Rational(1, 2), 2: Rational(1, 4), 3: Rational(1, 4)}
    mean = sum(k * v for k, v in deg.items())
    exc = {k - 1: k * v / mean for k, v in deg.items()}
    card = {2: Rational(1, 2), 3: Rational(1, 2)}
    cm = sum(k * v for k, v in card.items())
    cexc = {k - 1: k * v / cm for k, v in card.items()}
    phi[0][1], phibar[0][1] = finite_pgf(deg), finite_pgf(exc)
    g[1][0], gbar[1][0] = finite_pgf(card), finite_pgf(cexc)
    return Chygraph(phi, phibar, g, gbar)


def test_er_graph_reproduces_borel_below_and_above_threshold():
    G = hypergraph_giant(graph=True)
    for c in (0.5, 1.5):
        D = ComponentDistribution(G, {'k': c, 'p': 1, 'q': 1})
        P = D.distribution(300, radius=0.99)
        assert np.max(np.abs(P - borel(c, 300))) < 1e-11
        # the mass of the finite components is 1 - S, up to the truncation
        # (the Borel tail at c = 1.5 decays as 0.91^s)
        assert abs(P.sum() - D.finite_fraction()) < 1e-9


def test_numerical_matches_exact_series():
    M = _finite_two_layer()
    S = symbolic_series(M, 6)
    J = ComponentDistribution(M).joint(8)
    assert S, "the exact series produced no terms"
    for mon, c in S.items():
        assert abs(float(c) - J[mon]) < 1e-12, mon


def test_good_inversion_matches_exact_series():
    M = _finite_two_layer()
    S = symbolic_series(M, 4)
    # good_by_layer counts complexes beyond the root, so (a, b) in the series
    # (a nodes including the root, b hyperedges) is (a - 1, b) here
    for mon in ((2, 1), (3, 1)):
        assert good_by_layer(M, (mon[0] - 1, mon[1])) == S[mon], mon


def test_joint_marginal_is_the_univariate_distribution():
    M = household_epidemic({3: 1})
    subs = {'k': 2, 'p_H': 0.5, 'T': 0.2}
    D = ComponentDistribution(M, subs)
    P = D.distribution(40, weights=[1, 0, 0], radius=0.9)
    J = D.joint((40, 30), layers=(0, 1), radius=0.9)
    assert np.max(np.abs(J.sum(axis=1) - P)) < 1e-9
    # an outbreak of s individuals touches at most s households and at
    # least ceil(s / 3) of them, in households of size three; the floor is
    # the absolute accuracy of the transform, 1e-16 / r^s
    for s in range(1, 41):
        row = J[s]
        assert np.abs(row[s + 1:]).sum() < 1e-8
        assert np.abs(row[:(s + 2) // 3]).sum() < 1e-8


def test_household_distribution_mass_is_one_minus_S():
    M = household_epidemic({2: Rational(1, 2), 4: Rational(1, 2)})
    for T in (0.1, 0.4):
        subs = {'k': 2, 'p_H': 0.5, 'T': T}
        D = ComponentDistribution(M, subs)
        P = D.distribution(300, radius=0.99)
        S = M.node_fraction(subs)
        assert abs(P.sum() - (1 - S)) < 1e-6
