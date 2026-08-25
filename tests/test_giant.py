"""Tests for the giant component extension.

The key structural claim is that the published threshold tensor ``A`` is the
Jacobian at ``Q = 1`` of the non-linear map whose fixed point gives the order
parameter, so the two calculations share one object and cannot drift apart.
"""

import math

from sympy import symbols, zeros, simplify, expand, nsolve
import pytest

from chygraph.percolation import PercolationMatrix, GraphWithTriangles
from chygraph.giant import (
    GiantComponent, hypergraph_giant, multiplex_hypergraph_giant,
    graph_with_triangles_giant,
)


# ---------------------------------------------------------------------------
# A is the Jacobian of the non-linear map
# ---------------------------------------------------------------------------

def test_hypergraph_jacobian_is_published_A():
    p, q, k, c = symbols('p q k c')
    kk, KK, ss, SS = (zeros(2, 2) for _ in range(4))
    kk[0, 1] = KK[0, 1] = p * k
    ss[1, 0] = SS[1, 0] = q * c
    assert simplify(hypergraph_giant().A() - PercolationMatrix(kk, KK, ss, SS).A) == zeros(8, 8)


def test_hypergraph_theta_matches_paper():
    p, q, k, c = symbols('p q k c')
    assert simplify(hypergraph_giant().theta() - (c * k * p * q - 1)) == 0


def test_graph_with_triangles_theta_matches_paper():
    """Eq. (33): theta = q k_L + 2q(1 + q - q^2) k_T - 1 for Poisson degrees."""
    q, kL, kT = symbols('q k_L k_T')
    got = expand(simplify(graph_with_triangles_giant().theta()))
    want = expand(q * kL + 2 * q * (1 + q - q**2) * kT - 1)
    assert expand(got - want) == 0
    assert expand(got - expand(simplify(GraphWithTriangles(poisson=True).A.theta()))) == 0


def test_multiplex_theta_is_additive_for_poisson():
    """Eq. (22): theta = sum_l q_l <k>_l <c>_l - 1."""
    MH = multiplex_hypergraph_giant(number_of_types=2)
    k01, k02, c10, c20, q1, q2 = symbols('k_01 k_02 c_10 c_20 q_1 q_2')
    want = q1 * k01 * c10 + q2 * k02 * c20 - 1
    assert expand(simplify(MH.theta() - want)) == 0


# ---------------------------------------------------------------------------
# The order parameter itself
# ---------------------------------------------------------------------------

def test_er_graph_recovers_S_equals_one_minus_exp():
    G = hypergraph_giant(graph=True)
    for k in (1.5, 2.0, 3.0, 5.0):
        S = G.node_fraction({'k': k, 'p': 1, 'q': 1})
        assert abs(S - (1 - math.exp(-k * S))) < 1e-9


def test_no_giant_component_below_threshold():
    G = hypergraph_giant()
    # k*c*p*q = 0.8 < 1
    assert G.node_fraction({'k': 2, 'c': 2, 'p': 1, 'q': 0.2}) < 1e-6
    assert G.node_fraction({'k': 2, 'c': 2, 'p': 1, 'q': 0.4}) > 0.1


def test_site_percolation_bounds_node_fraction_by_p():
    """Absent nodes are never in the giant component, so S <= p."""
    G = hypergraph_giant()
    for p in (0.3, 0.6, 0.9):
        assert G.node_fraction({'k': 6, 'c': 4, 'p': p, 'q': 1}) <= p + 1e-12


def test_triangle_model_matches_reference_values():
    """Values independently confirmed by Monte Carlo on configuration-model
    graphs with n = 3e5 (agreement to within one MC standard deviation)."""
    T = graph_with_triangles_giant()
    for (kL, kT, q), want in {
        (1.0, 0.5, 1.0): 0.67400,
        (1.0, 0.5, 0.7): 0.47775,
        (0.0, 1.0, 1.0): 0.54924,
        (0.5, 1.5, 0.8): 0.81832,
        (3.0, 2.0, 0.35): 0.85963,
    }.items():
        got = T.node_fraction({'k_L': kL, 'k_T': kT, 'q': q})
        assert abs(got - want) < 1e-4, (kL, kT, q, got)


# ---------------------------------------------------------------------------
# Critical amplitude:  S ~ B * Lambda
# ---------------------------------------------------------------------------

def test_er_critical_amplitude_is_four():
    """S ~ 2(k-1) and Lambda = sqrt(k)-1 ~ (k-1)/2, so B = 4."""
    assert abs(hypergraph_giant(graph=True).amplitude({'k': 1, 'p': 1, 'q': 1}) - 4) < 1e-9


def test_amplitude_predicts_the_approach_to_the_threshold():
    G = hypergraph_giant()
    B = G.amplitude({'k': 2, 'c': 2, 'p': 1, 'q': 0.25})
    q = 0.2501
    S = G.node_fraction({'k': 2, 'c': 2, 'p': 1, 'q': q}, tol=1e-16, maxiter=4_000_000)
    assert abs(S / (math.sqrt(4 * q) - 1) - B) < 1e-3
