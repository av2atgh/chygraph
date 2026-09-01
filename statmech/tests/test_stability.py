"""WP1 checks.  Every one of these is a statement the README makes."""

import itertools

import pytest
import sympy as sp

import statmech as cs
from percolation.percolation import PercolationMatrix
from statmech.stability import StabilityMatrix

t = sp.Symbol('t', positive=True)
beta, J = sp.Symbol('beta', positive=True), sp.Symbol('J', positive=True)


def _tables(L, tag):
    """Generic symbolic moment tables."""
    m = [sp.Matrix(L, L, lambda i, j: sp.Symbol(f'{tag}{n}_{i}{j}'))
         for n in range(4)]
    return m


# ---------------------------------------------------------------------------
# Prediction 1: unit weights reproduce chygraph identically
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('L', [1, 2, 3])
def test_unit_weights_are_the_percolation_matrix(L):
    k, K, s, S = _tables(L, 'a')
    assert sp.simplify(StabilityMatrix(k, K, s, S).A
                       - PercolationMatrix(k, K, s, S).A) == sp.zeros(2 * L * L)


def test_graph_percolation_matches_chygraph():
    from percolation import HypergraphPercolation
    p = sp.Symbol('p')
    theirs = HypergraphPercolation(graph=True).A.A.subs(
        {sp.Symbol('q'): p, sp.Symbol('p'): 1})
    assert sp.simplify(cs.graph_percolation(p).A - theirs) == sp.zeros(8, 8)


def test_graph_percolation_is_molloy_reed():
    p = sp.Symbol('p')
    assert sp.factor(cs.graph_percolation(p).theta()) == sp.Symbol('K') * p - 1


# ---------------------------------------------------------------------------
# Known Ising results on a configuration-model graph
# ---------------------------------------------------------------------------

def test_edge_cavity_derivative():
    assert sp.simplify(cs.ising_edge_derivative() - sp.tanh(beta * J)) == 0


def test_graph_ising_ferromagnetic_line():
    """<kbar> tanh(beta J) = 1."""
    th = sp.factor(sp.simplify(cs.graph_ising().theta()))
    assert sp.simplify(th - (sp.Symbol('K') * sp.tanh(beta * J) - 1)) == 0


def test_graph_ising_at_line():
    """<kbar> tanh^2(beta J) = 1."""
    th = sp.factor(sp.simplify(cs.graph_ising(squared=True).theta()))
    assert sp.simplify(th - (sp.Symbol('K') * sp.tanh(beta * J)**2 - 1)) == 0


# ---------------------------------------------------------------------------
# The complex solved exactly inside itself
# ---------------------------------------------------------------------------

def test_triangle_cavity_derivative_closed_form():
    """u'_triangle = t / (1 - t + t^2), t = tanh(beta J)."""
    got = cs.in_tanh(cs.ising_triangle_derivative())
    assert sp.simplify(got - t / (1 - t + t**2)) == 0


def test_triangle_derivative_by_finite_difference():
    """Independent check of the symbolic derivative against the definition."""
    energy = cs.ising_clique(3, 1, sp.Rational(7, 10))
    h1, h2 = sp.symbols('h1 h2')
    u = cs.emitted_field(3, energy, (h1, h2))
    eps = sp.Rational(1, 10**6)
    fd = ((u.subs({h1: eps, h2: 0}) - u.subs({h1: -eps, h2: 0})) / (2 * eps))
    exact = cs.ising_triangle_derivative(1, sp.Rational(7, 10))
    assert abs(complex(fd.evalf()) - complex(exact.evalf())) < 1e-9


def test_triangle_transmits_more_than_an_edge():
    """u'_triangle > t on 0 < t < 1: clustering helps ferromagnetic order."""
    u_T = t / (1 - t + t**2)
    for v in (0.1, 0.3, 0.5, 0.7, 0.9, 0.99):
        assert float(u_T.subs(t, v)) > v


def test_triangle_is_two_edges_at_leading_order():
    """A triangle looks like two independent edges as t -> 0."""
    u_T = t / (1 - t + t**2)
    assert sp.series(u_T, t, 0, 2).removeO() == t


# ---------------------------------------------------------------------------
# Graph with links and triangles: the new number
# ---------------------------------------------------------------------------

def _triangle_theta():
    return sp.factor(sp.simplify(
        cs.in_tanh(sp.expand(cs.graph_with_triangles_ising().theta()))))


def test_triangle_model_is_the_branching_determinant():
    """theta = 0 is det(I - B) = 0 for the 2x2 layer branching matrix B,
    with the triangle's u_T in place of t."""
    K_L, K_T, k_L, k_T = sp.symbols('K_L K_T k_L k_T')
    u_T = t / (1 - t + t**2)
    B = sp.Matrix([[K_L * t, 2 * k_T * u_T],
                   [k_L * t, 2 * K_T * u_T]])
    expected = -(sp.eye(2) - B).det()
    assert sp.simplify(_triangle_theta() - sp.simplify(expected)) == 0


def test_triangle_model_reduces_to_plain_graph():
    K_L, K_T, k_L, k_T = sp.symbols('K_L K_T k_L k_T')
    got = _triangle_theta().subs({K_T: 0, k_T: 0})
    assert sp.simplify(got - (K_L * t - 1)) == 0


def test_pure_triangle_graph():
    """No links: 2 <kbar_T> u_T = 1."""
    K_L, K_T, k_L, k_T = sp.symbols('K_L K_T k_L k_T')
    got = _triangle_theta().subs({K_L: 0, k_L: 0})
    assert sp.simplify(got - (2 * K_T * t / (1 - t + t**2) - 1)) == 0


def test_clustering_raises_the_critical_temperature():
    """At fixed mean degree, moving edges into triangles lowers t_c.

    A node in <k_T> triangles has 2<k_T> neighbours; hold the total excess
    degree fixed and compare all-links against all-triangles (Poisson, K = k).
    """
    K_L, K_T, k_L, k_T = sp.symbols('K_L K_T k_L k_T')
    th = _triangle_theta()
    links = sp.solve(th.subs({K_L: 6, k_L: 6, K_T: 0, k_T: 0}), t)
    tris = sp.solve(th.subs({K_L: 0, k_L: 0, K_T: 3, k_T: 3}), t)
    tc_links = min(float(x) for x in links if x.is_real and 0 < x < 1)
    tc_tris = min(float(x) for x in tris if x.is_real and 0 < x < 1)
    assert tc_tris < tc_links
