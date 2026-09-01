"""Critical temperature on an arbitrary chygraph."""

import numpy as np
import pytest
import sympy as sp

import chygraph_statmech as cs
from chygraph_statmech.ising import (branching_matrix, clique_derivative,
                                     critical_coupling, critical_temperature,
                                     perron_root)

t = sp.Symbol('t', positive=True)
beta, J = sp.symbols('beta J', positive=True)


def _root_in_unit(expr):
    return min(float(x) for x in sp.solve(expr, t) if 0 < float(x) < 1)


# ---------------------------------------------------------------------------
# Against the closed forms
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('c', [2, 3, 4, 5])
@pytest.mark.parametrize('bJ', [0.2, 0.6, 1.2])
def test_clique_derivative_matches_the_symbolic_one(c, bJ):
    exact = float(cs.cavity_derivative(c, cs.ising_clique(c)).subs(
        {beta: bJ, J: 1}))
    assert clique_derivative(c, bJ) == pytest.approx(exact, abs=1e-7)


@pytest.mark.parametrize('k', [3.0, 6.0, 12.0])
def test_graph_critical_coupling(k):
    """tanh(beta J_c) = 1 / <kbar>."""
    assert critical_coupling([2], [k]) == pytest.approx(
        np.arctanh(1 / k), rel=1e-9)


@pytest.mark.parametrize('k', [3.0, 6.0])
def test_triangle_critical_coupling(k):
    """2 k t / (1 - t + t^2) = 1, the WP1 closed form."""
    want = np.arctanh(_root_in_unit(sp.Eq(2 * k * t / (1 - t + t**2), 1)))
    assert critical_coupling([3], [k]) == pytest.approx(want, rel=1e-9)


def test_two_layer_matches_the_branching_determinant():
    """k_L t + 2 k_T u_T = 1 for Poisson layers, from WP1."""
    kL, kT = 4.0, 2.0
    want = np.arctanh(_root_in_unit(
        sp.Eq(kL * t + 2 * kT * t / (1 - t + t**2), 1)))
    assert critical_coupling([2, 3], [kL, kT]) == pytest.approx(want, rel=1e-9)


# ---------------------------------------------------------------------------
# Structure
# ---------------------------------------------------------------------------

def test_perron_root_is_one_at_the_transition():
    for spec in (([2], [6.0]), ([3], [3.0]), ([4], [2.0]), ([2, 3], [3.0, 1.0])):
        bc = critical_coupling(*spec)
        assert perron_root(*spec, bc) == pytest.approx(1.0, abs=1e-9)


def test_ordering_is_easier_with_more_neighbours():
    prev = np.inf
    for k in (2.0, 4.0, 8.0, 16.0):
        bc = critical_coupling([2], [k])
        assert bc < prev
        prev = bc


@pytest.mark.parametrize('c', [3, 4, 5, 6])
def test_cliques_order_before_the_matched_graph(c):
    """At matched neighbour count, a clique layer orders at weaker coupling.

    The WP1 triangle result, now for any cardinality: solving the complex
    exactly beats treating its members as independent edges.
    """
    k = 2.0
    clique = critical_coupling([c], [k])
    matched = critical_coupling([2], [k * (c - 1)])
    assert clique < matched


def test_at_line_is_below_the_ferromagnetic_one():
    """u'^2 < u', so the spin-glass instability needs stronger coupling."""
    for spec in (([2], [6.0]), ([3], [3.0]), ([4], [2.0])):
        assert critical_coupling(*spec, squared=True) > critical_coupling(*spec)


def test_critical_temperature_is_the_reciprocal():
    assert critical_temperature([2], [6.0]) == pytest.approx(
        1 / critical_coupling([2], [6.0]), rel=1e-12)


def test_no_transition_is_reported_not_guessed():
    with pytest.raises(ValueError, match='no transition'):
        critical_coupling([2], [0.5])          # <kbar> < 1, never orders


def test_excess_defaults_to_the_mean():
    B1 = branching_matrix([2, 3], [3.0, 1.0], 0.2)
    B2 = branching_matrix([2, 3], [3.0, 1.0], 0.2, excess=[3.0, 1.0])
    assert np.allclose(B1, B2)
