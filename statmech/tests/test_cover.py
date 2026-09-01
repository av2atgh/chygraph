"""Minimum vertex cover of the graph a chygraph induces."""

import numpy as np
import pytest
from scipy.special import lambertw

import chygraph_statmech.cover as cv
from chygraph_statmech import antimonotone as am
import chygraph_statmech.core as cp


@pytest.mark.parametrize('k', [0.5, 1.0, 2.0, 2.7, 5.0, 10.0])
def test_graph_reduces_to_weigt_hartmann(k):
    """c = 2 is the graph: x_c = 1 - (2W + W^2)/(2k)."""
    W = lambertw(k).real
    assert cv.poisson([2], [k]).cover_size() == pytest.approx(
        1 - (2 * W + W**2) / (2 * k), abs=1e-12)


def test_map_is_anti_monotone():
    for spec in (([2], [4.0]), ([3], [2.0]), ([2, 3], [2.0, 1.0])):
        m = cv.poisson(*spec)
        assert am.is_anti_monotone(m.F, m.L)


def test_cover_grows_with_density():
    prev = 0.0
    for k in (0.5, 1.0, 2.0, 4.0, 8.0):
        x = cv.poisson([2], [k]).cover_size()
        assert x > prev
        prev = x


def test_bracket_contains_the_reported_value():
    for spec in (([2], [2.0]), ([3], [1.0]), ([2, 3], [1.0, 0.5])):
        m = cv.poisson(*spec)
        lo, hi = m.cover_bracket()
        assert lo <= m.cover_size() <= hi


# ---------------------------------------------------------------------------
# Where it is exact, and where it is not
# ---------------------------------------------------------------------------

def test_graph_is_certified_and_cliques_are_not():
    assert cv.poisson([2], [2.0]).certified()
    for c in (3, 4, 5):
        assert not cv.poisson([c], [1.0]).certified()
    assert not cv.poisson([2, 3], [2.0, 1.0]).certified()


@pytest.mark.parametrize('c', [3, 4, 5])
def test_lack_of_certification_is_the_core_result(c):
    """`certified` is not a separate assumption: it is exactly the statement
    that a cardinality >= 3 layer has no core-free branch, so leaf removal
    leaves an extensive core and never proves a cover minimal."""
    assert (not cv.poisson([c], [1.0]).certified()) == (
        not cp.clique_network(c, 1.0).has_core_free_branch())


def test_the_triangle_counterexample():
    """One triangle per node and nothing else -- isolated triangles.

    Every vertex sits at z = 0, so the graph degeneracy rule counts half of
    them and returns 1/2.  The truth is 2/3: a triangle's maximum independent
    set is one vertex, so two of three must be covered.  The bracket is honest
    but vacuous.  This is the number behind `certified` being False at
    cardinality >= 3.

    Note the pgf: `x0` is exactly one complex, not Poisson(1).
    """
    from sympy import Symbol
    m = cv.CliqueCover([3], Symbol('x0'))
    assert m.cover_size() == pytest.approx(0.5, abs=1e-9)
    lo, hi = m.cover_bracket()
    assert lo <= 2 / 3 <= hi
    assert abs(m.cover_size() - 2 / 3) > 0.15


def test_isolated_edges_are_exact():
    """The same construction at c = 2: one edge per node, cover exactly half.
    Here the degeneracy rule is right, which is the contrast."""
    from sympy import Symbol
    assert cv.CliqueCover([2], Symbol('x0')).cover_size() == pytest.approx(
        0.5, abs=1e-9)
