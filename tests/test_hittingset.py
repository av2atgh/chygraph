"""WP3 checks: minimum hitting set on a hypergraph."""

import numpy as np
import pytest
from scipy.special import lambertw

from chygraph_statmech import antimonotone as am
from chygraph_statmech import hittingset as hs

E = np.e


def _matched(cs, m, sign, s):
    """Two node classes, equal weight, identical marginals, correlated across
    cardinality layers.  ``s = 0`` is the uncorrelated Poisson ensemble."""
    def f(t):
        if sign > 0:
            a = [t * mi * (1 + s) for mi in m]
            b = [t * mi * (1 - s) for mi in m]
        else:
            a = [t * m[0] * (1 + s), t * m[1] * (1 - s)]
            b = [t * m[0] * (1 - s), t * m[1] * (1 + s)]
        return hs.HittingSet(cs, hs.two_class_phi(a, b, 0.5))
    return f


# ---------------------------------------------------------------------------
# Reduction to the graph case
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('k', [0.5, 1.0, 2.0, 2.7, 5.0, 10.0])
def test_cardinality_two_is_weigt_hartmann(k):
    """c = 2 is vertex cover: x_c = 1 - (2W + W^2)/(2k)."""
    W = lambertw(k).real
    assert hs.poisson([2], [k]).cover_size() == pytest.approx(
        1 - (2 * W + W**2) / (2 * k), rel=1e-11)


def test_graph_rsb_point_is_e():
    """Bauer-Golinelli core percolation at c = e -- the Erdos-Renyi control in
    ~/av2atg/computational_complexity."""
    assert hs.rsb_point([2], [1.0]) == pytest.approx(E, rel=1e-12)


# ---------------------------------------------------------------------------
# The new closed form
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('c', [2, 3, 4, 5, 10, 20])
def test_fixed_cardinality_rsb_is_e_over_c_minus_one(c):
    """k_RSB = e/(c-1): counted in neighbours, cardinality does not move it."""
    k = hs.rsb_point([c], [1.0])
    assert k == pytest.approx(E / (c - 1), rel=1e-12)
    assert k * (c - 1) == pytest.approx(E, rel=1e-12)


def test_rs_below_and_broken_above():
    for c in (2, 3, 5):
        k = E / (c - 1)
        assert not hs.poisson([c], [0.9 * k]).is_unstable()
        assert hs.poisson([c], [1.1 * k]).is_unstable()


# ---------------------------------------------------------------------------
# The anti-monotone solver
# ---------------------------------------------------------------------------

def test_map_is_anti_monotone():
    for cs, m in (([2], [4.0]), ([3], [2.0]), ([2, 6], [3.0, 1.0])):
        model = hs.poisson(cs, m)
        assert am.is_anti_monotone(model.F, model.L)


def test_bracket_agrees_with_the_jacobian():
    """Two independent RS criteria: does F o F close, and is |lambda| < 1."""
    for c in (2, 3, 5):
        k = E / (c - 1)
        for scale in (0.5, 0.9, 0.99, 1.01, 1.2, 2.0):
            model = hs.poisson([c], [scale * k])
            assert model.is_replica_symmetric() == (not model.is_unstable())


def test_bracket_splits_into_a_two_cycle_past_rsb():
    """a < b with F(a) = b: the period-2 orbit that VW03 sees as
    non-convergence."""
    model = hs.poisson([2], [6.0])
    a, b = am.bracket(model.F, 1)
    assert b[0] - a[0] > 1e-3
    assert model.F(a) == pytest.approx(b, abs=1e-9)
    assert model.F(b) == pytest.approx(a, abs=1e-9)


def test_fixed_point_found_past_the_instability():
    model = hs.poisson([2], [6.0])
    sigma = model.solve()
    assert model.is_unstable(sigma)
    assert model.F(sigma) == pytest.approx(sigma, abs=1e-12)


# ---------------------------------------------------------------------------
# The physics question
# ---------------------------------------------------------------------------

def test_cardinality_heterogeneity_postpones_rsb():
    """At fixed mean excess cardinality, mixing cardinalities pushes RSB up.

    Same direction as degree heterogeneity in VW03: heterogeneity keeps the
    problem easy longer.  The AND-structure of a hyperedge does not offset it.
    """
    homo = hs.rsb_point([3], [1.0]) * 2.0                  # <cbar> = 2
    hetero = hs.rsb_point([2, 6], [0.75, 0.25]) * 2.0      # <cbar> = 2
    assert homo == pytest.approx(E, rel=1e-12)
    assert hetero > 2.0 * homo


def test_uncorrelated_limit_is_the_poisson_ensemble():
    """s = 0 collapses both correlated families onto the independent one."""
    cs, m = [2, 6], [0.75, 0.25]
    base = hs.rsb_point(cs, m)
    for sign in (+1, -1):
        assert hs.rsb_scale(_matched(cs, m, sign, 0.0)) == pytest.approx(
            base, rel=1e-9)


def test_negative_correlation_brings_rsb_forward():
    """Anti-correlating participation in small and large hyperedges makes
    hitting set harder -- an axis with no counterpart in VW03."""
    cs, m = [2, 6], [0.75, 0.25]
    prev = hs.rsb_scale(_matched(cs, m, -1, 0.0))
    for s in (0.4, 0.6, 0.8):
        cur = hs.rsb_scale(_matched(cs, m, -1, s))
        assert cur < prev
        prev = cur


def test_dilution_confound_is_visible():
    """The naive 'positive correlation' construction empties a class.

    At large spread the positive family's low class approaches zero
    hyperdegree, so a shift in its RSB point is dilution, not correlation.
    isolated_fraction() makes that visible instead of silent.
    """
    cs, m = [2, 6], [0.75, 0.25]
    lo = _matched(cs, m, +1, 0.2)(hs.rsb_scale(_matched(cs, m, +1, 0.2)))
    hi = _matched(cs, m, +1, 0.95)(hs.rsb_scale(_matched(cs, m, +1, 0.95)))
    assert lo.isolated_fraction() < 0.15
    assert hi.isolated_fraction() > 0.35
