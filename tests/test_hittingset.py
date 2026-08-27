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


# ---------------------------------------------------------------------------
# Known limitations, pinned against the literature
# ---------------------------------------------------------------------------

def test_cover_size_is_certified_only_at_cardinality_two():
    assert hs.poisson([2], [3.0]).certified()
    for c in (3, 4, 6):
        assert not hs.poisson([c], [1.0]).certified()


def test_disjoint_triples_expose_the_degeneracy_rule():
    """One 3-hyperedge per vertex and nothing else: no interaction between
    complexes, so replica symmetry is trivially valid.  One of every three
    vertices must be taken, so the truth is 1/3.  The graph degeneracy rule
    returns 1/2.  The bracket is honest; the midpoint is not."""
    from sympy import Symbol
    m = hs.HittingSet([3], Symbol('x0'))
    assert m.cover_size() == pytest.approx(0.5, abs=1e-9)
    lo, hi = m.cover_bracket()
    assert lo <= 1 / 3 <= hi
    assert abs(m.cover_size() - 1 / 3) > 0.15


def test_against_mezard_tarzia_regular_benchmark():
    """Mezard & Tarzia, PRE 76, 041124 (2007), Fig. 5: a random regular
    hypergraph with L = 4 tests per variable and K = 6 variables per test has
    rho_cov = 0.178.  The hard-field ansatz here returns 0.252, and its own
    stability test already refuses to certify that value."""
    from sympy import Symbol
    m = hs.HittingSet([6], Symbol('x0') ** 4)
    assert m.cover_size() == pytest.approx(0.252, abs=2e-3)
    assert not m.certified()
    assert m.is_unstable()
    lo, hi = m.cover_bracket()
    assert lo <= 0.178 <= hi


def test_cardinality_two_regular_breaks_for_every_degree():
    """Mezard & Tarzia: 'for K = 2 the solution of VC exhibits higher order RSB
    for every value of L'.  The hard-field criterion agrees from L = 3 up."""
    from sympy import Symbol
    x = Symbol('x0')
    for L in (3, 4, 6, 10):
        assert hs.HittingSet([2], x ** L).is_unstable()


def test_mixed_cardinality_derivation_in_the_text():
    """The worked example given in the manuscript for the three-to-one mixture
    of c=2 and c=6: sigma = 0.3777, <k> = 3.415, <k> nu = 6.83."""
    import numpy as np
    from scipy.optimize import brentq
    kL, kT, cs = 0.75, 0.25, (2, 6)

    def cond(scale):
        ks = (scale * kL, scale * kT)
        s = brentq(lambda x: np.exp(-sum(k * x**(c - 1)
                                         for k, c in zip(ks, cs))) - x, 1e-12, 1.0)
        return sum(k * (c - 1) * s**(c - 1) for k, c in zip(ks, cs)) - 1, s

    scale = brentq(lambda x: cond(x)[0], 1e-6, 50.0)
    _, sigma = cond(scale)
    nu = (kL * 1 + kT * 5) / (kL + kT)
    assert sigma == pytest.approx(0.3777, abs=5e-4)
    assert scale == pytest.approx(3.415, abs=5e-3)
    assert scale * nu == pytest.approx(6.83, abs=1e-2)


def test_fixed_cardinality_derivation_in_the_text():
    """<k> = e/(c-1) follows from y = 1/(c-1) and <k> = y exp((c-1)y)."""
    import numpy as np
    for c in (2, 3, 5, 10):
        y = 1.0 / (c - 1)
        assert y * np.exp((c - 1) * y) == pytest.approx(np.e / (c - 1), rel=1e-12)
        assert hs.rsb_point([c], [1.0]) == pytest.approx(np.e / (c - 1), rel=1e-10)
