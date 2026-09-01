"""WP4 checks.  The expensive one lives in examples/wp4_validates_wp1.py."""

import numpy as np
import pytest
import sympy as sp

import statmech as cs
from statmech.population import CavityPopulation

beta, J = sp.symbols('beta J', positive=True)


# ---------------------------------------------------------------------------
# The bridge to WP1: same complex, symbolic vs numeric
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('bJ', [0.2, 0.5, 1.0, 2.0])
@pytest.mark.parametrize('c,symbolic', [
    (2, cs.ising_edge_derivative), (3, cs.ising_triangle_derivative)])
def test_emitted_derivative_matches_cavity(c, symbolic, bJ):
    """d(emitted)/dh at zero field is exactly WP1's u'.

    WP4's numeric enumeration and WP1's sympy enumeration are independent
    implementations of the same complex, so this pins one against the other.
    """
    p = CavityPopulation([c], [1.0], bJ)
    eps = 1e-6
    h = np.zeros((2, c - 1))
    h[0, 0], h[1, 0] = eps, -eps
    out = p.emitted(h, 0)
    numeric = (out[0] - out[1]) / (2 * eps)
    exact = float(symbolic().subs({beta: bJ, J: 1}))
    assert numeric == pytest.approx(exact, abs=1e-8)


@pytest.mark.parametrize('bJ', [0.3, 1.0])
def test_edge_emission_is_the_textbook_form(bJ):
    """c = 2 must be atanh(tanh(bJ) tanh(h)) at every field, not just at 0."""
    p = CavityPopulation([2], [1.0], bJ)
    h = np.array([[-2.0], [-0.5], [0.0], [0.7], [3.0]])
    want = np.arctanh(np.tanh(bJ) * np.tanh(h[:, 0]))
    assert p.emitted(h, 0) == pytest.approx(want, abs=1e-12)


def test_triangle_emission_exceeds_two_edges_at_zero_field():
    """The WP1 claim, in the field representation: a triangle transmits more
    per neighbour than an independent edge."""
    for bJ in (0.2, 0.5, 1.0):
        eps = 1e-6
        edge = CavityPopulation([2], [1.0], bJ)
        tri = CavityPopulation([3], [1.0], bJ)
        h2 = np.array([[eps]])
        h3 = np.array([[eps, 0.0]])
        assert tri.emitted(h3, 0)[0] > edge.emitted(h2, 0)[0]


# ---------------------------------------------------------------------------
# The map itself
# ---------------------------------------------------------------------------

def test_zero_field_is_a_fixed_point():
    """The paramagnet is exact at zero external field, which is why the runs
    have to be initialised magnetised."""
    for c in (2, 3):
        p = CavityPopulation([c], [4.0], 0.5, size=2000, seed=0)
        p.initialise(0.0)
        p.sweep()
        assert np.abs(p.P[0]).max() == pytest.approx(0.0, abs=1e-12)


def test_convolution_preserves_mean():
    """The chy-degree step is a Poisson-compound sum: its mean must be
    <k> times the mean of the down-message."""
    p = CavityPopulation([2], [3.0], 0.5, size=200_000, seed=3)
    p.initialise(1.0)
    p.Q[0] = np.full(p.size, 0.25)
    got = float(np.mean(p._sum_over_complexes(200_000)))
    assert got == pytest.approx(3.0 * 0.25, rel=0.02)


# ---------------------------------------------------------------------------
# The transition, cheaply
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('spec,bJc', [
    (([2], [6.0]), np.arctanh(1 / 6.0)),
    (([3], [3.0]), np.arctanh((7 - np.sqrt(45)) / 2)),
])
def test_orders_below_and_not_above(spec, bJc):
    """WP1's threshold, bracketed rather than bisected: no magnetisation well
    above T_c, clear magnetisation well below."""
    hot = CavityPopulation(*spec, 0.7 * bJc, size=40_000, seed=1)
    cold = CavityPopulation(*spec, 1.5 * bJc, size=40_000, seed=1)
    assert hot.run(200).magnetisation() == pytest.approx(0.0, abs=1e-3)
    assert cold.run(200).magnetisation() > 0.4


def test_triangles_order_before_the_matched_graph():
    """At six neighbours either way, the triangle network magnetises at a
    coupling where the graph is still paramagnetic.  This is the 14.5% T_c
    gain of WP1, as a strict inequality that needs no bisection."""
    bJ = np.arctanh(1 / 6.0) * 0.93        # between the two thresholds
    graph = CavityPopulation([2], [6.0], bJ, size=60_000, seed=1)
    tri = CavityPopulation([3], [3.0], bJ, size=60_000, seed=1)
    assert graph.run(250).magnetisation() == pytest.approx(0.0, abs=1e-3)
    assert tri.run(250).magnetisation() > 0.05
