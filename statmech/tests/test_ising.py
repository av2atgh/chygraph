"""Critical temperature on an arbitrary chygraph."""

import numpy as np
import pytest
import sympy as sp

import statmech as cs
from statmech.ising import (branching_matrix, clique_derivative,
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


# ---------------------------------------------------------------------------
# The order parameter and the response
# ---------------------------------------------------------------------------

from statmech.ising import (emitted_common, emitted_common_derivative,
                            interior_cumulants, magnetisation,
                            magnetisation_amplitude, susceptibility,
                            susceptibility_amplitude)


@pytest.mark.parametrize('c', [2, 3, 4, 5])
@pytest.mark.parametrize('bJ', [0.1, 0.4, 0.9])
def test_a1_is_the_transmission(c, bJ):
    """a1 = (c-1) u': the amplitude expansion starts where the threshold ends."""
    a1, _ = interior_cumulants(c, bJ)
    assert a1 == pytest.approx((c - 1) * clique_derivative(c, bJ), rel=1e-12)


@pytest.mark.parametrize('bJ', [0.2, 0.7, 1.3])
def test_a3_closed_form_on_an_edge(bJ):
    """ubar = atanh(t tanh h) on a link, so a3 = t(t^2-1)/3."""
    t = np.tanh(bJ)
    _, a3 = interior_cumulants(2, bJ)
    assert a3 == pytest.approx(t * (t * t - 1) / 3, abs=1e-13)


@pytest.mark.parametrize('c', [2, 3, 4])
@pytest.mark.parametrize('bJ', [0.15, 0.6])
def test_cumulants_are_the_taylor_coefficients(c, bJ):
    """a1 and a3 against finite differences of the enumeration itself."""
    a1, a3 = interior_cumulants(c, bJ)
    e = 1e-3
    d1 = (emitted_common(c, bJ, e) - emitted_common(c, bJ, -e)) / (2 * e)
    d3 = (emitted_common(c, bJ, 2 * e) - 2 * emitted_common(c, bJ, e)
          + 2 * emitted_common(c, bJ, -e) - emitted_common(c, bJ, -2 * e))
    assert float(d1) == pytest.approx(a1, abs=1e-5)
    assert float(d3) / (2 * e ** 3 * 6) == pytest.approx(a3, abs=1e-4)


@pytest.mark.parametrize('c', [2, 3, 4])
@pytest.mark.parametrize('h', [0.0, 0.3, 1.5])
def test_emitted_derivative_is_exact(c, h):
    e = 1e-5
    fd = (emitted_common(c, 0.45, h + e) - emitted_common(c, 0.45, h - e)) / (2 * e)
    assert emitted_common_derivative(c, 0.45, h) == pytest.approx(float(fd),
                                                                 abs=1e-8)


@pytest.mark.parametrize('k', [3, 4, 6])
@pytest.mark.parametrize('bJ', [0.05, 0.15, 0.25])
def test_susceptibility_bethe_lattice(k, bJ):
    """chi = (1+t)/(1-(k-1)t) on a regular graph -- the textbook result."""
    t = np.tanh(bJ)
    got = susceptibility([2], [float(k)], bJ, excess=[k - 1.0])
    assert got == pytest.approx((1 + t) / (1 - (k - 1) * t), rel=1e-11)


@pytest.mark.parametrize('k', [2.0, 3.0, 5.0])
@pytest.mark.parametrize('bJ', [0.05, 0.15])
def test_susceptibility_poisson_graph(k, bJ):
    """chi = 1/(1 - <k> t), a Poisson degree being its own excess."""
    t = np.tanh(bJ)
    assert susceptibility([2], [k], bJ) == pytest.approx(1 / (1 - k * t),
                                                         rel=1e-11)


@pytest.mark.parametrize('spec', [([2], [4.0], [3.0]), ([3], [2.0], [1.0]),
                                  ([2, 3], [2.0, 1.0], [1.0, 0.0]),
                                  ([2, 3], [3.0, 1.0], None)])
def test_susceptibility_diverges_at_the_transition(spec):
    """1/chi -> 0 exactly where det(I - B) = 0.  One condition, two readings."""
    cards, k, kbar = spec
    bc = critical_coupling(cards, k, excess=kbar)
    prev = np.inf
    for eps in (1e-3, 1e-4, 1e-5):
        inv = 1.0 / susceptibility(cards, k, bc * (1 - eps), excess=kbar)
        assert 0 < inv < prev
        prev = inv
    assert prev < 1e-4


@pytest.mark.parametrize('spec', [([2], [4.0]), ([3], [2.0]), ([2, 3], [2.0, 1.0])])
def test_susceptibility_amplitude_is_the_limit(spec):
    """chi (T/T_c - 1) -> C, so gamma = 1."""
    cards, k = spec
    kbar = [x - 1.0 for x in k]
    C = susceptibility_amplitude(cards, k, excess=kbar)
    Tc = 1.0 / critical_coupling(cards, k, excess=kbar)
    got = susceptibility(cards, k, 1.0 / (Tc * (1 + 1e-4)), excess=kbar) * 1e-4
    assert got == pytest.approx(C, rel=2e-3)


@pytest.mark.parametrize('spec', [([2], [4.0]), ([3], [2.0]), ([2], [6.0]),
                                  ([4], [2.0]), ([2, 3], [2.0, 1.0])])
def test_magnetisation_amplitude_is_the_limit(spec):
    """m / sqrt(1 - T/T_c) -> A, so beta = 1/2."""
    cards, k = spec
    A = magnetisation_amplitude(cards, k)
    Tc = 1.0 / critical_coupling(cards, k, excess=[x - 1.0 for x in k])
    for eps, rel in ((1e-4, 3e-3), (1e-6, 2e-4)):
        got = magnetisation(cards, k, 1.0 / (Tc * (1 - eps))) / np.sqrt(eps)
        assert got == pytest.approx(A, rel=rel)


@pytest.mark.parametrize('k', [3, 4, 6, 10])
def test_graph_amplitude_closed_form(k):
    """A = k/(k-1) sqrt(3(k-1)J/T_c) on a k-regular graph."""
    Tc = 1.0 / critical_coupling([2], [float(k)], excess=[k - 1.0])
    want = k / (k - 1) * np.sqrt(3 * (k - 1) / Tc)
    assert magnetisation_amplitude([2], [float(k)]) == pytest.approx(want,
                                                                     rel=1e-8)


def test_amplitude_tends_to_the_curie_weiss_value():
    """A -> sqrt(3) as the degree grows: mean field, from outside the book."""
    prev = np.inf
    for k in (10.0, 100.0, 1000.0):
        A = magnetisation_amplitude([2], [k])
        assert A > np.sqrt(3) and A < prev
        prev = A
    assert prev == pytest.approx(np.sqrt(3), abs=2e-3)


def test_magnetisation_vanishes_above_the_transition():
    for cards, k in (([2], [4.0]), ([3], [2.0]), ([2, 3], [2.0, 1.0])):
        Tc = 1.0 / critical_coupling(cards, k, excess=[x - 1.0 for x in k])
        assert magnetisation(cards, k, 1.0 / (Tc * 1.02)) == 0.0
        assert magnetisation(cards, k, 1.0 / (Tc * 0.98)) > 0.0


def test_magnetisation_matches_the_population():
    """The scalar closure against a population of samples, regular chy-degree.

    A regular chygraph collapses the population onto one value, so this checks
    the excess bookkeeping of both routes rather than the physics of either.
    """
    for cards, k in (([2], [4.0]), ([3], [2.0])):
        g = cs.Chygraph(cards, k, regular=True)
        Tc = g.critical_temperature()
        for frac in (0.6, 0.9):
            bJ = 1.0 / (Tc * frac)
            pop = g.population(bJ, size=40_000, seed=1).run(300).magnetisation()
            assert g.magnetisation(bJ) == pytest.approx(pop, abs=1e-9)


@pytest.mark.parametrize('spec', [([2], [4.0]), ([3], [2.0]), ([4], [2.0]),
                                  ([2, 3], [2.0, 1.0])])
@pytest.mark.parametrize('frac', [1.05, 1.3, 2.0])
def test_susceptibility_is_the_field_derivative_of_the_magnetisation(spec, frac):
    """chi = dm/dB, the two halves of Sec. 9.4 against each other.

    One is a linearisation in closed form and the other a finite difference of
    the full closure at a small field; nothing but the definition connects them.
    """
    cards, k = spec
    kbar = [x - 1.0 for x in k]
    bJ = critical_coupling(cards, k, excess=kbar) / frac
    e = 1e-5
    fd = (magnetisation(cards, k, bJ, field=e)
          - magnetisation(cards, k, bJ, field=-e)) / (2 * e)
    assert susceptibility(cards, k, bJ, excess=kbar) == pytest.approx(fd,
                                                                     rel=1e-5)


def test_magnetisation_follows_the_field_it_is_given():
    """m rises with B either side of the transition, and is odd above it.

    Below the transition it is not: started from saturation the solver returns
    the branch continuously connected to it, which at a negative field is the
    metastable state a sweep would follow, not the equilibrium one.
    """
    cards, k = [2], [4.0]
    Tc = 1.0 / critical_coupling(cards, k, excess=[3.0])
    for frac in (0.8, 1.2):
        bJ = 1.0 / (Tc * frac)
        row = [magnetisation(cards, k, bJ, field=B)
               for B in (0.0, 0.01, 0.05, 0.2)]
        assert all(b > a for a, b in zip(row, row[1:]))
    above = 1.0 / (Tc * 1.2)
    assert magnetisation(cards, k, above, field=-0.05) == pytest.approx(
        -magnetisation(cards, k, above, field=0.05), abs=1e-9)
    below = 1.0 / (Tc * 0.8)
    assert magnetisation(cards, k, below, field=-0.05) > 0.0


def test_clustering_raises_the_amplitude_it_lowers_the_temperature():
    """At degree four: T_c falls by 13.9 per cent, A rises by 20 per cent."""
    graph = cs.Chygraph([2], [4.0], regular=True)
    triangles = cs.Chygraph([3], [2.0], regular=True)
    assert triangles.critical_temperature() < graph.critical_temperature()
    assert triangles.magnetisation_amplitude() > graph.magnetisation_amplitude()
    assert triangles.susceptibility_amplitude() > graph.susceptibility_amplitude()


@pytest.mark.parametrize('frac', [1.1, 1.5, 2.5])
def test_susceptibility_takes_a_cardinality_distribution(frac):
    """A mixed layer against the same layer split, and against its mean.

    The two routes of Sec. 8.3, one order up: splitting a mixed-cardinality
    layer into one layer per cardinality with chy-degrees in the ratio c p_c
    must give the same chi, and substituting the mean cardinality must not.
    """
    mixed = cs.Chygraph([{3: 0.5, 5: 0.5}], [4.0])
    split = cs.Chygraph([3, 5], [1.5, 2.5])
    naive = cs.Chygraph([4], [4.0])
    bJ = mixed.critical_coupling() / frac
    assert mixed.susceptibility(bJ) == pytest.approx(split.susceptibility(bJ),
                                                     rel=1e-12)
    assert abs(naive.susceptibility(bJ) - mixed.susceptibility(bJ)) > 0.05


def test_mixed_cardinality_is_refused_for_the_closure():
    """A layer of two sizes emits two fields, so the closure does not apply."""
    with pytest.raises(NotImplementedError, match='cardinalit'):
        cs.Chygraph([{3: 0.5, 5: 0.5}], [4.0], regular=True).magnetisation(0.1)


def test_scalar_closure_is_refused_off_the_regular_ensemble():
    g = cs.Chygraph([2], [4.0])
    with pytest.raises(NotImplementedError, match='regular'):
        g.magnetisation(0.5)
