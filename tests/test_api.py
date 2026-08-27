"""The consolidated Chygraph class.

Every method must agree with the module it delegates to; the point of the class
is one handle, not a second implementation.
"""

import numpy as np
import pytest
from sympy import Symbol

import chygraph_statmech.core as _core
import chygraph_statmech.cover as _cover
import chygraph_statmech.freeenergy as _fe
import chygraph_statmech.hittingset as _hs
import chygraph_statmech.ising as _ising
import chygraph_statmech.simplicial as _simp
from chygraph_statmech import Chygraph
from chygraph_statmech.cavity import ising_edge_derivative

E = np.e


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

def test_poisson_excess_defaults_to_the_mean():
    g = Chygraph([2, 3], [4.0, 2.0])
    assert np.allclose(g.kbar, g.k)


def test_regular_excess_is_one_less():
    g = Chygraph([6], [4], regular=True)
    assert np.allclose(g.kbar, [3.0])


def test_shape_and_cardinality_are_checked():
    with pytest.raises(ValueError, match='one cardinality'):
        Chygraph([2, 3], [1.0])
    with pytest.raises(ValueError, match='at least 2'):
        Chygraph([1], [1.0])


# ---------------------------------------------------------------------------
# Sec. IV, Ising
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('spec', [([2], [6.0]), ([3], [3.0]), ([2, 3], [4.0, 2.0])])
def test_critical_coupling_delegates(spec):
    g = Chygraph(*spec)
    assert g.critical_coupling() == pytest.approx(
        _ising.critical_coupling(g.c, g.k, excess=g.kbar), rel=1e-12)


def test_graph_critical_coupling_is_the_textbook_value():
    assert Chygraph([2], [6.0]).critical_coupling() == pytest.approx(
        np.arctanh(1 / 6.0), rel=1e-9)


def test_critical_temperature_is_the_reciprocal():
    g = Chygraph([2], [6.0])
    assert g.critical_temperature() == pytest.approx(
        1 / g.critical_coupling(), rel=1e-12)


def test_branching_matrix_matches_the_module():
    g = Chygraph([2, 3], [4.0, 2.0])
    assert np.allclose(g.branching_matrix(0.2),
                       _ising.branching_matrix(g.c, g.k, 0.2, excess=g.kbar))


def test_u_prime_selects_the_interaction():
    g = Chygraph([3], [2.0])
    assert g.u_prime(0, 0.8) == pytest.approx(
        _ising.clique_derivative(3, 0.8), abs=1e-12)
    assert g.u_prime(0, 0.8, 'simplicial') == pytest.approx(
        _simp.uprime(3, 0.8), abs=1e-12)
    assert Chygraph([2], [1.0]).u_prime(0, 0.5) == pytest.approx(
        float(ising_edge_derivative().subs(
            {Symbol('beta', positive=True): 0.5,
             Symbol('J', positive=True): 1})), abs=1e-7)


def test_simplicial_needs_its_own_spinodal():
    with pytest.raises(NotImplementedError, match='simplicial'):
        Chygraph([3], [4.0]).critical_coupling(interaction='simplicial')


def test_paramagnetic_free_energy_delegates():
    g = Chygraph([2, 3], [4.0, 2.0])
    assert g.paramagnetic_free_energy(0.1) == pytest.approx(
        _fe.paramagnetic(g.c, g.k, 0.1), abs=1e-14)


# ---------------------------------------------------------------------------
# Sec. V, hard core
# ---------------------------------------------------------------------------

def test_hitting_set_delegates():
    g = Chygraph([2], [3.0])
    assert g.hitting_set().cover_size() == pytest.approx(
        _hs.poisson([2], [3.0]).cover_size(), abs=1e-12)


def test_clique_cover_delegates():
    g = Chygraph([2], [3.0])
    assert g.clique_cover().cover_size() == pytest.approx(
        _cover.poisson([2], [3.0]).cover_size(), abs=1e-12)


def test_hitting_set_bp_carries_the_regular_flag():
    g = Chygraph([6], [4], regular=True)
    m = g.hitting_set_bp(mu=60.0, size=40_000, seed=1, damping=0.5).run(500)
    assert m.regular
    assert m.density() == pytest.approx(1 / 6, abs=1e-4)


# ---------------------------------------------------------------------------
# Sec. VI, core percolation
# ---------------------------------------------------------------------------

def test_core_delegates_and_recovers_e():
    g = Chygraph([2], [1.0])
    assert g.core().core_fraction() == pytest.approx(
        _core.clique_network(2, 1.0).core_fraction(), abs=1e-9)
    assert _core.core_threshold(
        lambda k: Chygraph([2], [k]).core()) == pytest.approx(E, rel=1e-9)


def test_from_samples_uses_the_measured_ensemble():
    rng = np.random.default_rng(0)
    K = rng.poisson(1.0, (200_000, 1)).astype(float)
    g = Chygraph.from_samples([3], K)
    assert g.k[0] == pytest.approx(1.0, rel=0.02)
    assert g.core_from_samples().core_fraction() == pytest.approx(
        1 - np.exp(-1.0), rel=0.02)


def test_from_samples_is_required_for_the_empirical_core():
    with pytest.raises(ValueError, match='from_samples'):
        Chygraph([3], [1.0]).core_from_samples()


# ---------------------------------------------------------------------------
# Sec. IV C and VIII
# ---------------------------------------------------------------------------

def test_simplicial_delegates():
    g = Chygraph([2, 16], [4, 4], regular=True)
    M = g.simplicial([0.7, 0.3])
    assert isinstance(M, _simp.SimplicialChygraph)
    assert M.spinodal() == pytest.approx(1.0108, abs=1e-3)


def test_regions_and_overlap_profile():
    g = Chygraph([3], [1.0])
    rg = g.regions([[0, 1, 2], [2, 3, 4]])
    assert rg.is_bethe() and rg.counting_is_valid()
    assert Chygraph.overlap_profile([[0, 1, 2], [1, 2, 3]])['shared_2plus'] == 1.0


def test_emitted_field_is_reachable_for_any_interior():
    """Eq. (3) is exposed so a new interaction needs no new class."""
    from sympy import Symbol as S
    h = [S('h1')]
    u = Chygraph.emitted_field(2, lambda s: 0.5 * s[0] * s[1], h)
    assert u.subs({h[0]: 0}) == 0


# ---------------------------------------------------------------------------
# The methods added so every module has a manuscript entry
# ---------------------------------------------------------------------------

def test_fixed_point_stability_gives_the_exchange_relation():
    """rho(J(Q*)) = 1 - Lambda + O(Lambda^2), Sec. II D."""
    import sympy as sp
    from chygraph import hypergraph_giant
    M = hypergraph_giant()
    S = Chygraph([2], [1.0]).fixed_point_stability(M)
    k, c, p, q = sp.symbols('k c p q')
    for pv, tol in ((0.145, 2e-2), (0.140, 3e-3)):
        sub = {k: 3, c: 3, p: pv, q: 0.8}
        lam = float(M.Lambda().subs(sub))
        assert abs((1 - S.spectral_radius(sub)) / lam - 1) < tol
    assert S.monotonicity({k: 3, c: 3, p: 0.5, q: 0.8}) == 1


def test_vertex_cover_correlated_reproduces_weigt_hartmann_shape():
    """Sec. V B: x_c rises with r, and gamma=2.5 lies below gamma=3.0."""
    prev = {2.5: -1.0, 3.0: -1.0}
    for r in (0.0, 0.3, 0.6):
        xs = {}
        for gamma in (2.5, 3.0):
            xs[gamma], _ = Chygraph.vertex_cover_correlated(gamma, r, dmax=300)
            assert xs[gamma] > prev[gamma]
            prev[gamma] = xs[gamma]
        assert xs[2.5] < xs[3.0]
    assert not Chygraph.vertex_cover_correlated(2.5, 0.0, dmax=300)[1]


def test_excess_cardinality_matches_the_definition():
    """<c^2>/<c> - 1 on an explicit complex list (Table III)."""
    assert Chygraph.excess_cardinality([[0, 1], [2, 3]]) == pytest.approx(1.0)
    assert Chygraph.excess_cardinality([[0, 1, 2], [2, 3]]) == pytest.approx(
        (9 + 4) / 2 / ((3 + 2) / 2) - 1)
