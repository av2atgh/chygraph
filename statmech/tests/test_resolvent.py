"""The resolvent as a message (Chapter 18).

Pinned against the explicit inverse inside a complex, exact inversion on
incidence trees, the closed forms of the clean tree, and the line-graph
identity that ties the (0,2) cactus to the degree-three tree.
"""

import numpy as np
import pytest

from statmech import Chygraph
from statmech import resolvent as R


# ---------------------------------------------------------------------------
# The down step
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('m', [1, 2, 3, 5, 12])
def test_transmit_is_the_schur_complement(m):
    """hop^2 1^T M^{-1} 1 with M = diag(m_J) - hop A_J, A_J the clique."""
    rng = np.random.default_rng(m)
    hop = 0.7
    m_J = rng.normal(size=m) + 1j * rng.uniform(0.1, 1.0, m)
    eta_J = rng.uniform(0.5, 2.0, m)
    M = np.diag(m_J) - hop * (np.ones((m, m)) - np.eye(m))
    one = np.ones(m)
    v = np.linalg.solve(M, one)
    out = R.transmit(m_J, hop, eta_J)
    assert out.sigma == pytest.approx(hop * hop * one @ v, rel=1e-12)
    assert np.allclose(out.v, v, rtol=1e-12)
    assert out.eta == pytest.approx(hop * one @ np.linalg.solve(M, eta_J), rel=1e-12)


def test_transmit_reduces_to_the_tree_and_the_triangle():
    """Eq. (18.9) at one other member, Eq. (18.6) at two."""
    hop, m1, m2 = 1.0, 0.8 - 0.3j, -1.7 + 0.2j
    assert R.transmit([m1], hop).sigma == pytest.approx(hop * hop / m1)
    want = hop * hop * (m1 + m2 + 2 * hop) / (m1 * m2 - hop * hop)
    assert R.transmit([m1, m2], hop).sigma == pytest.approx(want)


def test_transmit_is_vectorised():
    rng = np.random.default_rng(0)
    m_J = rng.normal(size=(50, 3)) + 0.5j
    batched = R.transmit(m_J, 1.0).sigma
    single = np.array([R.transmit(row, 1.0).sigma for row in m_J])
    assert np.allclose(batched, single)


def test_linearised_imaginary_is_the_first_order_map():
    """d Im Sigma_{a->i} / d Im Sigma_{j->a} = hop^2 v_j^2 with real parts frozen."""
    rng = np.random.default_rng(3)
    hop, m_J = 0.9, rng.normal(size=4)
    v = R.transmit(m_J, hop).v
    im = rng.uniform(0, 1e-6, 4)
    exact = R.transmit(m_J - 1j * im, hop).sigma.imag      # Sigma_{j->a} gains +i im
    assert R.linearised_imaginary(v, im, hop) == pytest.approx(exact, rel=1e-4)


# ---------------------------------------------------------------------------
# Instances: exact on an incidence tree
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('degrees,depth', [([3, 0], 5), ([2, 1], 4), ([1, 2], 3), ([0, 2], 3)])
def test_cavity_is_exact_on_an_incidence_tree(degrees, depth):
    rng = np.random.default_rng(1)
    n, cx = R.incidence_tree([2, 3], degrees, depth)
    eps = rng.uniform(-1.5, 1.5, n)
    H = R.hamiltonian(n, cx, eps)
    z = 0.7 - 0.2j
    Gex = np.diag(np.linalg.inv(H - z * np.eye(n)))
    G, _, _ = R.cavity_instance(n, cx, eps, z, sweeps=2 * depth + 5)
    assert np.abs(G - Gex).max() < 1e-12


def test_landscape_is_exact_on_an_incidence_tree():
    rng = np.random.default_rng(2)
    n, cx = R.incidence_tree([2, 3], [1, 1], 4)
    eps = rng.uniform(-1.5, 1.5, n)
    H = R.hamiltonian(n, cx, eps).real
    emin = np.linalg.eigvalsh(H).min() - 0.05
    uex = np.linalg.solve(H - emin * np.eye(n), np.ones(n))
    u = R.landscape_instance(n, cx, eps, emin, sweeps=15)
    assert np.abs(u / uex - 1).max() < 1e-10


def test_regular_instance_has_the_right_degrees():
    rng = np.random.default_rng(5)
    n = 300
    cx = R.regular_instance(n, [2, 3], [2, 1], rng)
    edges = sum(1 for c in cx if len(c) == 2)
    tris = sum(1 for c in cx if len(c) == 3)
    assert edges == n and tris == n // 3 * 1 or abs(edges - n) <= 3
    assert all(len(set(c)) == len(c) for c in cx)


# ---------------------------------------------------------------------------
# The clean ensemble
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('degrees,want', [
    ([3, 0], -2 * np.sqrt(2)),          # tree of degree three
    ([4, 0], -2 * np.sqrt(3)),          # tree of degree four
    ([0, 2], -(2 * np.sqrt(2) + 1)),    # line graph of the degree-three tree
])
def test_band_bottom_closed_forms(degrees, want):
    assert R.band_bottom([2, 3], degrees) == pytest.approx(want, abs=2e-4)


def test_spectral_bottom_switches_to_the_isolated_eigenvalue():
    """E_min = e_min - W/2 above W_min = 2(degree + e_min), the uniform vector below."""
    emin = R.band_bottom([2, 3], [3, 0])
    Wmin = 2 * (3 + emin)
    assert R.spectral_bottom([2, 3], [3, 0], 0.5 * Wmin) == pytest.approx(-3.0)
    assert R.spectral_bottom([2, 3], [3, 0], 2 * Wmin) == pytest.approx(emin - Wmin)


def test_clean_tree_matches_kesten_mckay():
    K = 2
    for E in (0.0, 1.0, 2.5):
        _, Gii = R.clean_tree(K, E + 1e-9j)
        assert Gii.imag / np.pi == pytest.approx(R.kesten_mckay(K, E), rel=1e-6)
    assert R.kesten_mckay(K, 3.0) == 0.0
    # normalised to one state per site
    E = np.linspace(-2 * np.sqrt(K), 2 * np.sqrt(K), 200001)
    assert np.trapz(R.kesten_mckay(K, E), E) == pytest.approx(1.0, abs=1e-4)


def test_clean_tree_is_the_uniform_fixed_point_below_the_band():
    """Below the band the recursion's real fixed point is the real branch of Eq. (18.3)."""
    K, E = 2, -3.5
    sig = R.uniform_fixed_point([2, 3], [K + 1, 0], E)
    g, _ = R.clean_tree(K, E)
    assert sig[0] == pytest.approx(g.real, rel=1e-8)


# ---------------------------------------------------------------------------
# The incidence tree's branching, Table 18.1
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('degrees,want', [
    ([3, 0], 2.0), ([1, 1], np.sqrt(2)), ([4, 0], 3.0),
    ([2, 1], (1 + np.sqrt(17)) / 2), ([0, 2], 2.0)])
def test_incidence_branching(degrees, want):
    B = R.incidence_branching([2, 3], degrees)
    assert np.max(np.abs(np.linalg.eigvals(B))) == pytest.approx(want, rel=1e-12)


# ---------------------------------------------------------------------------
# The handle
# ---------------------------------------------------------------------------

def test_handle_delegates():
    g = Chygraph([2, 3], [1, 1], regular=True)
    assert np.allclose(g.incidence_branching(), R.incidence_branching([2, 3], [1, 1]))
    assert g.band_bottom() == pytest.approx(R.band_bottom([2, 3], [1, 1]))
    assert g.spectral_bottom(3.0) == pytest.approx(R.spectral_bottom([2, 3], [1, 1], 3.0))
    cx = g.instance(120, seed=0)
    assert all(len(c) in (2, 3) for c in cx)
    with pytest.raises(NotImplementedError):
        Chygraph([2], [3.0]).instance(10)
