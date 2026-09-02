"""Core percolation as a chygraph fixed point (TODO item 1, step 2)."""

import numpy as np
import pytest

import statmech.core as cp
from statmech.hittingset import poisson_phi

E = np.e


# ---------------------------------------------------------------------------
# The graph case, against published values
# ---------------------------------------------------------------------------

def test_er_threshold_is_e():
    """Bauer & Golinelli: the leaf-removal core appears at c = e."""
    assert cp.core_threshold(lambda k: cp.graph(mean=k)) == pytest.approx(
        E, rel=1e-9)


@pytest.mark.parametrize('c', [0.5, 1.0, 2.0, 2.5, 2.7])
def test_no_core_below_threshold(c):
    m = cp.graph(mean=c)
    assert m.state()[2][0] == pytest.approx(0.0, abs=1e-9)
    assert m.core_fraction() == pytest.approx(0.0, abs=1e-9)
    assert not m.has_core()


@pytest.mark.parametrize('c', [3.0, 3.5, 4.0, 6.0, 10.0])
def test_core_above_threshold(c):
    m = cp.graph(mean=c)
    assert m.has_core()
    assert 0.0 < m.core_fraction() < 1.0


def test_core_fraction_matches_simulation():
    """Values checked against pure leaf removal on n = 4x10^5 Erdos-Renyi
    graphs (probe/validate_core.py); agreement is at the 1e-3 finite-size
    level, so these are pinned a little looser than that."""
    expected = {3.0: 0.3136, 3.5: 0.6099, 4.0: 0.7700, 6.0: 0.9657,
                10.0: 0.9990}
    for c, sim in expected.items():
        assert cp.graph(mean=c).core_fraction() == pytest.approx(sim, abs=3e-3)


# ---------------------------------------------------------------------------
# The structure of the map
# ---------------------------------------------------------------------------

def test_map_is_monotone():
    """Order-preserving, unlike the hitting-set map: solvable by iteration
    upward from zero, exactly as percolation.giant.Chygraph.solve does."""
    for m in (cp.graph(mean=4.0), cp.clique_network(3, 1.0),
              cp.clique_network(4, 0.5)):
        assert (m.jacobian() >= -1e-9).all()


def test_jacobian_has_zero_diagonal_block():
    """lambda depends on delta only; that is the bipartite structure of WP2."""
    m = cp.graph(mean=4.0)
    J = m.jacobian()
    assert J[0, 0] == pytest.approx(0.0, abs=1e-6)


def test_core_free_branch_is_stable_below_and_unstable_above():
    for c in (1.0, 2.0, 2.7):
        assert cp.graph(mean=c).core_free_spectral() < 1.0
    for c in (2.8, 3.5, 6.0):
        assert cp.graph(mean=c).core_free_spectral() > 1.0


@pytest.mark.parametrize('split', [[0.5, 0.5], [0.9, 0.1],
                                   [0.4, 0.35, 0.25], [0.2] * 5])
def test_threshold_is_total_chy_degree_e_however_the_layers_split(split):
    """Poisson layers of cardinality two: rho(A) = 1 at sum_l kappa_l = e.

    A_ml = kappa_l lambda has rank one there, so only the total matters.  The
    single-layer derivative this used to take instead returned e / (number of
    layers) worth of stability and put the threshold at 14.8 for two layers.
    """
    def build(total):
        return cp.CorePercolation([2] * len(split),
                                  poisson_phi([total * f for f in split]))

    assert cp.core_threshold(build) == pytest.approx(E, rel=1e-8)


def test_layers_that_differ_in_shape_move_the_threshold():
    """Pooling is exact for Poisson layers and not in general.

    A matching layer plus a Poisson layer cores earlier than the single graph
    with the same pooled degree distribution: multiplexing correlates the
    degrees at the ends of a matching link, and pooling destroys that.
    """
    import sympy as sp
    from scipy.optimize import brentq
    from statmech.hittingset import layer_symbols

    q = 0.6
    x, (y,) = layer_symbols(2), layer_symbols(1)
    layered = lambda m: cp.CorePercolation(
        [2, 2], ((1 - q) + q * x[0]) * sp.exp(m * (x[1] - 1)))
    pooled = lambda m: cp.CorePercolation(
        [2], ((1 - q) + q * y) * sp.exp(m * (y - 1)))
    f = lambda build: brentq(lambda m: build(m).core_free_spectral() - 1.0,
                             0.5, 4.0, xtol=1e-10)
    assert q + f(layered) == pytest.approx(2.6196, abs=1e-3)
    assert q + f(pooled) == pytest.approx(2.6994, abs=1e-3)


def test_three_state_message_sums_to_one():
    for m in (cp.graph(mean=4.0), cp.clique_network(3, 0.7)):
        lam, delta, gamma = m.state()
        assert (lam + delta + gamma) == pytest.approx(np.ones(m.L), abs=1e-12)
        assert (gamma >= -1e-12).all()


# ---------------------------------------------------------------------------
# The result that matters for prediction 4
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('c', [3, 4, 5])
def test_cardinality_three_or_more_has_no_core_free_branch(c):
    """Every member of a c-clique has degree c-1 >= 2, so leaf removal can
    never reach it: the complex is a core by itself."""
    assert not cp.clique_network(c, 1.0).has_core_free_branch()


@pytest.mark.parametrize('c', [3, 4, 5])
@pytest.mark.parametrize('k', [0.01, 0.05, 0.2, 1.0])
def test_cored_at_every_density(c, k):
    """No threshold at all, however sparse."""
    m = cp.clique_network(c, k)
    assert m.has_core()
    assert m.core_fraction() > 0.0


@pytest.mark.parametrize('k', [0.05, 0.1, 0.5, 1.0, 2.0])
def test_sparse_clique_core_is_membership_probability(k):
    """At low density the core is exactly the nodes in at least one complex,
    1 - exp(-k) for Poisson memberships -- the whole complex survives."""
    for c in (3, 4):
        assert cp.clique_network(c, k).core_fraction() == pytest.approx(
            1 - np.exp(-k), rel=2e-2)


def test_complexes_beat_the_degree_matched_graph():
    """The prediction-4 comparison, analytically.

    A node in `k` triangles has `2k` neighbours.  Give an ordinary graph the
    same mean degree and it has no core at all until 2k > e, while the
    triangle network is cored at every k.
    """
    for k in (0.05, 0.2, 0.5, 1.0):
        tri = cp.clique_network(3, k).core_fraction()
        matched = cp.graph(mean=2 * k).core_fraction()
        assert tri > 0.0
        if 2 * k < E:
            assert matched == pytest.approx(0.0, abs=1e-9)
        assert tri > matched


def test_a_triangle_is_not_automatically_core():
    """Leaf removal deletes a leaf's *neighbour*, so a complex member can go.

    Triangle abc with a pendant on each vertex: d-a fires, removing a; then
    b and c are left with degree 1 and the triangle unravels completely.  The
    'cardinality >= 3' statement is about the core-free branch not existing,
    not about complexes surviving whole.
    """
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path.home() / 'av2atg' / 'computational_complexity' / 'code'))
    import leafremoval as lr
    from hrg import to_csr
    s = np.array([0, 1, 2, 0, 1, 2])
    d = np.array([1, 2, 0, 3, 4, 5])
    assert lr.core(*to_csr(6, s, d))[0] == 0


def test_core_positive_but_tiny_with_a_sparse_triangle_layer():
    """No core-free branch, yet the core can be 1e-4.  Both halves matter."""
    m = cp.CorePercolation([2, 3], poisson_phi([1.0, 0.001]))
    assert not m.has_core_free_branch()
    assert 0.0 < m.core_fraction() < 1e-3


def test_core_threshold_is_the_vertex_cover_instability():
    """On a graph the two are one point, not two.

    Eq. (19) puts the hard-field vertex-cover instability at <k>(c-1) = e, so at
    <k> = e for c = 2; Bauer-Golinelli put the core-percolation threshold at the
    same place.  The paper computes both, and they must agree.
    """
    import statmech.hittingset as _hs
    core_thr = cp.core_threshold(lambda k: cp.graph(mean=k))
    hard_thr = _hs.rsb_point([2], [1.0])
    assert core_thr == pytest.approx(E, rel=1e-9)
    assert hard_thr == pytest.approx(E, rel=1e-9)
    assert core_thr == pytest.approx(hard_thr, rel=1e-8)


@pytest.mark.parametrize('cs', [[2], [3], [4], [2, 3]])
def test_certification_is_the_core_free_branch(cs):
    """The cover size of Eq. (18) can be believed exactly when Eq. (21) admits
    gamma = 0.  Not two conditions that happen to coincide -- one condition."""
    import statmech.cover as _cv
    means = [1.0] * len(cs)
    assert (_cv.poisson(cs, means).certified()
            == cp.CorePercolation(cs, poisson_phi(means)).has_core_free_branch())
