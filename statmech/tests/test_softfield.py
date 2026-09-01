"""Hitting set with the O(1) cavity fields kept.

Pinned against Mezard & Tarzia, Phys. Rev. E 76, 041124 (2007), and against
the cases the hard-field ansatz of :mod:`statmech.hittingset` gets
wrong.
"""

import numpy as np
import pytest
from scipy.special import lambertw
from sympy import Symbol

import statmech.hittingset as hs
from statmech.softfield import (HittingSetBP, regular_density,
                                         regular_entropy, regular_field)

MU = 60.0


def _run(cs, deg, regular=False, mu=MU, size=40_000, sweeps=500, seed=1):
    return HittingSetBP(cs, deg, regular=regular, mu=mu, size=size,
                        seed=seed, damping=0.5).run(sweeps)


# ---------------------------------------------------------------------------
# Mezard-Tarzia, regular hypergraphs
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('L,K', [(1, 3), (2, 6), (4, 6), (3, 4), (2, 3)])
def test_regular_field_matches_mezard_tarzia(L, K):
    """h_RS = -mu/L - ((L-1)/L) ln(K-1), their Eq. (11)."""
    m = _run([K], [L], regular=True)
    assert m.P[0].mean() == pytest.approx(regular_field(L, K, MU), abs=1e-4)


@pytest.mark.parametrize('L,K', [(1, 3), (2, 6), (4, 6), (3, 4), (2, 3)])
def test_regular_density_is_one_over_cardinality(L, K):
    """rho = 1/K: every complex holds exactly one member of the set."""
    m = _run([K], [L], regular=True)
    assert m.density() == pytest.approx(regular_density(K), abs=1e-4)


def test_entropy_signs_match_the_paper():
    """MT report positive entropy at (L,K) = (2,6) and negative at (6,12)."""
    assert regular_entropy(2, 6) > 0
    assert regular_entropy(6, 12) < 0
    assert regular_entropy(2, 6) == pytest.approx(
        (5 * np.log(5) - 4 * np.log(6)) / 6, rel=1e-12)


# ---------------------------------------------------------------------------
# Reduction to the graph, where hard and soft fields must agree
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('k', [0.5, 1.0, 2.0, 2.5])
def test_poisson_graph_returns_weigt_hartmann(k):
    """K = 2 is vertex cover, where the hard-field limit is already exact, so
    the soft treatment has to reproduce it -- and with it the exact leaf-removal
    cover size, since the core is empty below k = e."""
    W = lambertw(k).real
    want = 1 - (2 * W + W**2) / (2 * k)
    got = _run([2], [k], mu=50.0, size=200_000, sweeps=600).density()
    assert got == pytest.approx(want, abs=3e-3)
    assert hs.poisson([2], [k]).cover_size() == pytest.approx(want, abs=1e-9)


# ---------------------------------------------------------------------------
# The cases the hard-field ansatz gets wrong
# ---------------------------------------------------------------------------

def test_disjoint_triples_are_now_right():
    """One 3-hyperedge per vertex: the truth is 1/3.  The hard-field rule
    returns 1/2; keeping the O(1) fields returns 1/3."""
    soft = _run([3], [1], regular=True).density()
    hard = 1 - hs.HittingSet([3], Symbol('x0')).cover_size()
    assert soft == pytest.approx(1 / 3, abs=1e-4)
    assert abs(hard - 1 / 3) > 0.15


def test_regular_benchmark_moves_the_right_way():
    """MT Fig. 5, L = 4 K = 6: rho_cov = 0.178 by 1RSB.  The RS answer is
    1/K = 0.1667, below it as an RS answer should be when the RS entropy is
    negative; the hard-field answer is 0.252, above it and further away."""
    soft = _run([6], [4], regular=True).density()
    assert soft == pytest.approx(1 / 6, abs=1e-4)
    assert regular_entropy(4, 6) < 0
    assert abs(soft - 0.178) < abs(0.252 - 0.178)


# ---------------------------------------------------------------------------
# The bug this module was written around
# ---------------------------------------------------------------------------

def test_damping_replaces_population_and_does_not_average():
    """Averaging field values contracts the distribution and moves the fixed
    point; replacing a subset does not.  With averaging the K = 2 Poisson case
    returns ~0.172 where Weigt-Hartmann and exact leaf removal give 0.272."""
    m = _run([2], [1.0], mu=50.0, size=100_000, sweeps=400)
    old = np.full(m.size, -1.0)
    new = np.full(m.size, 1.0)
    mixed = m._mix(old, new)
    assert set(np.unique(mixed)) <= {-1.0, 1.0}      # no intermediate values
    assert 0.3 < float((mixed == 1.0).mean()) < 0.7  # about half replaced


# ---------------------------------------------------------------------------
# The validity criterion
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('L,K', [(1, 3), (2, 6), (4, 6), (3, 4), (2, 4)])
def test_bethe_entropy_matches_mezard_tarzia(L, K):
    """The general Bethe form reproduces their Eq. (13) on regular ensembles.

    Exact there because symmetry removes the sampling: every field is equal.
    """
    m = _run([K], [L], regular=True, size=60_000)
    assert m.entropy() == pytest.approx(regular_entropy(L, K), abs=1e-5)


def test_entropy_is_sufficient_not_necessary():
    """A negative entropy proves the RS answer wrong; a positive one does not
    prove it right.  Vertex cover on Erdos-Renyi breaks replica symmetry at
    mean degree e, yet the entropy is still positive at 1 and at 3.
    """
    m = _run([2], [1.0], mu=20.0, size=100_000, sweeps=300)
    s, err = m.entropy_averaged(keep=100)
    assert s > 0 and err < 0.02
    assert hs.poisson([2], [3.0]).is_unstable()      # RS broken at 3 > e


def test_entropy_averaging_reduces_the_error():
    m = _run([2], [1.0], mu=20.0, size=100_000, sweeps=300)
    _, err = m.entropy_averaged(keep=100)
    assert err < 0.05
