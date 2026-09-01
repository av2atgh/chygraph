"""WP2 checks: stability of the non-trivial fixed point."""

import numpy as np
import pytest
import sympy as sp
from scipy.special import lambertw

from percolation import hypergraph_giant
from statmech import FixedPointStability
from statmech import vertexcover as vc

k, c, p, q = sp.symbols('k c p q')


def _sub(pv):
    return {k: 3, c: 3, p: pv, q: 0.8}


@pytest.fixture(scope='module')
def stab():
    return FixedPointStability(hypergraph_giant())


# ---------------------------------------------------------------------------
# Percolation: exchange of stability
# ---------------------------------------------------------------------------

def test_trivial_root_is_one_plus_Lambda(stab):
    """J at Q = 1 reproduces the WP1 threshold diagnostic."""
    M = stab.model
    for pv in (0.2, 0.5, 0.9):
        assert stab.trivial_perron_root(_sub(pv)) == pytest.approx(
            1 + float(M.Lambda().subs(_sub(pv))), rel=1e-9)


def test_physical_fixed_point_is_stable(stab):
    """rho(J(Q*)) < 1 wherever a giant component exists."""
    for pv in (0.2, 0.3, 0.5, 0.7, 0.9):
        assert stab.is_stable(_sub(pv))


def test_stability_is_monotone_in_occupation(stab):
    """Deeper into the percolating phase is more stable."""
    rho = [stab.spectral_radius(_sub(pv)) for pv in (0.2, 0.3, 0.5, 0.7, 0.9)]
    assert all(a > b for a, b in zip(rho, rho[1:]))


def test_rho_approaches_one_at_threshold(stab):
    """Q* -> 1 and rho -> 1 as Lambda -> 0+: the two fixed points exchange.

    Threshold is c k p q = 1, i.e. p = 1/7.2 at these parameters.
    """
    M = stab.model
    prev = 0.0
    for pv in (0.20, 0.16, 0.15, 0.145, 0.140, 0.1392):
        assert float(M.Lambda().subs(_sub(pv))) > 0
        rho = stab.spectral_radius(_sub(pv))
        assert rho > prev
        prev = rho
    assert prev > 0.998


def test_stability_falls_at_the_rate_the_threshold_rises(stab):
    """1 - rho(J(Q*)) = Lambda + O(Lambda^2).

    The trivial fixed point's Perron root is 1 + Lambda (WP1); the physical
    one's spectral radius is 1 - Lambda to leading order.  The pair exchange
    stability symmetrically through the transition.
    """
    M = stab.model
    for pv, tol in ((0.145, 2e-2), (0.140, 3e-3), (0.1392, 1.5e-3)):
        lam = float(M.Lambda().subs(_sub(pv)))
        rho = stab.spectral_radius(_sub(pv))
        assert abs((1 - rho) / lam - 1) < tol


def test_percolation_map_is_monotone(stab):
    """J >= 0 entrywise: generating functions are order-preserving."""
    assert stab.monotonicity(_sub(0.5)) == 1


def test_core_is_bipartite(stab):
    """Spectrum symmetric under lambda -> -lambda, so the eigenvalue sign
    cannot diagnose monotonicity and monotonicity() reads entry signs."""
    assert stab.is_bipartite_core(_sub(0.5))


# ---------------------------------------------------------------------------
# Vertex cover: the anti-monotone branch
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('mean', [1.0, 2.0, 3.0, 5.0, 10.0])
def test_cover_size_matches_weigt_hartmann(mean):
    """Uncorrelated Poisson: x_c = 1 - (2W + W^2)/(2c), W = LambertW(c)."""
    pd = vc.poisson(mean, 200)
    e, qq, _ = vc.excess(pd)
    W = lambertw(mean).real
    assert vc.cover_size(pd, vc.solve(0.0, e, qq)) == pytest.approx(
        1 - (2 * W + W**2) / (2 * mean), rel=1e-9)


def test_map_is_anti_monotone():
    """Every entry of J is <= 0, unlike any generating-function map."""
    e, qq, _ = vc.excess(vc.scale_free(2.5, 200))
    pi = vc.solve(0.5, e, qq)
    a = vc.jacobian_spectrum(0.5, e, qq, pi)
    J = -0.5 * np.diag(a) - 0.5 * np.outer(np.ones_like(a), qq * a)
    assert (J <= 1e-15).all()
    assert vc.leading_eigenvalue(0.5, e, qq, pi).real < 0


def test_secular_criterion_matches_eigensolve():
    """The rank-one shortcut agrees with a dense eigenvalue computation."""
    e, qq, _ = vc.excess(vc.scale_free(2.5, 300))
    for r in (0.0, 0.2, 0.5, 0.7, 0.75, 0.9):
        pi = vc.solve(r, e, qq)
        assert vc.is_unstable(r, e, qq, pi) == (
            abs(vc.leading_eigenvalue(r, e, qq, pi)) > 1)


def test_fixed_point_found_past_the_instability():
    """The bisection reaches the fixed point where plain iteration diverges."""
    e, qq, _ = vc.excess(vc.scale_free(2.5, 300))
    r = 0.95
    pi = vc.solve(r, e, qq)
    assert vc.is_unstable(r, e, qq, pi)                       # RS broken
    resid = r * (1 - pi) ** e + (1 - r) * (qq * (1 - pi) ** e).sum() - pi
    assert np.abs(resid).max() < 1e-12                        # yet solved
    # plain iteration from any start does not settle there
    x = np.zeros_like(pi)
    for _ in range(500):
        x = r * (1 - x) ** e + (1 - r) * (qq * (1 - x) ** e).sum()
    assert np.abs(x - pi).max() > 1e-3


def test_rs_breaks_inside_the_unit_interval():
    """VW03: 'the RS solution breaks at a certain value of r that depends on
    gamma'.  Both sit in (0,1), and the heavier tail breaks earlier."""
    r25, r30 = vc.rsb_point(2.5), vc.rsb_point(3.0)
    assert 0.0 < r25 < 1.0 and 0.0 < r30 < 1.0
    assert r25 < r30


def test_uncorrelated_scale_free_is_replica_symmetric():
    """VW03: uncorrelated power-law networks are simple."""
    for gamma in (2.5, 3.0):
        e, qq, _ = vc.excess(vc.scale_free(gamma, 400))
        assert not vc.is_unstable(0.0, e, qq, vc.solve(0.0, e, qq))


def test_cover_grows_with_assortativity():
    """VW03 Fig. 1: x_c increases with r, and gamma=2.5 lies below gamma=3.0."""
    prev = {2.5: -1.0, 3.0: -1.0}
    for r in (0.0, 0.2, 0.4, 0.6):
        xs = {}
        for gamma in (2.5, 3.0):
            pd = vc.scale_free(gamma, 400)
            e, qq, _ = vc.excess(pd)
            xs[gamma] = vc.cover_size(pd, vc.solve(r, e, qq))
            assert xs[gamma] > prev[gamma]
            prev[gamma] = xs[gamma]
        assert xs[2.5] < xs[3.0]


def test_exchange_residual_is_quadratic(stab):
    """1 - rho - Lambda divided by Lambda^2 converges, which is what makes the
    O(Lambda^2) in the exchange relation a statement about the order rather
    than a bound at one point (referee minor 6)."""
    M = stab.model
    ratios = []
    for pv in (0.15, 0.145, 0.140, 0.1392):
        lam = float(M.Lambda().subs(_sub(pv)))
        ratios.append((1 - stab.spectral_radius(_sub(pv)) - lam) / lam**2)
    assert all(r < 0 for r in ratios)
    assert abs(ratios[-1] - ratios[-2]) < abs(ratios[1] - ratios[0])   # converging
    assert ratios[-1] == pytest.approx(-0.61, abs=0.02)
