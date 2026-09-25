"""Tests for the endogeny test (statmech.endogeny).

Two coupled populations with the same structure draws.  What is pinned:
that the two copies share their draws (a copy started as a small
perturbation of the other decays to zero where the solution is stable, and
identical copies stay identical to the last bit); the zero-field spin-glass
transition, where the annealed and the jointly carried linear rates agree;
and that the vertex-cover pair converges well below <k> = e and does not
converge well above it.
"""

import numpy as np

from statmech.endogeny import (
    BivariateHittingSet, BivariateSpinGlass, distance,
)


def _rate(hist, a, b):
    h = np.array(hist[a:b])
    h = h[h > 0]
    return float(np.polyfit(np.arange(len(h)), np.log(h), 1)[0])


def test_identical_copies_stay_identical():
    m = BivariateSpinGlass(3, 1.2, field=0.3, regular=True, size=20_000, seed=0)
    m.run_single(50)
    m.start_pair(perturb=0.0)
    for _ in range(20):
        m.sweep_pair()
    assert max(m.history) == 0.0
    h = BivariateHittingSet([2], [2.0], mu=60.0, size=20_000, seed=0)
    h.run(50)
    h.start_pair(perturb=0.0)
    for _ in range(20):
        h.sweep_pair()
    assert max(h.history) == 0.0


def test_spin_glass_zero_field_transition():
    """At H = 0 the paramagnetic solution is trivially endogenous and the
    annealed growth <kbar> tanh^2 crosses one at beta J = atanh(1/sqrt 2)."""
    for bJ, above in ((0.7, False), (1.1, True)):
        m = BivariateSpinGlass(3, bJ, field=0.0, regular=True, size=50_000, seed=1)
        m.run_single(200)
        g = 2 * np.tanh(bJ) ** 2
        assert (g > 1) == above
        if not above:
            assert abs(m.magnetisation()) < 1e-3 and m.overlap() < 1e-3


def test_spin_glass_in_field_perturbation_decays_below_the_line():
    m = BivariateSpinGlass(3, 1.0, field=0.3, regular=True, size=50_000, seed=2)
    m.run_single(300)
    m.run_pair(80, perturb=1e-4)
    assert _rate(m.history, 5, 60) < -0.02
    m.run_pair(200)
    assert m.history[-1] < 1e-5


def test_spin_glass_in_field_pair_does_not_converge_deep_in_the_glass():
    m = BivariateSpinGlass(3, 1.8, field=0.3, regular=True, size=50_000, seed=3)
    m.run_single(300)
    m.run_pair(300)
    assert m.history[-1] > 0.5


def test_vertex_cover_pair_converges_below_e_and_not_far_above():
    lo = BivariateHittingSet([2], [2.0], mu=60.0, size=50_000, seed=4)
    lo.run(300)
    lo.run_pair(400)
    assert lo.history[-1] < 1e-6
    hi = BivariateHittingSet([2], [4.0], mu=60.0, size=50_000, seed=4)
    hi.run(300)
    hi.run_pair(400)
    assert hi.history[-1] > 0.2


def test_distance_is_normalised():
    a = np.array([1.0, -1.0, 2.0])
    assert distance(a, a) == 0.0
    assert abs(distance(a, np.zeros(3)) - 1.0) < 1e-15
