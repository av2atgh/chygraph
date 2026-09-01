"""WP6 checks: the Bethe free energy on a chygraph."""

import numpy as np
import pytest

from statmech.freeenergy import (BetheFreeEnergy, free_energy_gap,
                                          graph_paramagnetic, paramagnetic)
from statmech.population import CavityPopulation
from statmech.region import RegionGraph

T_TRI = (7 - np.sqrt(45)) / 2
SYSTEMS = {'graph': ([2], [6.0], np.arctanh(1 / 6.0)),
           'triangles': ([3], [3.0], np.arctanh(T_TRI))}


# ---------------------------------------------------------------------------
# Closed forms
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('mean', [1.0, 3.0, 6.0])
@pytest.mark.parametrize('bJ', [0.05, 0.3, 1.0])
def test_graph_paramagnetic_is_the_textbook_form(mean, bJ):
    """ln 2 + (c/2) ln cosh(beta J)."""
    assert paramagnetic([2], [mean], bJ) == pytest.approx(
        graph_paramagnetic(mean, bJ), abs=1e-14)


@pytest.mark.parametrize('spec', [([2], [4.0]), ([3], [2.0]), ([2, 3], [2.0, 1.0])])
def test_free_spins_at_infinite_temperature(spec):
    """Every complex decouples at beta J = 0, leaving ln 2 per spin."""
    assert paramagnetic(*spec, 0.0) == pytest.approx(np.log(2.0), abs=1e-14)


# ---------------------------------------------------------------------------
# The population estimator
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('spec', [([2], [6.0]), ([3], [3.0]), ([2, 3], [2.0, 1.0])])
@pytest.mark.parametrize('bJ', [0.05, 0.12])
def test_zero_field_estimator_is_exact(spec, bJ):
    """With the ``ln 2`` piece taken exactly rather than sampled, the
    paramagnetic value comes out to machine precision.  Sampling it instead
    costs ~5e-3, which is larger than the free-energy differences at issue."""
    p = CavityPopulation(*spec, bJ, size=20_000, seed=1)
    p.initialise(0.0)
    for _ in range(20):
        p.sweep()
    assert BetheFreeEnergy(p).minus_beta_f() == pytest.approx(
        paramagnetic(*spec, bJ), abs=1e-12)


# ---------------------------------------------------------------------------
# The transition, located thermodynamically
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('name', list(SYSTEMS))
def test_no_gap_above_the_transition(name):
    cs, ms, bc = SYSTEMS[name]
    for f in (0.8, 0.95):
        assert free_energy_gap(cs, ms, f * bc, sweeps=200,
                               size=40_000, seed=1) == pytest.approx(0.0, abs=1e-9)


@pytest.mark.parametrize('name', list(SYSTEMS))
def test_ordered_branch_wins_below_the_transition(name):
    """-beta f is larger for the ordered branch, i.e. its free energy is lower.

    This locates T_c using neither WP1's linearisation nor WP4's order
    parameter, and it agrees with both.
    """
    cs, ms, bc = SYSTEMS[name]
    prev = 0.0
    for f in (1.05, 1.2, 1.5):
        gap = free_energy_gap(cs, ms, f * bc, sweeps=200, size=40_000, seed=1)
        assert gap > prev
        prev = gap
    assert prev > 1e-2


def test_triangles_gain_free_energy_before_the_matched_graph():
    """At six neighbours either way, the triangle network has already ordered
    at a coupling where the graph's ordered branch does not exist."""
    bJ = np.arctanh(1 / 6.0) * 0.93
    assert free_energy_gap([2], [6.0], bJ, sweeps=200, size=40_000,
                           seed=1) == pytest.approx(0.0, abs=1e-9)
    assert free_energy_gap([3], [3.0], bJ, sweeps=200, size=40_000, seed=1) > 1e-5


# ---------------------------------------------------------------------------
# WP5 and WP6 are the same counting
# ---------------------------------------------------------------------------

def test_free_energy_weights_are_the_mobius_counting_numbers():
    """The ``1`` per complex and ``1 - k_v`` per node that the free energy sums
    are exactly the region-graph counting numbers, in the treelike case."""
    complexes = [[0, 1, 2], [2, 3, 4], [4, 5, 6]]
    rg = RegionGraph(complexes)
    assert rg.is_bethe()
    k = {v: sum(v in a for a in complexes) for v in range(7)}
    for a in complexes:
        assert rg.counting[frozenset(a)] == 1
    for v in range(7):
        assert rg.counting[frozenset({v})] == 1 - k[v]
