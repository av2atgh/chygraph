"""WP5 checks: complexes as regions."""

import pytest

from chygraph_statmech.region import RegionGraph, overlap_profile


def test_counting_numbers_are_valid():
    """Every node counted exactly once, whatever the overlap."""
    for cx in ([[0, 1, 2], [2, 3, 4]], [[0, 1, 2], [1, 2, 3]],
               [[0, 1, 2, 3], [2, 3, 4, 5], [0, 3, 5]]):
        assert RegionGraph(cx).counting_is_valid()


def test_single_node_overlaps_reduce_to_bethe():
    """The reduction that makes Kikuchi the right generalisation: when
    complexes meet in at most one node the counting is 1 and 1 - k_v, which is
    exactly the Bethe counting WP1 and WP4 assume."""
    rg = RegionGraph([[0, 1, 2], [2, 3, 4], [4, 5, 6]])
    assert rg.is_bethe()
    assert rg.bethe_error() == 0
    assert rg.counting[frozenset({2})] == -1        # 1 - 2 complexes
    assert rg.counting[frozenset({0})] == 0         # 1 - 1 complex


def test_edge_sharing_breaks_bethe():
    """Two triangles on a shared edge: Kikuchi subtracts the edge, Bethe
    subtracts its two endpoints.  Both count nodes correctly; only one
    subtracts the correlation."""
    rg = RegionGraph([[0, 1, 2], [1, 2, 3]])
    assert not rg.is_bethe()
    assert rg.counting_is_valid()
    assert rg.counting[frozenset({1, 2})] == -1
    assert rg.bethe_error() > 0


def test_overlap_profile_detects_treelike():
    tree = overlap_profile([[0, 1, 2], [2, 3, 4], [4, 5, 6]])
    assert tree['treelike'] and tree['shared_2plus'] == 0.0
    assert tree['edge_cover_mean'] == pytest.approx(1.0)
    loopy = overlap_profile([[0, 1, 2], [1, 2, 3]])
    assert not loopy['treelike'] and loopy['shared_2plus'] == 1.0
    assert loopy['edge_cover_max'] == 2


def test_disjoint_complexes_have_trivial_counting():
    rg = RegionGraph([[0, 1], [2, 3]])
    assert rg.is_bethe() and rg.bethe_error() == 0
    assert all(rg.counting[frozenset({v})] == 0 for v in range(4))


def test_kikuchi_beats_bethe_on_two_triangles():
    """Enumerate ln Z for four spins, two triangles sharing an edge, and check
    the Mobius counting is closer than the Bethe counting at every coupling.

    Referee-requested: the correction is demonstrated, not asserted.
    """
    import itertools
    import numpy as np
    cx = [[0, 1, 2], [1, 2, 3]]
    rg = RegionGraph(cx)
    edges = sorted({tuple(sorted((a[i], a[j]))) for a in cx
                    for i in range(len(a)) for j in range(i + 1, len(a))})

    def lnZ(bJ, region):
        R = sorted(region)
        e = [(i, j) for i, j in edges if i in R and j in R]
        tot = 0.0
        for s in itertools.product((1, -1), repeat=len(R)):
            d = dict(zip(R, s))
            tot += np.exp(bJ * sum(d[i] * d[j] for i, j in e))
        return np.log(tot)

    for bJ in (0.2, 0.5, 1.0, 2.0):
        exact = lnZ(bJ, range(4))
        kik = sum(c * lnZ(bJ, R) for R, c in rg.counting.items() if c)
        bet = sum(c * lnZ(bJ, R) for R, c in rg.bethe_counting().items() if c)
        assert abs(kik - exact) < abs(bet - exact)      # Kikuchi is closer
        assert bet > exact                              # Bethe overestimates
