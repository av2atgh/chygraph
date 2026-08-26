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
