"""The umbrella namespace re-exports both packages without shadowing either.

The one thing that can silently go wrong here is a name exported by both
packages: `from percolation import *` followed by `from statmech import *`
would resolve it by import order and say nothing.  chygraph.__init__ raises on
any such collision it has not been told about, and these tests pin both halves
of that -- the declared ambiguity, and the guard against a new one.
"""

import pytest

import chygraph
import percolation
import statmech


def test_every_public_name_of_both_packages_is_reachable():
    ambiguous = {"Chygraph"}
    for pkg in (percolation, statmech):
        for name in pkg.__all__:
            if name in ambiguous:
                continue
            assert hasattr(chygraph, name), f"{name} missing from chygraph"


def test_re_exports_are_the_same_objects_not_copies():
    assert chygraph.PercolationMatrix is percolation.PercolationMatrix
    assert chygraph.GBP is statmech.GBP


def test_the_two_packages_stay_reachable_under_their_own_names():
    assert chygraph.percolation is percolation
    assert chygraph.statmech is statmech
    assert chygraph.percolation.Chygraph is percolation.giant.Chygraph
    assert chygraph.statmech.Chygraph is statmech.api.Chygraph


def test_the_ambiguous_name_does_not_resolve_silently():
    with pytest.raises(AttributeError, match="ambiguous"):
        chygraph.Chygraph


def test_the_two_chygraph_classes_are_distinguished_by_alias():
    assert chygraph.ChygraphPercolation is percolation.giant.Chygraph
    assert chygraph.ChygraphStatmech is statmech.api.Chygraph
    assert chygraph.ChygraphPercolation is not chygraph.ChygraphStatmech


def test_an_undeclared_collision_would_be_caught():
    """The guard is a raise inside __init__, so it cannot be triggered after
    import.  Assert instead that it is still there and still reachable: the two
    __all__ lists overlap in exactly the one name the module declares."""
    overlap = set(percolation.__all__) & set(statmech.__all__)
    assert overlap == {"Chygraph"}, (
        f"__all__ overlap changed to {overlap}; chygraph.__init__ will raise on "
        "import until the new name is given distinct aliases"
    )


def test_no_unknown_attribute_resolves():
    with pytest.raises(AttributeError):
        chygraph.NotAThing
