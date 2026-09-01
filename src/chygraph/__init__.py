"""chygraph -- percolation and statistical mechanics on complex hypergraphs.

One namespace over the two packages this repository holds:

    percolation/   the self-consistency map, threshold tensor, critical
                   amplitude and correlated layers  (Part II of the book)
    statmech/      the cavity recursion with a general interior, the branching
                   matrix, and the Ising, hitting-set, vertex-cover, colouring,
                   satisfiability and region-graph solvers  (Parts III and IV)

Every public name of both is re-exported here, so

    from chygraph import PercolationMatrix, CriticalAmplitude, GBP

works without knowing which package a name comes from.  The two packages remain
reachable under their own names, ``chygraph.percolation`` and
``chygraph.statmech``, and installing or importing either directly is
unaffected -- this module adds a namespace and takes nothing away.

One name is deliberately not re-exported.  Both packages define a class called
``Chygraph`` and they are different objects: ``percolation.giant.Chygraph`` is
the non-linear self-consistency map whose fixed point gives the giant component,
while ``statmech.api.Chygraph`` is the handle that carries every calculation the
statistical mechanics performs on a chygraph.  Re-exporting both would make one
silently shadow the other, so ``chygraph.Chygraph`` raises instead, and the two
are available as ``ChygraphPercolation`` and ``ChygraphStatmech`` -- or, in
full, as ``chygraph.percolation.Chygraph`` and ``chygraph.statmech.Chygraph``.
"""

import percolation as percolation
import statmech as statmech

from percolation.giant import Chygraph as ChygraphPercolation
from statmech.api import Chygraph as ChygraphStatmech

_AMBIGUOUS = {"Chygraph"}

# Re-export both packages' public API, minus the ambiguous name.  Anything that
# collides and is *not* listed in _AMBIGUOUS is a bug rather than a decision, so
# it is caught here rather than resolved by import order.
_seen: dict[str, str] = {}
__all__ = ["percolation", "statmech", "ChygraphPercolation", "ChygraphStatmech"]

for _pkg, _label in ((percolation, "percolation"), (statmech, "statmech")):
    for _name in _pkg.__all__:
        if _name in _AMBIGUOUS:
            continue
        if _name in _seen:
            raise ImportError(
                f"chygraph: {_name!r} is exported by both {_seen[_name]} and "
                f"{_label} and is not a declared ambiguity. Add it to "
                f"_AMBIGUOUS with distinct aliases, or rename it in one package."
            )
        _seen[_name] = _label
        globals()[_name] = getattr(_pkg, _name)
        __all__.append(_name)

del _pkg, _label, _name, _seen


def __getattr__(name):                                    # PEP 562
    if name in _AMBIGUOUS:
        raise AttributeError(
            "chygraph.Chygraph is ambiguous: percolation and statmech each "
            "define a different Chygraph. Use chygraph.ChygraphPercolation "
            "(percolation.giant.Chygraph, the self-consistency map) or "
            "chygraph.ChygraphStatmech (statmech.api.Chygraph, the calculation "
            "handle); chygraph.percolation.Chygraph and "
            "chygraph.statmech.Chygraph also work."
        )
    raise AttributeError(f"module 'chygraph' has no attribute {name!r}")
