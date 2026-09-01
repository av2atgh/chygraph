"""Minimum vertex cover of the graph underlying a chygraph.

:mod:`statmech.hittingset` covers every *hyperedge*: a hyperedge is
satisfied as soon as one member is taken, so at most ``c-1`` members may be
left out.  Vertex cover of the graph a chygraph induces is the other end of the
same family.  Two members of a clique are adjacent, so leaving both out leaves
an edge uncovered: **at most one member per complex may be left out**, and the
uncovered set is exactly an independent set of the underlying graph.

The cavity in the ``mu -> inf`` limit differs from the hyperedge case in one
place, and it is the place the complex has to be solved rather than factorised.
Taking node ``i`` blocks every other member of every complex containing it --
but a complex could only ever have contributed *one* member anyway, so the cost
of taking ``i`` is the **maximum** over the complex, not the sum:

    z_{i->a} = 1 - sum_{b != a} max(0, max_{j in b, j != i} z_{j->b}).

With ``sigma_l`` the chance a node reached through a layer-``l`` complex has
``z = 1``, a complex blocks unless none of its other ``c-1`` members has
``z = 1``, so

    sigma_m = Phibar^(m)(tau),      tau_l = (1 - sigma_l)^{c_l - 1},

against ``tau_l = sigma_l^{c_l-1}`` and an outer complement for hitting set.
Both collapse to the same thing at ``c = 2``, which is the graph, and both are
order-reversing, so :mod:`statmech.antimonotone` solves them.

**Where this is exact, and where it is not.**  At ``c = 2`` the cover size
reduces to Weigt-Hartmann and the solution is exact wherever leaf removal
certifies it, i.e. below the core threshold ``<kappabar> = e``.  At ``c >= 3``
it is *never* certified, and the formalism says so itself: a chygraph carrying
any layer of cardinality three or more has no core-free branch
(:mod:`statmech.core`), so the leaf-removal core is always extensive
and replica symmetry is never established.  The degeneracy rule inherited from
the graph -- a node at ``z = 0`` counted as half taken -- is the visible
symptom.  On a single triangle every vertex sits at ``z = 0``, the rule gives
``1/2``, and the truth is ``1/3``.
"""

import numpy as np

from statmech import antimonotone as am
from statmech.hittingset import (layer_symbols, poisson_phi,  # noqa: F401
                                          two_class_phi)
from sympy import diff, lambdify


class CliqueCover:
    """Replica-symmetric minimum vertex cover of a chygraph's induced graph.

    Args:
        cardinalities: clique size per complex layer; ``[2]`` is a graph.
        phi: joint chy-degree generating function, sympy expression in
            :func:`~statmech.hittingset.layer_symbols`.
    """

    def __init__(self, cardinalities, phi):
        self.c = np.asarray(cardinalities, dtype=int)
        self.L = L = len(self.c)
        self.x = layer_symbols(L)
        one = {s: 1 for s in self.x}
        d = [diff(phi, s) for s in self.x]
        self.means = np.array([float(di.subs(one)) for di in d])
        if (self.means <= 0).any():
            raise ValueError("every layer needs a positive mean chy-degree")
        self._phi = lambdify(self.x, phi, 'numpy')
        self._phibar = [lambdify(self.x, di / di.subs(one), 'numpy')
                        for di in d]

    # -- the map ------------------------------------------------------------

    def tau(self, sigma):
        """``(1 - sigma)^{c-1}``: the complex leaves the node free."""
        return (1.0 - np.asarray(sigma, float)) ** (self.c - 1)

    def F(self, sigma):
        t = self.tau(sigma)
        return np.array([float(f(*t)) for f in self._phibar])

    def solve(self, **kw):
        return am.fixed_point(self.F, self.L, **kw)

    def is_replica_symmetric(self, **kw):
        return am.is_replica_symmetric(self.F, self.L, **kw)

    # -- observables --------------------------------------------------------

    def cover_size(self, sigma=None, degeneracy=0.5):
        """Fraction of vertices in a minimum cover.

        ``z > 0`` leaves a vertex out of the cover and ``z < 0`` puts it in;
        ``z = 0`` is degenerate and is counted with weight ``degeneracy``.  The
        default ``1/2`` is the graph rule of Vazquez & Weigt and is exact at
        ``c = 2``; at ``c >= 3`` it is an approximation (see the module
        docstring), and :meth:`cover_bracket` is the honest statement.
        """
        sigma = self.solve() if sigma is None else sigma
        t = self.tau(sigma)
        free = float(self._phi(*t))                       # P(z = 1)
        marginal = float((self.means * (1.0 - t) * sigma).sum())  # P(z = 0)
        return 1.0 - free - degeneracy * marginal

    def cover_bracket(self, sigma=None):
        """``(lower, upper)`` from counting the degenerate vertices either way.

        Every ``z = 0`` vertex is in some minimum cover and out of another, so
        the truth lies between counting them all in and all out.  At ``c = 2``
        the bracket is narrow and the midpoint is exact; at ``c >= 3`` the
        bracket is what can be asserted.
        """
        return (self.cover_size(sigma, 1.0), self.cover_size(sigma, 0.0))

    def certified(self):
        """Whether the replica-symmetric value can be trusted.

        False whenever any layer has cardinality >= 3, because then
        :mod:`statmech.core` has no core-free branch: leaf removal
        leaves an extensive core at every density and never certifies a
        minimum cover.
        """
        return bool((self.c == 2).all())


def poisson(cardinalities, means):
    return CliqueCover(cardinalities, poisson_phi(means))
