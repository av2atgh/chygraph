"""Stability of the *non-trivial* fixed point (WP2).

``percolation.giant.Chygraph.jacobian`` computes ``J = dF/dQ`` and substitutes
``Q = 1``.  That is the trivial fixed point, and its Perron root is
``lambda = 1 + Lambda``: the threshold diagnostic of WP1.  It says where the
giant component appears, and nothing about the solution once it has.

This module evaluates the same symbolic Jacobian at the fixed point ``Q*``
actually reached, which answers a different question — is the solution the map
converges to locally stable?  For percolation the answer is always yes and the
spectral radius is a consistency check.  For maps outside the chygraph class it
is the diagnostic that matters:

* Vazquez & Weigt use exactly this instability as their RSB detector — "an
  instability prevents the program from convergence, and thus provides a precise
  tool to detect RSB".
* The period-doubling route of triadic percolation (``chygraph/TODO.md`` item 3)
  is governed by the same derivative at the same fixed point.

Both are loss of monotonicity in the message map.  What distinguishes them is
the sign of the *entries* of ``J``, not the sign of its leading eigenvalue: a
map built from generating functions has ``J >= 0`` and loses stability through
``+1``, an anti-monotone one has ``J <= 0`` and loses it through ``-1``, a
period-2 orbit.

The eigenvalue sign will not do the job because the coupled core of a chygraph
is *bipartite* in the up/down message role -- a node message feeds a complex
message feeds a node message -- so ``J`` is a non-negative matrix of period 2
and Perron-Frobenius puts eigenvalues at both ``+rho`` and ``-rho``.  The
spectral radius is the quantity that decides stability in either case;
:meth:`monotonicity` reads the entry signs, which are unambiguous.

Nothing here modifies ``chygraph``; ``Chygraph`` is used as given.
"""

import numpy as np
from sympy import Matrix, lambdify


class FixedPointStability:
    """Jacobian of a chygraph map at the fixed point it actually reaches.

    Args:
        model: a ``percolation.giant.Chygraph`` (or ``GiantComponent``).
    """

    def __init__(self, model):
        self.model = model
        self._J = None
        self._compiled = {}

    def symbolic_jacobian(self):
        """``J(Q)``, *not* evaluated at ``Q = 1``."""
        if self._J is None:
            self._J = Matrix(self.model.F()).jacobian(self.model.Q)
        return self._J

    def _compile(self, params):
        key = tuple(str(s) for s in params)
        if key not in self._compiled:
            args = list(self.model.Q) + list(params)
            self._compiled[key] = lambdify(args, self.symbolic_jacobian(), 'numpy')
        return self._compiled[key]

    def at(self, Q, subs=None):
        """``J`` as a numpy array at an explicit ``Q``."""
        subs = dict(subs or {})
        params = sorted(subs, key=str)
        f = self._compile(params)
        return np.array(f(*(list(Q) + [float(subs[k]) for k in params])),
                        dtype=float)

    def jacobian(self, subs=None, **kw):
        """``J`` at the fixed point reached by ``model.solve``."""
        return self.at(self.model.solve(subs, **kw), subs)

    def eigenvalues(self, subs=None, **kw):
        return np.linalg.eigvals(self.jacobian(subs, **kw))

    def perron_root(self, subs=None, **kw):
        """Eigenvalue of largest real part at the fixed point."""
        w = self.eigenvalues(subs, **kw)
        return complex(w[np.argmax(w.real)])

    def spectral_radius(self, subs=None, **kw):
        """``rho(J(Q*))``.  This is what decides stability, for either sign.

        For percolation it is ``< 1`` strictly above threshold and rises to
        exactly ``1`` at it, where ``Q* -> 1`` and the two fixed points
        exchange stability with the trivial one.
        """
        return float(np.max(np.abs(self.eigenvalues(subs, **kw))))

    def monotonicity(self, subs=None, **kw):
        """``+1`` if every entry of ``J`` is ``>= 0``, ``-1`` if all ``<= 0``,
        ``0`` if mixed.

        The unambiguous form of the question the eigenvalue sign cannot answer.
        Generating functions give ``+1``; the ``mu -> inf`` vertex-cover map of
        :mod:`statmech.vertexcover` gives ``-1``.
        """
        J = self.jacobian(subs, **kw)
        if (J >= -1e-12).all():
            return 1
        return -1 if (J <= 1e-12).all() else 0

    def is_stable(self, subs=None, tol=1e-9, **kw):
        return self.spectral_radius(subs, **kw) < 1 - tol

    def is_bipartite_core(self, subs=None, tol=1e-9, **kw):
        """True when the spectrum is symmetric under ``lambda -> -lambda``.

        Holds for every chygraph whose core alternates between up and down
        message roles, and is why :meth:`monotonicity` exists.
        """
        w = np.sort(np.abs(self.eigenvalues(subs, **kw)))
        return bool(np.max(np.abs(w - np.sort(np.abs(-w)))) < tol)

    def trivial_perron_root(self, subs=None):
        """Perron root at ``Q = 1``, for comparison; equals ``1 + Lambda``."""
        w = np.linalg.eigvals(self.at([1.0] * self.model.n, subs))
        return float(np.max(w.real))
