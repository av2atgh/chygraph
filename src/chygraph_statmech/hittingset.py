"""Minimal hitting set on a hypergraph (WP3).

A hitting set (hypergraph vertex cover, transversal) is a set of vertices
meeting every hyperedge.  Finding a minimum one is NP-hard, and it is the
natural higher-order lift of the vertex-cover problem of Vazquez & Weigt.

Hard-core lattice gas, following Weigt & Hartmann: ``x_i = 1`` when ``i`` is
*outside* the cover, so the complement of a hitting set is a set containing no
hyperedge whole, and the constraint on a hyperedge ``e`` is

    e^{w} = 1 - prod_{i in e} x_i

which is ``1 - x_i x_j`` at cardinality 2 — the graph case, exactly.  Minimum
hitting sets are maximal packings, reached as ``mu -> inf``.

In that limit the fields collapse onto integers, ``h = mu z`` with
``z in {1, 0, -1, ...}``, and the factor-graph cavity recursion becomes

    z_i = 1 - #{ hyperedges of i whose every *other* member has z = 1 }.

Writing ``sigma`` for the probability that a node reached along a hyperedge has
``z = 1``, and separating hyperedges into layers by cardinality so that a
node's participation across cardinalities can be correlated,

    sigma_m = Phibar^(m)( 1 - sigma_1^{c_1-1}, ..., 1 - sigma_L^{c_L-1} )

with ``Phi`` the joint generating function of the hyperdegrees ``(k_1..k_L)``
and ``Phibar^(m) = (dPhi/dx_m) / (dPhi/dx_m|_1)`` its inclusion-biased excess,
exactly as in ``chygraph``'s ``JointChygraph``.  Layer refinement by cardinality
is the device of ``chygraph.applications.correlated_cardinality_hypergraph``.

The map is order-reversing, so :mod:`chygraph_statmech.antimonotone` solves it.
At ``L = 1``, ``c = 2`` and Poisson hyperdegree it reduces to
``sigma = exp(-k sigma)`` with instability at ``k sigma = 1``, i.e. ``k = e``:
the Bauer-Golinelli core-percolation point, which is the Erdos-Renyi control in
``~/av2atg/computational_complexity``.
"""

import numpy as np
from sympy import Symbol, diff, lambdify, prod, exp as sexp
from scipy.optimize import brentq

from chygraph_statmech import antimonotone as am


class HittingSet:
    """Replica-symmetric solution of minimum hitting set.

    Args:
        cardinalities: length-``L`` list of hyperedge cardinalities, one per
            layer.
        phi: joint hyperdegree generating function, a sympy expression in the
            symbols returned by :func:`layer_symbols`.  Factorising ``phi``
            means uncorrelated layers.
    """

    def __init__(self, cardinalities, phi):
        self.c = np.asarray(cardinalities, dtype=int)
        self.L = L = len(self.c)
        self.x = layer_symbols(L)
        self.phi = phi
        one = {s: 1 for s in self.x}
        d = [diff(phi, s) for s in self.x]
        self.means = np.array([float(di.subs(one)) for di in d])
        if (self.means <= 0).any():
            raise ValueError("every layer needs a positive mean hyperdegree")
        self._phibar = [lambdify(self.x, di / di.subs(one), 'numpy')
                        for di in d]
        self._phi = lambdify(self.x, phi, 'numpy')
        self._dphibar = [[lambdify(self.x, diff(di / di.subs(one), s), 'numpy')
                          for s in self.x] for di in d]

    # -- the map ------------------------------------------------------------

    def tau(self, sigma):
        """``tau_l``: probability a cardinality-``c_l`` hyperedge constrains."""
        return np.asarray(sigma, float) ** (self.c - 1)

    def F(self, sigma):
        """One step of the message map; order-reversing in ``sigma``."""
        y = 1.0 - self.tau(sigma)
        return np.array([float(f(*y)) for f in self._phibar])

    def solve(self, **kw):
        """The fixed point, stable or not."""
        return am.fixed_point(self.F, self.L, **kw)

    # -- diagnostics --------------------------------------------------------

    def jacobian(self, sigma):
        """``d sigma_m / d sigma_n``; every entry is non-positive."""
        s = np.asarray(sigma, float)
        y = 1.0 - self.tau(s)
        dtau = (self.c - 1) * s ** (self.c - 2)
        return np.array([[-float(self._dphibar[m][n](*y)) * dtau[n]
                          for n in range(self.L)] for m in range(self.L)])

    def leading_eigenvalue(self, sigma=None):
        sigma = self.solve() if sigma is None else sigma
        w = np.linalg.eigvals(self.jacobian(sigma))
        return complex(w[np.argmax(np.abs(w))])

    def is_unstable(self, sigma=None):
        """``|lambda| > 1``: replica symmetry is broken."""
        return bool(abs(self.leading_eigenvalue(sigma)) > 1.0)

    def is_replica_symmetric(self, **kw):
        """The independent check: does the ``F o F`` bracket close?"""
        return am.is_replica_symmetric(self.F, self.L, **kw)

    # -- observables --------------------------------------------------------

    def isolated_fraction(self):
        """``Phi(0,...,0)``: the fraction of vertices in no hyperedge at all.

        A diagnostic, not an observable.  Mixture ensembles built to correlate
        layer participation can silently push one class toward zero
        hyperdegree, in which case a shift in the RSB point is dilution rather
        than correlation.  Report this alongside any such comparison.
        """
        return float(self._phi(*np.zeros(self.L)))

    def cover_size(self, sigma=None):
        """``x_c``, the fraction of vertices in a minimum hitting set.

        ``z > 0`` fixes a vertex outside the cover and ``z < 0`` inside;
        ``z = 0`` is degenerate and contributes one half, which is the
        ``(d-2)/2`` term of Vazquez & Weigt Eq. (17) in the graph case.
        """
        sigma = self.solve() if sigma is None else sigma
        tau = self.tau(sigma)
        return float(1.0 - self._phi(*(1.0 - tau))
                     - 0.5 * float((self.means * tau * sigma).sum()))


# ---------------------------------------------------------------------------
# Generating functions
# ---------------------------------------------------------------------------

def layer_symbols(L):
    return [Symbol(f'x{l}') for l in range(L)]


def poisson_phi(means):
    """Independent Poisson hyperdegrees: ``prod_l exp(k_l (x_l - 1))``."""
    x = layer_symbols(len(means))
    return prod([sexp(k * (xi - 1)) for k, xi in zip(means, x)])


def two_class_phi(means_a, means_b, weight):
    """A mixture of two node classes with different layer preferences.

    ``weight`` of nodes draw Poisson hyperdegrees with means ``means_a``, the
    rest with ``means_b``.  Correlates participation across cardinalities while
    leaving each marginal mean at ``weight*a + (1-weight)*b``; the independent
    comparison is :func:`poisson_phi` of those marginals.
    """
    return (weight * poisson_phi(means_a)
            + (1 - weight) * poisson_phi(means_b))


# ---------------------------------------------------------------------------
# Convenience
# ---------------------------------------------------------------------------

def poisson(cardinalities, means):
    """Uncorrelated Poisson hyperdegrees over the given cardinality layers."""
    return HittingSet(cardinalities, poisson_phi(means))


def rsb_point(cardinalities, direction, lo=1e-6, hi=50.0):
    """Scale ``direction`` until replica symmetry breaks.

    ``direction`` is a vector of layer mean hyperdegrees; the returned scalar
    ``t`` is where ``t * direction`` first becomes unstable.

    ``hi`` is capped well below where ``exp(k(x-1))`` overflows in double
    precision: sympy factors it as ``exp(-k) exp(kx)``, so mean hyperdegrees
    past a few hundred lose the cancellation.  Every RSB point here is O(1), so
    this costs nothing; raise ``hi`` deliberately if that ever stops being true.
    """
    d = np.asarray(direction, float)

    return rsb_scale(lambda t: poisson(cardinalities, t * d), lo, hi)


def rsb_scale(factory, lo=1e-6, hi=50.0):
    """Where a one-parameter family of hypergraphs breaks replica symmetry.

    ``factory(t)`` returns a :class:`HittingSet`; the returned ``t`` is where
    ``|lambda|`` first reaches 1.  This is the general form of
    :func:`rsb_point`, which is the uncorrelated-Poisson special case.
    """
    def gap(t):
        return abs(factory(t).leading_eigenvalue()) - 1.0

    if gap(lo) * gap(hi) > 0:
        raise ValueError(
            f"no sign change of |lambda| - 1 on [{lo}, {hi}]: "
            f"gap({lo}) = {gap(lo):.3g}, gap({hi}) = {gap(hi):.3g}")
    return brentq(gap, lo, hi, xtol=1e-13, rtol=8.9e-16)
