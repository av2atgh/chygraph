"""Leaf-removal (Karp-Sipser) core percolation as a chygraph fixed point.

Step 2 of TODO item 1.  The leaf-removal core is the order parameter
``~/av2atg/computational_complexity`` measures: repeatedly delete a degree-1
vertex together with its neighbour, and the residue -- every vertex of degree
>= 2 -- is the core.  It is empty below a threshold and extensive above, at
``c = e`` for Erdos-Renyi (Bauer & Golinelli).  It is *not* the giant component,
so nothing in ``chygraph`` computes it.

**It is a chygraph fixed point, with a three-state message.**  Percolation gets
away with a scalar because "is this branch connected to the giant component" has
two answers.  Leaf removal has three, and they are all needed: run the process
on the cavity branch at ``i`` with complex ``a`` deleted, and ``i`` ends as

    L  leaf-ready -- survives with no surviving neighbour outside ``a``, so in
       the full graph its only neighbours are inside ``a``;
    D  deleted    -- removed during the cavity process;
    C  core-side  -- survives with at least one surviving neighbour outside.

Write ``(lambda, delta, gamma)`` for their probabilities on a directed
inclusion, one triple per layer.  The node step is the chy-degree step of any
chygraph, a generating function of the messages coming down:

    lambda_m = Phibar^(m)( zeta ),        delta_m = 1 - Phibar^(m)( 1 - kappa )

where ``zeta_l`` is the chance a layer-``l`` complex gives the node no surviving
neighbour and ``kappa_l`` the chance it forces the node out.  Those two come
from solving leaf removal *inside* the complex, which is the intra-complex step
and the part percolation never has to do.

**Inside a clique of cardinality c.**  Let ``m`` be the number of its other
``c-1`` members that are live (not ``D``); each is live with probability
``1 - delta``.  Every live member has clique-degree exactly ``m``, so:

* ``m = 0``  -- the complex gives ``i`` no neighbour at all;
* ``m = 1``  -- the one live member has clique-degree 1, so it is a leaf if it
  has no outside neighbour either (state ``L``), and then ``i`` is deleted as
  that leaf's partner.  If it is ``C`` it has an outside neighbour, degree 2,
  and gives ``i`` one surviving neighbour;
* ``m >= 2`` -- every live member has degree >= 2, nothing inside the complex is
  a leaf, and ``i`` gets ``m`` surviving neighbours.

So a complex kills only through the single-live-member-with-no-outside-edge
configuration, and

    zeta  = delta^{c-1},        kappa = (c-1) delta^{c-2} lambda.

At ``c = 2`` these are ``delta`` and ``lambda``, the map collapses to
``lambda = Gbar(delta)``, ``delta = 1 - Gbar(1 - lambda)``, and the threshold is
``Gbar'(1 - lambda) = 1``, which is ``c = e`` on Erdos-Renyi.

The map is order-*preserving* in ``(lambda, delta)``, so it is solved by
monotone iteration from zero exactly as ``chygraph.giant.Chygraph.solve`` does,
and unlike the hitting-set map of :mod:`chygraph_statmech.hittingset`.  Its
Jacobian has zero diagonal blocks and a spectrum symmetric about zero -- the
bipartite core structure of WP2.
"""

import numpy as np
from sympy import diff, lambdify

from chygraph_statmech.hittingset import layer_symbols, poisson_phi, two_class_phi  # noqa: F401


class CorePercolation:
    """Leaf-removal core of a graph presented as a chygraph of cliques.

    Args:
        cardinalities: length-``L`` list of clique sizes, one per layer.
            ``[2]`` is an ordinary graph.
        phi: joint chy-degree generating function -- the number of layer-``l``
            complexes a node belongs to -- as a sympy expression in
            :func:`~chygraph_statmech.hittingset.layer_symbols`.
    """

    def __init__(self, cardinalities, phi):
        self.c = np.asarray(cardinalities, dtype=int)
        if (self.c < 2).any():
            raise ValueError("cardinalities must be at least 2")
        self.L = L = len(self.c)
        self.x = layer_symbols(L)
        one = {s: 1 for s in self.x}
        d = [diff(phi, s) for s in self.x]
        self.means = np.array([float(di.subs(one)) for di in d])
        if (self.means <= 0).any():
            raise ValueError("every layer needs a positive mean chy-degree")
        self._phi = lambdify(self.x, phi, 'numpy')
        self._dphi = [lambdify(self.x, di, 'numpy') for di in d]
        self._phibar = [lambdify(self.x, di / di.subs(one), 'numpy')
                        for di in d]

    # -- the intra-complex step --------------------------------------------

    def zeta(self, delta):
        """``delta^{c-1}``: the complex gives the node no surviving neighbour."""
        return np.asarray(delta, float) ** (self.c - 1)

    def kappa(self, lam, delta):
        """``(c-1) delta^{c-2} lambda``: the complex forces the node out."""
        d = np.asarray(delta, float)
        return (self.c - 1) * d ** (self.c - 2) * np.asarray(lam, float)

    # -- the map ------------------------------------------------------------

    def F(self, state):
        """One step.  ``state`` is ``[lambda_0..lambda_{L-1}, delta_0..]``."""
        lam, delta = state[:self.L], state[self.L:]
        z, k = self.zeta(delta), self.kappa(lam, delta)
        lam_new = np.array([float(f(*z)) for f in self._phibar])
        del_new = np.array([1.0 - float(f(*(1.0 - k))) for f in self._phibar])
        return np.concatenate((lam_new, del_new))

    def solve(self, iters=200000, tol=1e-14):
        """Smallest fixed point, by monotone iteration from zero.

        The map has non-negative partial derivatives, so this converges upward
        to the physical solution, exactly as for a percolation chygraph.
        """
        s = np.zeros(2 * self.L)
        for _ in range(iters):
            t = self.F(s)
            if np.abs(t - s).max() < tol:
                return t
            s = t
        return s

    def state(self, s=None):
        """``(lambda, delta, gamma)`` at the fixed point."""
        s = self.solve() if s is None else s
        lam, delta = s[:self.L], s[self.L:]
        return lam, delta, 1.0 - lam - delta

    # -- the order parameter ------------------------------------------------

    def core_fraction(self, s=None):
        """Fraction of nodes in the leaf-removal core.

        A node is in the core when no complex forces it out and at least two
        surviving neighbours remain.  Let ``N_l(y)`` generate the number of
        surviving neighbours one layer-``l`` complex contributes, with the
        killing configuration excluded:

            N_l(y) = (delta_l + (1-delta_l) y)^{c_l-1} - (c_l-1) delta_l^{c_l-2} lambda_l y

        Then ``Phi(N(y))`` generates the total, and the core fraction is what is
        left after dropping the ``y^0`` and ``y^1`` terms.
        """
        lam, delta, gamma = self.state(s)
        N1 = 1.0 - self.kappa(lam, delta)                  # N_l(1)
        N0 = self.zeta(delta)                              # N_l(0)
        dN0 = (self.c - 1) * delta ** (self.c - 2) * gamma  # N_l'(0)
        no_kill = float(self._phi(*N1))
        none = float(self._phi(*N0))
        one = sum(float(self._dphi[l](*N0)) * dN0[l] for l in range(self.L))
        return max(0.0, no_kill - none - one)

    # -- the threshold ------------------------------------------------------

    def jacobian(self, s=None):
        """``dF/dstate`` at the fixed point, by central differences.

        Zero diagonal blocks: ``lambda`` depends only on ``delta`` through
        ``zeta``, and ``delta`` only on ``lambda`` and ``delta`` through
        ``kappa``.  On the core-free branch the spectrum is symmetric about
        zero, the WP2 bipartite structure.
        """
        s = self.solve() if s is None else s
        n, h = 2 * self.L, 1e-7
        J = np.zeros((n, n))
        for j in range(n):
            a, b = s.copy(), s.copy()
            a[j] += h
            b[j] = max(0.0, b[j] - h)
            J[:, j] = (self.F(a) - self.F(b)) / (a[j] - b[j])
        return J

    def spectral_radius(self, s=None):
        return float(np.max(np.abs(np.linalg.eigvals(self.jacobian(s)))))

    def has_core(self, tol=1e-9, s=None):
        """True when ``gamma > 0`` in some layer: the core is extensive."""
        return bool(np.max(self.state(s)[2]) > tol)


# ---------------------------------------------------------------------------
# Convenience
# ---------------------------------------------------------------------------

    # -- the core-free branch ----------------------------------------------

    def has_core_free_branch(self):
        """Whether a solution with ``gamma = 0`` can exist at all.

        Setting ``gamma_l = 0`` means ``delta_l = 1 - lambda_l``, and the two
        halves of the map agree only if ``zeta_l = 1 - kappa_l`` for every
        layer, i.e. with ``u = 1 - lambda``

            u^{c-1} = 1 - (c-1) u^{c-2} (1 - u).

        At ``c = 2`` this is ``u = u``, an identity: the branch always exists.
        At ``c >= 3`` it rearranges to ``(1-u)^{c-1} = 0``, so it holds only at
        ``lambda = 0``, which is not a fixed point of a chygraph with positive
        chy-degree.

        **A complex of cardinality three or more is a core by itself.**  Every
        member has clique-degree ``c-1 >= 2``, so leaf removal can never reach
        it, and the core is extensive at *any* density.  Nothing analogous
        happens on a graph, where the core needs ``Gbar'(1-lambda) > 1``.
        """
        return bool((self.c == 2).all())

    def core_free_lambda(self):
        """``lambda`` on the ``gamma = 0`` branch, from ``lambda = Phibar(1-lambda)``."""
        if not self.has_core_free_branch():
            raise ValueError("no core-free branch: some cardinality is >= 3")
        from scipy.optimize import brentq
        f = self._phibar[0]
        return brentq(lambda x: float(f(*(np.full(self.L, 1.0 - x)))) - x,
                      1e-15, 1.0, xtol=1e-15, rtol=8.9e-16)

    def core_free_spectral(self):
        """``Phibar'(1-lambda)`` on the core-free branch; the core appears at 1.

        The Jacobian there has zero diagonal blocks and eigenvalues
        ``+-Phibar'(1-lambda)``, the bipartite structure of WP2, so this single
        number decides stability.  On Erdos-Renyi it is ``k lambda`` with
        ``lambda = exp(-k lambda)``, giving the threshold ``k = e`` exactly.
        """
        lam = self.core_free_lambda()
        u = np.full(self.L, 1.0 - lam)
        h = 1e-6
        up, dn = u.copy(), u.copy()
        up[0] += h
        dn[0] -= h
        return (float(self._phibar[0](*up)) - float(self._phibar[0](*dn))) / (2 * h)


def graph(degree_pgf=None, mean=None):
    """Ordinary graph: one layer of cardinality 2.

    With ``mean`` given, Poisson degree -- the Erdos-Renyi case whose core
    appears at ``mean = e``.
    """
    if degree_pgf is None:
        degree_pgf = poisson_phi([mean])
    return CorePercolation([2], degree_pgf)


def clique_network(cardinality, mean_membership):
    """Poisson number of cliques per node, all of one cardinality."""
    return CorePercolation([cardinality], poisson_phi([mean_membership]))


def core_threshold(factory, lo=1e-6, hi=50.0, iters=80):
    """Where the core first becomes extensive.

    Bisects on the stability of the core-free branch, which is exact, rather
    than on ``gamma > tol``, which is limited by how continuously ``gamma``
    rises out of the transition.  Raises when the family has no core-free
    branch -- cardinality three or more, where the core exists at every
    density and there is no threshold to find.
    """
    from scipy.optimize import brentq
    if not factory(0.5 * (lo + hi)).has_core_free_branch():
        raise ValueError("no threshold: cardinality >= 3 is always cored")
    return brentq(lambda t: factory(t).core_free_spectral() - 1.0,
                  lo, hi, xtol=1e-14, rtol=8.9e-16)
