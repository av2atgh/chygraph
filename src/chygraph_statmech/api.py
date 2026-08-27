"""One object for one chygraph, with a method per calculation.

The modules of this package each solve one problem and each grew from one
section of the paper.  :class:`Chygraph` is the single handle on all of them: it
holds the *structure* --- the complex layers, their cardinalities, and the
chy-degree distribution --- and its methods are the calculations, delegating to
those modules rather than reimplementing them.

    >>> from chygraph_statmech import Chygraph
    >>> g = Chygraph([2, 3], [4.0, 2.0])          # links and triangles
    >>> g.critical_coupling()                      # Ising T_c, Sec. IV
    >>> g.core().core_fraction()                   # leaf-removal core, Sec. VI
    >>> g.hitting_set_bp(mu=60).run().density()    # hitting set, Sec. V E

A chygraph is specified by ``cardinalities`` and ``degrees``; ``excess``
defaults to ``degrees`` (Poisson chy-degree, its own excess) and is
``degrees - 1`` when ``regular=True``.  Where a calculation needs the full
generating function rather than its moments, ``phi`` may be given as a sympy
expression in :func:`~chygraph_statmech.hittingset.layer_symbols`; otherwise the
independent-Poisson form is built from ``degrees``.

Which method gives which result is listed in the paper; the docstrings name the
section each belongs to.
"""

import numpy as np

from chygraph_statmech import antimonotone as _am
from chygraph_statmech import cavity as _cavity
from chygraph_statmech import core as _core
from chygraph_statmech import cover as _cover
from chygraph_statmech import freeenergy as _fe
from chygraph_statmech import hittingset as _hs
from chygraph_statmech import ising as _ising
from chygraph_statmech import population as _pop
from chygraph_statmech import region as _region
from chygraph_statmech import simplicial as _simp
from chygraph_statmech import softfield as _soft
from chygraph_statmech import stability as _stab


class Chygraph:
    """A chygraph, and every calculation the paper performs on one.

    Args:
        cardinalities: ``c_l``, the size of a layer-``l`` complex.
        degrees: ``<kappa>_l``, the mean number of layer-``l`` complexes a node
            belongs to.
        excess: ``<kappabar>_l``.  Defaults to ``degrees`` (Poisson) or to
            ``degrees - 1`` when ``regular``.
        regular: fixed rather than Poisson chy-degree.
        phi: joint chy-degree generating function, if the layers are correlated
            or non-Poisson.  A sympy expression in ``layer_symbols(L)``.
    """

    def __init__(self, cardinalities, degrees, excess=None, regular=False,
                 phi=None):
        self.c = np.asarray(cardinalities, dtype=int)
        self.k = np.asarray(degrees, dtype=float)
        if self.c.shape != self.k.shape:
            raise ValueError("one cardinality and one chy-degree per layer")
        if (self.c < 2).any():
            raise ValueError("cardinalities must be at least 2")
        self.L = len(self.c)
        self.regular = bool(regular)
        if excess is not None:
            self.kbar = np.asarray(excess, dtype=float)
        elif regular:
            self.kbar = np.maximum(self.k - 1.0, 0.0)
        else:
            self.kbar = self.k.copy()
        self._phi = phi

    def __repr__(self):
        kind = 'regular' if self.regular else 'Poisson'
        return (f"Chygraph(cardinalities={list(self.c)}, "
                f"degrees={list(self.k)}, {kind})")

    # -- the generating function -------------------------------------------

    def phi(self):
        """Joint chy-degree generating function, as given or independent-Poisson."""
        if self._phi is None:
            return _hs.poisson_phi(self.k)
        return self._phi

    # ======================================================================
    # Sec. II-III: structure
    # ======================================================================

    def regions(self, complexes):
        """Region graph and Mobius counting numbers for an explicit complex
        list (Sec. VIII).  Returns :class:`~chygraph_statmech.region.RegionGraph`."""
        return _region.RegionGraph(complexes)

    @staticmethod
    def overlap_profile(complexes):
        """How far an explicit complex family is from treelike (Sec. VIII)."""
        return _region.overlap_profile(complexes)

    # ======================================================================
    # Sec. IV: the Ising model
    # ======================================================================

    def u_prime(self, layer, beta_J, interaction='clique'):
        """``u'``, the cavity derivative of a layer's complex (Sec. IV A).

        ``interaction='clique'`` couples every pair inside the complex;
        ``'simplicial'`` is the unanimity rule of Sec. IV C.
        """
        q = int(self.c[layer])
        if interaction == 'simplicial':
            return _simp.uprime(q, beta_J)
        return _ising.clique_derivative(q, beta_J)

    def branching_matrix(self, beta_J, interaction='clique', squared=False):
        """``B_{lm}`` of Eq. (5) (Sec. III).  ``squared`` gives the AT line."""
        bj = np.broadcast_to(np.asarray(beta_J, float), (self.L,))
        u = np.array([self.u_prime(m, bj[m], interaction)
                      for m in range(self.L)])
        if squared:
            u = u ** 2
        B = np.empty((self.L, self.L))
        for l in range(self.L):
            for m in range(self.L):
                deg = self.kbar[m] if m == l else self.k[m]
                B[l, m] = deg * (self.c[m] - 1) * u[m]
        return B

    def critical_coupling(self, interaction='clique', squared=False, **kw):
        """``beta J`` where the Perron root of ``B`` reaches 1 (Sec. IV B).

        The transition when it is continuous; the spinodal when it is not ---
        see :meth:`transition`.
        """
        if interaction != 'clique':
            raise NotImplementedError(
                "use simplicial().spinodal() for the unanimity interaction")
        return _ising.critical_coupling(self.c, self.k, excess=self.kbar,
                                        squared=squared, **kw)

    def critical_temperature(self, J=1.0, **kw):
        """``T_c / J`` for the continuous case (Sec. IV B)."""
        return 1.0 / (J * self.critical_coupling(**kw))

    def stability_tensor(self, k, K, s, S, wkappa=None, ws=None):
        """The symbolic ``2L^2`` tensor of Sec. III, reweighted (Eq. 4).

        Takes the four moment tables directly, as
        :class:`~chygraph_statmech.stability.StabilityMatrix` does.
        """
        return _stab.StabilityMatrix(k, K, s, S, wkappa=wkappa, ws=ws)

    def population(self, beta_J, **kw):
        """Field distributions by population dynamics (Sec. IV D).

        Returns :class:`~chygraph_statmech.population.CavityPopulation`.
        """
        return _pop.CavityPopulation(self.c, self.k, beta_J, **kw)

    def free_energy(self, beta_J, sweeps=250, **kw):
        """Bethe free energy of the Ising model (Sec. VII).

        Returns :class:`~chygraph_statmech.freeenergy.BetheFreeEnergy` on a
        converged population.
        """
        return _fe.BetheFreeEnergy(self.population(beta_J, **kw).run(sweeps))

    def paramagnetic_free_energy(self, beta_J):
        """``-beta f`` at zero cavity field, in closed form (Sec. VII)."""
        return _fe.paramagnetic(self.c, self.k, beta_J)

    # ======================================================================
    # Sec. IV C: the simplicial (unanimity) interaction
    # ======================================================================

    def simplicial(self, couplings):
        """The simplicial Ising model on this chygraph (Sec. IV C).

        Returns :class:`~chygraph_statmech.simplicial.SimplicialChygraph`,
        which carries ``spinodal``, ``coexistence``, ``transition`` and
        ``minus_beta_f``.  Fixed chy-degree is assumed, as in Ref. [SLG26].
        """
        return _simp.SimplicialChygraph(self.c, self.k, couplings)

    # ======================================================================
    # Sec. V: hitting set and vertex cover
    # ======================================================================

    def hitting_set(self):
        """Hard-field hitting set, Eq. (13) (Sec. V A).

        Thresholds are exact; the cover size only at cardinality two --- see
        :meth:`hitting_set_bp`.
        """
        return _hs.HittingSet(self.c, self.phi())

    def hitting_set_bp(self, mu=60.0, **kw):
        """Hitting set with the ``O(1)`` fields kept, Eq. (16) (Sec. V E).

        Returns :class:`~chygraph_statmech.softfield.HittingSetBP`; call
        ``.run()`` before reading ``density`` or ``entropy``.
        """
        return _soft.HittingSetBP(self.c, self.k, regular=self.regular, mu=mu,
                                  **kw)

    def clique_cover(self):
        """Vertex cover of the induced graph, Eq. (15) (Sec. V C)."""
        return _cover.CliqueCover(self.c, self.phi())

    # ======================================================================
    # Sec. VI: core percolation
    # ======================================================================

    def core(self):
        """Leaf-removal core as a chygraph fixed point, Eq. (17) (Sec. VI).

        Returns :class:`~chygraph_statmech.core.CorePercolation`.
        """
        return _core.CorePercolation(self.c, self.phi())

    @classmethod
    def from_samples(cls, cardinalities, memberships):
        """Build from a *measured* chy-degree ensemble (Sec. VIII).

        ``memberships`` is ``(n, L)``: how many layer-``l`` complexes each
        sampled node belongs to.  Only :meth:`core` uses the empirical
        generating function; the moment-level methods use its means.
        """
        K = np.asarray(memberships, dtype=float)
        obj = cls(cardinalities, K.mean(axis=0))
        obj._samples = K
        return obj

    def core_from_samples(self):
        """:meth:`core` on the measured ensemble, if built by
        :meth:`from_samples`."""
        if not hasattr(self, '_samples'):
            raise ValueError("build with Chygraph.from_samples first")
        return _core.CorePercolation.from_samples(self.c, self._samples)

    def fixed_point_stability(self, model):
        """Jacobian at the fixed point reached, not at ``Q = 1`` (Sec.~II D).

        Takes a ``chygraph.giant.Chygraph`` and returns
        :class:`~chygraph_statmech.fixedpoint.FixedPointStability`, whose
        ``spectral_radius`` obeys ``rho = 1 - Lambda + O(Lambda^2)`` and whose
        ``monotonicity`` reads entry signs rather than the eigenvalue sign.
        """
        from chygraph_statmech import fixedpoint as _fp
        return _fp.FixedPointStability(model)

    @staticmethod
    def vertex_cover_correlated(gamma, r, dmax=800):
        """Vertex cover on a scale-free graph with the degree correlations of
        Ref. [VW03] Eq. (18) (Sec. V B).  Returns ``(cover_size, unstable)``."""
        from chygraph_statmech import vertexcover as _vc
        p = _vc.scale_free(gamma, dmax)
        e, q, _ = _vc.excess(p)
        pi = _vc.solve(r, e, q)
        return _vc.cover_size(p, pi), _vc.is_unstable(r, e, q, pi)

    @staticmethod
    def excess_cardinality(complexes):
        """``<c^2>/<c> - 1`` of a complex family (Sec. VIII B).

        The quantity the threshold tensor needs to be finite, measured on an
        explicit clique list rather than assumed.
        """
        import numpy as _np
        from collections import Counter as _C
        cnt = _C(len(a) for a in complexes)
        c = _np.array(sorted(cnt), float)
        w = _np.array([cnt[int(i)] for i in c], float)
        w /= w.sum()
        m1, m2 = (w * c).sum(), (w * c * c).sum()
        return float(m2 / m1 - 1.0)

    # ======================================================================
    # solvers shared across sections
    # ======================================================================

    @staticmethod
    def solve_anti_monotone(F, n, **kw):
        """Fixed point of an order-reversing map by ``F o F`` bracketing
        (Sec. V A)."""
        return _am.fixed_point(F, n, **kw)

    @staticmethod
    def emitted_field(c, minus_betaH, hs):
        """Eq. (3): the field a complex emits, for any interior Hamiltonian."""
        return _cavity.emitted_field(c, minus_betaH, hs)
