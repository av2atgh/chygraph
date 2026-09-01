"""One object for one chygraph, with a method per calculation.

The modules of this package each solve one problem and each grew from one
section of the paper.  :class:`Chygraph` is the single handle on all of them: it
holds the *structure* --- the complex layers, their cardinalities, and the
chy-degree distribution --- and its methods are the calculations, delegating to
those modules rather than reimplementing them.

    >>> from statmech import Chygraph
    >>> g = Chygraph([2, 3], [4.0, 2.0])          # links and triangles
    >>> g.critical_coupling()                      # Ising T_c, Sec. IV
    >>> g.core().core_fraction()                   # leaf-removal core, Sec. VI
    >>> g.hitting_set_bp(mu=60).run().density()    # hitting set, Sec. V E

A chygraph is specified by ``cardinalities`` and ``degrees``; ``excess``
defaults to ``degrees`` (Poisson chy-degree, its own excess) and is
``degrees - 1`` when ``regular=True``.  Where a calculation needs the full
generating function rather than its moments, ``phi`` may be given as a sympy
expression in :func:`~statmech.hittingset.layer_symbols`; otherwise the
independent-Poisson form is built from ``degrees``.

Which method gives which result is listed in the paper; the docstrings name the
section each belongs to.
"""

import numpy as np

from statmech import antimonotone as _am
from statmech import cavity as _cavity
from statmech import core as _core
from statmech import cover as _cover
from statmech import freeenergy as _fe
from statmech import hittingset as _hs
from statmech import ising as _ising
from statmech import population as _pop
from statmech import region as _region
from statmech import simplicial as _simp
from statmech import softfield as _soft
from statmech import stability as _stab


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
        # A layer's cardinality may be a single integer or a distribution given
        # as {c: weight}.  The general branching matrix needs the size-biased
        # average over that distribution; a single integer is the special case
        # where it collapses to c - 1.
        self.cdist = [self._as_dist(c) for c in cardinalities]
        self.c = np.array([self._mean_c(d) if len(d) > 1
                           else next(iter(d)) for d in self.cdist])
        if len(self.cdist) == 1 and not isinstance(cardinalities, (list, tuple,
                                                                   np.ndarray)):
            raise ValueError("cardinalities must be a sequence, one per layer")
        self.k = np.asarray(degrees, dtype=float)
        if len(self.cdist) != self.k.size:
            raise ValueError("one cardinality and one chy-degree per layer")
        if any(min(d) < 2 for d in self.cdist):
            raise ValueError("cardinalities must be at least 2")
        self.L = len(self.cdist)
        self.regular = bool(regular)
        if excess is not None:
            self.kbar = np.asarray(excess, dtype=float)
        elif regular:
            self.kbar = np.maximum(self.k - 1.0, 0.0)
        else:
            self.kbar = self.k.copy()
        self._phi = phi

    @staticmethod
    def _as_dist(c):
        """Normalise a layer's cardinality to ``{c: probability}``."""
        if isinstance(c, dict):
            tot = float(sum(c.values()))
            if tot <= 0:
                raise ValueError("cardinality weights must be positive")
            return {int(k): float(v) / tot for k, v in c.items() if v > 0}
        return {int(c): 1.0}

    @staticmethod
    def _mean_c(d):
        return sum(c * p for c, p in d.items())

    def size_biased(self, layer, f):
        """``sum_c (c p_c / <c>) f(c)``, the average over a complex reached by
        following one of its inclusions.

        Arriving at a complex from a member samples it in proportion to its
        cardinality, so every per-complex quantity in Eq.~(8) is weighted this
        way.  For a single-cardinality layer it is just ``f(c)``.
        """
        d = self.cdist[layer]
        cbar = self._mean_c(d)
        return sum(c * p / cbar * f(c) for c, p in d.items())

    def excess_cardinality_layer(self, layer):
        """``<sbar>_m = <c^2>_m/<c>_m - 1``: the mean number of *other* members
        seen on arrival.  Equals ``c - 1`` when the layer has one cardinality."""
        return self.size_biased(layer, lambda c: c - 1.0)

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
        list (Sec. VII).  Returns :class:`~statmech.region.RegionGraph`."""
        return _region.RegionGraph(complexes)

    @staticmethod
    def overlap_profile(complexes):
        """How far an explicit complex family is from treelike (Sec. VII)."""
        return _region.overlap_profile(complexes)

    @staticmethod
    def region_gbp(complexes, beta_J, damping=0.5, field=0.0, **kw):
        """Generalised belief propagation on the region graph (Sec. VII).

        Takes an explicit list of complexes and the Ising coupling of Eq. (16),
        builds the closed region family and its Mobius counting numbers, and
        returns a :class:`~statmech.gbp.GBP` ready to ``run()``.  Each
        pair inside a complex is counted once however many complexes contain
        it, so overlapping complexes do not double the coupling on a shared
        bond.

        This is the instance-level calculation; every other method on this
        class works at the level of an ensemble.
        """
        from statmech import gbp as _gbp
        rg = _region.RegionGraph(complexes)
        factors = _gbp.ising_factors(_gbp.clique_edges(complexes), beta_J,
                                     field=field)
        return _gbp.GBP(rg, factors, damping=damping, **kw)

    # ======================================================================
    # Sec. IV: the Ising model
    # ======================================================================

    def _u_of_c(self, c, beta_J, interaction):
        if interaction == 'simplicial':
            return _simp.uprime(int(c), beta_J)
        return _ising.clique_derivative(int(c), beta_J)

    def u_prime(self, layer, beta_J, interaction='clique'):
        """``u'``, the cavity derivative per neighbour (Sec. IV A).

        ``interaction='clique'`` couples every pair inside the complex;
        ``'simplicial'`` is the unanimity rule of Sec. IV D.  For a layer with a
        distribution of cardinalities this returns the size-biased average of
        ``u'(c)``; the transmission that enters Eq.~(8) is
        :meth:`transmission`, which weights by ``c-1`` as well.
        """
        return self.size_biased(layer,
                                lambda c: self._u_of_c(c, beta_J, interaction))

    def transmission(self, layer, beta_J, interaction='clique', squared=False):
        """``<sbar u'>_m``, the entry Eq.~(8) actually needs.

        The size-biased average of ``(c-1) u'(c)`` over the layer's cardinality
        distribution, which is *not* ``<sbar> <u'>`` unless ``u'`` is common to
        the layer.  Collapses to ``(c-1) u'`` for a single cardinality.
        """
        def f(c):
            u = self._u_of_c(c, beta_J, interaction)
            return (c - 1.0) * (u * u if squared else u)
        return self.size_biased(layer, f)

    def branching_matrix(self, beta_J, interaction='clique', squared=False):
        """``B_{lm}`` of Eq.~(8) (Sec. II C).  ``squared`` gives the AT line.

        General in the layer's cardinality distribution: the per-complex factor
        is the size-biased ``<sbar u'>_m`` of :meth:`transmission`, which
        reduces to ``(c_m - 1) u'_m`` when the layer carries one cardinality.
        """
        bj = np.broadcast_to(np.asarray(beta_J, float), (self.L,))
        w = np.array([self.transmission(m, bj[m], interaction, squared)
                      for m in range(self.L)])
        B = np.empty((self.L, self.L))
        for l in range(self.L):
            for m in range(self.L):
                deg = self.kbar[m] if m == l else self.k[m]
                B[l, m] = deg * w[m]
        return B

    def critical_coupling(self, interaction='clique', squared=False,
                          lo=1e-9, hi=20.0):
        """``beta J`` where the Perron root of ``B`` reaches 1 (Sec. IV B).

        Computed from :meth:`branching_matrix`, so a layer carrying a
        distribution of cardinalities is handled by the size-biased
        ``<sbar u'>`` and not by substituting a mean cardinality.

        The transition when it is continuous; the spinodal when it is not.
        """
        if interaction != 'clique':
            raise NotImplementedError(
                "use simplicial().spinodal() for the unanimity interaction")

        def gap(bj):
            B = self.branching_matrix(bj, interaction, squared)
            return float(np.max(np.abs(np.linalg.eigvals(B)))) - 1.0

        if gap(lo) * gap(hi) > 0:
            raise ValueError(
                f"no transition on beta J in [{lo}, {hi}]: "
                f"Perron root runs {gap(lo)+1:.4g} to {gap(hi)+1:.4g}")
        from scipy.optimize import brentq
        return brentq(gap, lo, hi, xtol=1e-14, rtol=8.9e-16)

    def critical_temperature(self, J=1.0, **kw):
        """``T_c / J`` for the continuous case (Sec. IV B)."""
        return 1.0 / (J * self.critical_coupling(**kw))

    def stability_tensor(self, k, K, s, S, wkappa=None, ws=None):
        """The symbolic ``2L^2`` tensor of Sec. II C, reweighted (Eq. 7).

        Takes the four moment tables directly, as
        :class:`~statmech.stability.StabilityMatrix` does.
        """
        return _stab.StabilityMatrix(k, K, s, S, wkappa=wkappa, ws=ws)

    def population(self, beta_J, **kw):
        """Field distributions by population dynamics (Sec. IV C).

        Returns :class:`~statmech.population.CavityPopulation`.
        """
        return _pop.CavityPopulation(self.c, self.k, beta_J, **kw)

    def free_energy(self, beta_J, sweeps=250, **kw):
        """Bethe free energy of the Ising model, Eq. (13) (Sec. II E).

        Returns :class:`~statmech.freeenergy.BetheFreeEnergy` on a
        converged population.
        """
        return _fe.BetheFreeEnergy(self.population(beta_J, **kw).run(sweeps))

    def paramagnetic_free_energy(self, beta_J):
        """``-beta f`` at zero cavity field, Eq. (14) (Sec. II E)."""
        return _fe.paramagnetic(self.c, self.k, beta_J)

    # ======================================================================
    # Sec. IV D: the simplicial (unanimity) interaction
    # ======================================================================

    def simplicial(self, couplings):
        """The simplicial Ising model on this chygraph (Sec. IV D).

        Returns :class:`~statmech.simplicial.SimplicialChygraph`,
        which carries ``spinodal``, ``coexistence``, ``transition`` and
        ``minus_beta_f``.  Fixed chy-degree is assumed, as in Ref. [SLG26].
        """
        return _simp.SimplicialChygraph(self.c, self.k, couplings)

    # ======================================================================
    # Sec. V: hitting set and vertex cover
    # ======================================================================

    def hitting_set(self):
        """Hard-field hitting set, Eq. (26) (Sec. V A).

        Thresholds are exact; the cover size only at cardinality two --- see
        :meth:`hitting_set_bp`.
        """
        return _hs.HittingSet(self.c, self.phi())

    def hitting_set_bp(self, mu=60.0, **kw):
        """Hitting set with the ``O(1)`` fields kept, Eq. (31) (Sec. V E).

        Returns :class:`~statmech.softfield.HittingSetBP`; call
        ``.run()`` before reading ``density`` or ``entropy``.
        """
        return _soft.HittingSetBP(self.c, self.k, regular=self.regular, mu=mu,
                                  **kw)

    def clique_cover(self):
        """Vertex cover of the induced graph, Eq. (30) (Sec. V C)."""
        return _cover.CliqueCover(self.c, self.phi())

    # ======================================================================
    # Sec. VI: core percolation
    # ======================================================================

    def core(self):
        """Leaf-removal core as a chygraph fixed point, Eq. (35) (Sec. VI).

        Returns :class:`~statmech.core.CorePercolation`.
        """
        return _core.CorePercolation(self.c, self.phi())

    @classmethod
    def from_samples(cls, cardinalities, memberships):
        """Build from a *measured* chy-degree ensemble (Sec. VII).

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

        Takes a ``percolation.giant.Chygraph`` and returns
        :class:`~statmech.fixedpoint.FixedPointStability`, whose
        ``spectral_radius`` obeys ``rho = 1 - Lambda + O(Lambda^2)`` and whose
        ``monotonicity`` reads entry signs rather than the eigenvalue sign.
        """
        from statmech import fixedpoint as _fp
        return _fp.FixedPointStability(model)

    @staticmethod
    def vertex_cover_correlated(gamma, r, dmax=800):
        """Vertex cover on a scale-free graph with the degree correlations of
        Ref. [VW03] Eq. (18) (Sec. V B).  Returns ``(cover_size, unstable)``."""
        from statmech import vertexcover as _vc
        p = _vc.scale_free(gamma, dmax)
        e, q, _ = _vc.excess(p)
        pi = _vc.solve(r, e, q)
        return _vc.cover_size(p, pi), _vc.is_unstable(r, e, q, pi)

    @staticmethod
    def excess_cardinality(complexes):
        """``<c^2>/<c> - 1`` of a complex family (Sec. VII A).

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
        """Eq. (5): the field a complex emits, for any interior Hamiltonian."""
        return _cavity.emitted_field(c, minus_betaH, hs)
