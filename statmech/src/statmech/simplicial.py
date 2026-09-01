"""The simplicial Ising model as a chygraph.

Son, Lee and Goh [Commun. Phys. (2026), arXiv:2411.19080] study an Ising model
whose hyperedges lower the energy only when *all* their members agree,

    H = - sum_e J_{|e|} delta_{{S_i}, i in e} - H sum_i S_i,

the higher-order Kronecker delta being 1 on the two unanimous states of a
hyperedge and 0 on the other ``2^q - 2``.  It produces continuous,
discontinuous, mixed-order and double transitions depending on the cardinalities
present.

Nothing about the chygraph formalism has to change to take it.  The
intra-complex step is a sum over the interior of a complex given its members'
cavity fields, and the interior is whatever Hamiltonian is put there; the
chy-degree step never sees it.  Their Bethe-Peierls equation is that sum with
the simplicial energy substituted, on a chygraph with two layers of fixed
chy-degree.

For a simplicial complex the sum closes.  With ``a = beta J_q`` and all other
members at cavity field ``h``,

    Z(S_0 = +-1) = (2 cosh h)^{q-1} + (e^a - 1) e^{+- h (q-1)},
    u = (1/2) ln[ Z(+) / Z(-) ],

so the emitted field is elementary at any cardinality, where the all-pairs
clique of Sec. IV needs an enumeration.  Differentiating at zero field,

    u' = (e^a - 1) / (2^{q-1} + e^a - 1)      (per neighbour),

which at ``q = 2`` is ``tanh(beta J / 2)`` -- the simplicial rule on a pair is
an ordinary Ising bond of half the coupling, since
``delta_{S_0 S_1} = (1 + S_0 S_1)/2``.  The multiplicity ``q-1`` is supplied
once, by the branching matrix of Sec. II C.

**One caveat, and it is the important one.**  The branching matrix of
Sec. II C locates a *linear* instability of the paramagnet.  Where the transition
is discontinuous -- which for this model is most of the interesting range --
that is the spinodal and not the transition.  The ordered branch has to be
found by iterating from a magnetised start and the true temperature by
comparing free energies, exactly as in Sec. II E.  ``det(I - B) = 0`` is still
the right object; it just answers a different question.
"""

import numpy as np


def _log2cosh(x):
    ax = np.abs(x)
    return ax + np.log1p(np.exp(-2.0 * ax))


def emitted(q, beta_J, h):
    """Field emitted onto one member of a simplicial ``q``-complex.

    ``h`` is the cavity field common to the other ``q-1`` members.
    """
    q = int(q)
    h = np.asarray(h, dtype=float)
    bulk = (q - 1) * _log2cosh(h)
    amp = np.log(np.expm1(beta_J))          # ln(e^{beta J} - 1)
    lo = np.logaddexp(bulk, amp + h * (q - 1))
    hi = np.logaddexp(bulk, amp - h * (q - 1))
    return 0.5 * (lo - hi)


def uprime(q, beta_J):
    """``du/dh_j`` at zero field for ONE other member,

        u' = (e^a - 1) / (2^{q-1} + e^a - 1),

    per neighbour, matching :func:`statmech.ising.clique_derivative`
    so that the multiplicity ``c-1`` is supplied once, by the branching matrix,
    and not twice.  Differentiating instead with respect to a field common to
    all ``q-1`` other members multiplies this by ``q-1``; that is the form which
    enters a symmetric fixed point directly, and conflating the two double-counts
    the multiplicity.  The two coincide at ``q = 2``, where ``u' = tanh(beta J/2)``.
    """
    # In log space: e^{beta J} - 1 overflows for beta J past ~700, which a
    # spinodal search reaches whenever it probes low temperature.
    a = float(beta_J)
    lnA = a + np.log1p(-np.exp(-a)) if a > 1e-8 else np.log(np.expm1(a))
    return float(np.exp(lnA - np.logaddexp((q - 1) * np.log(2.0), lnA)))


# ---------------------------------------------------------------------------
# The chygraph: layers of fixed chy-degree
# ---------------------------------------------------------------------------

class SimplicialChygraph:
    """Simplicial Ising on a Bethe hyperlattice of several cardinalities.

    Args:
        cardinalities: ``q_l`` per layer.
        degrees: ``k_l``, the number of layer-``l`` complexes each node is in.
        couplings: ``J_l`` per layer.
    """

    def __init__(self, cardinalities, degrees, couplings):
        self.q = np.asarray(cardinalities, dtype=int)
        self.k = np.asarray(degrees, dtype=float)
        self.J = np.asarray(couplings, dtype=float)
        self.L = len(self.q)

    # -- the self-consistency ----------------------------------------------

    def solve(self, T, u0=None, iters=20000, tol=1e-13, damping=0.5):
        """Fixed point of ``u_l = emitted(q_l, J_l/T, h - u_l)``, ``h = sum k u``.

        ``u0`` seeds the iteration: a large value finds the ordered branch on
        cooling, zero the paramagnetic one.  Where the two differ the
        transition is discontinuous and the pair brackets the hysteresis loop.
        """
        beta = 1.0 / T
        u = np.zeros(self.L) if u0 is None else np.full(self.L, float(u0))
        for _ in range(iters):
            h = float((self.k * u).sum())
            new = np.array([emitted(self.q[l], beta * self.J[l], h - u[l])
                            for l in range(self.L)])
            new = damping * new + (1 - damping) * u
            if np.abs(new - u).max() < tol:
                return new
            u = new
        return u

    # -- observables --------------------------------------------------------

    def magnetisation(self, T, u0=None, **kw):
        u = self.solve(T, u0, **kw)
        return float(np.tanh((self.k * u).sum())), u

    def components(self, T, u0=None, **kw):
        """``m_l`` of Eq. (10) of Ref. [SLG]: the share of ``m`` each layer
        carries, with ``m = sum_l m_l`` by the tanh addition formula."""
        u = self.solve(T, u0, **kw)
        t = np.tanh(self.k * u)
        denom = 1.0
        for i in range(self.L):
            for j in range(i + 1, self.L):
                denom += t[i] * t[j]
            # exact for two layers; the general product form is not needed here
        return t / denom, u

    # -- the free energy ----------------------------------------------------

    def minus_beta_f(self, T, u0=None, **kw):
        """``-beta f`` per node on the branch reached from ``u0``.

        The Bethe form of Sec. II E with the simplicial interior:

            Z_i = 2 cosh h,
            Z_a = (2 cosh h_cav)^q + (e^{beta J} - 1) 2 cosh(q h_cav),

        the complex term being the sum over its interior, the second piece
        carrying the two unanimous states, with ``h_cav = h - u_l`` the cavity
        field of each member.  Comparing this between branches is what locates a
        discontinuous transition; the spinodal does not.
        """
        u = self.solve(T, u0, **kw)
        beta = 1.0 / T
        h = float((self.k * u).sum())
        out = float((1.0 - self.k.sum()) * _log2cosh(h))
        for l in range(self.L):
            hc = h - u[l]
            q = int(self.q[l])
            bulk = q * _log2cosh(hc)
            # _log2cosh already carries the factor 2 of 2 cosh(q h_cav)
            amp = np.log(np.expm1(beta * self.J[l])) + _log2cosh(q * hc)
            out += (self.k[l] / q) * float(np.logaddexp(bulk, amp))
        return out

    def branch_gap(self, T, u_hi=8.0, u_lo=1e-8, **kw):
        """``m`` from a magnetised start minus ``m`` from a near-zero one.

        Non-zero exactly where two solutions coexist, which is the signature of
        a discontinuous transition.  Outside the window the two starts converge
        on the same fixed point and this is zero to machine precision -- which
        is why a free-energy comparison must be bracketed *inside* the window
        and not across it, where the difference is rounding noise of random
        sign.

        One false positive to know about: just below a *continuous* transition
        the near-zero start converges slowly, so a finite iteration budget can
        leave the two branches apart when they are really the same one.  That is
        critical slowing down, not coexistence.  Keep the scan range away from a
        continuous transition, or raise ``iters``.
        """
        a, _ = self.magnetisation(T, u0=u_hi, **kw)
        b, _ = self.magnetisation(T, u0=u_lo, **kw)
        return a - b

    def coexistence(self, lo, hi, n=120, tol=1e-4, iters=40, **kw):
        """``(T_star, T_starstar)``, the limits of metastability.

        Scans ``[lo, hi]`` for a temperature where two branches coexist, then
        bisects each edge.  Returns ``(nan, nan)`` if the scan finds none, which
        means the transition is continuous on this interval or the grid is too
        coarse to catch the window.
        """
        grid = np.linspace(lo, hi, n)
        inside = [T for T in grid if abs(self.branch_gap(T, **kw)) > tol]
        if not inside:
            return float('nan'), float('nan')
        mid = inside[len(inside) // 2]
        edges = []
        for a, b in ((mid, lo), (mid, hi)):        # a inside, b outside
            x, y = a, b
            for _ in range(iters):
                m = 0.5 * (x + y)
                x, y = (m, y) if abs(self.branch_gap(m, **kw)) > tol else (x, m)
            edges.append(0.5 * (x + y))
        return min(edges), max(edges)

    def transition(self, lo=None, hi=None, scan=None, iters=100, **kw):
        """First-order transition temperature: where the free energies cross.

        With ``scan=(a, b)`` the coexistence window is found first and used as
        the bracket, which is the only place the comparison is meaningful.
        Returns ``nan`` when no window exists.
        """
        if scan is not None:
            lo, hi = self.coexistence(*scan, **kw)
        if lo is None or hi is None or not np.isfinite(lo) or not np.isfinite(hi):
            return float('nan')
        eps = 1e-9
        lo, hi = lo + eps, hi - eps

        def gap(T):
            return (self.minus_beta_f(T, u0=8.0) - self.minus_beta_f(T, u0=1e-8))

        if gap(lo) * gap(hi) > 0:
            return float('nan')
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if gap(mid) * gap(hi) < 0 else (lo, mid)
        return 0.5 * (lo + hi)

    # -- the linear instability --------------------------------------------

    def branching_matrix(self, T):
        """``B_{lm}`` of Sec. II C with the simplicial ``u'``.

        Fixed chy-degree, so the excess is ``k_l - 1``.
        """
        beta = 1.0 / T
        up = np.array([uprime(self.q[m], beta * self.J[m])
                       for m in range(self.L)])
        B = np.empty((self.L, self.L))
        for l in range(self.L):
            for m in range(self.L):
                deg = self.k[m] - 1 if m == l else self.k[m]
                B[l, m] = deg * (self.q[m] - 1) * up[m]
        return B

    def spinodal(self, lo=1e-3, hi=50.0, iters=200):
        """Temperature where the paramagnet loses linear stability.

        The transition temperature only when the transition is continuous;
        otherwise the lower limit of metastability.
        """
        f = lambda T: np.max(np.abs(np.linalg.eigvals(
            self.branching_matrix(T)))) - 1.0
        if f(lo) < 0 or f(hi) > 0:
            return float('nan')
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
        return 0.5 * (lo + hi)
