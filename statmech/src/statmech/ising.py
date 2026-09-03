"""Critical temperature of the Ising model on an arbitrary chygraph.

WP1 built the stability tensor and applied it to two hand-written cases, a
graph and a graph with triangles.  This closes the general case.

Arrive at a node through a layer-``l`` complex.  It belongs to
``<kappabar>_l`` further layer-``l`` complexes and ``<kappa>_m`` layer-``m``
complexes for ``m != l``; each layer-``m`` complex offers ``c_m - 1`` other
members, and each of those receives the field with derivative ``u'_m``.  So the
branching matrix on the ``L`` complex layers is

    B_{lm} = [ <kappabar>_m if m == l else <kappa>_m ] (c_m - 1) u'_m

and the transition is at ``det(I - B) = 0``, i.e. Perron root 1.  For the
ferromagnetic line ``u'_m`` is the cavity derivative of
:mod:`statmech.cavity`; squaring it gives the de Almeida-Thouless line
instead, exactly as in WP1.

``u'`` has no useful closed form beyond a triangle --

    c = 2:  t
    c = 3:  t / (1 - t + t^2)
    c = 4:  (t + t^3) / (1 - 2t + 3t^2)

with ``t = tanh(beta J)`` -- but it is exact by enumeration for any cardinality
small enough to enumerate, which is the same condition the rest of the package
already lives under.  So the answer to "can this give a critical temperature"
is yes, for any chygraph whose complexes can be enumerated: the linear algebra
is ``L x L`` however large the complexes are.
"""

import numpy as np
from scipy.optimize import brentq


def clique_derivative(c, beta_J):
    """``u'`` at zero field for a ``c``-clique, numerically.

    The same enumeration as :func:`statmech.cavity.cavity_derivative`,
    evaluated rather than kept symbolic so it can be swept over temperature.
    """
    m = c - 1
    bits = ((np.arange(2 ** m)[:, None] >> np.arange(m)) & 1)
    cfg = 1.0 - 2.0 * bits
    S1 = cfg.sum(axis=1)
    base = beta_J * (S1 ** 2 - m) / 2.0
    # u = (1/2)[ ln Z(+) - ln Z(-) ]; differentiating in h_1 at h = 0 gives the
    # difference of <sigma_1> under the two weightings.  Analytic, so no
    # finite-difference error enters the critical coupling.
    out = []
    for sign in (+1.0, -1.0):
        w = base + sign * beta_J * S1
        w = w - w.max()
        e = np.exp(w)
        out.append(float((e * cfg[:, 0]).sum() / e.sum())) 
    return 0.5 * (out[0] - out[1])


def branching_matrix(cardinalities, means, beta_J, excess=None, squared=False):
    """``B_{lm}``, the layer-to-layer branching matrix.

    Args:
        cardinalities: ``c_l`` per complex layer.
        means: ``<kappa>_l``, mean chy-degree to layer ``l``.
        excess: ``<kappabar>_l``; defaults to ``means`` (Poisson).
        squared: use ``u'^2`` for the de Almeida-Thouless line.
    """
    c = np.asarray(cardinalities, dtype=int)
    k = np.asarray(means, dtype=float)
    K = k if excess is None else np.asarray(excess, dtype=float)
    L = len(c)
    u = np.array([clique_derivative(int(c[m]), beta_J) for m in range(L)])
    if squared:
        u = u ** 2
    B = np.empty((L, L))
    for l in range(L):
        for m in range(L):
            deg = K[m] if m == l else k[m]
            B[l, m] = deg * (c[m] - 1) * u[m]
    return B


def perron_root(cardinalities, means, beta_J, excess=None, squared=False):
    B = branching_matrix(cardinalities, means, beta_J, excess, squared)
    return float(np.max(np.abs(np.linalg.eigvals(B))))


def critical_coupling(cardinalities, means, excess=None, squared=False,
                      lo=1e-9, hi=20.0):
    """``beta J`` at which the Perron root of ``B`` reaches 1.

    ``squared=True`` returns the de Almeida-Thouless line instead of the
    ferromagnetic one.  Raises if the bracket contains no transition, which
    happens when the chygraph is too sparse to order at any temperature.
    """
    def gap(bj):
        return perron_root(cardinalities, means, bj, excess, squared) - 1.0

    if gap(lo) * gap(hi) > 0:
        raise ValueError(
            f"no transition on beta J in [{lo}, {hi}]: "
            f"Perron root runs {gap(lo)+1:.4g} to {gap(hi)+1:.4g}")
    return brentq(gap, lo, hi, xtol=1e-14, rtol=8.9e-16)


def critical_temperature(cardinalities, means, excess=None, J=1.0):
    """``T_c / J``, the reciprocal of :func:`critical_coupling`."""
    return 1.0 / (J * critical_coupling(cardinalities, means, excess))


# ---------------------------------------------------------------------------
# The order parameter and the response
# ---------------------------------------------------------------------------
#
# ``critical_coupling`` above linearises the recursion and returns where the
# paramagnetic fixed point loses stability.  What follows carries the same
# expansion one order further, which is the step Ch. 5 takes for percolation:
# the threshold needs the first derivative of the interior sum, the critical
# amplitude needs the third, and the order parameter away from the transition
# needs the whole field distribution.
#
# Two objects do all of it.  ``interior_cumulants`` returns the first and third
# Taylor coefficients of the field a complex emits when every incoming field is
# the same number, and ``emitted_common`` is that field itself.


def emitted_common(c, beta_J, h):
    """``ubar(h)``: field emitted when all ``c-1`` incoming fields equal ``h``.

    The enumeration of :func:`statmech.cavity.emitted_field` with its arguments
    tied together, which is the only case a *regular* chygraph ever needs: every
    message on a layer is then the same number by symmetry.  Vectorised in ``h``.
    """
    from scipy.special import logsumexp
    m = int(c) - 1
    bits = ((np.arange(2 ** m)[:, None] >> np.arange(m)) & 1)
    cfg = 1.0 - 2.0 * bits
    S1 = cfg.sum(axis=1)
    # -beta H = bJ (S2 + sigma_0 S1) with S2 = (S1^2 - m)/2 the interior pair
    # sum of the other members; sigma_0 = +1 here and Z(-1) is Z(+1) at -h.
    base = beta_J * (S1 ** 2 - m) / 2.0 + beta_J * S1
    h = np.asarray(h, dtype=float)
    g = lambda x: logsumexp(base + np.asarray(x)[..., None] * S1, axis=-1)
    return 0.5 * (g(h) - g(-h))


def emitted_common_derivative(c, beta_J, h):
    """``d ubar / dh``, exactly.

    ``ubar(h) = [g(h) - g(-h)]/2`` with ``g`` the cumulant generating function
    of the interior magnetisation ``S1``, so the derivative is
    ``[<S1>_h + <S1>_{-h}]/2`` and needs no finite difference.
    """
    m = int(c) - 1
    bits = ((np.arange(2 ** m)[:, None] >> np.arange(m)) & 1)
    cfg = 1.0 - 2.0 * bits
    S1 = cfg.sum(axis=1)
    base = beta_J * (S1 ** 2 - m) / 2.0 + beta_J * S1

    def mean(x):
        w = base + x * S1
        w = np.exp(w - w.max())
        return float((w * S1).sum() / w.sum())

    return 0.5 * (mean(h) + mean(-h))


def interior_cumulants(c, beta_J):
    """``(a1, a3)`` in ``ubar(h) = a1 h + a3 h^3 + O(h^5)``.

    Both are cumulants of the interior magnetisation.  Writing ``S1`` for the
    total spin of the other ``c-1`` members under the Gibbs weight of the
    isolated complex with ``sigma_0 = +1``,

        a1 = <S1>,      a3 = kappa_3(S1) / 6,

    because ``ubar(h) = [g(h) - g(-h)]/2`` with ``g = ln Z(+1)`` the cumulant
    generating function of ``S1``.  So ``a1 = (c-1) u'`` --- the transmission
    that Eq. (8.8) already needs --- and the amplitude below needs one further
    odd cumulant of the same interior and nothing else.
    """
    m = int(c) - 1
    bits = ((np.arange(2 ** m)[:, None] >> np.arange(m)) & 1)
    cfg = 1.0 - 2.0 * bits
    S1 = cfg.sum(axis=1)
    w = beta_J * (S1 ** 2 - m) / 2.0 + beta_J * S1
    w = np.exp(w - w.max())
    w = w / w.sum()
    m1 = float((w * S1).sum())
    m2 = float((w * S1 ** 2).sum())
    m3 = float((w * S1 ** 3).sum())
    return m1, (m3 - 3.0 * m2 * m1 + 2.0 * m1 ** 3) / 6.0


def _w_vector(cardinalities, means, beta_J):
    """``w_l = <kappa>_l (c_l - 1) u'_l``: what one node's layer-l complexes
    deliver to it, with the ordinary chy-degree rather than the excess."""
    c = np.asarray(cardinalities, dtype=int)
    k = np.asarray(means, dtype=float)
    return np.array([k[l] * (c[l] - 1) * clique_derivative(int(c[l]), beta_J)
                     for l in range(len(c))])


def susceptibility(cardinalities, means, beta_J, excess=None):
    """``chi = dm/dB`` at ``B -> 0`` in the paramagnetic phase.

    At small external field every message is ``O(B)``, so the down step is
    ``u = u' sum_j h_j`` and the up step adds one ``B`` at the node.  Averaging,
    the field sent into a layer-``m`` complex obeys

        hhat_m = B + sum_l B_{ml} hhat_l,

    with ``B`` the *same* branching matrix that gives the transition.  Hence
    ``hhat = B (I - B)^{-1} 1``, and the full field at a node adds the ordinary
    rather than the excess chy-degree, giving

        chi = 1 + w . (I - B)^{-1} 1.

    It diverges exactly where ``det(I - B) = 0``, which is the transition: the
    threshold condition and the divergence of the response are one statement.

    Valid for any chy-degree distribution, since only the linearisation enters.
    """
    B = branching_matrix(cardinalities, means, beta_J, excess)
    L = B.shape[0]
    w = _w_vector(cardinalities, means, beta_J)
    return 1.0 + float(w @ np.linalg.solve(np.eye(L) - B, np.ones(L)))


def _perron_pair(B):
    """Right and left Perron vectors of ``B``, normalised so ``l . r = 1``."""
    ev, R = np.linalg.eig(B)
    i = int(np.argmax(ev.real))
    r = np.real(R[:, i])
    r = r / np.sign(r[int(np.argmax(np.abs(r)))])
    ev2, Lv = np.linalg.eig(B.T)
    j = int(np.argmax(ev2.real))
    l = np.real(Lv[:, j])
    return r, l / float(l @ r)


def _dlambda_dT(cardinalities, degrees, Tc, excess=None):
    """``|d lambda / dT|`` at ``T_c``, by a central difference in temperature."""
    def lam(T):
        B = branching_matrix(cardinalities, degrees, 1.0 / T, excess)
        return float(np.max(np.abs(np.linalg.eigvals(B))))
    e = 1e-5 * Tc
    return abs((lam(Tc + e) - lam(Tc - e)) / (2.0 * e))


def susceptibility_amplitude(cardinalities, degrees, excess=None):
    """``C`` in ``chi ~ C / (T/T_c - 1)`` as ``T -> T_c`` from above.

    Near the transition ``(I - B)^{-1}`` is dominated by the Perron projector
    ``r l^T / (1 - lambda)``, so ``chi ~ (w . r)(l . 1) / (1 - lambda)`` and
    ``1 - lambda`` vanishes linearly in ``T - T_c``.  The exponent is therefore
    ``gamma = 1`` whenever ``d lambda / dT`` is finite and non-zero, which it is
    for any enumerable interior.
    """
    Tc = 1.0 / critical_coupling(cardinalities, degrees, excess)
    B = branching_matrix(cardinalities, degrees, 1.0 / Tc, excess)
    r, l = _perron_pair(B)
    w = _w_vector(cardinalities, degrees, 1.0 / Tc)
    return float((w @ r) * (l.sum()) / (Tc * _dlambda_dT(
        cardinalities, degrees, Tc, excess)))


def _excess_matrix(degrees):
    """``N_{ml} = kappa_l - delta_{lm}``: how many layer-``l`` complexes a node
    reached through a layer-``m`` one still has to send to."""
    k = _regular_degrees(degrees)
    return k[None, :] - np.eye(len(k))


def _regular_degrees(degrees):
    """The chy-degrees of a regular chygraph, checked to be whole numbers.

    The scalar closure below is the statement that every message on a layer is
    the same number, which needs every node to lie in exactly ``kappa_l``
    layer-``l`` complexes.  A fractional chy-degree describes a mixture, whose
    messages carry a distribution, and the closure says nothing about it.
    """
    k = np.asarray(degrees, dtype=float)
    if np.any(np.abs(k - np.rint(k)) > 1e-12):
        raise ValueError(
            "a regular chygraph needs whole-number chy-degrees; a mixture "
            "carries a distribution of fields, for which the route is "
            "statmech.population.CavityPopulation")
    return k


def magnetisation(cardinalities, degrees, beta_J, field=0.0, tol=1e-13,
                  steps=200):
    """``m`` on a *regular* chygraph, exactly.

    The cavity construction gives the single-site partition function.  Cutting
    node ``i`` out leaves its complexes independent, each contributing a factor
    ``exp(u_{b->i} sigma_i)``, so

        Z_i = sum_sigma exp(H_i sigma) = 2 cosh H_i,
        H_i = B + sum_{b in i} u_{b->i},

    and the magnetisation is the first derivative of that free energy in the
    field, ``m_i = d ln Z_i / dB = tanh H_i``.

    On a regular chygraph every node lies in exactly ``kappa_l`` layer-``l``
    complexes, so by symmetry every message on a layer is one number and the
    recursion closes on ``L`` scalars:

        h_m = B + sum_l (kappa_l - delta_{lm}) ubar_l(h_l),
        m   = tanh( B + sum_l kappa_l ubar_l(h_l) ).

    No population is needed.  Solved by Newton from above rather than by
    iteration: ``ubar`` is increasing and concave on ``h >= 0``, so Newton
    started above the largest root descends monotonically onto it, and does so
    quadratically where plain iteration slows to the critical rate and stops
    early with an answer that is too large.  At ``field=0`` above the transition
    the only fixed point is ``h = 0`` and this returns ``0.0``.

    Started from saturation, Newton descends onto the branch continuously
    connected to it.  Below the transition with ``field < 0`` that branch is the
    metastable one -- the state a field sweep follows rather than the
    equilibrium state -- so ``m`` is odd in the field above the transition and
    hysteretic below it.  Nothing in Sec. 9.4 needs the second case: the
    magnetisation is taken at zero field and the susceptibility at ``B -> 0``
    above the transition.
    """
    c = np.asarray(cardinalities, dtype=int)
    k = np.asarray(degrees, dtype=float)
    N, L = _excess_matrix(k), len(c)
    B = float(field)

    def emitted(h):
        return np.array([float(emitted_common(int(c[l]), beta_J, h[l]))
                         for l in range(L)])

    def slope(h):
        return np.array([emitted_common_derivative(int(c[l]), beta_J, h[l])
                         for l in range(L)])

    # ubar saturates at (c-1) beta_J, so N @ that bounds every root from above
    h = N @ ((c - 1.0) * beta_J) + abs(B) + 1.0
    for _ in range(steps):
        r = h - B - N @ emitted(h)
        if np.max(np.abs(r)) < tol:
            break
        step = np.linalg.solve(np.eye(L) - N * slope(h)[None, :], r)
        hn = h - step
        if B >= 0:
            hn = np.maximum(hn, 0.0)
        if np.max(np.abs(hn - h)) < tol:
            h = hn
            break
        h = hn
    if B == 0.0 and np.max(np.abs(h)) < 1e-9:
        return 0.0
    return float(np.tanh(B + k @ emitted(h)))


def magnetisation_amplitude(cardinalities, degrees):
    """``A`` in ``m ~ A (1 - T/T_c)^(1/2)`` on a regular chygraph.

    Third order in the same expansion the threshold linearises.  Writing
    ``h = eps r`` along the Perron direction and projecting on ``l``,

        eps^2 = (lambda - 1) / D,    D = - l . N diag(a3) r^3,

    and ``m = (w . r) eps`` to leading order, so with ``lambda - 1`` vanishing
    linearly in ``T_c - T`` the exponent is ``beta = 1/2`` and

        A = (w . r) sqrt( T_c |d lambda / dT| / D ).

    Only ``a3`` is new: the threshold already needs ``a1 = (c-1) u'``.
    """
    c = np.asarray(cardinalities, dtype=int)
    k = np.asarray(degrees, dtype=float)
    kbar = np.maximum(k - 1.0, 0.0)
    bc = critical_coupling(cardinalities, k, excess=kbar)
    Tc, N = 1.0 / bc, _excess_matrix(k)
    a1, a3 = zip(*(interior_cumulants(int(x), bc) for x in c))
    a1, a3 = np.asarray(a1), np.asarray(a3)
    B = N * a1[None, :]
    r, l = _perron_pair(B)
    D = -float(l @ (N @ (a3 * r ** 3)))
    if D <= 0:
        raise ValueError(
            "no continuous branch at the reported T_c: either the cubic term "
            "does not saturate, so the transition is not second order, or the "
            "structure has no transition at all and critical_coupling has "
            "returned the edge of its bracket (a two-regular graph does this)")
    w = k * a1
    return float((w @ r) * np.sqrt(Tc * _dlambda_dT(c, k, Tc, kbar) / D))
