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
