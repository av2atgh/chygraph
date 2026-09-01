"""Minimal vertex cover on a correlated random graph (WP2 demonstration).

Vazquez & Weigt, Phys. Rev. E 67, 027101 (2003), Eqs. (16)-(18).  In the
``mu -> inf`` limit the field distribution collapses onto integers and the
replica-symmetric solution reduces to one scalar per excess-degree class,

    pi_d = sum_{d'} p(d'|d) (1 - pi_{d'})^{d'}                        Eq. (16)

with ``pi_d`` the probability that an edge entering a vertex of degree ``d+1``
carries a constraint.  No field distribution survives, so this is the same shape
``chygraph.giant.Chygraph`` solves — with one difference that matters.  A
generating function is order-*preserving*, so the percolation map is monotone
and ``Chygraph.solve`` can iterate upward from ``Q = 0``.  Eq. (16) is
order-*reversing*: raising any ``pi`` lowers every ``pi``.  Plain iteration
oscillates, and past the instability it never converges — which is precisely
what Vazquez & Weigt use as their RSB detector.

This module solves it anyway.  Under the correlation model of Eq. (18) the map
closes on a single scalar, so bisecting on that scalar reaches the fixed point
whether or not it is stable, and the stability question is then asked separately
of the Jacobian there.  That separation is what WP2 buys: a fixed point that is
found, and a diagnosis of it that is computed, rather than one instrument doing
both by failing to converge.
"""

import numpy as np


# ---------------------------------------------------------------------------
# Ensembles
# ---------------------------------------------------------------------------

def scale_free(gamma, dmax):
    """``p_d ~ d^-gamma`` for ``d = 1..dmax``; no isolated vertices, as in VW03."""
    d = np.arange(dmax + 1, dtype=float)
    p = np.zeros(dmax + 1)
    p[1:] = d[1:] ** (-gamma)
    return p / p.sum()


def poisson(c, dmax):
    from scipy.special import gammaln
    d = np.arange(dmax + 1)
    p = np.exp(-c + d * np.log(c) - gammaln(d + 1))
    return p / p.sum()


def excess(p):
    """Degree distribution -> excess degrees ``e`` and ``q_e = (e+1)p_{e+1}/c``."""
    d = np.arange(len(p))
    c = float((d * p).sum())
    e = np.arange(len(p) - 1)
    return e, (e + 1) * p[e + 1] / c, c


# ---------------------------------------------------------------------------
# The map
# ---------------------------------------------------------------------------

def solve(r, e, q, iters=80):
    """Fixed point of Eq. (16) under the correlations of Eq. (18).

    ``p(d'|d) = r delta_{dd'} + (1-r) q_{d'}`` reduces Eq. (16) to

        pi_e = r (1 - pi_e)^e + (1 - r) A,   A = sum_e q_e (1 - pi_e)^e

    so the whole map closes on the scalar ``A``.  Both the inner and outer
    relations are monotone decreasing, so nested bisection converges to the
    unique fixed point regardless of its stability.
    """
    def pi_of_A(A):
        lo, hi = np.zeros_like(e, float), np.ones_like(e, float)
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            f = r * (1 - mid) ** e + (1 - r) * A - mid
            lo, hi = np.where(f > 0, mid, lo), np.where(f > 0, hi, mid)
        return 0.5 * (lo + hi)

    lo, hi = 0.0, 1.0
    for _ in range(iters):
        A = 0.5 * (lo + hi)
        if (q * (1 - pi_of_A(A)) ** e).sum() - A > 0:
            lo = A
        else:
            hi = A
    return pi_of_A(0.5 * (lo + hi))


def cover_size(p, pi):
    """VW03 Eq. (17): the fraction of vertices in a minimal cover.

    Reduces to the Weigt-Hartmann ``1 - (2W + W^2)/(2c)`` for uncorrelated
    Poisson graphs, with ``W = LambertW(c)``; that identity is the test.
    """
    d = np.arange(len(p))
    dd, pp, pm = d[1:], p[1:], pi[d[1:] - 1]
    return 1.0 - (p[0] + (pp * (1 - pm) ** (dd - 1)
                          * (1 + (dd - 2) / 2 * pm)).sum())


# ---------------------------------------------------------------------------
# Stability of that fixed point
# ---------------------------------------------------------------------------

def jacobian_spectrum(r, e, q, pi):
    """Leading eigenvalue of ``J = d pi / d pi`` at the fixed point.

    ``J_{d d'} = -[r delta + (1-r) q_{d'}] a_{d'}`` with
    ``a_e = e (1 - pi_e)^{e-1}`` — diagonal plus rank one, so the spectrum
    follows from a secular equation and no eigensolve is needed:

        eigenvalues  lambda = -r a_e,  and the roots of
        (1 - r) sum_e q_e a_e / (-lambda - r a_e) = 1.

    Every entry of ``J`` is non-positive: the map is anti-monotone, the leading
    eigenvalue is negative, and the instability at ``-1`` is a period-2 orbit,
    not an exchange of fixed points.
    """
    a = np.where(e > 0, e * (1 - pi) ** np.maximum(e - 1, 0), 0.0)
    return a


def is_unstable(r, e, q, pi):
    """``|lambda| > 1``: RS is broken.  VW03's convergence failure, computed."""
    a = jacobian_spectrum(r, e, q, pi)
    if r * a.max() >= 1.0:
        return True
    return bool((1 - r) * (q * a / (1 - r * a)).sum() > 1.0)


def leading_eigenvalue(r, e, q, pi):
    """The signed leading eigenvalue, by dense eigensolve.

    Slower than :func:`is_unstable` and used to cross-check it.
    """
    a = jacobian_spectrum(r, e, q, pi)
    J = -r * np.diag(a) - (1 - r) * np.outer(np.ones_like(a), q * a)
    w = np.linalg.eigvals(J)
    return complex(w[np.argmax(np.abs(w))])


def rsb_point(gamma, dmax=800, iters=50):
    """The ``r`` at which RS breaks, by bisection."""
    e, q, _ = excess(scale_free(gamma, dmax))
    lo, hi = 0.0, 1.0
    for _ in range(iters):
        r = 0.5 * (lo + hi)
        lo, hi = ((lo, r) if is_unstable(r, e, q, solve(r, e, q)) else (r, hi))
    return 0.5 * (lo + hi)
