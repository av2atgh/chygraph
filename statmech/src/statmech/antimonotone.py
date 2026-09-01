"""Solving order-reversing message maps (WP3).

``chygraph.giant.Chygraph.solve`` iterates upward from ``Q = 0`` and relies on
the map being order-*preserving*: a probability generating function has
non-negative coefficients, so raising any message raises every message, and the
iteration converges monotonically to the smallest fixed point.

Hard-core problems invert that.  In the ``mu -> inf`` limit of a vertex cover or
a hitting set, raising any message *lowers* every message: covering more of a
hyperedge's members makes the remaining member less constrained.  Plain
iteration then oscillates, and past the replica-symmetry-breaking point it never
settles — which is exactly what Vazquez & Weigt use as their RSB detector.

For an order-reversing ``F`` on ``[0,1]^n``, ``G = F o F`` is order-preserving.
Iterating ``G`` from ``0`` increases to its least fixed point ``a``, from ``1``
decreases to its greatest ``b``, and ``F(a) = b``, ``F(b) = a``.  Two outcomes:

* ``a == b``   — one fixed point, stable, RS.
* ``a < b``    — a period-2 orbit; the fixed point of ``F`` lies strictly
                 between ``a`` and ``b`` and is unstable.  RS is broken.

So the bracket both *finds* the solution and *diagnoses* it, and it turns
"the program fails to converge" into two numbers that say by how much.
"""

import numpy as np
from scipy.optimize import brentq, root


def bracket(F, n, iters=20000, tol=1e-14):
    """Least and greatest fixed points of ``F o F`` on ``[0,1]^n``.

    Returns ``(a, b)`` with ``a <= b`` elementwise.  They coincide exactly when
    ``F`` has a unique, stable fixed point.
    """
    a, b = np.zeros(n), np.ones(n)
    for _ in range(iters):
        na, nb = F(F(a)), F(F(b))
        if max(np.abs(na - a).max(), np.abs(nb - b).max()) < tol:
            a, b = na, nb
            break
        a, b = na, nb
    return np.minimum(a, b), np.maximum(a, b)


def is_replica_symmetric(F, n, tol=1e-9, **kw):
    """True when the ``F o F`` bracket closes: one fixed point, stable."""
    a, b = bracket(F, n, **kw)
    return bool(np.max(b - a) < tol)


def fixed_point(F, n, tol=1e-13, **kw):
    """The fixed point of ``F``, stable or not.

    Inside the bracket ``[a, b]`` from :func:`bracket`, which contains it in
    either case.  For ``n = 1`` this is a sign change of ``F(x) - x`` and
    Brent's method applies; above that, a root solve started from the bracket
    midpoint.
    """
    a, b = bracket(F, n, **kw)
    if n == 1:
        lo, hi = float(a[0]), float(b[0])
        if hi - lo < tol:
            return np.array([0.5 * (lo + hi)])
        return np.array([brentq(lambda x: F(np.array([x]))[0] - x, lo, hi,
                                xtol=1e-15, rtol=8.9e-16)])
    sol = root(lambda x: F(x) - x, 0.5 * (a + b), tol=tol)
    return np.clip(sol.x, 0.0, 1.0)


def is_anti_monotone(F, n, trials=8, seed=0):
    """Check numerically that ``F`` reverses the order on ``[0,1]^n``."""
    rng = np.random.default_rng(seed)
    for _ in range(trials):
        x = rng.random(n)
        y = np.minimum(1.0, x + rng.random(n) * (1 - x))
        if not (F(y) <= F(x) + 1e-12).all():
            return False
    return True
