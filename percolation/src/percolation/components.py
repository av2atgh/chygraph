"""Finite-component size distribution of a random chygraph, by layer.

The giant-component map of :mod:`percolation.giant`,

    Q^{ml}_-  =  prod_k Phi^l_k(Q^{lk}_-)  prod_k [Gbar^l_k]_{k=m}(Q^{lk}_+)
    Q^{ml}_+  =  prod_k [Phibar^l_k]_{k=m}(Q^{lk}_-)  prod_k G^l_k(Q^{lk}_+)

returns one number per directed inclusion type: the probability that the
component reached is finite.  Mark every layer-``l`` complex reached with a
variable ``x_l`` and the same map becomes an equation between *generating
functions* of the finite component,

    H^{ml}_i(x)  =  x_l  F_{(i,m,l)}(H),          P^l(x) = x_l * root_l(H),

whose coefficient of ``prod_l x_l^{n_l}`` in ``P^l`` is the probability that a
random layer-``l`` complex lies in a finite component holding exactly ``n_k``
layer-``k`` complexes, for every ``k`` at once.  At ``x = 1`` this is the map
of Chapter 4 and ``P^l(1)`` is ``1 - S_l``; the coefficients are what Sec. 21.3
of the book calls the full distribution, of which Chapter 5 extracted the
leading singular term.

Two routes to the coefficients are implemented, and they are each other's
check.

``ComponentDistribution``
    Numerical.  The marked map is iterated on a circle ``|x_l| = r < 1`` in the
    complex plane, where iteration from zero converges to the power-series
    fixed point, and the coefficients are read off by the inverse discrete
    Fourier transform (Cauchy's integral on a grid).  Costs one map
    evaluation per grid point per iteration, vectorised; a distribution to
    size 400 on a 2-layer chygraph takes well under a second.

``good_coefficient``, ``symbolic_series``
    Exact, for chygraphs whose generating functions are polynomials.  The
    first is Good's multivariate Lagrange inversion (Proc. Cambridge Philos.
    Soc. 56, 367, 1960): for the system ``w_j = x_j phi_j(w)``,

        [x^n] f(w) = [w^n]  f(w)  prod_j phi_j(w)^{n_j}
                              det( delta_ij - (w_j / phi_i) d phi_i / d w_j ),

    one marking variable per unknown, summed over the unknowns that share a
    layer's marking variable.  The second iterates the marked map in exact
    rational arithmetic, truncated at a total degree.

What the marking assumes.  ``x_l`` multiplies the message of every layer-``l``
complex *reached*, so the model's generating functions must fold occupation
into whether a complex is reached (bond-type thinning inside ``G`` and
``Gbar``, as ``household_epidemic`` and ``graph_with_triangles_giant`` do),
not into whether a reached complex exists (site-type thinning of ``Phi`` by
``thin``).  With site thinning an absent complex would be counted; pass
``p = 1`` or thin the intra-complex functions instead.
"""

import math

import numpy as np
from sympy import Matrix, Poly, Rational, lambdify, symbols, sympify, eye

from percolation.giant import Chygraph


# ---------------------------------------------------------------------------
# Numerical route: Cauchy integral on a grid
# ---------------------------------------------------------------------------

class ComponentDistribution:
    """Finite-component size distribution of a :class:`percolation.giant.Chygraph`.

    Args:
        model: the chygraph map.
        subs: numerical values for the model's free symbols.
    """

    def __init__(self, model, subs=None):
        if not isinstance(model, Chygraph):
            raise TypeError("model must be a percolation.giant.Chygraph")
        self.model = model
        self.L = model.L
        self.n = model.n
        subs = dict(subs or {})
        self.params = sorted(subs, key=str)
        self.pv = [float(subs[k]) for k in self.params]
        args = list(model.Q) + [sympify(k) for k in self.params]
        self._f = lambdify(args, model.F(), 'numpy')
        self._root = lambdify(args, model.root(model.Q), 'numpy')
        # the layer each unknown belongs to, for the marking
        self._layer = [l for i in range(2) for m in range(self.L)
                       for l in range(self.L)]
        self.iterations = None

    # -- the marked map -----------------------------------------------------

    def _iterate(self, xl, tol=1e-14, maxiter=20000):
        """Fixed point of ``H = x_l F(H)`` at marking values ``xl[l]`` (arrays)."""
        shape = np.broadcast(*xl).shape
        H = [np.zeros(shape, dtype=complex) for _ in range(self.n)]
        for it in range(1, maxiter + 1):
            out = self._f(*(H + self.pv))
            Hn = [np.broadcast_to(np.asarray(v, dtype=complex), shape)
                  * xl[self._layer[j]] for j, v in enumerate(out)]
            err = max(float(np.max(np.abs(a - b))) for a, b in zip(Hn, H))
            H = Hn
            if err < tol:
                self.iterations = it
                return H
        self.iterations = maxiter
        return H

    def marked_root(self, xl, **kw):
        """``P^l(x)`` for every layer, at marking values ``xl``."""
        H = self._iterate(list(xl), **kw)
        shape = np.broadcast(*xl).shape
        vals = self._root(*(H + self.pv))
        out = []
        occ = self.model.root_occupation
        for l, v in enumerate(vals):
            v = np.broadcast_to(np.asarray(v, dtype=complex), shape)
            if occ is None:
                out.append(xl[l] * v)
            else:
                pi = float(occ[l])
                # root() returned 1 - pi + pi * val; an absent root has size 0
                val = (v - (1 - pi)) / pi if pi else v * 0
                out.append(1 - pi + pi * xl[l] * val)
        return out

    # -- distributions ------------------------------------------------------

    @staticmethod
    def _grid(nmax, radius, points):
        N = points or 1 << (2 * nmax + 1).bit_length()
        r = radius if radius is not None else 10 ** (-4.0 / max(nmax, 1))
        s = np.arange(N)
        x = r * np.exp(2j * np.pi * s / N)
        return N, r, x

    def distribution(self, nmax, weights=None, layer=0, radius=None,
                     points=None, **kw):
        """``P(size = s)`` for ``s = 0..nmax``, the size being
        ``sum_l weights[l] * (number of layer-l complexes)`` in the finite
        component of a random layer-``layer`` complex.

        The default weight counts layer 0 only (atoms).  ``sum`` of the result
        over all ``s`` is ``1 - S_layer`` up to the truncation.
        """
        w = [1] + [0] * (self.L - 1) if weights is None else list(weights)
        N, r, x = self._grid(nmax, radius, points)
        xl = [x ** int(w[l]) if w[l] else np.ones_like(x) for l in range(self.L)]
        vals = self.marked_root(xl, **kw)[layer]
        # vals_k = sum_j c_j r^j e^{2 pi i k j / N}, so the forward transform
        # divided by N returns c_s r^s (aliased by r^N, which the radius and
        # the grid size keep negligible)
        coef = np.fft.fft(vals).real / N / r ** np.arange(N)
        return coef[:nmax + 1]

    def joint(self, nmax, layers=(0, 1), layer=0, radius=None, points=None,
              **kw):
        """Joint distribution of the numbers of complexes of two layers.

        Returns an array ``J[a, b] = P(n_{layers[0]} = a, n_{layers[1]} = b)``
        for ``a, b <= nmax`` (``nmax`` may be a pair).
        """
        na, nb = (nmax, nmax) if isinstance(nmax, int) else nmax
        Na, ra, xa = self._grid(na, radius, points)
        Nb, rb, xb = self._grid(nb, radius, points)
        XA, XB = np.meshgrid(xa, xb, indexing='ij')
        xl = [np.ones_like(XA) for _ in range(self.L)]
        xl[layers[0]] = XA
        xl[layers[1]] = XB
        vals = self.marked_root(xl, **kw)[layer]
        coef = np.fft.fft2(vals).real / (Na * Nb)
        coef /= np.outer(ra ** np.arange(Na), rb ** np.arange(Nb))
        return coef[:na + 1, :nb + 1]

    def finite_fraction(self, layer=0, **kw):
        """``1 - S_layer`` at ``x = 1``: the mass the distribution must carry."""
        one = [np.ones(1, dtype=complex) for _ in range(self.L)]
        return float(self.marked_root(one, **kw)[layer].real[0])


# ---------------------------------------------------------------------------
# Exact routes, for polynomial generating functions
# ---------------------------------------------------------------------------

def _marked_symbols(model):
    xs = symbols(f'x0:{model.L}')
    return list(xs)


def symbolic_series(model, degree, subs=None, layer=0):
    """Exact coefficients of ``P^layer(x)`` up to total degree ``degree``.

    Iterates the marked map in exact arithmetic, dropping every monomial of
    total degree above ``degree``; after ``degree`` iterations the kept
    coefficients are exact when every layer is marked (each iteration fixes
    one more degree).  Returns ``{(n_0, ..., n_{L-1}): probability}``.
    """
    subs = {sympify(k): Rational(v) if not isinstance(v, str) else sympify(v)
            for k, v in (subs or {}).items()}
    xs = _marked_symbols(model)
    F = [sympify(e).subs(subs) for e in model.F()]
    layer_of = [l for i in range(2) for m in range(model.L) for l in range(model.L)]

    def truncate(expr):
        p = Poly(sympify(expr).expand(), *xs)
        return sum(c * math.prod(x ** e for x, e in zip(xs, mon))
                   for mon, c in p.terms() if sum(mon) <= degree)

    H = [sympify(0)] * model.n
    for _ in range(degree + 1):
        H = [truncate(xs[layer_of[j]] * F[j].subs(dict(zip(model.Q, H))))
             for j in range(model.n)]
    root = sympify(model.root(model.Q)[layer]).subs(subs)
    P = truncate(xs[layer] * root.subs(dict(zip(model.Q, H))))
    return {mon: c for mon, c in Poly(P, *xs).terms()}


def good_coefficient(model, n_by_unknown, subs=None, layer=0):
    """Good's inversion: the coefficient of ``prod_j x_j^{n_j}`` in
    ``P^layer``, with one marking variable per unknown ``j`` (flat index of
    :meth:`Chygraph.index`) and the root's own marking counted separately.

    ``n_by_unknown`` is the multi-index over the ``2 L^2`` unknowns.  The
    root's factor ``x_layer`` is not part of the system, so the returned
    number is the probability that the finite component of a random
    layer-``layer`` complex holds, *beyond the root itself*, exactly ``n_j``
    complexes reached through inclusions of type ``j``.
    """
    subs = {sympify(k): Rational(v) if not isinstance(v, str) else sympify(v)
            for k, v in (subs or {}).items()}
    w = list(model.Q)
    phi = [sympify(e).subs(subs) for e in model.F()]
    f = sympify(model.root(model.Q)[layer]).subs(subs)
    n = model.n
    M = eye(n)
    for i in range(n):
        for j in range(n):
            M[i, j] = M[i, j] - w[j] * phi[i].diff(w[j]) / phi[i]
    expr = f * math.prod(phi[j] ** int(n_by_unknown[j]) for j in range(n)) * M.det()
    # the determinant is a rational function of the unknowns (1/phi_i), so
    # expand it as a series in a scaling variable t before reading the
    # coefficient: w_j -> t w_j, keep total degree <= sum(n)
    t = symbols('t_good')
    total = int(sum(n_by_unknown))
    expr = sympify(expr).subs({wj: t * wj for wj in w})
    expr = expr.series(t, 0, total + 1).removeO().subs(t, 1).expand()
    p = Poly(expr, *w)
    return p.coeff_monomial(math.prod(wj ** int(k) for wj, k in zip(w, n_by_unknown)))


def good_by_layer(model, n_by_layer, subs=None, layer=0):
    """Sum :func:`good_coefficient` over the multi-indices that put ``n_l``
    complexes in layer ``l`` (the root excluded), giving the coefficient of
    ``prod_l x_l^{n_l}`` in ``P^layer / x_layer``.  For small ``n`` only."""
    layer_of = [l for i in range(2) for m in range(model.L) for l in range(model.L)]
    unknowns_of = {l: [j for j in range(model.n) if layer_of[j] == l]
                   for l in range(model.L)}

    def compositions(total, parts):
        if parts == 1:
            yield (total,)
            return
        for first in range(total + 1):
            for rest in compositions(total - first, parts - 1):
                yield (first,) + rest

    total = sympify(0)
    per_layer = [list(compositions(int(n_by_layer[l]), len(unknowns_of[l])))
                 for l in range(model.L)]

    def product(lists, acc=()):
        if not lists:
            yield acc
            return
        for item in lists[0]:
            yield from product(lists[1:], acc + (item,))

    for choice in product(per_layer):
        nvec = [0] * model.n
        for l in range(model.L):
            for j, k in zip(unknowns_of[l], choice[l]):
                nvec[j] = k
        total += good_coefficient(model, nvec, subs, layer)
    return total


def borel(mean, nmax):
    """Component-size distribution of a random node of an Erdos-Renyi graph,
    ``P(s) = e^{-cs} (cs)^{s-1} / s!``, the Borel distribution; the anchor."""
    out = np.zeros(nmax + 1)
    for s in range(1, nmax + 1):
        out[s] = math.exp(-mean * s + (s - 1) * math.log(mean * s)
                          - math.lgamma(s + 1))
    return out


__all__ = ["ComponentDistribution", "symbolic_series", "good_coefficient",
           "good_by_layer", "borel"]
