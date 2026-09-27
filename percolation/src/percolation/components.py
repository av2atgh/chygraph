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


# ---------------------------------------------------------------------------
# Sec. 5.5's joint laws measured on a graph's clique chygraph (Sec. 21.4)
# ---------------------------------------------------------------------------

from collections import Counter  # noqa: E402
from sympy import Rational  # noqa: E402


def _clique_classes(cliques, bins):
    """Class index of every clique from ``bins``, a list of (lo, hi)."""
    def cls(c):
        for i, (lo, hi) in enumerate(bins):
            if lo <= c <= hi:
                return i
        raise ValueError(c)
    return [cls(len(c)) for c in cliques]


def default_bins(g, cliques=None):
    """One class per cardinality up to four and one for the rest, unless the
    cliques run large, when they are binned by size."""
    import networkx as nx
    top = max(len(c) for c in (cliques or nx.find_cliques(g)))
    if top <= 8:
        return [(2, 2), (3, 3), (4, 4), (5, top)] if top >= 5 else [(c, c) for c in range(2, top + 1)]
    return [(2, 2), (3, 3), (4, 4), (5, 7), (8, 14), (15, 24), (25, 49), (50, 99), (100, top)][: (3 + (top >= 5) + (top >= 8) + (top >= 15) + (top >= 25) + (top >= 50) + (top >= 100))]


def joint_clique_model(g, bins=None, thinned=False, types=None, complexes=None):
    """Sec. 5.5's joint construction on the clique chygraph of ``g``.

    Atoms are split into layers by chy-degree (``types``, a list of
    ``(lo, hi)`` bins; the default is one layer for chy-degree one and one
    for the rest) and cliques into layers by cardinality (``bins``).  Each
    atom layer carries the joint distribution of the numbers of cliques of
    each cardinality class containing an atom of that type, and each clique
    layer the joint distribution of the numbers of members of each atom
    type, both measured, so that a node's chy-degree and the cardinalities
    of its cliques are correlated as in the data and so are the chy-degrees
    of the members of one clique -- an isolated complex is a clique all of
    whose members are of the chy-degree-one type.  ``types=[(1, 10**9)]``
    keeps the atoms in one layer, the node-level law alone;
    ``thinned=True`` replaces both measured laws by independent draws, the
    ensemble of :func:`clique_model`, which the marked map must reproduce.
    """
    import networkx as nx
    from percolation.joint import JointChygraph
    cliques = ([c for c in nx.find_cliques(g) if len(c) >= 2] if complexes is None
               else [tuple(c) for c in complexes if len(c) >= 2])
    bins = bins or default_bins(g, cliques)
    types = types or [(1, 1), (2, 10 ** 9)]
    A, L = len(types), len(bins)
    cls = _clique_classes(cliques, bins)
    n = g.number_of_nodes()
    kappa = Counter()
    for c in cliques:
        for v in c:
            kappa[v] += 1

    def atype(v):
        for i, (lo, hi) in enumerate(types):
            if lo <= kappa[v] <= hi:
                return i
        raise ValueError(kappa[v])

    vec = {v: [0] * L for v in g}
    for c, k in zip(cliques, cls):
        for v in c:
            vec[v][k] += 1
    # atom layers 0..A-1, clique layers A..A+L-1
    Phi = [None] * (A + L)
    for t in range(A):
        nodes = [v for v in g if atype(v) == t]
        count = Counter(tuple(vec[v]) for v in nodes)
        tot = len(nodes)
        if thinned:
            # global inclusion fractions by class: no correlation kept
            incl = [sum(vec[v][l] for v in g) for l in range(L)]
            w = [Rational(m, sum(incl)) for m in incl]
            kap = Counter(sum(vec[v]) for v in nodes)

            def phi(x, w=w, kap=kap, tot=tot):
                z = sum(w[l] * x[A + l] for l in range(L))
                return sum(Rational(m, tot) * z ** k for k, m in kap.items())
        else:
            def phi(x, count=count, tot=tot):
                out = 0
                for v, m in count.items():
                    term = Rational(m, tot)
                    for l in range(L):
                        term = term * x[A + l] ** v[l]
                    out = out + term
                return out
        Phi[t] = phi
    G = [None] * (A + L)
    for l in range(L):
        members = [Counter(atype(v) for v in c) for c, k in zip(cliques, cls) if k == l]
        tot = len(members)
        if thinned:
            # global type fractions of inclusions: no correlation kept
            incl = [sum(1 for c in cliques for v in c if atype(v) == t) for t in range(A)]
            w = [Rational(x, sum(incl)) for x in incl]
            card = Counter(sum(m.values()) for m in members)

            def gg(y, w=w, card=card, tot=tot):
                z = sum(w[t] * y[t] for t in range(A))
                return sum(Rational(m, tot) * z ** c for c, m in card.items())
        else:
            count = Counter(tuple(m[t] for t in range(A)) for m in members)

            def gg(y, count=count, tot=tot):
                out = 0
                for v, m in count.items():
                    term = Rational(m, tot)
                    for t in range(A):
                        term = term * y[t] ** v[t]
                    out = out + term
                return out
        G[A + l] = gg
    model = JointChygraph(Phi=Phi, G=G)
    model.atom_layers = A
    model.atom_weights = [sum(1 for v in g if atype(v) == t) / n for t in range(A)]
    return model


def atom_distribution(D, nmax, radius=0.9):
    """``P(s)`` over atoms for a model with several atom layers: mark every
    atom layer and mix the root laws by the node fractions."""
    A = getattr(D.model, 'atom_layers', 1)
    wts = getattr(D.model, 'atom_weights', [1.0])
    mark = [1] * A + [0] * (D.L - A)
    return sum(w * D.distribution(nmax, weights=mark, layer=t, radius=radius)
               for t, w in enumerate(wts))


def atom_finite_fraction(D):
    A = getattr(D.model, 'atom_layers', 1)
    wts = getattr(D.model, 'atom_weights', [1.0])
    return sum(w * D.finite_fraction(layer=t) for t, w in enumerate(wts))


def atom_distribution(D, nmax, radius=0.9):
    """``P(s)`` over atoms for a model with several atom layers: mark every
    atom layer and mix the root laws by the node fractions."""
    A = getattr(D.model, 'atom_layers', 1)
    wts = getattr(D.model, 'atom_weights', [1.0])
    mark = [1] * A + [0] * (D.L - A)
    return sum(w * D.distribution(nmax, weights=mark, layer=t, radius=radius)
               for t, w in enumerate(wts))


def atom_finite_fraction(D):
    A = getattr(D.model, 'atom_layers', 1)
    wts = getattr(D.model, 'atom_weights', [1.0])
    return sum(w * D.finite_fraction(layer=t) for t, w in enumerate(wts))


def interactome_panel(ax, nmax=16):
    import networkx as nx
    from real_chygraphs import load
    out = {}
    for (name, stem), style in zip(INTERACTOMES, (('o', DARK), ('s', MID))):
        g = load(stem)
        n = g.number_of_nodes()
        comps = sorted((len(c) for c in nx.connected_components(g)), reverse=True)
        lcc = comps[0] / n
        real = np.zeros(nmax + 1)
        for s in comps[1:]:
            if s <= nmax:
                real[s] += s / n
        Pg = ComponentDistribution(graph_model(g), {'p': 1, 'q': 1}).distribution(nmax, radius=0.9)
        Sg = 1 - ComponentDistribution(graph_model(g), {'p': 1, 'q': 1}).finite_fraction()
        Dc = ComponentDistribution(clique_model(g))
        Pc = Dc.distribution(nmax, radius=0.9)
        Sc = 1 - Dc.finite_fraction()
        Dj = ComponentDistribution(joint_clique_model(g, complexes=merged_family(g)))
        Pj = atom_distribution(Dj, nmax)
        out[name] = dict(n=n, lcc=lcc, Sg=Sg, Sc=Sc, real=real, Pg=Pg, Pc=Pc,
                         Pj=Pj, Sj=1 - atom_finite_fraction(Dj), small=sum(comps[1:]))
        print(f'  {name}: n = {n}, largest component {lcc:.3f} of the nodes; '
              f'graph model S = {Sg:.3f}, clique chygraph S = {Sc:.3f}; '
              f'{sum(comps[1:])} nodes in {len(comps) - 1} finite components')
        for s in (1, 2, 3, 4, 6, 8):
            print(f'    P({s}): real {real[s]:.4f}  graph {Pg[s]:.4f}  clique {Pc[s]:.4f}')
        ss = np.arange(1, nmax + 1)
        ax.plot(ss, Pg[1:], ':', color=style[1], lw=1.0)
        ax.plot(ss, Pc[1:], '-', color=style[1], lw=1.0)
        ax.plot(ss, Pj[1:], '--', color=style[1], lw=1.0)
        mask = real[1:] > 0
        ax.plot(ss[mask], real[1:][mask], style[0], color=style[1], ms=3.2,
                mfc='white', mew=0.8, label=name)
    ax.plot([], [], ':', color=DARK, lw=1.0, label='degree distribution')
    ax.plot([], [], '-', color=DARK, lw=1.0, label='clique chygraph')
    ax.plot([], [], '--', color=DARK, lw=1.0, label='merged, joint laws')
    ax.set_yscale('log')
    ax.set_xlabel('finite component size $s$', fontsize=8)
    ax.set_ylabel('fraction of nodes', fontsize=8)
    ax.set_ylim(1e-3, 1)
    ax.set_xlim(0.5, 16)
    ax.set_xticks(range(2, 17, 2))
    ax.legend(fontsize=6.5, frameon=False)
    _tidy(ax)
    return out


def table(nmax=30):
    """Table 21.1: real against predicted components, nine networks."""
    import networkx as nx
    from real_chygraphs import load
    rows = []
    for name, stem in TABLE:
        g = load(stem)
        n = g.number_of_nodes()
        comps = sorted((len(c) for c in nx.connected_components(g)), reverse=True)
        real = np.zeros(nmax + 1)
        for s in comps[1:]:
            if s <= nmax:
                real[s] += s / n
        Dg = ComponentDistribution(graph_model(g), {'p': 1, 'q': 1})
        Pg = Dg.distribution(nmax, radius=0.9)
        Dc = ComponentDistribution(clique_model(g))
        Pc = Dc.distribution(nmax, radius=0.9)
        rows.append((name, n, comps[0] / n, 1 - Dg.finite_fraction(),
                     1 - Dc.finite_fraction(), real[2], Pg[2], Pc[2],
                     real[3], Pg[3], Pc[3]))
        print(f'  {name:20s} n={n:5d} largest {comps[0] / n:.3f} | S graph '
              f'{1 - Dg.finite_fraction():.3f} clique {1 - Dc.finite_fraction():.3f} | '
              f'P2 {real[2]:.3f} {Pg[2]:.3f} {Pc[2]:.3f} | P3 {real[3]:.3f} {Pg[3]:.3f} {Pc[3]:.3f}')
    with open(OUT / 'tab-components.tex', 'w') as f:
        f.write('% generated by figs/components.py; do not edit\n')
        f.write('\\begin{tabular}{lrrrrrrrr}\n\\hline\\hline\n')
        f.write('network & $n$ & largest & \\multicolumn{2}{c}{$S$} '
                '& \\multicolumn{2}{c}{$P(2)$} & \\multicolumn{2}{c}{$P(3)$}\\\\\n')
        f.write(' & & real & graph & clique & real & clique & real & clique\\\\\n\\hline\n')
        for r in rows:
            name, n, lcc, Sg, Sc, r2, g2, c2, r3, g3, c3 = r
            vals = ' & '.join(f'{v:.3f}' for v in (lcc, Sg, Sc, r2, c2, r3, c3))
            f.write(f'{name} & {n} & {vals}\\\\\n')
        f.write('\\hline\\hline\n\\end{tabular}\n')
    print('  wrote tab-components.tex')


