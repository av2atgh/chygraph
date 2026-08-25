"""Chygraphs with statistically dependent layers.

:mod:`chygraph.giant` takes the chy-degree generating function of a layer-``l``
complex to factorise over target layers,
``Phi^l(x_0, ..., x_{L-1}) = prod_k Phi^l_k(x_k)``, and likewise for the
intra-complex component sizes.  That covers every mapping in

    Alexei Vazquez, "Percolation in higher order networks via mapping to
    chygraphs", https://doi.org/10.1093/comnet/cnae047

but not constructions in which a complex's participation in different layers is
correlated -- a node whose number of triangles is tied to its number of links,
a hyperedge whose cardinalities in two vertex layers move together.

This module drops the factorisation.  A layer is described by a single
*multivariate* generating function, and the excess generating functions are no
longer independent inputs: they are derivatives of the joint one,

    Phibar^{l,(m)}(x) = (dPhi^l/dx_m)(x) / (dPhi^l/dx_m)(1),

and similarly ``Gbar^{l,(m)}``.  This derivation is the configuration-model
size-biasing relation and is correct whenever a complex is reached by sampling
one of its inclusions uniformly, which is always the case for the chy-degrees.
It is *not* correct when the fine-grained structure inside a complex makes the
entry vertex matter -- bond percolation inside a motif, where ``G^l`` is a
component size rather than a cardinality.  Those excess functions must be
supplied explicitly through ``Gbar``; see the triangle enumeration in
Ref. above.  Two consequences follow.

The map becomes

    Q^{ml}_- = Phi^l(Q^{l.}_-)          Gbar^{l,(m)}(Q^{l.}_+),
    Q^{ml}_+ = Phibar^{l,(m)}(Q^{l.}_-) G^l(Q^{l.}_+),
    P^l      = Phi^l(Q^{l.}_-)          G^l(Q^{l.}_+),

which reduces term by term to the factorised map when the generating functions
factorise.

Its Jacobian generalises the published threshold tensor.  Define the
*inclusion-biased* moment matrices

    <kappabar^(m)>_{lk} = <kappa_lm (kappa_lk - delta_mk)> / <kappa_lm>
                        = (d^2 Phi^l / dx_m dx_k)(1) / <kappa>_{lm},
    <sbar^(m)>_{lk}     = <s_lm (s_lk - delta_mk)> / <s_lm>,

i.e. the expected chy-degree to layer ``k`` of a complex reached through an
inclusion into layer ``m``.  Then

    (A_--)^{ml}_{nk} = delta_mn delta_lk - <kappa>_{nk} delta_nl,
    (A_-+)^{ml}_{nk} = -<sbar^(m)>_{nk} delta_nl,
    (A_+-)^{ml}_{nk} = -<kappabar^(m)>_{nk} delta_nl,
    (A_++)^{ml}_{nk} = delta_mn delta_lk - <s>_{nk} delta_nl,

a single expression per block in place of the ``delta_mk`` / ``(1 - delta_mk)``
split of the published tensor.  Under independence
``<kappa_m kappa_k> = <kappa_m><kappa_k>`` for ``k != m`` and the two forms
coincide, so the published tensor is the independent-layer special case.  When
the layers are correlated the off-diagonal entries must carry the joint second
moment, and both the threshold and the order parameter change.
"""

from sympy import Symbol, diff, simplify, sympify

from chygraph.giant import GiantComponent


class JointGiantComponent(GiantComponent):
    """A chygraph whose layers may be statistically dependent.

    Args:
        Phi: length-``L`` list.  ``Phi[l]`` is a callable taking a length-``L``
            sequence and returning the joint generating function of the
            chy-degree vector of a layer-``l`` complex.  ``None`` means the
            layer is included in nothing.
        G: length-``L`` list.  ``G[l]`` is a callable taking a length-``L``
            sequence and returning the joint generating function of the
            intra-complex component size vector of a layer-``l`` complex.
            ``None`` means the layer is an atom.
        Phibar: optional ``L x L`` table overriding the derived excess
            chy-degree generating functions, ``Phibar[l][m]``.
        Gbar: optional ``L x L`` table overriding the derived excess
            intra-complex generating functions, ``Gbar[l][m]``.  Required
            whenever ``G[l]`` generates a component size rather than a
            cardinality, because size-biasing by differentiation is then the
            wrong relation.
        occupation: optional length-``L`` list of occupation probabilities.
            ``occupation[l]`` multiplies the whole of layer ``l``'s contribution,
            ``1 - pi_l + pi_l Phi^l G^l``: an absent complex cannot be left in
            any direction.  It must be applied here and not folded into the
            generating functions, because the excess functions are obtained by
            differentiating them and size-biasing would condition the complex to
            be present.

    The excess generating functions are derived, not supplied.
    """

    def __init__(self, Phi, G, Phibar=None, Gbar=None, occupation=None,
                 root_occupation=None):
        L = len(Phi)
        if len(G) != L:
            raise ValueError("Phi and G must have the same length")
        self.occupation = [1] * L if occupation is None else list(occupation)
        self.root_occupation = root_occupation
        if len(self.occupation) != L:
            raise ValueError("occupation must have one entry per layer")
        self.L = L
        self.n = 2 * L * L
        self.Phi, self.G = Phi, G
        self.Q = [Symbol(f'Q{i}_{m}_{l}')
                  for i in range(2) for m in range(L) for l in range(L)]
        self._F = None
        self._Pexpr = None
        self._J = None

        self._x = [Symbol(f'_x{k}') for k in range(L)]
        self._y = [Symbol(f'_y{k}') for k in range(L)]
        self._PHI, self._PHIBAR = self._build(Phi, self._x, Phibar)
        self._G, self._GBAR = self._build(G, self._y, Gbar)

    # -- generating functions and their derived excess counterparts ---------

    @staticmethod
    def _build(table, vars_, override=None):
        L = len(vars_)
        at_one = {v: 1 for v in vars_}
        base, excess = [], []
        for l in range(L):
            if table[l] is None:
                base.append(None)
                excess.append([None] * L)
                continue
            expr = sympify(table[l](vars_))
            base.append(expr)
            row = []
            for m in range(L):
                if override is not None and override[l][m] is not None:
                    row.append(sympify(override[l][m](vars_)))
                    continue
                d = diff(expr, vars_[m])
                norm = d.subs(at_one)
                # no inclusions into layer m: that state is unreachable and the
                # excess is 0/0.  Fall back to the unbiased function, which is
                # the convention of the published threshold tensor
                row.append(expr if simplify(norm) == 0 else d / norm)
            excess.append(row)
        return base, excess

    def _eval(self, expr, vars_, values):
        if expr is None:
            return sympify(1)
        return expr.subs(dict(zip(vars_, values)), simultaneous=True)

    # -- the map ------------------------------------------------------------

    def apply(self, Qv):
        L = self.L
        out = [None] * self.n
        for i in range(2):
            for m in range(L):
                for l in range(L):
                    xs = [Qv[self.index(0, l, k)] for k in range(L)]
                    ys = [Qv[self.index(1, l, k)] for k in range(L)]
                    ph = self._PHIBAR[l][m] if i == 1 else self._PHI[l]
                    gg = self._GBAR[l][m] if i == 0 else self._G[l]
                    val = (self._eval(ph, self._x, xs)
                           * self._eval(gg, self._y, ys))
                    pi = self.occupation[l]
                    out[self.index(i, m, l)] = 1 - pi + pi * val
        return out

    def root(self, Qv):
        L = self.L
        out = []
        for l in range(L):
            xs = [Qv[self.index(0, l, k)] for k in range(L)]
            ys = [Qv[self.index(1, l, k)] for k in range(L)]
            val = (self._eval(self._PHI[l], self._x, xs)
                   * self._eval(self._G[l], self._y, ys))
            pi = self.occupation[l]
            val = 1 - pi + pi * val
            if self.root_occupation is not None:
                r = self.root_occupation[l]
                val = 1 - r + r * val
            out.append(val)
        return out

    # -- moment matrices ----------------------------------------------------

    def _moment(self, base, vars_, l, k):
        if base[l] is None:
            return sympify(0)
        return diff(base[l], vars_[k]).subs({v: 1 for v in vars_})

    def kappa(self, l, k):
        """``<kappa>_{lk}``."""
        return self._moment(self._PHI, self._x, l, k)

    def s(self, l, k):
        """``<s>_{lk}``."""
        return self._moment(self._G, self._y, l, k)

    def kappa_bar(self, m, l, k):
        """``<kappabar^(m)>_{lk}``: expected chy-degree to layer ``k`` of a
        layer-``l`` complex reached through an inclusion into layer ``m``."""
        e = self._PHIBAR[l][m]
        if e is None:
            return sympify(0)
        return diff(e, self._x[k]).subs({v: 1 for v in self._x})

    def s_bar(self, m, l, k):
        """``<sbar^(m)>_{lk}``, the intra-complex counterpart."""
        e = self._GBAR[l][m]
        if e is None:
            return sympify(0)
        return diff(e, self._y[k]).subs({v: 1 for v in self._y})

    def layers_independent(self):
        """``True`` when every joint second moment factorises, i.e. when the
        published threshold tensor applies unchanged."""
        L = self.L
        for l in range(L):
            for m in range(L):
                for k in range(L):
                    if k == m:
                        continue
                    if self._PHIBAR[l][m] is not None:
                        if simplify(self.kappa_bar(m, l, k) - self.kappa(l, k)) != 0:
                            return False
                    if self._GBAR[l][m] is not None:
                        if simplify(self.s_bar(m, l, k) - self.s(l, k)) != 0:
                            return False
        return True


__all__ = ["JointGiantComponent"]
