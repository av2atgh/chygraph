"""Giant component fraction in chygraphs.

The percolation threshold in ``percolation.percolation`` is obtained from the
tensor ``A`` of Vazquez, Phys. Rev. E 107, 024316 (2023), Eq. (34), whose
entries are the first moments <k>, <K>, <s>, <S>.  That tensor is exactly the
Jacobian, at the trivial fixed point, of the nonlinear self-consistency map of
Eqs. (25)-(28) of the same reference:

    A = I - J,      J_ab = d F_a / d Q_b |_{Q=1}

so ``theta = -det(A)`` and ``Lambda = max eig(-A) = max eig(J) - 1`` are
diagnostics of the *linearisation*.  The order parameter itself is the
non-trivial fixed point of ``F``, which requires the probability generating
functions rather than their first derivatives.

Unknowns.  ``Q[i][m][l]`` is the probability that a layer-``l`` complex, reached
from a layer-``m`` complex, is not connected to the giant component through any
route other than the one used to reach it.  ``i = 0`` ("-") means the complex was
reached from below (from a complex it includes), ``i = 1`` ("+") from above (from
a complex that includes it).  The flat index is ``i*L*L + m*L + l``, matching
``(vec A)_{mL+l, nL+k}`` of the threshold calculation.

Inputs.  Four ``L x L`` tables of *univariate* generating functions:

    phi[l][k](x)     number of layer-k complexes that include a layer-l complex
    phibar[l][k](x)  the excess (size-biased) counterpart
    g[l][k](y)       number of layer-k complexes in the intra-complex component
                     of a layer-l complex
    gbar[l][k](y)    the excess counterpart

``None`` means "identically 1" (no such coupling).  Occupation probabilities are
folded into the generating functions by thinning, e.g. a node present with
probability ``p`` gives ``phi(x) -> 1 - p + p*phi(x)``; this reproduces the
``<k>_{0l} = p<k>`` convention of the threshold calculation and makes
``1 - P^0`` the fraction of *all* nodes in the giant component.

Reference:
    Alexei Vazquez, "Percolation in higher order networks via mapping to
    chygraphs", https://doi.org/10.1093/comnet/cnae047
"""

from sympy import Symbol, Matrix, eye, exp, lambdify, sympify

from percolation.percolation import PercolationMatrix


# ---------------------------------------------------------------------------
# Generating function helpers
# ---------------------------------------------------------------------------

def poisson_pgf(mean):
    """PGF of a Poisson variable; equal to its own excess PGF."""
    return lambda x: exp(mean * (x - 1))


def thin(pgf, p):
    """Thin a PGF: the whole object is absent with probability ``1 - p``."""
    if pgf is None:
        return None
    return lambda x: 1 - p + p * pgf(x)


def finite_pgf(probs):
    """PGF of a finite distribution given as ``{value: probability}``."""
    return lambda x: sum(pr * x**int(v) for v, pr in probs.items())


def moment_pgf(mean, second):
    """Quadratic surrogate ``1 + <x>(t-1) + <x(x-1)>(t-1)^2/2``.

    Only the first two derivatives at ``t = 1`` enter the critical amplitude
    (:mod:`percolation.amplitude`), so this surrogate makes the amplitude of an
    *arbitrary* distribution exactly computable from its first two factorial
    moments.  It is not a probability generating function and must not be used
    for the fixed point itself.
    """
    return lambda x: 1 + mean * (x - 1) + second * (x - 1)**2 / 2


def constant_pgf(c):
    """PGF of a variable equal to ``c`` with probability one."""
    return lambda x: x**int(c)


# ---------------------------------------------------------------------------
# Core map
# ---------------------------------------------------------------------------

class Chygraph:
    """Non-linear self-consistency map whose fixed point gives the giant component.

    Args:
        phi, phibar, g, gbar: ``L x L`` nested lists of callables (or ``None``).
        root_occupation: optional length-``L`` list.  ``root_occupation[l]``
            multiplies layer ``l``'s contribution to ``P^l`` only, not to the
            map.  It is needed when a complex's existence is already implied by
            the connectivity rule and so must be reinstated at the root: under
            AND-logic percolation a node reached through an active hyperedge is
            present by construction, but a *randomly chosen* node need not be.
            With it, ``S_l`` is the fraction of all layer-``l`` complexes.
    """

    def __init__(self, phi, phibar, g, gbar, root_occupation=None,
                 perron_probe=None):
        self.L = L = len(phi)
        self.root_occupation = root_occupation
        self.perron_probe = perron_probe
        self.phi, self.phibar, self.g, self.gbar = phi, phibar, g, gbar
        self.n = 2 * L * L
        self.Q = [Symbol(f'Q{i}_{m}_{l}')
                  for i in range(2) for m in range(L) for l in range(L)]
        self._F = None
        self._Pexpr = None
        self._J = None

    # -- index bookkeeping --------------------------------------------------

    def index(self, i, m, l):
        """Flat index of ``Q[i][m][l]``, matching the vec2 layout of ``A``."""
        return i * self.L * self.L + m * self.L + l

    # -- the map ------------------------------------------------------------

    def _pgf(self, table, fallback, l, k):
        f = table[l][k]
        return f if f is not None else fallback[l][k]

    def apply(self, Qv):
        """Evaluate the map ``F`` at the vector ``Qv``."""
        L = self.L
        out = [None] * self.n
        for i in range(2):
            for m in range(L):
                for l in range(L):
                    val = sympify(1)
                    for k in range(L):
                        if self.phi[l][k] is None:
                            continue
                        f = (self._pgf(self.phibar, self.phi, l, k)
                             if (i == 1 and k == m) else self.phi[l][k])
                        val = val * f(Qv[self.index(0, l, k)])
                    for k in range(L):
                        if self.g[l][k] is None:
                            continue
                        f = (self._pgf(self.gbar, self.g, l, k)
                             if (i == 0 and k == m) else self.g[l][k])
                        val = val * f(Qv[self.index(1, l, k)])
                    out[self.index(i, m, l)] = val
        return out

    def root(self, Qv):
        """``P^l``: probability a randomly chosen layer-l complex is not in the
        giant component.  ``1 - P^0`` is the fraction of nodes."""
        L = self.L
        out = []
        for l in range(L):
            val = sympify(1)
            for k in range(L):
                if self.phi[l][k] is not None:
                    val = val * self.phi[l][k](Qv[self.index(0, l, k)])
            for k in range(L):
                if self.g[l][k] is not None:
                    val = val * self.g[l][k](Qv[self.index(1, l, k)])
            if self.root_occupation is not None:
                pi = self.root_occupation[l]
                val = 1 - pi + pi * val
            out.append(val)
        return out

    # -- symbolic objects ---------------------------------------------------

    def F(self):
        if self._F is None:
            self._F = self.apply(self.Q)
        return self._F

    def P(self):
        if self._Pexpr is None:
            self._Pexpr = self.root(self.Q)
        return self._Pexpr

    def jacobian(self):
        """``J = dF/dQ`` at ``Q = 1``."""
        if self._J is None:
            J = Matrix(self.F()).jacobian(self.Q)
            self._J = J.subs({q: 1 for q in self.Q})
        return self._J

    def A(self):
        """``A = I - J``: the threshold tensor of the published calculation."""
        return eye(self.n) - self.jacobian()

    def theta(self):
        return -self.A().det()

    # -- numerics -----------------------------------------------------------

    def _compile(self, params):
        key = tuple(str(s) for s in params)
        cache = getattr(self, '_compiled', None)
        if cache is None:
            cache = self._compiled = {}
        if key not in cache:
            args = list(self.Q) + list(params)
            cache[key] = (lambdify(args, self.F(), 'math'),
                          lambdify(args, self.P(), 'math'))
        return cache[key]

    def solve(self, subs=None, tol=1e-13, maxiter=100000):
        """Smallest fixed point, reached by monotone iteration from ``Q = 0``.

        The map has non-negative coefficients, so iterating upward from 0
        converges monotonically to the smallest (physical) fixed point.
        """
        subs = dict(subs or {})
        params = sorted(subs, key=str)
        f, _ = self._compile(params)
        pv = [float(subs[k]) for k in params]
        Q = [0.0] * self.n
        for _ in range(maxiter):
            Qn = f(*(Q + pv))
            if max(abs(a - b) for a, b in zip(Qn, Q)) < tol:
                return Qn
            Q = Qn
        return Q

    def fractions(self, subs=None, **kw):
        """``S_l = 1 - P^l`` for every layer."""
        subs = dict(subs or {})
        params = sorted(subs, key=str)
        _, p = self._compile(params)
        pv = [float(subs[k]) for k in params]
        Q = self.solve(subs, **kw)
        return [1.0 - v for v in p(*(Q + pv))]

    def node_fraction(self, subs=None, **kw):
        """Fraction of layer-0 complexes (atoms/nodes) in the giant component."""
        return self.fractions(subs, **kw)[0]

    # -- critical behaviour -------------------------------------------------

    # -- the linearisation and its next order ------------------------------
    #
    # Implemented in percolation.amplitude and exposed here so that one object
    # answers all three questions about a chygraph: where the transition is,
    # how steeply the order parameter rises out of it, and what it is away from
    # the threshold.  The import is deferred to avoid a circular one.

    @property
    def _amplitude(self):
        if getattr(self, '_amp_cache', None) is None:
            from percolation.amplitude import CriticalAmplitude
            self._amp_cache = CriticalAmplitude(
                self, probe=getattr(self, 'perron_probe', None))
        return self._amp_cache

    def perron_root(self):
        """``lambda = 1 + Lambda``, the Perron root of the Jacobian."""
        return self._amplitude.perron_root()

    def Lambda(self):
        """``max eig(-vec2 A)``, the order parameter of the threshold."""
        return self._amplitude.Lambda()

    def core(self):
        """Indices of the coupled core of the map."""
        return self._amplitude.core

    def amplitude(self, layer=0):
        """``B`` in ``S_layer = B Lambda + O(Lambda^2)``, in closed form."""
        return self._amplitude.amplitude(layer)

    def amplitude_at_threshold(self, layer=0, solve_for=None):
        """``B`` reduced on the critical manifold ``Lambda = 0``."""
        return self._amplitude.amplitude_at_threshold(layer, solve_for)

    def curvature(self):
        """``C = l . M[r, r]``; ``B`` holds when it is finite and positive."""
        return self._amplitude.curvature()

    def is_continuous(self):
        """``False`` when the curvature vanishes identically."""
        return self._amplitude.is_continuous()

    def verify(self, subs=None, tol=1e-9):
        """Check the core reduction against the full index set."""
        return self._amplitude.verify(subs, tol)

    def amplitude_numeric(self, subs, layer=0):
        """``B`` at a point in parameter space, by numerical linear algebra.

        An independent evaluation of :meth:`amplitude`, using numpy eigenvectors
        on the full unreduced index set, kept as a cross-check of the symbolic
        route, which works on the reduced core.

        Expanding ``Q = 1 - x`` gives ``x = J x - (1/2) M[x, x] + ...``.  With
        ``r`` and ``l`` the right and left Perron vectors of ``J`` normalised by
        ``l.r = 1``, the leading solution is ``x = e r`` with
        ``e = 2 Lambda / (l . M[r, r])``, so

            S_layer = (grad P^layer . r) * e = B * Lambda,
            B = 2 (grad P^layer . r) / (l . M[r, r]).

        Only *second* derivatives of the generating functions enter, i.e. one
        moment order beyond what the threshold needs.  ``B`` diverges
        (``l . M[r, r] -> 0``) exactly when the transition stops being
        continuous.
        """
        import numpy as np

        subs = {sympify(k): v for k, v in subs.items()}
        J = np.array(self.jacobian().subs(subs).evalf(), dtype=float)

        w, V = np.linalg.eig(J)
        i = int(np.argmax(w.real))
        r = np.abs(V[:, i].real)
        wl, Vl = np.linalg.eig(J.T)
        j = int(np.argmax(wl.real))
        lv = np.abs(Vl[:, j].real)
        lv = lv / float(lv @ r)

        # l . M[r, r] without forming the full Hessian tensor
        Fs = self.F()
        one = {q: 1 for q in self.Q}
        tot = 0.0
        for a in range(self.n):
            if abs(lv[a]) < 1e-14:
                continue
            expr = Fs[a]
            acc = 0.0
            for b in range(self.n):
                if abs(r[b]) < 1e-14:
                    continue
                d1 = expr.diff(self.Q[b])
                for c in range(self.n):
                    if abs(r[c]) < 1e-14:
                        continue
                    d2 = float(d1.diff(self.Q[c]).subs(one).subs(subs).evalf())
                    acc += d2 * r[b] * r[c]
            tot += lv[a] * acc

        Pexpr = self.P()[layer]
        grad = np.array([float(Pexpr.diff(q).subs(one).subs(subs).evalf())
                         for q in self.Q])
        return float(2.0 * (grad @ r) / tot)


# ---------------------------------------------------------------------------
# Models, mirroring those in percolation.percolation
# ---------------------------------------------------------------------------

def _tables(L):
    return ([[None] * L for _ in range(L)], [[None] * L for _ in range(L)],
            [[None] * L for _ in range(L)], [[None] * L for _ in range(L)])


def hypergraph_giant(degree=None, excess_degree=None,
                     cardinality=None, excess_cardinality=None,
                     p=None, q=None, poisson=True, graph=False):
    """Site-bond percolation on a hypergraph (layer 0 nodes, layer 1 hyperedges).

    With ``poisson=True`` the degree and cardinality PGFs are built from the
    means ``k`` and ``c``, so the input is exactly the ``<k>``/``<c>`` of the
    threshold calculation.
    """
    from sympy import symbols
    k, c = symbols('k c')
    p = symbols('p') if p is None else p
    q = symbols('q') if q is None else q

    if poisson:
        degree = degree or poisson_pgf(k)
        excess_degree = excess_degree or degree
    if graph:
        cardinality = cardinality or constant_pgf(2)
        excess_cardinality = excess_cardinality or constant_pgf(1)
    elif poisson:
        cardinality = cardinality or poisson_pgf(c)
        excess_cardinality = excess_cardinality or cardinality

    phi, phibar, g, gbar = _tables(2)
    phi[0][1] = thin(degree, p)
    phibar[0][1] = thin(excess_degree, p)
    g[1][0] = thin(cardinality, q)
    gbar[1][0] = thin(excess_cardinality, q)
    return Chygraph(phi, phibar, g, gbar)


def multiplex_hypergraph_giant(number_of_types=2, poisson=True, graph=False):
    """Multiplex hypergraph: layer 0 nodes, layers 1..L hyperedge types."""
    from sympy import symbols
    L = number_of_types + 1
    phi, phibar, g, gbar = _tables(L)
    for l in range(1, L):
        kl, cl, ql = symbols(f'k_0{l} c_{l}0 q_{l}')
        deg = poisson_pgf(kl) if poisson else None
        if deg is None:
            raise NotImplementedError("supply explicit PGFs for non-Poisson")
        phi[0][l] = deg
        phibar[0][l] = deg
        card = constant_pgf(2) if graph else poisson_pgf(cl)
        excess = constant_pgf(1) if graph else card
        g[l][0] = thin(card, ql)
        gbar[l][0] = thin(excess, ql)
    return Chygraph(phi, phibar, g, gbar)


def graph_with_triangles_giant(poisson=True):
    """Bond percolation on a graph with links (layer 1) and triangles (layer 2).

    ``gbar[2][0]`` is the full excess-component distribution of Fig. 3 of the
    chygraph percolation paper, of which ``<S>_20 = 2q(1 + q - q^2)`` is the
    mean.
    """
    from sympy import symbols
    q, kL, kT = symbols('q k_L k_T')
    phi, phibar, g, gbar = _tables(3)

    deg_L = poisson_pgf(kL)
    deg_T = poisson_pgf(kT)
    if not poisson:
        raise NotImplementedError("supply explicit PGFs for non-Poisson")
    phi[0][1], phibar[0][1] = deg_L, deg_L
    phi[0][2], phibar[0][2] = deg_T, deg_T

    # a link: the far endpoint is reachable with probability q
    g[1][0] = lambda y: (1 - q) * y + q * y**2
    gbar[1][0] = lambda y: (1 - q) + q * y

    # a triangle, entered at one corner: excess component size 0, 1 or 2
    tri = {2: 3 * q**2 - 2 * q**3,
           1: 2 * q * (1 - q)**2,
           0: (1 - q)**3 + q * (1 - q)**2}
    gbar[2][0] = finite_pgf(tri)
    g[2][0] = lambda y: y * finite_pgf(tri)(y)
    return Chygraph(phi, phibar, g, gbar)


# The class was called GiantComponent before it also carried the threshold and
# the amplitude; it is now simply a chygraph.
GiantComponent = Chygraph


__all__ = [
    "Chygraph",
    "GiantComponent",
    "poisson_pgf", "thin", "finite_pgf", "constant_pgf", "moment_pgf",
    "hypergraph_giant", "multiplex_hypergraph_giant", "graph_with_triangles_giant",
]
