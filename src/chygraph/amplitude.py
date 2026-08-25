"""Symbolic critical amplitude of the chygraph percolation transition.

The threshold calculation of Vazquez, Phys. Rev. E 107, 024316 (2023) returns
``theta = -det(A)`` and ``Lambda = max eig(-A)`` from the tensor ``A``, whose
entries are the first moments ``<k>``, ``<K>``, ``<s>``, ``<S>``.  ``A = I - J``
is the Jacobian of the non-linear self-consistency map ``F`` at its trivial
fixed point (see :mod:`chygraph.giant`).  Carrying the expansion of ``F`` one
order further gives the *amplitude* of the order parameter, still in closed
form.

Writing ``Q = 1 - x`` and expanding,

    x = J x - (1/2) M[x, x] + O(x^3),      M_abc = d^2 F_a / dQ_b dQ_c |_{Q=1}

Let ``lambda = 1 + Lambda`` be the Perron root of ``J`` and let ``r``, ``l`` be
its right and left Perron vectors normalised by ``l . r = 1``.  Setting
``x = e r`` and projecting on ``l`` kills the linear term and leaves
``e = 2 Lambda / (l . M[r, r])``, so for every layer

    S_layer = B * Lambda + O(Lambda^2),
    B = 2 (grad P^layer . r) / (l . M[r, r]).

Only *second* derivatives of the same generating functions whose *first*
derivatives build ``A`` enter ``B``.  This gives a three-tier hierarchy:

    threshold          <- first moments
    critical amplitude <- second moments
    full S(parameters) <- the generating functions themselves

The denominator ``C = l . M[r, r]`` is the curvature of the map along the
critical direction.  It is a sum of non-negative terms, and ``B`` describes the
transition only when ``C`` is *finite and positive*.

Finiteness is a condition on second factorial moments of the excess
distributions: for a graph ``C`` is proportional to
``<kbar(kbar-1)> = <k(k-1)(k-2)>/<k>``, which diverges for a scale-free degree
distribution with exponent ``gamma <= 4``.  There ``B -> 0``, the linear branch
has vanishing amplitude and ``beta != 1`` -- the familiar ``1/(gamma-3)`` of the
heterogeneous mean field for ``3 < gamma < 4``.  The formula signals its own
breakdown, but must not be applied outside the finite-moment range.

``C = 0`` identically means the map is affine along the critical direction and
the linearisation cannot locate the transition at all.

Reference:
    Alexei Vazquez, "Percolation in higher order networks via mapping to
    chygraphs", https://doi.org/10.1093/comnet/cnae047
"""

from sympy import Eq, Matrix, Symbol, eye, nsimplify, simplify, solve, sympify, zeros

from chygraph.giant import GiantComponent


def _default_probe(free):
    """A generic numeric point, safe for both probabilities and means."""
    vals = {}
    for i, sym in enumerate(sorted(free, key=str)):
        vals[sym] = 0.23 + 0.61 * ((7 * i + 3) % 11) / 11.0
    return vals


class CriticalAmplitude:
    """Closed-form critical amplitude ``B`` in ``S ~ B * Lambda``.

    Args:
        model: a :class:`chygraph.giant.GiantComponent`.
        probe: optional ``{symbol: value}`` used only to pick the Perron branch
            among the symbolic eigenvalues.  Any generic point works.

    The heavy lifting is done on the *core* of the map: iteratively discarding
    indices whose row or column of ``J`` vanishes identically leaves the
    physically coupled unknowns, and the non-zero spectrum is unchanged.  For a
    hypergraph this is 8 -> 2 unknowns, for a graph with triangles 18 -> 4.
    """

    def __init__(self, model, probe=None):
        if not isinstance(model, GiantComponent):
            raise TypeError("model must be a GiantComponent")
        self.model = model
        self.n = model.n
        self._probe = probe
        self._cache = {}

    # -- reduction ----------------------------------------------------------

    @property
    def J(self):
        if 'J' not in self._cache:
            self._cache['J'] = self.model.jacobian()
        return self._cache['J']

    @property
    def core(self):
        """Indices surviving iterated removal of zero rows and columns."""
        if 'core' not in self._cache:
            J, keep = self.J, list(range(self.n))
            while True:
                rows = {a for a in keep if any(J[a, b] != 0 for b in keep)}
                cols = {b for b in keep if any(J[a, b] != 0 for a in keep)}
                new = sorted(rows & cols)
                if new == keep:
                    break
                keep = new
            self._cache['core'] = keep
        return self._cache['core']

    def labels(self):
        """Human-readable ``(sign, from_layer, at_layer)`` labels for the core."""
        L = self.model.L
        out = []
        for a in self.core:
            i, rem = divmod(a, L * L)
            m, l = divmod(rem, L)
            out.append(('-' if i == 0 else '+', m, l))
        return out

    @property
    def Jcore(self):
        if 'Jcore' not in self._cache:
            # rationalise: float coefficients make the symbolic nullspace of
            # (lambda I - J) numerically rank-deficient and the Perron vectors
            # are then not found
            Jc = self.J[self.core, self.core]
            self._cache['Jcore'] = Jc.applyfunc(
                lambda e: nsimplify(e, rational=True))
        return self._cache['Jcore']

    def probe(self):
        if self._probe is None:
            free = set()
            for e in self.Jcore:
                free |= sympify(e).free_symbols
            self._probe = _default_probe(free)
        return self._probe

    # -- Perron pair --------------------------------------------------------

    def perron_root(self):
        """``lambda = 1 + Lambda``, the Perron root of ``J``, in closed form."""
        if 'lam' not in self._cache:
            pr = self.probe()
            best, bestval = None, None
            for ev in self.Jcore.eigenvals():
                try:
                    v = complex(ev.subs(pr).evalf())
                except (TypeError, ValueError):
                    continue
                if abs(v.imag) > 1e-9:
                    continue
                if bestval is None or v.real > bestval:
                    best, bestval = ev, v.real
            if best is None:
                raise ValueError("no real eigenvalue found at the probe point")
            self._cache['lam'] = simplify(best)
        return self._cache['lam']

    def Lambda(self):
        """The order parameter of the threshold calculation, ``max eig(-A)``."""
        return self.perron_root() - 1

    def _perron_vectors(self):
        if 'rl' not in self._cache:
            lam = self.perron_root()
            k = self.core
            nc = len(k)
            Mr = lam * eye(nc) - self.Jcore
            ns_r = Mr.nullspace()
            ns_l = Mr.T.nullspace()
            if not ns_r or not ns_l:
                raise ValueError("Perron vectors not found; check the probe")
            rc, lc = ns_r[0], ns_l[0]
            norm = (lc.T * rc)[0]
            lc = lc / norm

            # extend r to the full index set: rows outside the core are
            # determined by back-substitution, r_a = (J r)_a / lambda
            r = zeros(self.n, 1)
            for j, a in enumerate(k):
                r[a] = rc[j]
            outside = [a for a in range(self.n) if a not in set(k)]
            for _ in range(len(outside) + 1):
                for a in outside:
                    r[a] = simplify(sum(self.J[a, b] * r[b]
                                        for b in range(self.n)) / lam)
            l = zeros(self.n, 1)
            for j, a in enumerate(k):
                l[a] = lc[j]
            self._cache['rl'] = (r, l)
        return self._cache['rl']

    def right_vector(self):
        return self._perron_vectors()[0]

    def left_vector(self):
        return self._perron_vectors()[1]

    # -- second order -------------------------------------------------------

    def curvature(self):
        """``C = l . M[r, r]``, the curvature of ``F`` along the critical direction.

        Computed as a directional derivative, ``M_a[r,r] = d^2/dt^2 F_a(1 + t r)``
        at ``t = 0``, which avoids ever forming the rank-3 tensor ``M``.
        """
        if 'C' not in self._cache:
            r, l = self._perron_vectors()
            t = Symbol('t')
            Ft = self.model.apply([1 + t * r[b] for b in range(self.n)])
            C = 0
            for a in range(self.n):
                if l[a] == 0:
                    continue
                C += l[a] * Ft[a].diff(t, 2).subs(t, 0)
            self._cache['C'] = simplify(C)
        return self._cache['C']

    def gradient(self, layer=0):
        """``grad P^layer . r``, again as a directional derivative."""
        r = self.right_vector()
        t = Symbol('t')
        Pt = self.model.root([1 + t * r[b] for b in range(self.n)])
        return simplify(Pt[layer].diff(t).subs(t, 0))

    def amplitude(self, layer=0):
        """``B`` in ``S_layer = B * Lambda + O(Lambda^2)``, in closed form."""
        return simplify(2 * self.gradient(layer) / self.curvature())

    def amplitude_at_threshold(self, layer=0, solve_for=None):
        """``B`` reduced on the critical manifold ``Lambda = 0``.

        The unconstrained ``B`` carries the Perron root explicitly; on the
        critical manifold ``lambda = 1`` and the expression collapses.  One
        parameter is eliminated using ``Lambda = 0``; pass ``solve_for`` to
        choose which, otherwise the first symbol the root depends on is used.
        """
        B, lam = self.amplitude(layer), self.perron_root()
        if solve_for is None:
            for sym in sorted(lam.free_symbols, key=str):
                sol = solve(Eq(lam, 1), sym)
                if sol:
                    solve_for = sym
                    break
            else:
                raise ValueError("could not eliminate a parameter; pass solve_for")
        sol = solve(Eq(lam, 1), solve_for)
        if not sol:
            raise ValueError(f"Lambda = 0 is not solvable for {solve_for}")
        return simplify(B.subs(solve_for, sol[0]))

    def is_continuous(self):
        """``False`` when the curvature vanishes identically: the expansion has
        no linear branch and the transition is not continuous, so ``Lambda = 0``
        no longer locates it."""
        return simplify(self.curvature()) != 0

    # -- checks -------------------------------------------------------------

    def verify(self, subs=None, tol=1e-9):
        """Numerically confirm the core reduction against the full space.

        Checks that ``J r = lambda r`` and ``l J = lambda l`` hold on all ``n``
        indices, and that ``lambda`` is the largest eigenvalue of the full ``J``.
        """
        import numpy as np
        pr = dict(self.probe())
        pr.update({sympify(k): v for k, v in (subs or {}).items()})

        J = np.array(self.J.subs(pr).evalf(), dtype=float)
        lam = float(self.perron_root().subs(pr).evalf())
        r = np.array(self.right_vector().subs(pr).evalf(), dtype=float).ravel()
        l = np.array(self.left_vector().subs(pr).evalf(), dtype=float).ravel()

        return {
            'right': float(np.max(np.abs(J @ r - lam * r))) < tol,
            'left': float(np.max(np.abs(l @ J - lam * l))) < tol,
            'normalised': abs(float(l @ r) - 1.0) < tol,
            'perron': abs(max(np.linalg.eigvals(J).real) - lam) < tol,
        }


__all__ = ["CriticalAmplitude"]
