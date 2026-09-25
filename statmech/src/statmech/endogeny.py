"""The endogeny test of Aldous and Bandyopadhyay, run as two coupled populations.

A recursive distributional equation X = g(xi; X_1, ..., X_N) has an
*endogenous* solution when X is a function of the innovations xi on the tree
below it and of nothing else.  Theorem 11(c) of Aldous and Bandyopadhyay
(Ann. Appl. Probab. 15, 1047, 2005) makes that computable: the solution mu is
endogenous if and only if the bivariate recursion

    (X, X') = ( g(xi; X_1, ..., X_N),  g(xi; X'_1, ..., X'_N) ),

started from independent copies mu (x) mu, converges to the diagonal X = X'.
Same innovations, independent leaves.  In cavity language the innovations
are the structure -- which complexes a node is in, which members a complex
has, which couplings -- and the two copies are two sets of cavity fields on
the same structure.  Chapter 20 of the book identifies endogeny with the
replica-symmetric assumption and its linearisation with the de
Almeida--Thouless line; this module runs the nonlinear test.

Two populations, one per model.  Each carries two copies of every message,
draws the structure once per update and applies it to both copies, and
reports the distance between the copies sweep by sweep.  The distance is
``mean |h - h'| / mean |h|``; the test is whether it goes to zero.

``BivariateSpinGlass``
    Ising with couplings +-J on a random graph with a given chy-degree
    distribution (Poisson or regular), in an external field.  The couplings
    are innovations, so the two copies see the same sign on every edge.  At
    zero field the paramagnetic solution is trivially endogenous and the
    spin-glass transition is at <kbar> tanh^2(beta J) = 1; in a field the
    replica-symmetric solution exists at every temperature and the de
    Almeida--Thouless line is where its linearised distance stops
    contracting.  ``at_growth`` evaluates that linear rate on the population.

``BivariateHittingSet``
    The soft-field population of :mod:`statmech.softfield`, doubled.  At
    cardinality two this is vertex cover, whose replica symmetry breaks at
    <k> = e; above two the book has two criteria that disagree, the
    hard-field line <k>(c-1) = e and the soft-field entropy, and the test is
    the third.
"""

import numpy as np

from statmech.softfield import HittingSetBP


def distance(a, b):
    """``mean |a - b| / mean |a|``, the normalised distance between copies."""
    return float(np.mean(np.abs(a - b)) / max(np.mean(np.abs(a)), 1e-300))


# ---------------------------------------------------------------------------
# Ising spin glass on a graph, as the anchor
# ---------------------------------------------------------------------------

class BivariateSpinGlass:
    """Two copies of the cavity-field population of an Ising model with
    couplings ``+-J`` on a random graph, in a field ``H``.

    Args:
        mean: mean degree (Poisson), or the degree when ``regular``.
        beta_J: ``beta J``.
        field: ``beta H``.
        p_plus: probability that a coupling is ``+J``.
    """

    def __init__(self, mean, beta_J, field=0.0, p_plus=0.5, regular=False,
                 size=100_000, seed=0):
        self.mean, self.t = float(mean), float(np.tanh(beta_J))
        self.field, self.p_plus, self.regular = float(field), float(p_plus), regular
        self.size, self.rng = size, np.random.default_rng(seed)
        self.h = None
        self.h2 = None
        self.history = []

    def _excess(self, n):
        if self.regular:
            return np.full(n, int(round(self.mean)) - 1)
        return self.rng.poisson(self.mean, n)

    def _sum_u(self, h, k, idx, sign, ends):
        # u = atanh( s t tanh h ) summed over the k drawn neighbours
        u = np.arctanh(sign * self.t * np.tanh(h[idx]))
        csum = np.concatenate(([0.0], np.cumsum(u)))
        return csum[ends] - csum[ends - k]

    def _draw(self):
        k = self._excess(self.size)
        total = int(k.sum())
        idx = self.rng.integers(0, self.size, total)
        sign = np.where(self.rng.random(total) < self.p_plus, 1.0, -1.0)
        ends = np.cumsum(k)
        return k, idx, sign, ends

    def run_single(self, sweeps=300):
        """Bring one copy to its fixed point from a random start."""
        self.h = self.rng.normal(0.0, 1.0, self.size) + self.field
        for _ in range(sweeps):
            k, idx, sign, ends = self._draw()
            self.h = self.field + self._sum_u(self.h, k, idx, sign, ends)
        return self

    def start_pair(self, perturb=None):
        """Second copy: an independent draw from the same law (a shuffle),
        or, with ``perturb``, the first copy plus noise of that size."""
        if perturb is None:
            self.h2 = self.rng.permutation(self.h)
        else:
            self.h2 = self.h + self.rng.normal(0.0, perturb, self.size)
        self.history = [distance(self.h, self.h2)]
        return self

    def sweep_pair(self):
        k, idx, sign, ends = self._draw()
        self.h = self.field + self._sum_u(self.h, k, idx, sign, ends)
        self.h2 = self.field + self._sum_u(self.h2, k, idx, sign, ends)
        self.history.append(distance(self.h, self.h2))
        return self.history[-1]

    def run_pair(self, sweeps=300, perturb=None):
        self.start_pair(perturb)
        for _ in range(sweeps):
            self.sweep_pair()
        return self

    def at_growth(self):
        """``<kbar> E[(du/dh)^2]`` on the current copy: the linearised growth
        rate of the mean-square distance, one where the de Almeida--Thouless
        line is crossed."""
        th = np.tanh(self.h)
        du = self.t * (1 - th ** 2) / (1 - self.t ** 2 * th ** 2)
        kbar = self.mean if not self.regular else self.mean - 1
        return float(kbar * np.mean(du ** 2))

    def magnetisation(self):
        return float(np.mean(np.tanh(self.h)))

    def overlap(self):
        return float(np.mean(np.tanh(self.h) ** 2))


# ---------------------------------------------------------------------------
# Hitting set with soft fields, doubled
# ---------------------------------------------------------------------------

class BivariateHittingSet(HittingSetBP):
    """:class:`statmech.softfield.HittingSetBP` with a second copy of every
    message, updated with the same structure draws and the same damping
    mask.  ``run`` brings the first copy to its fixed point; ``start_pair``
    makes the second an independent draw from the same law; ``sweep_pair``
    advances both and records their distance."""

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.P2 = None
        self.Q2 = None
        self.history = []

    # the structure draws, taken once and applied to both copies

    def _draws(self, n, exclude=None):
        out = []
        for l in range(self.L):
            d = self._draw(n, l, excess=(l == exclude))
            total = int(d.sum())
            idx = self.rng.integers(0, self.size, total) if total else None
            out.append((d, idx))
        return out

    @staticmethod
    def _sum_with(Q, draws):
        n = len(draws[0][0])
        out = np.zeros(n)
        for l, (d, idx) in enumerate(draws):
            if idx is None:
                continue
            draws_l = Q[l][idx]
            ends = np.cumsum(d)
            csum = np.concatenate(([0.0], np.cumsum(draws_l)))
            out += csum[ends] - csum[ends - d]
        return out

    def start_pair(self, perturb=None):
        if perturb is None:
            self.P2 = [self.rng.permutation(p) for p in self.P]
            self.Q2 = [self.rng.permutation(q) for q in self.Q]
        else:
            self.P2 = [p + self.rng.normal(0.0, perturb, self.size) for p in self.P]
            self.Q2 = [q + self.rng.normal(0.0, perturb, self.size) for q in self.Q]
        self.history = [self.pair_distance()]
        return self

    def pair_distance(self):
        return max(distance(p, q) for p, q in zip(self.P, self.P2))

    def sweep_pair(self):
        newQ, newQ2 = [], []
        for l in range(self.L):
            idx = self.rng.integers(0, self.size,
                                    (self.size, int(self.c[l]) - 1))
            newQ.append(self._emit(self.P[l][idx]))
            newQ2.append(self._emit(self.P2[l][idx]))
        keep = [self.rng.random(self.size) >= self.damping for _ in range(self.L)]
        self.Q = [np.where(k, o, q) for k, o, q in zip(keep, self.Q, newQ)]
        self.Q2 = [np.where(k, o, q) for k, o, q in zip(keep, self.Q2, newQ2)]
        newP, newP2 = [], []
        for l in range(self.L):
            draws = self._draws(self.size, exclude=l)
            newP.append(-self.mu + self._sum_with(self.Q, draws))
            newP2.append(-self.mu + self._sum_with(self.Q2, draws))
        keep = [self.rng.random(self.size) >= self.damping for _ in range(self.L)]
        self.P = [np.where(k, o, q) for k, o, q in zip(keep, self.P, newP)]
        self.P2 = [np.where(k, o, q) for k, o, q in zip(keep, self.P2, newP2)]
        self.history.append(self.pair_distance())
        return self.history[-1]

    def run_pair(self, sweeps=300, perturb=None):
        self.start_pair(perturb)
        for _ in range(sweeps):
            self.sweep_pair()
        return self


def endogeny_boundary(factory, lo, hi, sweeps=300, tol=1e-3, iters=12,
                      equilibrate=300, verbose=False):
    """Bisect a one-parameter family on whether the pair converges.

    ``factory(x)`` returns a bivariate population; the pair is judged to
    converge when its final distance is below ``tol``.  Returns the
    boundary ``x`` between a converging ``lo`` and a non-converging ``hi``.
    """
    def converges(x):
        m = factory(x)
        m.run(equilibrate) if hasattr(m, 'run') else m.run_single(equilibrate)
        m.run_pair(sweeps)
        d = m.history[-1]
        if verbose:
            print(f'    x = {x:.5f}: distance {d:.2e}')
        return d < tol

    if not converges(lo) or converges(hi):
        raise ValueError('bracket does not straddle the boundary')
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if converges(mid):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


__all__ = ['BivariateSpinGlass', 'BivariateHittingSet', 'endogeny_boundary',
           'distance']
