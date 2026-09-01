"""Population dynamics for chygraph cavity messages (WP4).

WP1 replaced the scalar percolation message by a *number*, the mean derivative
``<u'>`` of the field update, and kept everything else — so it linearises the
map and returns thresholds in closed form.  It never represents the field
distribution itself.

Here the message is the distribution.  ``P_l(h)`` is the law of the cavity field
a node sends up into a layer-``l`` complex, ``Q_l(u)`` the law of the field the
complex sends back down, and each is carried as a population of samples:

    h = sum over the node's *other* complexes of their down-fields    (up step)
    u = the field complex ``a`` emits onto ``i`` given the fields of its
        other ``c-1`` members                                       (down step)

The up step is the chy-degree step of any chygraph, with a convolution in place
of a generating-function product.  The down step is the exact solve *inside* the
complex from :mod:`chygraph_statmech.cavity`, evaluated numerically rather than
symbolically.  That is the substitution the README states as the whole content
of the extension: **argument scalar -> argument measure, product -> convolution.**

The point of this module is that it shares no code path with WP1.  WP1 gets
``T_c`` by symbolic linearisation; this solves the full non-linear stochastic
problem and reads ``T_c`` off the order parameter.  They agree, including on the
triangle result — which is WP1's one genuinely new number.
"""

import numpy as np
from scipy.special import logsumexp


class CavityPopulation:
    """Finite-temperature Ising on a chygraph, by population dynamics.

    Args:
        cardinalities: clique size of each complex layer.  ``[2]`` is an
            ordinary graph, ``[2, 3]`` a graph with links and triangles.
        means: Poisson mean chy-degree per layer -- how many layer-``l``
            complexes a node belongs to.  Poisson is its own excess, so the
            cavity and full distributions coincide.
        beta_J: inverse temperature times coupling, one per layer or a scalar.
        size: population size.
        seed: RNG seed.
    """

    def __init__(self, cardinalities, means, beta_J, size=100_000, seed=0):
        self.c = np.asarray(cardinalities, dtype=int)
        self.L = len(self.c)
        self.means = np.asarray(means, dtype=float)
        bj = np.asarray(beta_J, dtype=float)
        self.beta_J = np.full(self.L, float(bj)) if bj.ndim == 0 else bj
        self.size = size
        self.rng = np.random.default_rng(seed)
        # enumerate the +-1 configurations of the other c-1 members, per layer
        self._cfg = []
        for c in self.c:
            m = c - 1
            bits = ((np.arange(2 ** m)[:, None] >> np.arange(m)) & 1)
            self._cfg.append(1.0 - 2.0 * bits)          # (2^m, m)
        self.P = [None] * self.L                        # node -> complex
        self.Q = [None] * self.L                        # complex -> node

    # -- the intra-complex step --------------------------------------------

    def emitted(self, h, l):
        """Field emitted onto one member of a layer-``l`` clique.

        ``h`` is ``(n, c-1)``, the cavity fields of the other members.  Exact
        enumeration inside the complex:

            u = (1/2) ln [ Z(sigma_0 = +1) / Z(sigma_0 = -1) ]

        with the clique energy split as ``bJ (S2 + sigma_0 S1)`` where ``S1`` is
        the magnetisation of the other members and ``S2`` their internal pair
        sum, so the whole sum is a matrix product plus two log-sum-exps.
        """
        cfg, bJ, c = self._cfg[l], self.beta_J[l], int(self.c[l])
        S1 = cfg.sum(axis=1)                            # (M,)
        S2 = (S1 ** 2 - (c - 1)) / 2.0
        base = h @ cfg.T + bJ * S2                      # (n, M)
        return 0.5 * (logsumexp(base + bJ * S1, axis=1)
                      - logsumexp(base - bJ * S1, axis=1))

    # -- the chy-degree step ------------------------------------------------

    def _sum_over_complexes(self, n):
        """Sample ``sum_l sum_{d_l} u`` with ``d_l ~ Poisson(means[l])``.

        The convolution that replaces the generating-function product.
        """
        out = np.zeros(n)
        for k in range(self.L):
            d = self.rng.poisson(self.means[k], n)
            total = int(d.sum())
            if total == 0:
                continue
            draws = self.Q[k][self.rng.integers(0, self.size, total)]
            ends = np.cumsum(d)
            starts = ends - d
            csum = np.concatenate(([0.0], np.cumsum(draws)))
            out += csum[ends] - csum[starts]
        return out

    # -- iteration ----------------------------------------------------------

    def initialise(self, field=1.0):
        """Start fully magnetised, so the ferromagnetic branch is found if it
        exists.  Starting at zero field would sit on the paramagnetic fixed
        point forever, since it is exact at zero external field."""
        for l in range(self.L):
            self.P[l] = np.full(self.size, float(field))
            self.Q[l] = np.full(self.size, float(field))

    def sweep(self):
        """One parallel update of every population."""
        newQ = []
        for l in range(self.L):
            idx = self.rng.integers(0, self.size,
                                    (self.size, int(self.c[l]) - 1))
            newQ.append(self.emitted(self.P[l][idx], l))
        self.Q = newQ
        self.P = [self._sum_over_complexes(self.size) for _ in range(self.L)]

    def run(self, sweeps=300, field=1.0):
        self.initialise(field)
        for _ in range(sweeps):
            self.sweep()
        return self

    # -- observables --------------------------------------------------------

    def magnetisation(self, samples=None):
        """``<tanh h>`` over the *full* field, all complexes included."""
        n = samples or self.size
        return float(np.mean(np.tanh(self._sum_over_complexes(n))))


def critical_coupling(cardinalities, means, lo=1e-3, hi=3.0, iters=22,
                      sweeps=250, size=60_000, seed=0, tol=1e-3):
    """``beta J`` at which the magnetisation first survives, by bisection.

    Noisy by construction — population dynamics has sampling error and critical
    slowing down — so this is good to a percent or so, which is enough to test
    a closed form derived by an entirely different route.
    """
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        m = CavityPopulation(cardinalities, means, mid, size=size,
                             seed=seed).run(sweeps).magnetisation()
        lo, hi = ((lo, mid) if m > tol else (mid, hi))
    return 0.5 * (lo + hi)
